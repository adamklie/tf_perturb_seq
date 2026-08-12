#!/usr/bin/env python3
"""
nontargeting_trans_analysis.py

Investigates significant trans hits from non-targeting (NT) guides.

Non-targeting guides should, by definition, not perturb any gene. Any
significant trans associations they show reflect:
  - Residual false positives after BH correction
  - Systematic technical effects (cell stress, batch, etc.)
  - Genes that are genuinely sensitive to off-target effects

This script identifies which genes are hit by NT guides and asks whether
they overlap with the most frequently trans-regulated genes identified from
targeting guides.

Analyses:
  1. BH correction over the full NT guide trans universe
  2. Histogram: distribution of how many NT guides significantly hit each gene
  3. Barplot: top 30 genes most frequently hit by NT guides
  4. Overlap: are the top NT-hit genes also the top targeting-hit genes?

Data source:
  trans_per_guide_results.tsv.gz -- gene_id | guide_id | log2_fc | log2_fc_std | p_value
  Guide type determined by guide_id prefix:
    non-targeting guides: id starts with "non-targeting"
    targeting guides:     all others (excluding positive/negative control)

Usage:
    python nontargeting_trans_analysis.py \\
        --output_dir /path/to/output_dir
"""

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# --- Paths --------------------------------------------------------------------
RESULTS_DIR = (
    "/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/"
    "helen_output_poolabcd_prod_v9"
)
PIPELINE_OUTPUTS_DIR = (
    "/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/"
    "helen_output_poolabcd_prod_v9/pipeline_outputs"
)
# Calibrated trans results (empirical p-values from NT distribution):
#   element_id | element_symbol | element_label | tested_gene_id | tested_gene_symbol
#   | n_cells | log2fc | log2fc_se | is_cis | is_direct_target
#   | posterior_pval | empirical_pval | empirical_pval_adj
CALIB_TRANS_ELEM = os.path.join(
    PIPELINE_OUTPUTS_DIR, "gersbach_ihep_calibrated_trans_results.tsv"
)
# Symbol reference for gene ID -> symbol mapping
ENSG_TO_SYM = "/hpc/home/seg95/lab-storage/ref/ensembl_to_symbol.tsv"


# --- CLI ----------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--output_dir", required=True)
    p.add_argument("--fdr_thresh", type=float, default=0.05)
    p.add_argument("--top_n",      type=int,   default=30,
                   help="Number of top genes to show in barplot (default: 30)")
    p.add_argument("--top_n_targeting", type=int, default=100,
                   help="Top N targeting-hit genes to compare against for overlap "
                        "(default: 100)")
    return p.parse_args()


# --- Helpers ------------------------------------------------------------------
def build_ensg_to_symbol():
    """Load Ensembl ID -> gene symbol map from reference TSV."""
    ensg2sym = {}
    if not os.path.exists(ENSG_TO_SYM):
        print(f"  WARNING: symbol reference not found at {ENSG_TO_SYM}")
        return ensg2sym
    ref = pd.read_csv(ENSG_TO_SYM, sep="\t", header=None,
                      names=["versioned_id", "symbol"])
    ref["ensg"] = ref["versioned_id"].str.split(".").str[0]
    ensg2sym = dict(zip(ref["ensg"], ref["symbol"]))
    print(f"  Loaded {len(ensg2sym):,} Ensembl->symbol mappings")
    return ensg2sym


def bh_correct(pvals):
    """BH FDR correction over a 1-D array of p-values."""
    pvals  = np.asarray(pvals, dtype=float)
    finite = np.isfinite(pvals)
    padj   = np.full(len(pvals), np.nan)
    pf     = pvals[finite]
    m      = finite.sum()
    if m == 0:
        return padj
    try:
        from statsmodels.stats.multitest import multipletests
        _, adj, _, _ = multipletests(pf, method="fdr_bh")
    except ImportError:
        from scipy.stats import rankdata
        ranks = rankdata(pf, method="ordinal")
        adj   = np.minimum(
            np.minimum.accumulate((pf * m / ranks)[::-1])[::-1], 1.0
        )
    padj[finite] = adj
    return padj


def is_nontargeting(guide_id):
    """Return True if guide_id belongs to a non-targeting guide."""
    return str(guide_id).lower().startswith("non-targeting")


# --- Main ---------------------------------------------------------------------
def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Output directory : {args.output_dir}")
    print(f"FDR threshold    : {args.fdr_thresh}")

    # --- Load symbol map -------------------------------------------------------
    ensg2sym = build_ensg_to_symbol()

    # --- Load calibrated trans element file -----------------------------------
    print(f"\nLoading calibrated trans results: {CALIB_TRANS_ELEM}")
    df = pd.read_csv(CALIB_TRANS_ELEM, sep="\t")
    print(f"  {len(df):,} rows, cols: {list(df.columns)}")

    type_counts = df["element_label"].value_counts()
    print(f"  Element label distribution:")
    for t, n in type_counts.items():
        print(f"    {t}: {n:,} rows")

    # --- Split into NT and targeting ------------------------------------------
    # element_label values: "targeting", "non-targeting", "negative control",
    # "positive control". empirical_pval_adj is already BH-corrected.
    nt_df  = df[df["element_label"] == "non-targeting"].copy()
    tgt_df = df[df["element_label"] == "targeting"].copy()
    print(f"\n  NT rows:        {len(nt_df):,}")
    print(f"  Targeting rows: {len(tgt_df):,}")

    # empirical_pval_adj is pre-corrected -- no BH correction needed
    n_sig = (nt_df["empirical_pval_adj"] < args.fdr_thresh).sum()
    print(f"\n  {n_sig:,} significant NT trans pairs "
          f"(empirical_pval_adj < {args.fdr_thresh}, pre-corrected)")

    sig_nt = nt_df[nt_df["empirical_pval_adj"] < args.fdr_thresh].copy()

    if len(sig_nt) == 0:
        print("  No significant NT trans hits -- nothing to plot")
        return

    # Use tested_gene_symbol directly (available in calibrated file)
    sig_nt["gene_symbol"] = sig_nt["tested_gene_symbol"].fillna(
        sig_nt["tested_gene_id"].map(ensg2sym).fillna(sig_nt["tested_gene_id"])
    )

    # --- Analysis 1: How many NT elements hit each gene? ----------------------
    # Count: per tested_gene_id, how many distinct NT elements are significant
    hits_per_gene = (
        sig_nt.groupby(["tested_gene_id", "gene_symbol"])["element_id"]
        .nunique()
        .reset_index()
        .rename(columns={"element_id": "n_nt_elements_sig",
                         "tested_gene_id": "gene_id"})
        .sort_values("n_nt_elements_sig", ascending=False)
    )
    print(f"\n  Genes hit by >=1 NT element: {len(hits_per_gene):,}")
    print(f"  Top 10 genes by NT element hits:")
    print(hits_per_gene.head(10)[["gene_symbol", "n_nt_elements_sig"]].to_string(index=False))

    csv_path = os.path.join(args.output_dir, "nontargeting_trans_hits_per_gene.csv")
    hits_per_gene.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    # --- Analysis 2: How many genes does each NT element hit? -----------------
    hits_per_element = (
        sig_nt.groupby("element_id")["tested_gene_id"]
        .nunique()
        .reset_index()
        .rename(columns={"tested_gene_id": "n_genes_sig"})
        .sort_values("n_genes_sig", ascending=False)
    )
    n_nt_elem_total   = nt_df["element_id"].nunique()
    n_nt_elem_any_hit = hits_per_element["element_id"].nunique()
    print(f"\n  NT elements with >=1 significant trans hit: "
          f"{n_nt_elem_any_hit:,} / {n_nt_elem_total:,}")
    print(f"  Top 10 NT elements by n_genes_sig:")
    # Annotate element symbol
    elem_sym = dict(zip(df["element_id"], df["element_symbol"]))
    hits_per_element["element_symbol"] = hits_per_element["element_id"].map(elem_sym)
    print(hits_per_element.head(10)[["element_symbol", "n_genes_sig"]].to_string(index=False))

    # --- Plot 1: Histogram of genes-per-NT-element ----------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    counts = hits_per_element["n_genes_sig"]
    max_bin = min(int(counts.max()) + 2, 200)
    ax.hist(counts, bins=range(0, max_bin), color="#2ca02c",
            edgecolor="white", linewidth=0.4)
    ax.set_xlabel("Number of genes significantly hit per NT element", fontsize=11)
    ax.set_ylabel("Number of NT elements", fontsize=11)
    ax.set_title(
        f"Distribution of significant trans hits per non-targeting element\n"
        f"(calibrated empirical_pval_adj < {args.fdr_thresh}, pre-BH-corrected)\n"
        f"{n_nt_elem_any_hit:,} / {n_nt_elem_total:,} NT elements have >=1 hit",
        fontsize=11, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    hist_path = os.path.join(args.output_dir, "nontargeting_trans_hist_per_element.png")
    plt.savefig(hist_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\n  Saved: {hist_path}")

    # --- Plot 2: Top N genes most frequently hit by NT guides -----------------
    top_genes = hits_per_gene.head(args.top_n).sort_values(
        "n_nt_elements_sig", ascending=True
    )
    fig, ax = plt.subplots(figsize=(9, max(5, len(top_genes) * 0.35)))
    bars = ax.barh(
        top_genes["gene_symbol"], top_genes["n_nt_elements_sig"],
        color="#2ca02c", edgecolor="white", linewidth=0.4,
    )
    ax.set_xlabel("Number of NT elements significantly hitting this gene", fontsize=11)
    ax.set_title(
        f"Top {args.top_n} genes most frequently hit by non-targeting guides\n"
        f"(BH FDR < {args.fdr_thresh})",
        fontsize=11, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    # Annotate bars with count
    for bar, val in zip(bars, top_genes["n_nt_elements_sig"]):
        ax.text(val + 0.1, bar.get_y() + bar.get_height() / 2,
                str(int(val)), va="center", fontsize=8)
    plt.tight_layout()
    bar_path = os.path.join(args.output_dir, "nontargeting_trans_top_genes.png")
    plt.savefig(bar_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {bar_path}")

    # --- Analysis 3: Overlap with top targeting trans hits --------------------
    # Use calibrated targeting data -- empirical_pval_adj already BH-corrected
    print(f"\nComputing overlap with top {args.top_n_targeting} targeting trans hits ...")
    sig_tgt = tgt_df[tgt_df["empirical_pval_adj"] < args.fdr_thresh].copy()
    print(f"  {len(sig_tgt):,} significant targeting trans pairs "
          f"(empirical_pval_adj < {args.fdr_thresh})")

    tgt_hits_per_gene = (
        sig_tgt.groupby("tested_gene_id")["element_id"]
        .nunique()
        .reset_index()
        .rename(columns={"element_id":      "n_targeting_elements_sig",
                         "tested_gene_id":  "gene_id"})
        .sort_values("n_targeting_elements_sig", ascending=False)
    )
    tgt_hits_per_gene["gene_symbol"] = (
        sig_tgt.groupby("tested_gene_id")["tested_gene_symbol"]
        .first().reset_index()
        .rename(columns={"tested_gene_id": "gene_id",
                         "tested_gene_symbol": "gene_symbol"})
        .set_index("gene_id")["gene_symbol"]
        .reindex(tgt_hits_per_gene["gene_id"])
        .fillna(tgt_hits_per_gene["gene_id"].map(ensg2sym))
        .values
    )

    top_tgt_genes = set(tgt_hits_per_gene.head(args.top_n_targeting)["gene_id"])
    top_nt_genes  = set(hits_per_gene["gene_id"])

    overlap = top_nt_genes & top_tgt_genes
    print(f"  NT-hit genes: {len(top_nt_genes):,}")
    print(f"  Top {args.top_n_targeting} targeting-hit genes: {len(top_tgt_genes):,}")
    print(f"  Overlap: {len(overlap):,} genes")

    if overlap:
        overlap_syms = sorted(
            ensg2sym.get(g, g) for g in overlap
        )
        print(f"  Overlapping genes: {', '.join(overlap_syms)}")

    # Save overlap table
    overlap_df = (
        hits_per_gene[hits_per_gene["gene_id"].isin(overlap)]
        .merge(
            tgt_hits_per_gene[["gene_id", "n_targeting_elements_sig"]],
            on="gene_id", how="left",
        )
        .sort_values("n_nt_elements_sig", ascending=False)
    )
    overlap_csv = os.path.join(args.output_dir, "nontargeting_vs_targeting_overlap.csv")
    overlap_df.to_csv(overlap_csv, index=False)
    print(f"  Saved overlap table: {overlap_csv}")

    # --- Plot 3: Overlap scatter -----------------------------------------------
    # x = n NT guides hitting gene, y = n targeting guides hitting gene
    merged = hits_per_gene.merge(
        tgt_hits_per_gene[["gene_id", "n_targeting_elements_sig"]],
        on="gene_id", how="left",
    ).fillna(0)

    fig, ax = plt.subplots(figsize=(8, 7))
    in_overlap = merged["gene_id"].isin(overlap)

    ax.scatter(
        merged.loc[~in_overlap, "n_nt_elements_sig"],
        merged.loc[~in_overlap, "n_targeting_elements_sig"],
        color="#aaaaaa", s=20, alpha=0.6, linewidths=0,
        rasterized=True, label="NT-hit genes (not in top targeting)",
    )
    ax.scatter(
        merged.loc[in_overlap, "n_nt_elements_sig"],
        merged.loc[in_overlap, "n_targeting_elements_sig"],
        color="#d62728", s=40, alpha=0.9, linewidths=0,
        rasterized=True,
        label=f"Overlap with top {args.top_n_targeting} targeting-hit genes",
        zorder=5,
    )
    # Label top overlapping genes
    for _, row in merged[in_overlap].nlargest(15, "n_nt_elements_sig").iterrows():
        ax.text(
            row["n_nt_elements_sig"] + 0.1, row["n_targeting_elements_sig"],
            ensg2sym.get(row["gene_id"], row["gene_id"]),
            fontsize=7, va="center",
        )

    ax.set_xlabel("NT guides significantly hitting gene", fontsize=11)
    ax.set_ylabel(f"Targeting elements significantly hitting gene", fontsize=11)
    ax.set_title(
        f"NT trans hits vs targeting trans hits per gene\n"
        f"(BH FDR < {args.fdr_thresh})",
        fontsize=11, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=9, frameon=False)
    plt.tight_layout()
    scatter_path = os.path.join(args.output_dir, "nontargeting_vs_targeting_scatter.png")
    plt.savefig(scatter_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {scatter_path}")

    print("\n=== Complete ===")


if __name__ == "__main__":
    main()