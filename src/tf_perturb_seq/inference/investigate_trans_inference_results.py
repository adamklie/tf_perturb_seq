#!/usr/bin/env python3
"""
investigate_trans_inference_results.py

Investigates trans inference results from SCEPTRE (and optionally Perturbo)
on a TF perturb-seq screen.

Trans definition:
  element files: gene_id != intended_target_name
  guide files:   guide target symbol (parsed from guide_id prefix) != mapped gene symbol

File format notes:
  trans_per_guide_results.tsv.gz (merged, Perturbo):
      gene_id | guide_id | log2_fc | log2_fc_std | p_value
  trans_per_element_results.tsv.gz (merged, Perturbo):
      gene_id | intended_target_name | intended_target_chr | intended_target_start
              | intended_target_end | log2_fc | log2_fc_std | p_value
  inference/sceptre_per_guide_output.tsv.gz (--sceptre_only):
      gene_id | guide_id | p_value | log2_fc
  inference/sceptre_per_element_output.tsv.gz (--sceptre_only):
      gene_id | intended_target_name | p_value | log2_fc

  Calibrated element files (*_calibrated_trans_results.tsv) are auto-detected
  from inference_dir by suffix and used in preference to raw element files.
  Guide-level analysis is skipped if no per-guide file is found.

Usage:
    python investigate_trans_inference_results.py \\
        --pipeline_outputs_dir /path \\
        --inference_dir /path/to/calibrated_outs \\
        --output_dir /path/to/output
"""

import argparse
import glob
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

warnings.filterwarnings("ignore")

# --- Default paths ------------------------------------------------------------
_DEFAULT_RESULTS_DIR          = (
    "/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/"
    "helen_output_poolabcd_prod_v9"
)
_DEFAULT_PIPELINE_OUTPUTS_DIR = _DEFAULT_RESULTS_DIR + "/pipeline_outputs"
_DEFAULT_INFERENCE_DIR        = _DEFAULT_RESULTS_DIR + "/inference"
_DEFAULT_GUIDE_META           = (
    "/hpc/group/gersbachlab/seg95/crispr-pipeline-personal/"
    "example-data/production_guide_metadata_v6.tsv"
)
_DEFAULT_MUDATA_NAME          = "inference_mudata.h5mu"

FDR_THRESH_DEFAULT = 0.05
TOP_N_DEFAULT      = 20


def _find_calibrated_file(directory: str, suffix: str) -> str | None:
    """Auto-detect a calibrated result file by suffix. Returns first match or None."""
    matches = glob.glob(os.path.join(directory, f"*{suffix}"))
    if not matches:
        return None
    if len(matches) > 1:
        print(f"  WARNING: multiple files match *{suffix} in {directory}; "
              f"using first: {os.path.basename(matches[0])}")
    return matches[0]


def _build_paths(args):
    global PIPELINE_OUTPUTS_DIR, INFERENCE_DIR
    global MERGED_TRANS_GUIDE, MERGED_TRANS_ELEM
    global CALIB_TRANS_ELEM
    global SCEPTRE_GUIDE, SCEPTRE_ELEM
    global MUDATA_PATH
    global GUIDE_METADATA_PATH

    PIPELINE_OUTPUTS_DIR = args.pipeline_outputs_dir
    INFERENCE_DIR        = args.inference_dir
    GUIDE_METADATA_PATH  = args.guide_metadata

    MERGED_TRANS_GUIDE = os.path.join(PIPELINE_OUTPUTS_DIR,
                                      "perturbo_trans_per_guide_output.tsv.gz")
    MERGED_TRANS_ELEM  = os.path.join(PIPELINE_OUTPUTS_DIR,
                                      "perturbo_trans_per_element_output.tsv.gz")
    SCEPTRE_GUIDE      = os.path.join(INFERENCE_DIR,
                                      "sceptre_per_guide_output.tsv.gz")
    SCEPTRE_ELEM       = os.path.join(INFERENCE_DIR,
                                      "sceptre_per_element_output.tsv.gz")
    MUDATA_PATH        = os.path.join(PIPELINE_OUTPUTS_DIR, args.mudata_name)

    # Auto-detect calibrated trans file from inference_dir
    CALIB_TRANS_ELEM = _find_calibrated_file(INFERENCE_DIR,
                                             "_calibrated_trans_results.tsv")
    if CALIB_TRANS_ELEM:
        print(f"  Auto-detected calibrated trans file: "
              f"{os.path.basename(CALIB_TRANS_ELEM)}")
    else:
        print(f"  WARNING: no *_calibrated_trans_results.tsv found in {INFERENCE_DIR}")

    # --calib_prefix is accepted for backwards compatibility but ignored
    if getattr(args, "calib_prefix", None):
        print("  NOTE: --calib_prefix is deprecated; calibrated file is now auto-detected.")


# --- CLI ----------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--output_dir", default="./figures",
                   help="Directory for output figures and tables (default: ./figures)")
    p.add_argument("--sceptre_only", action="store_true", default=False,
                   help="Use standalone SCEPTRE files instead of calibrated Perturbo.")
    p.add_argument("--fdr_thresh", type=float, default=FDR_THRESH_DEFAULT,
                   help=f"BH FDR threshold (default: {FDR_THRESH_DEFAULT})")
    p.add_argument("--top_n", type=int, default=TOP_N_DEFAULT,
                   help=f"Number of top trans hits to show (default: {TOP_N_DEFAULT})")
    p.add_argument("--keep_nontargeting", action="store_true", default=False,
                   help="Include non-targeting control guides/elements.")
    p.add_argument("--pipeline_outputs_dir", default=_DEFAULT_PIPELINE_OUTPUTS_DIR,
                   help="Directory containing per-guide .tsv.gz and MuData files.")
    p.add_argument("--inference_dir", default=_DEFAULT_INFERENCE_DIR,
                   help="Directory containing SCEPTRE files and/or calibrated TSVs "
                        "(auto-detected by suffix).")
    p.add_argument("--mudata_name", default=_DEFAULT_MUDATA_NAME,
                   help=f"MuData filename inside pipeline_outputs_dir "
                        f"(default: {_DEFAULT_MUDATA_NAME})")
    p.add_argument("--guide_metadata", default=_DEFAULT_GUIDE_META,
                   help="Path to guide metadata TSV.")
    # Kept for backwards compatibility; ignored
    p.add_argument("--calib_prefix", default=None,
                   help="[DEPRECATED] Calibrated file is now auto-detected; ignored.")
    return p.parse_args()


# --- Helpers ------------------------------------------------------------------
def load_tsv(path, label):
    print(f"  Loading {label}:\n    {path}")
    df = pd.read_csv(path, sep="\t", compression="gzip")
    n_raw = len(df)
    df = df.drop_duplicates()
    print(f"    {n_raw:,} rows -> {len(df):,} after drop_duplicates  |  "
          f"cols: {list(df.columns)}")
    return df


def parse_target_from_guide_id(guide_id):
    return str(guide_id).split("#")[0]


def build_ensg_to_symbol_mygene(ensg_ids, verbose=True) -> dict:
    ref_path = os.path.expanduser("~/lab-storage/ref/ensembl_to_symbol.tsv")
    if not os.path.exists(ref_path):
        print(f"  WARNING: {ref_path} not found -- labels will be Ensembl IDs")
        return {}
    ref = pd.read_csv(ref_path, sep="\t", header=None,
                      names=["ensembl_versioned", "symbol"])
    ref["ensembl_id"] = ref["ensembl_versioned"].str.split(".").str[0]
    ref = ref.drop_duplicates("ensembl_id")
    full_map = dict(zip(ref["ensembl_id"], ref["symbol"]))
    query = set(str(x) for x in ensg_ids if str(x).startswith("ENSG"))
    mapping = {k: v for k, v in full_map.items() if k in query}
    if verbose:
        n_mapped   = len(mapping)
        n_unmapped = len(query) - n_mapped
        print(f"  Mapped {n_mapped:,} / {len(query):,} Ensembl IDs from {ref_path} "
              f"({n_unmapped:,} unmapped -> will use Ensembl ID as label)")
    return mapping


def build_ensg_to_symbol(guide_df):
    ensg_ids = guide_df["gene_id"].unique().tolist()
    return build_ensg_to_symbol_mygene(ensg_ids)


def bh_correct(df, p_cols):
    from scipy.stats import rankdata

    def _bh(pvals):
        n      = len(pvals)
        finite = np.isfinite(pvals)
        padj   = np.full(n, np.nan)
        pf     = pvals[finite]
        m      = finite.sum()
        ranks  = rankdata(pf, method="ordinal")
        raw    = pf * m / ranks
        padj_f = np.minimum.accumulate(raw[::-1])[::-1]
        padj_f = np.minimum(padj_f, 1.0)
        padj[finite] = padj_f
        return padj

    try:
        from statsmodels.stats.multitest import multipletests
        def _correct(pvals):
            finite = np.isfinite(pvals)
            out    = np.full(len(pvals), np.nan)
            if finite.sum() > 0:
                _, padj, _, _ = multipletests(pvals[finite], method="fdr_bh")
                out[finite]   = padj
            return out
    except ImportError:
        print("    statsmodels not found -- using built-in BH implementation")
        _correct = _bh

    df = df.copy()
    for col in p_cols:
        if col not in df.columns:
            continue
        bh_col = col + "_bh"
        df[bh_col] = _correct(df[col].values.astype(float))
        n_sig_raw = (df[col]    < 0.05).sum()
        n_sig_bh  = (df[bh_col] < 0.05).sum()
        print(f"    {col}: {n_sig_raw:,} raw p<0.05  ->  {n_sig_bh:,} "
              f"BH-adjusted FDR<0.05  (over {len(df):,} tests)")
    return df


def split_trans_guide(df):
    ensg2sym = build_ensg_to_symbol(df)
    df = df.copy()
    df["target_symbol"] = df["guide_id"].apply(parse_target_from_guide_id)
    df["gene_symbol"]   = df["gene_id"].map(ensg2sym).fillna(df["gene_id"])
    is_cis  = df["target_symbol"].str.upper() == df["gene_symbol"].str.upper()
    return df[~is_cis].copy()


def split_trans_element(df):
    is_cis = df["gene_id"] == df["intended_target_name"]
    return df[~is_cis].copy()


def is_nontargeting(series):
    def _all_nontargeting(val):
        if pd.isna(val):
            return False
        parts = str(val).split("&")
        return all(p.strip().lower().startswith("non-targeting") for p in parts)
    return series.map(_all_nontargeting)


def filter_nontargeting(df, cols, label):
    mask = pd.Series(False, index=df.index)
    for col in cols:
        if col in df.columns:
            mask |= is_nontargeting(df[col])
    n_removed = mask.sum()
    df_filtered = df[~mask].copy()
    print(f"  [{label}] Removed {n_removed:,} non-targeting rows "
          f"({100.0 * n_removed / max(len(df), 1):.1f}%)  "
          f"-> {len(df_filtered):,} rows remaining")
    return df_filtered


def print_trans_stats(label, summary, fdr_thresh, level="guide"):
    avg_tested = summary["n_trans_tested"].mean()
    avg_sig    = summary["n_trans_sig"].mean()
    pct_sig    = summary["pct_trans_sig"].mean()
    pct_any    = (summary["n_trans_sig"] > 0).mean() * 100
    print(f"\n  [{label}] BH FDR < {fdr_thresh}  (per-{level})")
    print(f"    Unique {level}s:                        {len(summary):,}")
    print(f"    Avg trans genes tested per {level}:    {avg_tested:.1f}")
    print(f"    Avg trans sig per {level}:             {avg_sig:.2f}")
    print(f"    Avg % of tested that are sig:          {pct_sig:.2f}%")
    print(f"    {level.capitalize()}s with >=1 trans sig hit:   {pct_any:.1f}%")
    print(f"\n    Top 10 {level}s by n_trans_sig:")
    print(summary.nlargest(10, "n_trans_sig").to_string(index=False))


# --- Analysis 1: Per-guide trans statistics -----------------------------------
def per_guide_trans(args, out_dir, guide_df, method_specs):
    """
    Per-guide trans statistics.
    Skipped automatically if guide_df is None (no per-guide file available).

    Intermediate TSV outputs (for R/ggplot):
        trans_per_guide_{method}.tsv  -- one row per guide with n_trans_tested,
                                         n_trans_sig, pct_trans_sig
        trans_per_guide_all_methods.tsv -- long-format combined
    """
    print("\n=== Analysis 1: Per-guide trans statistics ===")

    if guide_df is None:
        print("  Skipping: no per-guide file available for this dataset.")
        return {}

    trans = split_trans_guide(guide_df)
    print(f"  Trans guide-gene pairs: {len(trans):,}")

    all_summaries = {}
    for method, p_col, padj_col, lfc_col in method_specs:
        if padj_col not in trans.columns:
            print(f"  WARNING: {padj_col} not found, skipping {method}")
            continue

        df_m = trans[["guide_id", "target_symbol", "gene_id", "gene_symbol",
                       p_col, padj_col]].dropna(subset=[padj_col]).copy()
        df_m = (
            df_m.groupby(["guide_id", "target_symbol", "gene_id", "gene_symbol"],
                         as_index=False)
                .agg({p_col: "min", padj_col: "min"})
        )
        df_m["significant"] = df_m[padj_col] < args.fdr_thresh

        summary = (
            df_m.groupby(["guide_id", "target_symbol"])
                .agg(
                    n_trans_tested=("gene_id",     "nunique"),
                    n_trans_sig   =("significant", "sum"),
                )
                .reset_index()
        )
        summary["pct_trans_sig"] = 100.0 * summary["n_trans_sig"] / summary["n_trans_tested"]
        summary["method"]        = method
        summary["fdr_thresh"]    = args.fdr_thresh
        all_summaries[method]    = summary

        print_trans_stats(method, summary, args.fdr_thresh, level="guide")

        tsv_path = os.path.join(out_dir, f"trans_per_guide_{method.lower()}.tsv")
        summary.to_csv(tsv_path, sep="\t", index=False)
        print(f"  Saved R-ready table: {tsv_path}")

    if all_summaries:
        combined = pd.concat(all_summaries.values(), ignore_index=True)
        combined_path = os.path.join(out_dir, "trans_per_guide_all_methods.tsv")
        combined.to_csv(combined_path, sep="\t", index=False)
        print(f"  Saved R-ready combined table: {combined_path}")

    _plot_trans_distribution(
        all_summaries,
        id_col="guide_id",
        level="guide",
        fdr_thresh=args.fdr_thresh,
        out_path=os.path.join(out_dir, "trans_per_guide_distribution.png"),
    )
    return all_summaries


# --- Analysis 2: Per-element trans statistics ---------------------------------
def per_element_trans(args, out_dir, elem_df, method_specs_elem):
    """
    Per-element trans statistics.

    Intermediate TSV outputs (for R/ggplot):
        trans_per_element_{method}.tsv  -- one row per element
        trans_per_element_all_methods.tsv -- long-format combined
        trans_significant_all_methods.tsv -- all significant trans pairs
    """
    print("\n=== Analysis 2: Per-element trans statistics ===")

    trans = split_trans_element(elem_df)
    print(f"  Trans element-gene pairs: {len(trans):,}")

    trans = trans.copy()
    if "_ensg2sym" in elem_df.attrs:
        ensg2sym = elem_df.attrs["_ensg2sym"]
        trans["target_symbol"] = (
            trans["intended_target_name"].map(ensg2sym)
                                         .fillna(trans["intended_target_name"])
        )
    else:
        trans["target_symbol"] = trans["intended_target_name"]

    all_summaries = {}
    all_sig_frames = []

    for method, p_col, padj_col, lfc_col in method_specs_elem:
        if padj_col not in trans.columns:
            print(f"  WARNING: {padj_col} not found, skipping {method}")
            continue

        df_m = trans[["intended_target_name", "target_symbol", "gene_id",
                       p_col, padj_col]].dropna(subset=[padj_col]).copy()
        df_m = (
            df_m.groupby(["intended_target_name", "target_symbol", "gene_id"],
                         as_index=False)
                .agg({p_col: "min", padj_col: "min"})
        )
        df_m["significant"] = df_m[padj_col] < args.fdr_thresh

        summary = (
            df_m.groupby(["intended_target_name", "target_symbol"])
                .agg(
                    n_trans_tested=("gene_id",     "nunique"),
                    n_trans_sig   =("significant", "sum"),
                )
                .reset_index()
        )
        summary["pct_trans_sig"] = 100.0 * summary["n_trans_sig"] / summary["n_trans_tested"]
        summary["method"]        = method
        summary["fdr_thresh"]    = args.fdr_thresh
        all_summaries[method]    = summary

        print_trans_stats(method, summary, args.fdr_thresh, level="element")

        tsv_path = os.path.join(out_dir, f"trans_per_element_{method.lower()}.tsv")
        summary.to_csv(tsv_path, sep="\t", index=False)
        print(f"  Saved R-ready table: {tsv_path}")

        # Full significant pairs for R
        sig_all = df_m[df_m["significant"]].copy()
        sig_all["neg_log10_fdr"] = -np.log10(sig_all[padj_col].clip(lower=1e-300))
        sig_all["method"] = method

        # Attach gene_symbol if available
        if "_ensg2sym" in elem_df.attrs:
            sig_all["gene_symbol"] = sig_all["gene_id"].map(elem_df.attrs["_ensg2sym"]).fillna(sig_all["gene_id"])
        elif "gene_symbol" in trans.columns:
            sig_all = sig_all.merge(
                trans[["gene_id", "gene_symbol"]].drop_duplicates("gene_id"),
                on="gene_id", how="left"
            )

        # Include extra calibrated columns where present
        base_cols  = ["method", "intended_target_name", "target_symbol",
                      "gene_id", p_col, padj_col, "neg_log10_fdr"]
        extra_cols = [c for c in ["gene_symbol", "n_cells", "log2fc_se",
                                  "posterior_pval", "is_cis", "is_direct_target"]
                      if c in sig_all.columns]
        sig_out = sig_all[[c for c in base_cols + extra_cols if c in sig_all.columns]].copy()
        sig_out = sig_out.rename(columns={p_col: "p_value", padj_col: "padj",
                                          lfc_col: "log2fc"})
        all_sig_frames.append(sig_out)

    if all_summaries:
        combined = pd.concat(all_summaries.values(), ignore_index=True)
        combined_path = os.path.join(out_dir, "trans_per_element_all_methods.tsv")
        combined.to_csv(combined_path, sep="\t", index=False)
        print(f"  Saved R-ready combined table: {combined_path}")

    if all_sig_frames:
        all_sig = pd.concat(all_sig_frames, ignore_index=True)
        sig_path = os.path.join(out_dir, "trans_significant_all_methods.tsv")
        all_sig.to_csv(sig_path, sep="\t", index=False)
        print(f"  Saved R-ready significant pairs: {sig_path}")

    _plot_trans_distribution(
        all_summaries,
        id_col="intended_target_name",
        level="element",
        fdr_thresh=args.fdr_thresh,
        out_path=os.path.join(out_dir, "trans_per_element_distribution.png"),
    )
    return all_summaries


# --- Analysis 3: Top trans hits -----------------------------------------------
def top_trans_hits(args, out_dir, guide_df, elem_df, method_specs, method_specs_elem):
    """
    Top perturbations and top trans-regulated genes.

    Intermediate TSV outputs (for R/ggplot):
        trans_top_elements_{method}.tsv  -- top N elements by sig trans hit count
        trans_top_genes_{method}.tsv     -- top N trans genes by breadth
    """
    print("\n=== Analysis 3: Top trans hits ===")

    trans_elem = split_trans_element(elem_df)
    if "_ensg2sym" in elem_df.attrs:
        ensg2sym = elem_df.attrs["_ensg2sym"]
        trans_elem["target_symbol"] = (
            trans_elem["intended_target_name"].map(ensg2sym)
                                              .fillna(trans_elem["intended_target_name"])
        )
        trans_elem["gene_symbol"] = (
            trans_elem["gene_id"].map(ensg2sym).fillna(trans_elem["gene_id"])
        )
    else:
        trans_elem["target_symbol"] = trans_elem["intended_target_name"]
        trans_elem["gene_symbol"]   = trans_elem["gene_id"]

    for method, p_col, padj_col, lfc_col in method_specs_elem:
        if padj_col not in trans_elem.columns:
            continue

        df_m = trans_elem[["intended_target_name", "target_symbol",
                            "gene_id", "gene_symbol",
                            padj_col, lfc_col]].dropna(subset=[padj_col]).copy()
        df_m = (
            df_m.groupby(["intended_target_name", "target_symbol",
                          "gene_id", "gene_symbol"], as_index=False)
                .agg({padj_col: "min", lfc_col: "mean"})
        )
        df_m["significant"] = df_m[padj_col] < args.fdr_thresh
        sig = df_m[df_m["significant"]].copy()
        sig["neg_log10_fdr"] = -np.log10(sig[padj_col].clip(lower=1e-300))
        sig["method"] = method

        print(f"\n  [{method}] Total significant trans pairs: {len(sig):,}")

        # Top elements by number of sig trans genes
        top_elements = (
            sig.groupby(["intended_target_name", "target_symbol"])
               .size()
               .reset_index(name="n_trans_sig")
               .nlargest(args.top_n, "n_trans_sig")
        )
        top_elements["method"] = method
        print(f"  Top {args.top_n} targeting elements by trans sig hits:")
        print(top_elements.to_string(index=False))

        tsv_path = os.path.join(out_dir, f"trans_top_elements_{method.lower()}.tsv")
        top_elements.to_csv(tsv_path, sep="\t", index=False)
        print(f"  Saved R-ready table: {tsv_path}")

        # Top trans genes by breadth across targeting elements
        top_genes = (
            sig.groupby("gene_symbol")
               .agg(
                   n_targeting_elements=("intended_target_name", "nunique"),
                   median_lfc          =(lfc_col,                "median"),
                   min_fdr             =(padj_col,               "min"),
               )
               .reset_index()
               .nlargest(args.top_n, "n_targeting_elements")
        )
        top_genes["neg_log10_min_fdr"] = -np.log10(top_genes["min_fdr"].clip(lower=1e-300))
        top_genes["method"] = method
        print(f"\n  Top {args.top_n} trans-regulated genes:")
        print(top_genes.to_string(index=False))

        tsv_path = os.path.join(out_dir, f"trans_top_genes_{method.lower()}.tsv")
        top_genes.to_csv(tsv_path, sep="\t", index=False)
        print(f"  Saved R-ready table: {tsv_path}")

        # -- Plot --------------------------------------------------------------
        fig, axes = plt.subplots(1, 2, figsize=(16, max(5, args.top_n * 0.35)))
        fig.suptitle(
            f"Top trans hits -- {method}  (BH FDR < {args.fdr_thresh})",
            fontsize=13, fontweight="bold",
        )

        ax = axes[0]
        tp = top_elements.sort_values("n_trans_sig", ascending=True).reset_index(drop=True)
        ax.barh(range(len(tp)), tp["n_trans_sig"],
                color="#4C72B0", edgecolor="white", linewidth=0.4)
        ax.set_xlabel("Number of significant trans genes", fontsize=11)
        ax.set_title(f"Top {len(tp)} targeting elements\n"
                     f"ranked by significant trans genes", fontsize=11)
        ax.set_yticks(range(len(tp)))
        ax.set_yticklabels(tp["target_symbol"])
        ax.spines[["top", "right"]].set_visible(False)

        ax = axes[1]
        tg = top_genes.sort_values("n_targeting_elements", ascending=True).reset_index(drop=True)
        lfc_vals = tg["median_lfc"].values
        norm  = mcolors.TwoSlopeNorm(
            vmin=min(lfc_vals.min(), -0.1), vcenter=0, vmax=max(lfc_vals.max(), 0.1)
        )
        cmap   = plt.cm.RdBu_r
        colors = [cmap(norm(v)) for v in lfc_vals]
        ax.barh(range(len(tg)), tg["n_targeting_elements"],
                color=colors, edgecolor="white", linewidth=0.4)
        ax.set_xlabel("Number of targeting elements for which gene is significant", fontsize=11)
        ax.set_title(f"Top {len(tg)} trans-regulated genes\n"
                     f"color = median log2FC", fontsize=11)
        ax.set_yticks(range(len(tg)))
        ax.set_yticklabels(tg["gene_symbol"])
        ax.spines[["top", "right"]].set_visible(False)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.colorbar(sm, ax=ax,
                     label="Median log2FC across targeting elements",
                     shrink=0.6, pad=0.01)

        plt.tight_layout()
        out_path = os.path.join(out_dir, f"trans_top_hits_{method.lower()}.png")
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {out_path}")


# --- Shared plot helper -------------------------------------------------------
def _plot_trans_distribution(all_summaries, id_col, level, fdr_thresh, out_path):
    n_methods = len(all_summaries)
    if n_methods == 0:
        return

    fig, axes = plt.subplots(n_methods, 2,
                             figsize=(12, 4 * n_methods), squeeze=False)
    fig.suptitle(f"Trans results per {level}  (BH FDR < {fdr_thresh})",
                 fontsize=13, fontweight="bold")

    for row, (method, summary) in enumerate(all_summaries.items()):
        ax = axes[row][0]
        ax.hist(summary["n_trans_tested"], bins=40,
                color="#888888", edgecolor="white", linewidth=0.4)
        avg = summary["n_trans_tested"].mean()
        ax.axvline(avg, color="firebrick", linestyle="--", linewidth=1.5,
                   label=f"Mean = {avg:.1f}")
        ax.set_xlabel(f"Trans genes tested per {level}", fontsize=10)
        ax.set_ylabel("Count", fontsize=10)
        ax.set_title(f"{method}: tested", fontsize=11)
        ax.legend(fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)

        ax = axes[row][1]
        ax.hist(summary["n_trans_sig"], bins=40,
                color="#4C72B0", edgecolor="white", linewidth=0.4)
        avg_sig = summary["n_trans_sig"].mean()
        ax.axvline(avg_sig, color="firebrick", linestyle="--", linewidth=1.5,
                   label=f"Mean = {avg_sig:.2f}")
        ax.set_xlabel(f"Significant trans genes per {level}", fontsize=10)
        ax.set_ylabel("Count", fontsize=10)
        ax.set_title(f"{method}: significant (FDR < {fdr_thresh})", fontsize=11)
        ax.legend(fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


# --- Label-colored volcano plot -----------------------------------------------
LABEL_COLORS = {
    "targeting":          "#E07B00",
    "positive control":   "#CC0000",
    "negative control":   "#1155CC",
    "non-targeting":      "#888888",
}
LABEL_ORDER = ["targeting", "positive control", "negative control", "non-targeting"]


def load_guide_metadata() -> pd.DataFrame:
    path = GUIDE_METADATA_PATH
    if not os.path.exists(path):
        print(f"  WARNING: guide metadata not found at {path}")
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", low_memory=False)
    needed = {"guide_id", "type"}
    missing = needed - set(df.columns)
    if missing:
        print(f"  WARNING: guide metadata missing columns: {missing}")
        return pd.DataFrame()
    print(f"  Loaded guide metadata: {len(df):,} rows, "
          f"types: {df['type'].value_counts().to_dict()}")
    return df


def plot_label_volcano(df_results, p_col, lfc_col, level, method, out_path,
                       guide_meta, title_suffix="", p_thresh=0.01, fc_thresh=0.5,
                       out_dir=None):
    """
    Calibration volcano plot colored by guide/element type.
    Writes an intermediate TSV for R/ggplot if out_dir is provided.
    """
    if guide_meta.empty:
        print(f"  Skipping label volcano (no guide metadata)")
        return

    keep_cols = [col for col in
                 [p_col, lfc_col, "guide_id", "intended_target_name",
                  "gene_id", "type", "tested_gene_id"]
                 if col in df_results.columns]
    df = df_results[keep_cols].copy().dropna(subset=[p_col, lfc_col])

    if level == "guide":
        if "type" not in df.columns:
            if "guide_id" in df.columns:
                df = df.merge(guide_meta[["guide_id", "type"]],
                              on="guide_id", how="left")
            else:
                print(f"  WARNING: cannot merge type -- no guide_id col")
                return
    elif level == "element":
        if "type" not in df.columns:
            elem_labels = (
                guide_meta[["intended_target_name", "type"]]
                .dropna(subset=["intended_target_name"])
                .drop_duplicates("intended_target_name")
            )
            if "intended_target_name" in df.columns:
                df = df.merge(elem_labels, on="intended_target_name", how="left")
            elif "gene_id" in df.columns:
                df = df.merge(
                    elem_labels.rename(columns={"intended_target_name": "gene_id"}),
                    on="gene_id", how="left"
                )
            else:
                print(f"  WARNING: cannot merge type -- df cols: {list(df.columns)}")
                return
    else:
        print(f"  WARNING: unrecognised level={level}")
        return

    n_typed = df["type"].notna().sum()
    print(f"  Merged type: {n_typed:,} / {len(df):,} rows have a type label")

    df["-log10p"] = -np.log10(df[p_col].clip(lower=1e-350))
    df.loc[~np.isfinite(df["-log10p"]), "-log10p"] = 350

    sig  = df[p_col] < p_thresh
    up   = sig & (df[lfc_col] >  fc_thresh)
    down = sig & (df[lfc_col] < -fc_thresh)
    df["direction"] = np.where(up, "up", np.where(down, "down", "ns"))
    df["method"]    = method
    df["level"]     = level
    df["title_suffix"] = title_suffix

    # Intermediate TSV for R/ggplot
    if out_dir is not None:
        safe = method.lower().replace(" ", "_").replace("(", "").replace(")", "")
        tsv_path = os.path.join(out_dir, f"volcano_{title_suffix}_{level}_{safe}.tsv")
        df.rename(columns={p_col: "p_value", lfc_col: "log2fc"}).to_csv(
            tsv_path, sep="\t", index=False
        )
        print(f"  Saved R-ready volcano data: {tsv_path}")

    fig, ax = plt.subplots(figsize=(9, 7))

    unlabelled = df[df["type"].isna()]
    if len(unlabelled):
        ax.scatter(unlabelled[lfc_col], unlabelled["-log10p"],
                   color="#dddddd", s=6, alpha=0.25, linewidths=0,
                   rasterized=True, zorder=1)

    for lbl in LABEL_ORDER:
        sub = df[df["type"] == lbl]
        if len(sub) == 0:
            continue
        color  = LABEL_COLORS.get(lbl, "#aaaaaa")
        n_up   = (sub["direction"] == "up").sum()
        n_down = (sub["direction"] == "down").sum()
        n_ns   = (sub["direction"] == "ns").sum()
        ax.scatter(sub[lfc_col], sub["-log10p"],
                   color=color, s=18, alpha=0.75, linewidths=0,
                   rasterized=True, zorder=10,
                   label=f"{lbl}  (up={n_up}, down={n_down}, ns={n_ns})")

    ax.axhline(-np.log10(p_thresh), color="#666666",
               linewidth=0.9, linestyle="--", zorder=5, label=f"p = {p_thresh}")
    ax.axvline( fc_thresh, color="#666666", linewidth=0.9, linestyle="--", zorder=5)
    ax.axvline(-fc_thresh, color="#666666", linewidth=0.9, linestyle="--", zorder=5)
    ax.axvline(0,          color="#444444", linewidth=0.6, zorder=5)

    ax.set_xlabel("log2 Fold Change", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    present_types = sorted(df["type"].dropna().unique())
    absent_types  = [t for t in LABEL_ORDER if t not in present_types]
    absent_note   = f"  (no {', '.join(absent_types)})" if absent_types else ""
    ax.set_title(
        f"{method} {title_suffix} volcano  (per {level})\n"
        f"nominal p < {p_thresh}, |log2FC| > {fc_thresh}{absent_note}",
        fontsize=10, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper left", fontsize=8, frameon=False, bbox_to_anchor=(1.01, 1))

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


# --- Main ---------------------------------------------------------------------
def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    _build_paths(args)

    print(f"Output directory        : {args.output_dir}")
    print(f"Pipeline outputs dir    : {args.pipeline_outputs_dir}")
    print(f"Inference dir           : {args.inference_dir}")
    print(f"MuData filename         : {args.mudata_name}")
    print(f"BH FDR threshold        : {args.fdr_thresh}")
    print(f"Mode                    : {'SCEPTRE-only' if args.sceptre_only else 'Calibrated Perturbo (element) + raw guide'}")
    print(f"Non-targeting           : {'included' if args.keep_nontargeting else 'excluded (default)'}")

    # -- Load files ------------------------------------------------------------
    if args.sceptre_only:
        guide_df = None
        elem_df  = None
        method_specs_guide = []
        method_specs_elem  = []

        if os.path.exists(SCEPTRE_GUIDE):
            guide_df = load_tsv(SCEPTRE_GUIDE, "SCEPTRE per-guide")
            method_specs_guide = [("SCEPTRE", "p_value", "p_value_bh", "log2_fc")]
        else:
            print(f"  Skipping guide-level: {SCEPTRE_GUIDE} not found")

        if os.path.exists(SCEPTRE_ELEM):
            elem_df = load_tsv(SCEPTRE_ELEM, "SCEPTRE per-element")
            method_specs_elem = [("SCEPTRE", "p_value", "p_value_bh", "log2_fc")]
        else:
            print(f"  Skipping element-level: {SCEPTRE_ELEM} not found")

        sym_map_guide_df = guide_df if guide_df is not None else elem_df

    else:
        # Guide level: raw Perturbo file (no calibrated guide file exists)
        guide_df = None
        method_specs_guide = []
        if os.path.exists(MERGED_TRANS_GUIDE):
            print("\n  Loading Perturbo trans guide file ...")
            guide_df = load_tsv(MERGED_TRANS_GUIDE, "Merged trans per-guide")
            method_specs_guide = [("Perturbo", "p_value", "p_value_bh", "log2_fc")]
        else:
            print(f"  Skipping guide-level: {MERGED_TRANS_GUIDE} not found")

        # Element level: calibrated file (auto-detected)
        elem_df = None
        method_specs_elem = []
        if CALIB_TRANS_ELEM and os.path.exists(CALIB_TRANS_ELEM):
            print("\n  Loading calibrated Perturbo trans element file ...")
            elem_df = pd.read_csv(CALIB_TRANS_ELEM, sep="\t")
            n_raw   = len(elem_df)
            print(f"    {n_raw:,} rows, cols: {list(elem_df.columns)}")

            # Normalise column names to match the rest of the pipeline.
            # Use is_cis flag if available to filter trans rows explicitly.
            elem_df = elem_df.rename(columns={
                "element_id":         "intended_target_name",
                "element_symbol":     "target_symbol",
                "element_label":      "type",
                "tested_gene_id":     "gene_id",
                "tested_gene_symbol": "gene_symbol",
                "empirical_pval":     "p_value",
                "empirical_pval_adj": "p_value_bh",
                "log2fc":             "log2_fc",
            })

            # Strip compound IDs (e.g. "ENSG...|chr...")
            for col in ["intended_target_name", "target_symbol"]:
                if col in elem_df.columns:
                    elem_df[col] = elem_df[col].str.split("|").str[0]

            # Use is_cis pre-computed flag to keep only trans rows
            if "is_cis" in elem_df.columns:
                n_before = len(elem_df)
                elem_df = elem_df[elem_df["is_cis"] == False].copy()
                print(f"    is_cis==False filter: {n_before:,} -> {len(elem_df):,} trans rows")

            print(f"    Remapped cols: {list(elem_df.columns)}")
            method_specs_elem = [("Perturbo (calibrated)", "p_value", "p_value_bh", "log2_fc")]
        elif os.path.exists(MERGED_TRANS_ELEM):
            print("\n  Calibrated file not found; falling back to raw element file ...")
            elem_df = load_tsv(MERGED_TRANS_ELEM, "Merged trans per-element")
            method_specs_elem = [("Perturbo", "p_value", "p_value_bh", "log2_fc")]
        else:
            print("  Skipping element-level: no calibrated or raw element file found")

        sym_map_guide_df = guide_df if guide_df is not None else elem_df

    if sym_map_guide_df is None:
        print("ERROR: no input files found. Check --pipeline_outputs_dir and --inference_dir.")
        return

    # -- Filter non-targeting controls ----------------------------------------
    if not args.keep_nontargeting:
        print("\n  Filtering non-targeting controls ...")
        if guide_df is not None:
            guide_df = filter_nontargeting(guide_df, cols=["guide_id"], label="guide")
        if elem_df is not None:
            elem_df = filter_nontargeting(
                elem_df, cols=["intended_target_name"], label="element"
            )
    else:
        print("\n  Keeping non-targeting controls (--keep_nontargeting set)")

    # -- BH correction --------------------------------------------------------
    if guide_df is not None:
        guide_df = bh_correct(guide_df, ["p_value"])

    if elem_df is not None:
        if "p_value_bh" not in elem_df.columns:
            elem_df = bh_correct(elem_df, ["p_value"])
        else:
            print("  Skipping BH correction for element file "
                  "(empirical_pval_adj already present)")

    # -- Ensembl->symbol map --------------------------------------------------
    # Seed from the TSV reference file first, then layer in symbol columns
    # already present in the calibrated file (tested_gene_symbol covers the
    # ~5,500 genes in the tested set; element_symbol / target_symbol covers
    # the perturber genes, which may not be in the tested set at all).
    ensg2sym = build_ensg_to_symbol(sym_map_guide_df)

    if elem_df is not None:
        # tested genes: gene_id -> gene_symbol
        if "gene_symbol" in elem_df.columns:
            extra = dict(zip(elem_df["gene_id"], elem_df["gene_symbol"]))
            ensg2sym.update({k: v for k, v in extra.items()
                             if k not in ensg2sym and pd.notna(v) and not str(v).startswith("ENSG")})

        # perturber elements: intended_target_name -> target_symbol
        # This covers TF knockdown targets (e.g. ENSG00000141510 -> TP53) that
        # are not in the tested gene universe and therefore absent from the TSV.
        if "target_symbol" in elem_df.columns and "intended_target_name" in elem_df.columns:
            perturber_pairs = (
                elem_df[["intended_target_name", "target_symbol"]]
                .drop_duplicates("intended_target_name")
                .dropna()
            )
            for ensg, sym in zip(perturber_pairs["intended_target_name"],
                                 perturber_pairs["target_symbol"]):
                if (ensg not in ensg2sym
                        and str(ensg).startswith("ENSG")
                        and not str(sym).startswith("ENSG")):
                    ensg2sym[ensg] = sym

    print(f"\n  Built Ensembl->symbol map: {len(ensg2sym):,} entries")

    # Re-resolve any remaining raw ENSG IDs in target_symbol
    if elem_df is not None and "target_symbol" in elem_df.columns:
        is_ensg = elem_df["target_symbol"].str.startswith("ENSG", na=False)
        if is_ensg.any():
            elem_df.loc[is_ensg, "target_symbol"] = (
                elem_df.loc[is_ensg, "target_symbol"]
                       .map(ensg2sym)
                       .fillna(elem_df.loc[is_ensg, "target_symbol"])
            )
            n_still_ensg = elem_df["target_symbol"].str.startswith("ENSG", na=False).sum()
            print(f"  target_symbol: {n_still_ensg:,} IDs still unresolved after map")

    if elem_df is not None:
        elem_df.attrs["_ensg2sym"] = ensg2sym

    # -- Run analyses ----------------------------------------------------------
    per_guide_trans(args, args.output_dir, guide_df, method_specs_guide)

    if elem_df is not None:
        per_element_trans(args, args.output_dir, elem_df, method_specs_elem)
        top_trans_hits(args, args.output_dir,
                       guide_df, elem_df,
                       method_specs_guide, method_specs_elem)
    else:
        print("\n  Skipping Analyses 2 & 3: no element-level data available.")

    # --- Analysis 4: Label-colored volcano plots (trans) ----------------------
    print("\n=== Analysis 4: Label-colored volcano plots (trans) ===")
    guide_meta = load_guide_metadata()

    if args.sceptre_only:
        trans_file_specs = []
        if guide_df is not None:
            trans_file_specs.append(("SCEPTRE", "guide", guide_df, "p_value", "log2_fc"))
        if elem_df is not None:
            trans_file_specs.append(("SCEPTRE", "element", elem_df, "p_value", "log2_fc"))
    else:
        trans_file_specs = []
        if guide_df is not None:
            trans_file_specs.append(("Perturbo", "guide", guide_df, "p_value", "log2_fc"))
        if elem_df is not None:
            trans_file_specs.append(
                ("Perturbo (calibrated)", "element", elem_df, "p_value", "log2_fc")
            )

    for method, level, df_raw, p_col, lfc_col in trans_file_specs:
        safe = method.lower().replace(" ", "_").replace("(", "").replace(")", "")
        out_path = os.path.join(args.output_dir,
                                f"volcano_trans_{level}_{safe}.png")
        try:
            plot_label_volcano(
                df_results=df_raw,
                p_col=p_col, lfc_col=lfc_col,
                level=level, method=method,
                out_path=out_path,
                guide_meta=guide_meta,
                title_suffix="trans",
                p_thresh=0.01,
                fc_thresh=0.5,
                out_dir=args.output_dir,
            )
        except Exception as e:
            print(f"  WARNING: could not plot {method} {level} volcano: {e}")

    print("\n=== All analyses complete ===")
    print(f"Results saved to: {args.output_dir}/")


if __name__ == "__main__":
    main()