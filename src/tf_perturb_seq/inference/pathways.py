#!/usr/bin/env python3
"""
Lightweight per-element pathway enrichment analysis.

For each perturbed element, tests whether its significant DEGs are
overrepresented in gene sets (pathways) from a GMT file using Fisher's
exact test, then applies BH correction.

Input:
  - Calibrated trans results from calibrate.py (with empirical_pval_adj)
  - Gene set collection in GMT format (MSigDB, GO, KEGG, etc.)

Output:
  - {prefix}_calibrated_pathway_results.tsv

GMT format:
  Each line: pathway_name<tab>description<tab>gene1<tab>gene2<tab>...
  Genes are symbols (must match tested_gene_symbol in calibrated results).
"""

import argparse
import logging
import os
import sys
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)


# ---------------------------------------------------------------------
# GMT parsing
# ---------------------------------------------------------------------
def load_gmt(path: str) -> dict:
    """
    Load gene sets from a GMT file.

    Parameters
    ----------
    path : str
        Path to GMT file (plain text or .gz).

    Returns
    -------
    dict
        {pathway_id: {"name": str, "genes": set[str]}}
    """
    import gzip

    gene_sets = {}
    opener = gzip.open if path.endswith(".gz") else open

    with opener(path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            pathway_id = parts[0]
            description = parts[1]
            genes = set(g for g in parts[2:] if g)
            gene_sets[pathway_id] = {
                "name": description if description else pathway_id,
                "genes": genes,
            }

    logger.info(f"Loaded {len(gene_sets)} gene sets from {path}")

    # Summary stats
    sizes = [len(gs["genes"]) for gs in gene_sets.values()]
    logger.info(f"  Gene set sizes: min={min(sizes)}, median={np.median(sizes):.0f}, max={max(sizes)}")

    return gene_sets


# ---------------------------------------------------------------------
# Overrepresentation test
# ---------------------------------------------------------------------
def fisher_overrepresentation(
    deg_genes: set,
    pathway_genes: set,
    background_genes: set,
) -> tuple:
    """
    Fisher's exact test for overrepresentation.

    Tests whether deg_genes are enriched in pathway_genes relative to
    the background.

    Returns
    -------
    n_overlap : int
        DEGs in this pathway.
    fold_enrichment : float
        Observed / expected ratio.
    pvalue : float
        One-sided Fisher's exact test p-value (greater).
    """
    deg_in_pathway = len(deg_genes & pathway_genes)
    deg_not_in_pathway = len(deg_genes - pathway_genes)
    bg_in_pathway = len((pathway_genes & background_genes) - deg_genes)
    bg_not_in_pathway = len(background_genes - pathway_genes - deg_genes)

    table = [
        [deg_in_pathway, deg_not_in_pathway],
        [bg_in_pathway, bg_not_in_pathway],
    ]

    _, pvalue = fisher_exact(table, alternative="greater")

    # Fold enrichment
    n_bg = len(background_genes)
    n_pathway_in_bg = len(pathway_genes & background_genes)
    expected = len(deg_genes) * n_pathway_in_bg / n_bg if n_bg > 0 else 0
    fold = deg_in_pathway / expected if expected > 0 else np.inf

    return deg_in_pathway, fold, pvalue


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def run_pathway_analysis(
    calibrated_results_path: str,
    gmt_path: str,
    outdir: str,
    prefix: str,
    fdr_threshold: float = 0.10,
    min_genes: int = 3,
    max_pathway_size: int = 500,
    min_pathway_size: int = 10,
    gene_col: str = "tested_gene_symbol",
    pval_col: str = "empirical_pval_adj",
    element_id_col: str = "element_id",
    element_symbol_col: str = "element_symbol",
) -> None:
    """
    Run per-element pathway enrichment analysis.

    Parameters
    ----------
    calibrated_results_path : str
        Path to calibrated trans results TSV from calibrate.py.
    gmt_path : str
        Path to GMT file with gene sets.
    outdir : str
        Output directory.
    prefix : str
        Prefix for output filenames.
    fdr_threshold : float
        FDR cutoff to define significant DEGs from calibrated results.
    min_genes : int
        Minimum DEGs in a pathway to report the test.
    max_pathway_size : int
        Skip pathways larger than this (too broad to be informative).
    min_pathway_size : int
        Skip pathways smaller than this.
    gene_col : str
        Column with gene symbols in calibrated results.
    pval_col : str
        Column with adjusted p-values in calibrated results.
    element_id_col : str
        Column with element IDs.
    element_symbol_col : str
        Column with element symbols.
    """
    os.makedirs(outdir, exist_ok=True)

    # --- Load data ---
    logger.info(f"Loading calibrated results from {calibrated_results_path}")
    df = pd.read_csv(calibrated_results_path, sep="\t", low_memory=False)
    logger.info(f"  {len(df):,} rows")

    gene_sets = load_gmt(gmt_path)

    # --- Filter gene sets by size ---
    filtered_sets = {
        k: v for k, v in gene_sets.items()
        if min_pathway_size <= len(v["genes"]) <= max_pathway_size
    }
    n_dropped = len(gene_sets) - len(filtered_sets)
    if n_dropped > 0:
        logger.info(f"  Filtered to {len(filtered_sets)} gene sets (dropped {n_dropped} by size)")
    gene_sets = filtered_sets

    # --- Define background: all tested gene symbols ---
    background = set(df[gene_col].dropna().unique())
    logger.info(f"  Background: {len(background):,} genes")

    # --- Restrict gene sets to background ---
    for gs in gene_sets.values():
        gs["genes"] = gs["genes"] & background

    # Remove gene sets with no overlap
    gene_sets = {k: v for k, v in gene_sets.items() if len(v["genes"]) > 0}
    logger.info(f"  {len(gene_sets)} gene sets overlap with background")

    # --- Find significant DEGs per element ---
    sig_mask = df[pval_col].notna() & (df[pval_col] < fdr_threshold)
    sig_df = df[sig_mask].copy()
    logger.info(f"  {len(sig_df):,} significant tests at FDR < {fdr_threshold}")

    # Group DEGs by element
    element_degs = sig_df.groupby(element_id_col)[gene_col].apply(set).to_dict()

    # Also get element symbol mapping
    element_symbols = df.drop_duplicates(element_id_col).set_index(element_id_col)[element_symbol_col].to_dict()

    # Filter to elements with enough DEGs
    elements_to_test = {
        eid: genes for eid, genes in element_degs.items()
        if len(genes) >= min_genes
    }
    logger.info(f"  {len(elements_to_test)} elements with >= {min_genes} DEGs")

    if len(elements_to_test) == 0:
        logger.warning("No elements with enough DEGs for pathway analysis")
        # Write empty output
        out_path = os.path.join(outdir, f"{prefix}_calibrated_pathway_results.tsv")
        pd.DataFrame(columns=[
            "element_id", "element_symbol", "pathway_id", "pathway_name",
            "pathway_source", "n_degs_in_pathway", "n_degs_total",
            "pathway_size", "background_size", "fold_enrichment",
            "pvalue", "pvalue_adj",
        ]).to_csv(out_path, sep="\t", index=False)
        logger.info(f"Written (empty): {out_path}")
        return

    # --- Run enrichment per element ---
    logger.info(f"Running enrichment for {len(elements_to_test)} elements × {len(gene_sets)} gene sets")

    results = []
    for eid, deg_genes in elements_to_test.items():
        element_results = []
        for pathway_id, gs_info in gene_sets.items():
            n_overlap, fold, pval = fisher_overrepresentation(
                deg_genes, gs_info["genes"], background,
            )
            if n_overlap < 1:
                continue

            element_results.append({
                "element_id": eid,
                "element_symbol": element_symbols.get(eid, ""),
                "pathway_id": pathway_id,
                "pathway_name": gs_info["name"],
                "n_degs_in_pathway": n_overlap,
                "n_degs_total": len(deg_genes),
                "pathway_size": len(gs_info["genes"]),
                "background_size": len(background),
                "fold_enrichment": fold,
                "pvalue": pval,
            })

        # BH correct per element
        if element_results:
            elem_df = pd.DataFrame(element_results)
            _, pvals_adj, _, _ = multipletests(
                elem_df["pvalue"].values, alpha=0.05, method="fdr_bh",
            )
            elem_df["pvalue_adj"] = pvals_adj
            results.append(elem_df)

    if not results:
        logger.warning("No enrichment results found")
        return

    out_df = pd.concat(results, ignore_index=True)

    # Infer pathway_source from GMT filename
    gmt_basename = os.path.basename(gmt_path).split(".")[0]
    out_df["pathway_source"] = gmt_basename

    # Reorder columns
    col_order = [
        "element_id", "element_symbol", "pathway_id", "pathway_name",
        "pathway_source", "n_degs_in_pathway", "n_degs_total",
        "pathway_size", "background_size", "fold_enrichment",
        "pvalue", "pvalue_adj",
    ]
    out_df = out_df[col_order]

    # --- Write output ---
    out_path = os.path.join(outdir, f"{prefix}_calibrated_pathway_results.tsv")
    out_df.to_csv(out_path, sep="\t", index=False)
    logger.info(f"Written: {out_path} ({len(out_df):,} rows)")

    # Summary
    n_sig = (out_df["pvalue_adj"] < 0.05).sum()
    n_elements_with_hits = out_df.loc[out_df["pvalue_adj"] < 0.05, "element_id"].nunique()
    logger.info("=" * 50)
    logger.info("Pathway Analysis Summary:")
    logger.info(f"  Elements tested: {len(elements_to_test)}")
    logger.info(f"  Gene sets tested: {len(gene_sets)}")
    logger.info(f"  Total enrichment tests: {len(out_df):,}")
    logger.info(f"  Significant (FDR < 0.05): {n_sig:,}")
    logger.info(f"  Elements with significant pathways: {n_elements_with_hits}")
    logger.info("=" * 50)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Per-element pathway enrichment from calibrated PerTurbo results."
    )
    parser.add_argument(
        "--calibrated-results", required=True,
        help="Path to calibrated trans results TSV from calibrate.py",
    )
    parser.add_argument(
        "--gene-sets", required=True,
        help="Path to GMT file with gene sets (plain text or .gz)",
    )
    parser.add_argument(
        "--outdir", "-o", required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--prefix", required=True,
        help="Prefix for output filenames",
    )
    parser.add_argument(
        "--fdr-threshold", type=float, default=0.10,
        help="FDR cutoff to define DEGs from calibrated results (default: 0.10)",
    )
    parser.add_argument(
        "--min-genes", type=int, default=3,
        help="Minimum DEGs in a pathway to report (default: 3)",
    )
    parser.add_argument(
        "--max-pathway-size", type=int, default=500,
        help="Skip pathways larger than this (default: 500)",
    )
    parser.add_argument(
        "--min-pathway-size", type=int, default=10,
        help="Skip pathways smaller than this (default: 10)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_pathway_analysis(
        calibrated_results_path=args.calibrated_results,
        gmt_path=args.gene_sets,
        outdir=args.outdir,
        prefix=args.prefix,
        fdr_threshold=args.fdr_threshold,
        min_genes=args.min_genes,
        max_pathway_size=args.max_pathway_size,
        min_pathway_size=args.min_pathway_size,
    )
