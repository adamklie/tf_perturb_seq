#!/usr/bin/env python3
"""
Intended target inference QC metrics and visualizations.

This script evaluates how well guides knock down their intended target genes
using results from cis-regulatory analysis stored in MuData.

Data structures used:
  - mdata.uns["cis_per_guide_results"]: DataFrame with cis test results
      Columns: guide_id, gene_id, sceptre_log2_fc, sceptre_p_value,
               perturbo_log2_fc, perturbo_p_value, ...
  - guide.var: Per-guide metadata
      Key columns: guide_id, intended_target_name, gene_name, label
  - guide.layers["guide_assignment"]: Binary matrix (cells x guides)
  - gene.X or gene.layers[layer]: Gene expression matrix

Processing steps:
  1. Load cis_per_guide_results from mdata.uns
  2. Map gene_id to gene_name using guide.var["intended_target_name"] -> "gene_name"
  3. Extract target name from guide_id (e.g., "CD81#strong" -> "CD81")
  4. Filter to "intended target" tests where gene_name == target_name
  5. Compute knockdown metrics (strong knockdowns, significant tests)
  6. Generate volcano plots and summary tables

Metrics computed:
  - n_guides_tested: Number of guides with intended target test results
  - n_strong_knockdowns: Guides with log2FC <= threshold (e.g., 60% knockdown)
  - n_significant: Guides with p_value < threshold
  - frac_strong_knockdowns: Fraction of guides with strong knockdown
  - frac_significant: Fraction of significant tests

Outputs:
  - Summary metrics TSV
  - Volcano plot (log2FC vs -log10(p_value))
  - Intended target results TSV (all guides with their knockdown stats)
"""

import argparse
import logging
import os
import sys
from typing import Optional, Literal

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import mudata
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from mudata import MuData

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)

sns.set_style("white")


# ---------------------------------------------------------------------
# Utility functions for notebook use
# ---------------------------------------------------------------------
def get_perturbation_expression(
    mdata: MuData,
    guide_id: str,
    gene_id: str,
    guide_mod: str = "guide",
    gene_mod: str = "gene",
    gene_layer: Optional[str] = None,
    guide_assignment_layer: str = "guide_assignment",
    non_targeting_label: str = "non_targeting",
) -> pd.DataFrame:
    """
    Compare gene expression between cells with a specific guide vs non-targeting controls.

    This function identifies:
      1. "Perturbed" cells: cells assigned the specified guide
      2. "Control" cells: cells with ONLY non-targeting guides (no targeting guides)

    Then returns expression values for the specified gene in both groups.

    Parameters
    ----------
    mdata : MuData
        MuData object with guide and gene modalities.
    guide_id : str
        Guide identifier (must exist in guide.var_names).
    gene_id : str
        Gene identifier (must exist in gene.var_names).
    guide_mod : str
        Key for guide modality in mdata.mod.
    gene_mod : str
        Key for gene modality in mdata.mod.
    gene_layer : str, optional
        Layer in gene modality to use for expression. If None, uses .X.
    guide_assignment_layer : str
        Layer in guide modality containing binary assignment matrix.
    non_targeting_label : str
        Value in guide.var["label"] marking non-targeting guides.

    Returns
    -------
    pd.DataFrame
        Long-format DataFrame with columns:
          - "expression": Gene expression value
          - "group": Either the guide_id or "non_targeting"
        Index is barcode (cell ID).

    Example
    -------
    >>> expr_df = get_perturbation_expression(mdata, "CD81#strong", "ENSG00000110651")
    >>> sns.violinplot(data=expr_df, x="group", y="expression")
    """
    # Extract modalities
    guide: AnnData = mdata.mod[guide_mod]
    gene: AnnData = mdata.mod[gene_mod]

    # Validate inputs
    if guide_id not in guide.var_names:
        raise ValueError(f"guide_id '{guide_id}' not found in guide.var_names")
    if gene_id not in gene.var_names:
        raise ValueError(f"gene_id '{gene_id}' not found in gene.var_names")
    if "label" not in guide.var.columns:
        raise ValueError("guide.var must contain a 'label' column")

    # -----------------------------------------------------------------
    # Get gene expression for all cells
    # sc.get.obs_df returns DataFrame with index=cells, columns=genes
    # -----------------------------------------------------------------
    gene_expr_df = sc.get.obs_df(gene, keys=[gene_id], layer=gene_layer)
    gene_expr = gene_expr_df[gene_id]

    # -----------------------------------------------------------------
    # Identify cells with the specific guide (perturbed group)
    # guide_assignment_layer contains binary matrix: >0 means assigned
    # -----------------------------------------------------------------
    guide_assign = sc.get.obs_df(guide, keys=[guide_id], layer=guide_assignment_layer)
    cells_with_guide = guide_assign.index[guide_assign[guide_id] > 0].tolist()

    # -----------------------------------------------------------------
    # Identify cells with ONLY non-targeting guides (control group)
    # Step 1: Find cells with any non-targeting guide
    # Step 2: Find cells with any targeting guide
    # Step 3: Control = cells with NT guides but NO targeting guides
    # -----------------------------------------------------------------
    # Get non-targeting guide indices
    nt_guide_idx = guide.var.index[guide.var["label"] == non_targeting_label].tolist()
    if len(nt_guide_idx) == 0:
        raise ValueError(f"No guides with label '{non_targeting_label}' in guide.var['label']")

    # Cells with any non-targeting guide
    nt_assign = sc.get.obs_df(guide, keys=nt_guide_idx, layer=guide_assignment_layer)
    cells_with_nt = nt_assign.index[nt_assign.sum(axis=1) > 0].tolist()

    # Cells with any targeting (non-NT) guide
    targeting_guide_idx = guide.var.index[guide.var["label"] != non_targeting_label].tolist()
    targeting_assign = sc.get.obs_df(guide, keys=targeting_guide_idx, layer=guide_assignment_layer)
    cells_with_targeting = targeting_assign.index[targeting_assign.sum(axis=1) > 0].tolist()

    # Control cells: have NT guides but NO targeting guides
    cells_nt_only = list(set(cells_with_nt) - set(cells_with_targeting))

    # -----------------------------------------------------------------
    # Build output DataFrame
    # -----------------------------------------------------------------
    perturbed_expr = gene_expr.loc[cells_with_guide]
    control_expr = gene_expr.loc[cells_nt_only]

    expr_df = pd.DataFrame({
        "expression": pd.concat([perturbed_expr, control_expr]),
        "group": [guide_id] * len(perturbed_expr) + [non_targeting_label] * len(control_expr),
    })
    expr_df.index.name = "barcode"

    return expr_df


def build_intended_target_map(guide_var: pd.DataFrame) -> pd.Series:
    """
    Build a mapping from intended_target_name (gene ID) to gene_name (symbol).

    Parameters
    ----------
    guide_var : pd.DataFrame
        guide.var DataFrame with columns "intended_target_name" and "gene_name".

    Returns
    -------
    pd.Series
        Index = intended_target_name, values = gene_name.
        Duplicates are dropped (keeps first).
    """
    if "intended_target_name" not in guide_var.columns:
        raise ValueError("guide.var must have 'intended_target_name' column")
    if "gene_name" not in guide_var.columns:
        raise ValueError("guide.var must have 'gene_name' column")

    id_map = guide_var.set_index("intended_target_name")["gene_name"]
    id_map = id_map[~id_map.index.duplicated(keep="first")]
    return id_map


# ---------------------------------------------------------------------
# Core QC functions
# ---------------------------------------------------------------------
def load_cis_results(
    mdata: MuData,
    results_key: str = "cis_per_guide_results",
) -> pd.DataFrame:
    """
    Load cis-regulatory test results from MuData.uns.

    Parameters
    ----------
    mdata : MuData
        MuData object.
    results_key : str
        Key in mdata.uns containing cis results DataFrame.

    Returns
    -------
    pd.DataFrame
        Cis test results with columns like guide_id, gene_id, log2_fc, p_value.
    """
    if results_key not in mdata.uns:
        raise ValueError(
            f"Results key '{results_key}' not found in mdata.uns. "
            f"Available keys: {list(mdata.uns.keys())}"
        )
    return mdata.uns[results_key].copy()


def filter_to_intended_targets(
    cis_results: pd.DataFrame,
    guide_var: pd.DataFrame,
    log2fc_col: str = "perturbo_log2_fc",
    pvalue_col: str = "perturbo_p_value",
) -> pd.DataFrame:
    """
    Filter cis results to only intended target tests (guide targets its intended gene).

    Processing steps:
      1. Map gene_id to gene_name using guide.var
      2. Extract target_name from guide_id (e.g., "CD81#strong" -> "CD81")
      3. Keep rows where gene_name == target_name

    Parameters
    ----------
    cis_results : pd.DataFrame
        Raw cis results with guide_id, gene_id columns.
    guide_var : pd.DataFrame
        guide.var with intended_target_name and gene_name columns.
    log2fc_col : str
        Column name for log2 fold change.
    pvalue_col : str
        Column name for p-value.

    Returns
    -------
    pd.DataFrame
        Filtered results with added columns: gene_name, target_name.
    """
    # Build gene ID -> gene name mapping
    id_map = build_intended_target_map(guide_var)

    # Add gene_name column by mapping gene_id
    results = cis_results.copy()
    results["gene_name"] = results["gene_id"].map(id_map)

    # Extract target_name from guide_id (format: "TARGET#variant" -> "TARGET")
    results["target_name"] = results["guide_id"].str.split("#").str[0]

    # Filter to intended targets only
    intended = results[results["gene_name"] == results["target_name"]].copy()

    # Deduplicate by guide_id + gene_id
    intended = intended.drop_duplicates(subset=["guide_id", "gene_id"])

    # Merge with full guide metadata
    guide_meta = guide_var.reset_index(drop=True)
    intended = guide_meta.merge(
        intended[["guide_id", "gene_id", log2fc_col, pvalue_col]],
        on="guide_id",
        how="left",
    )

    return intended


def compute_knockdown_metrics(
    intended_results: pd.DataFrame,
    log2fc_col: str = "perturbo_log2_fc",
    pvalue_col: str = "perturbo_p_value",
    fc_threshold: float = 0.4,
    pval_threshold: float = 0.05,
) -> dict:
    """
    Compute summary metrics for intended target knockdown efficiency.

    Parameters
    ----------
    intended_results : pd.DataFrame
        Filtered cis results for intended targets.
    log2fc_col : str
        Column with log2 fold change values.
    pvalue_col : str
        Column with p-values.
    fc_threshold : float
        Fold change threshold for "strong" knockdown.
        Default 0.4 means 60% knockdown (log2(0.4) = -1.32).
    pval_threshold : float
        P-value threshold for significance.

    Returns
    -------
    dict
        Dictionary of metrics.
    """
    # Filter to rows with valid results
    valid = intended_results.dropna(subset=[log2fc_col, pvalue_col])

    n_guides_tested = len(valid)
    log2fc_threshold = np.log2(fc_threshold)

    # Strong knockdowns: log2FC <= threshold (negative = knockdown)
    strong_knockdowns = valid[valid[log2fc_col] <= log2fc_threshold]
    n_strong = len(strong_knockdowns)

    # Significant tests
    significant = valid[valid[pvalue_col] < pval_threshold]
    n_sig = len(significant)

    # Both strong AND significant
    strong_and_sig = valid[
        (valid[log2fc_col] <= log2fc_threshold) & (valid[pvalue_col] < pval_threshold)
    ]
    n_strong_and_sig = len(strong_and_sig)

    metrics = {
        "n_guides_total": len(intended_results),
        "n_guides_tested": n_guides_tested,
        "fc_threshold": fc_threshold,
        "log2fc_threshold": log2fc_threshold,
        "pval_threshold": pval_threshold,
        "n_strong_knockdowns": n_strong,
        "n_significant": n_sig,
        "n_strong_and_significant": n_strong_and_sig,
        "frac_strong_knockdowns": n_strong / n_guides_tested if n_guides_tested > 0 else 0,
        "frac_significant": n_sig / n_guides_tested if n_guides_tested > 0 else 0,
        "frac_strong_and_significant": n_strong_and_sig / n_guides_tested if n_guides_tested > 0 else 0,
        "median_log2fc": valid[log2fc_col].median() if n_guides_tested > 0 else np.nan,
        "mean_log2fc": valid[log2fc_col].mean() if n_guides_tested > 0 else np.nan,
    }

    return metrics


# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------
def plot_volcano(
    intended_results: pd.DataFrame,
    outdir: str,
    prefix: str = "intended_target",
    log2fc_col: str = "perturbo_log2_fc",
    pvalue_col: str = "perturbo_p_value",
    label_col: str = "label",
    pval_threshold: float = 0.05,
    n_label: int = 10,
) -> None:
    """
    Plot volcano plot of intended target knockdown results.

    X-axis: log2 fold change
    Y-axis: -log10(p-value)
    Points colored by guide label (e.g., positive_control, target, non_targeting)
    Top N most significant guides are labeled.

    Parameters
    ----------
    intended_results : pd.DataFrame
        Filtered results for intended targets.
    outdir : str
        Output directory.
    prefix : str
        Prefix for output filename.
    log2fc_col, pvalue_col : str
        Column names for fold change and p-value.
    label_col : str
        Column for coloring points.
    pval_threshold : float
        Significance threshold for horizontal line.
    n_label : int
        Number of top guides to label.
    """
    plot_data = intended_results.dropna(subset=[log2fc_col, pvalue_col]).copy()
    if len(plot_data) == 0:
        logger.warning("No valid data for volcano plot")
        return

    # Compute -log10(p-value), adding small epsilon to avoid log(0)
    plot_data["neg_log10_pval"] = -np.log10(plot_data[pvalue_col] + 1e-300)

    with sns.plotting_context("talk"):
        fig, ax = plt.subplots(figsize=(8, 6))

        # Scatter plot colored by label
        if label_col in plot_data.columns:
            plot_data[label_col] = plot_data[label_col].astype(str)
            sns.scatterplot(
                data=plot_data, x=log2fc_col, y="neg_log10_pval",
                hue=label_col, alpha=0.7, ax=ax,
            )
            ax.legend(title="Label", bbox_to_anchor=(1.02, 1), loc="upper left")
        else:
            ax.scatter(plot_data[log2fc_col], plot_data["neg_log10_pval"], alpha=0.7)

        # Reference lines
        ax.axhline(-np.log10(pval_threshold), color="red", linestyle="--", alpha=0.5)
        ax.axvline(0, color="grey", linestyle="--", alpha=0.5)

        # Label top N guides
        top_n = plot_data.nsmallest(n_label, pvalue_col)
        for _, row in top_n.iterrows():
            ax.text(
                row[log2fc_col], row["neg_log10_pval"],
                row["guide_id"], fontsize=8, alpha=0.8,
            )

        ax.set_xlabel("Log2 Fold Change")
        ax.set_ylabel("-Log10(p-value)")
        ax.set_title(f"Intended Target Knockdown (N={len(plot_data)} guides)")

        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_volcano.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


def plot_log2fc_distribution(
    intended_results: pd.DataFrame,
    outdir: str,
    prefix: str = "intended_target",
    log2fc_col: str = "perturbo_log2_fc",
    label_col: str = "label",
) -> None:
    """
    Plot distribution of log2 fold changes for intended targets.
    """
    plot_data = intended_results.dropna(subset=[log2fc_col]).copy()
    if len(plot_data) == 0:
        logger.warning("No valid data for log2FC distribution plot")
        return

    with sns.plotting_context("talk"):
        fig, ax = plt.subplots(figsize=(10, 5))

        if label_col in plot_data.columns:
            plot_data[label_col] = plot_data[label_col].astype(str)
            sns.histplot(
                data=plot_data, x=log2fc_col, hue=label_col,
                bins=50, ax=ax, element="step", stat="count",
            )
        else:
            sns.histplot(plot_data[log2fc_col], bins=50, ax=ax)

        # Reference lines
        ax.axvline(0, color="grey", linestyle="--", alpha=0.5)
        ax.axvline(np.log2(0.4), color="red", linestyle="--", alpha=0.5, label="60% KD")

        median_fc = plot_data[log2fc_col].median()
        ax.axvline(median_fc, color="blue", linestyle="--", alpha=0.5)
        ax.text(0.95, 0.95, f"Median: {median_fc:.2f}", transform=ax.transAxes,
                ha="right", va="top", fontsize=12)

        ax.set_xlabel("Log2 Fold Change")
        ax.set_ylabel("Count")
        ax.set_title("Distribution of Intended Target Log2 Fold Changes")

        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_log2fc_distribution.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def run_intended_target_qc(
    input_path: str,
    outdir: str,
    guide_mod_key: str = "guide",
    results_key: str = "cis_per_guide_results",
    log2fc_col: str = "perturbo_log2_fc",
    pvalue_col: str = "perturbo_p_value",
    fc_threshold: float = 0.4,
    pval_threshold: float = 0.05,
    prefix: str = "intended_target",
) -> None:
    """
    Run intended target inference QC.

    Parameters
    ----------
    input_path : str
        Path to MuData (.h5mu) file.
    outdir : str
        Output directory.
    guide_mod_key : str
        Modality key for guide data.
    results_key : str
        Key in mdata.uns for cis results.
    log2fc_col : str
        Column name for log2 fold change in results.
    pvalue_col : str
        Column name for p-value in results.
    fc_threshold : float
        Fold change threshold for strong knockdown (default 0.4 = 60% knockdown).
    pval_threshold : float
        P-value threshold for significance.
    prefix : str
        Prefix for output filenames.
    """
    os.makedirs(outdir, exist_ok=True)

    # Load data
    logger.info(f"Loading MuData from {input_path}")
    mdata = mudata.read_h5mu(input_path)

    if guide_mod_key not in mdata.mod:
        raise ValueError(f"Modality '{guide_mod_key}' not found in MuData")

    guide = mdata.mod[guide_mod_key]
    logger.info(f"Guide modality: {guide.n_obs} cells, {guide.n_vars} guides")

    # Load cis results
    cis_results = load_cis_results(mdata, results_key=results_key)
    logger.info(f"Loaded {len(cis_results)} cis test results from '{results_key}'")

    # Filter to intended targets
    intended = filter_to_intended_targets(
        cis_results, guide.var,
        log2fc_col=log2fc_col, pvalue_col=pvalue_col,
    )
    logger.info(f"Filtered to {len(intended)} intended target tests")

    # Compute metrics
    metrics = compute_knockdown_metrics(
        intended,
        log2fc_col=log2fc_col, pvalue_col=pvalue_col,
        fc_threshold=fc_threshold, pval_threshold=pval_threshold,
    )

    # Save metrics
    metrics_df = pd.DataFrame([metrics])
    metrics_path = os.path.join(outdir, f"{prefix}_metrics.tsv")
    metrics_df.to_csv(metrics_path, sep="\t", index=False)
    logger.info(f"Saved metrics to {metrics_path}")

    # Save full intended target results
    results_path = os.path.join(outdir, f"{prefix}_results.tsv")
    output_cols = ["guide_id", "gene_name", "label", log2fc_col, pvalue_col]
    output_cols = [c for c in output_cols if c in intended.columns]
    intended[output_cols].to_csv(results_path, sep="\t", index=False)
    logger.info(f"Saved results to {results_path}")

    # Generate plots
    plot_volcano(
        intended, outdir, prefix=prefix,
        log2fc_col=log2fc_col, pvalue_col=pvalue_col,
        pval_threshold=pval_threshold,
    )
    plot_log2fc_distribution(
        intended, outdir, prefix=prefix,
        log2fc_col=log2fc_col,
    )

    # Log summary
    logger.info("=" * 50)
    logger.info("Intended Target QC Summary:")
    logger.info(f"  Guides tested: {metrics['n_guides_tested']}")
    logger.info(f"  Strong knockdowns (>={int((1-fc_threshold)*100)}%): "
                f"{metrics['n_strong_knockdowns']} ({metrics['frac_strong_knockdowns']*100:.1f}%)")
    logger.info(f"  Significant (p<{pval_threshold}): "
                f"{metrics['n_significant']} ({metrics['frac_significant']*100:.1f}%)")
    logger.info(f"  Median log2FC: {metrics['median_log2fc']:.3f}")
    logger.info("=" * 50)

    logger.info("Intended target QC complete.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Intended target inference QC metrics and visualizations."
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Path to MuData (.h5mu) file."
    )
    parser.add_argument(
        "--outdir", "-o", required=True,
        help="Output directory."
    )
    parser.add_argument(
        "--guide-mod-key", default="guide",
        help="Modality key for guide data (default: 'guide')."
    )
    parser.add_argument(
        "--results-key", default="cis_per_guide_results",
        help="Key in mdata.uns for cis results (default: 'cis_per_guide_results')."
    )
    parser.add_argument(
        "--log2fc-col", default="perturbo_log2_fc",
        help="Column name for log2 fold change (default: 'perturbo_log2_fc')."
    )
    parser.add_argument(
        "--pvalue-col", default="perturbo_p_value",
        help="Column name for p-value (default: 'perturbo_p_value')."
    )
    parser.add_argument(
        "--fc-threshold", type=float, default=0.4,
        help="Fold change threshold for strong knockdown (default: 0.4 = 60%% KD)."
    )
    parser.add_argument(
        "--pval-threshold", type=float, default=0.05,
        help="P-value threshold for significance (default: 0.05)."
    )
    parser.add_argument(
        "--prefix", default="intended_target",
        help="Prefix for output filenames (default: 'intended_target')."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_intended_target_qc(
        input_path=args.input,
        outdir=args.outdir,
        guide_mod_key=args.guide_mod_key,
        results_key=args.results_key,
        log2fc_col=args.log2fc_col,
        pvalue_col=args.pvalue_col,
        fc_threshold=args.fc_threshold,
        pval_threshold=args.pval_threshold,
        prefix=args.prefix,
    )
