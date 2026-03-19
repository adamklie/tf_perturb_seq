#!/usr/bin/env python3
"""
Calibrate PerTurbo per-element inference results using an empirical null.

The IGVF CRISPR pipeline (PerTurbo) produces posterior p-values from a
negative binomial model. These are not well-calibrated for frequentist hit
calling. This script derives empirical p-values by comparing each targeting
element's test statistics against the distribution from non-targeting control
(NTC) elements.

Calibration methods:
  - ecdf (default): Rank-based empirical CDF. For each test, count the
    fraction of NTC statistics that are at least as extreme.
    p = (r + 1) / (B + 1)  where r = # NTC stats >= |observed|, B = # NTC stats
  - t-fit: Fit a t-distribution (location fixed at 0) to NTC z-values,
    then compute p-values from the fitted survival function.

Reference: https://github.com/pinellolab/PerTurbo/blob/main/src/perturbo/utils/utils.py

Inputs:
  - perturbo_trans_per_element_output.tsv.gz  (all trans tests)
  - perturbo_cis_per_element_output.tsv.gz    (cis-window tests, for annotation)
  - inference_mudata.h5mu                     (guide metadata, gene annotations)

Outputs:
  - {prefix}_calibrated_trans_results.tsv     (master table, all tests)
  - {prefix}_calibrated_direct_target_results.tsv
  - {prefix}_calibrated_cis_results.tsv
  - {prefix}_calibrated_trans_only_results.tsv
"""

import argparse
import logging
import os
import sys
from typing import Literal, Optional

import mudata
import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import issparse
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)


# ---------------------------------------------------------------------
# Empirical null calibration
# ---------------------------------------------------------------------
def empirical_pvals_ecdf(
    null_stats: np.ndarray,
    test_stats: np.ndarray,
    two_sided: bool = True,
) -> np.ndarray:
    """
    Compute empirical p-values via rank-based eCDF.

    For each test statistic, count the fraction of null statistics
    that are at least as extreme. Uses bias-corrected formula:
        p = (r + 1) / (B + 1)

    Parameters
    ----------
    null_stats : array
        Test statistics from non-targeting controls.
    test_stats : array
        Test statistics from targeting elements.
    two_sided : bool
        If True, use absolute values (test for any effect direction).

    Returns
    -------
    array
        Empirical p-values, same length as test_stats.
    """
    null = null_stats[np.isfinite(null_stats)].copy()
    test = test_stats.copy()

    if two_sided:
        null = np.abs(null)
        test = np.abs(test)

    null_sorted = np.sort(null)
    B = len(null_sorted)

    # For each test stat, count how many null stats are >= test stat
    # searchsorted with side="right" gives index where test would be inserted
    # so r = B - index = number of null values >= test value
    idx = np.searchsorted(null_sorted, test, side="right")
    r = B - idx

    pvals = (r + 1.0) / (B + 1.0)
    return pvals


def empirical_pvals_tfit(
    null_stats: np.ndarray,
    test_stats: np.ndarray,
    two_sided: bool = True,
    winsor: Optional[float] = 0.01,
) -> np.ndarray:
    """
    Compute empirical p-values by fitting a t-distribution to the null.

    Fits a t-distribution with location fixed at 0 to the NTC z-values,
    then computes p-values from the fitted survival function.

    Parameters
    ----------
    null_stats : array
        Test statistics (z-values) from NTC elements.
    test_stats : array
        Test statistics from targeting elements.
    two_sided : bool
        If True, compute two-sided p-values.
    winsor : float or None
        If not None, winsorize null at (winsor, 1-winsor) quantiles
        before fitting for robustness.

    Returns
    -------
    array
        Parametric empirical p-values.
    """
    null = null_stats[np.isfinite(null_stats)].copy()

    if winsor is not None:
        lo, hi = np.quantile(null, [winsor, 1.0 - winsor])
        null = np.clip(null, lo, hi)

    # Fit t-distribution with location fixed at 0
    df_hat, _, scale_hat = stats.t.fit(null, floc=0.0)
    logger.info(f"t-fit: df={df_hat:.2f}, scale={scale_hat:.4f}")

    # Standardize test statistics
    z_std = test_stats / scale_hat

    if two_sided:
        pvals = 2.0 * stats.t.sf(np.abs(z_std), df_hat)
    else:
        pvals = stats.t.sf(z_std, df_hat)

    return pvals


# ---------------------------------------------------------------------
# Data loading and annotation
# ---------------------------------------------------------------------
def load_trans_results(path: str) -> pd.DataFrame:
    """Load per-element trans results TSV."""
    logger.info(f"Loading trans results from {path}")
    df = pd.read_csv(path, sep="\t")
    logger.info(f"  {len(df):,} rows, columns: {list(df.columns)}")
    return df


def load_cis_results(path: str) -> pd.DataFrame:
    """Load per-element cis results TSV."""
    logger.info(f"Loading cis results from {path}")
    df = pd.read_csv(path, sep="\t")
    logger.info(f"  {len(df):,} rows, columns: {list(df.columns)}")
    return df


def load_guide_metadata(mdata_path: str, guide_mod_key: str = "guide"):
    """
    Load guide and gene metadata from MuData.

    Returns
    -------
    guide_var : pd.DataFrame
        Guide metadata with columns like intended_target_name, type, targeting, etc.
    gene_var : pd.DataFrame
        Gene metadata (from gene modality .var).
    n_cells_per_element : pd.Series
        Number of cells assigned to each element (intended_target_name).
    """
    logger.info(f"Loading MuData from {mdata_path}")
    mdata_obj = mudata.read_h5mu(mdata_path)

    guide = mdata_obj.mod[guide_mod_key]
    guide_var = guide.var.copy()

    # Compute n_cells per element from guide_assignment layer
    # Group by (intended_target_name, chr, start, end) to get per-promoter counts
    if "guide_assignment" in guide.layers:
        assignment = guide.layers["guide_assignment"]
        if issparse(assignment):
            cells_per_guide = np.array((assignment > 0).sum(axis=0)).flatten()
        else:
            cells_per_guide = (assignment > 0).sum(axis=0)
        guide_var["n_cells_assigned"] = cells_per_guide

        # Aggregate to element level using coordinates for uniqueness
        coord_cols = ["intended_target_name", "intended_target_chr",
                      "intended_target_start", "intended_target_end"]
        if all(c in guide_var.columns for c in coord_cols):
            n_cells_per_element = (
                guide_var.groupby(coord_cols, observed=True)["n_cells_assigned"]
                .sum()
            )
            logger.info(f"  n_cells computed per element (with coordinates)")
        else:
            # Fall back to gene-level aggregation
            n_cells_per_element = (
                guide_var.groupby("intended_target_name", observed=True)["n_cells_assigned"]
                .sum()
            )
            logger.info(f"  n_cells computed per gene (no coordinate columns in guide.var)")
    else:
        logger.warning("guide_assignment layer not found; n_cells will be NaN")
        n_cells_per_element = pd.Series(dtype=float)

    # Gene metadata for symbol lookup
    gene_var = mdata_obj.mod["gene"].var.copy() if "gene" in mdata_obj.mod else pd.DataFrame()

    logger.info(f"  {len(guide_var)} guides, {guide_var['intended_target_name'].nunique()} elements")
    if "type" in guide_var.columns:
        logger.info(f"  Element types: {guide_var.groupby('type', observed=True)['intended_target_name'].nunique().to_dict()}")

    return guide_var, gene_var, n_cells_per_element


def get_ntc_elements(guide_var: pd.DataFrame, non_targeting_label: str = "non_targeting") -> set:
    """
    Get the set of intended_target_name values that correspond to NTC elements.

    Identifies NTCs by the 'type' column (e.g., 'non_targeting') or the
    'targeting' column (False).
    """
    if "type" in guide_var.columns:
        # Normalize: replace hyphens/spaces with underscores for comparison
        normalized_type = guide_var["type"].str.lower().str.replace(r"[-\s]", "_", regex=True)
        normalized_label = non_targeting_label.lower().replace("-", "_").replace(" ", "_")
        ntc_mask = normalized_type == normalized_label
    elif "targeting" in guide_var.columns:
        ntc_mask = ~guide_var["targeting"].astype(bool)
    else:
        raise ValueError("Cannot identify NTC elements: no 'type' or 'targeting' column in guide.var")

    ntc_elements = set(guide_var.loc[ntc_mask, "intended_target_name"].unique())
    logger.info(f"  {len(ntc_elements)} NTC elements identified")
    return ntc_elements


def build_gene_symbol_map(gene_var: pd.DataFrame) -> dict:
    """Build gene_id → symbol mapping from gene.var."""
    if "symbol" in gene_var.columns:
        return gene_var["symbol"].to_dict()
    return {}


def build_element_symbol_map(guide_var: pd.DataFrame, gene_var: pd.DataFrame) -> dict:
    """Build intended_target_name → gene symbol mapping.

    Tries guide.var["gene_name"] first, then falls back to gene.var["symbol"]
    since intended_target_name is often an Ensembl gene ID.
    """
    if "gene_name" in guide_var.columns:
        mapping = guide_var.drop_duplicates("intended_target_name").set_index("intended_target_name")
        result = mapping["gene_name"].to_dict()
        if any(v and not v.startswith("ENSG") for v in result.values() if isinstance(v, str)):
            return result

    # Fall back to gene.var symbol lookup (intended_target_name = Ensembl ID)
    if "symbol" in gene_var.columns:
        return gene_var["symbol"].to_dict()
    return {}


def build_element_type_map(guide_var: pd.DataFrame) -> dict:
    """Build intended_target_name → type mapping from guide.var."""
    if "type" in guide_var.columns:
        mapping = guide_var.drop_duplicates("intended_target_name").set_index("intended_target_name")
        return mapping["type"].to_dict()
    return {}


# ---------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------
def annotate_cis(
    trans_df: pd.DataFrame,
    gene_var: pd.DataFrame,
    cis_window: int = 100_000,
) -> pd.Series:
    """
    Annotate which trans results are cis interactions based on coordinates.

    An element-gene pair is cis if the element and gene are on the same
    chromosome and within cis_window bp of each other.

    Parameters
    ----------
    trans_df : pd.DataFrame
        Trans results with intended_target_chr/start/end columns.
    gene_var : pd.DataFrame
        Gene metadata with chromosome/start/end for each gene_id.
    cis_window : int
        Maximum distance (bp) for a cis interaction.
    """
    if "intended_target_chr" not in trans_df.columns:
        logger.warning("No element coordinate columns; cannot annotate cis")
        return pd.Series(False, index=trans_df.index)

    # Build gene coordinate lookup from gene.var
    # Try common column names for gene coordinates
    gene_coords = None
    for chr_col, start_col, end_col in [
        ("gene_chr", "gene_start", "gene_end"),
        ("chromosome", "start", "end"),
        ("chr", "start", "end"),
        ("seqname", "start", "end"),
    ]:
        if chr_col in gene_var.columns and start_col in gene_var.columns:
            gene_coords = gene_var[[chr_col, start_col, end_col]].copy()
            gene_coords.columns = ["gene_chr", "gene_start", "gene_end"]
            break

    if gene_coords is None:
        logger.warning("No gene coordinate columns found in gene.var; cannot annotate cis")
        return pd.Series(False, index=trans_df.index)

    gene_coords["gene_id"] = gene_var.index.values
    gene_coords = gene_coords.reset_index(drop=True)
    logger.info(f"  Gene coordinates available for {len(gene_coords):,} genes")

    # Merge gene coords onto trans results
    merged = trans_df[["gene_id", "intended_target_chr", "intended_target_start", "intended_target_end"]].reset_index(drop=True).merge(
        gene_coords, on="gene_id", how="left",
    )

    # Cis = same chromosome AND within window
    same_chr = merged["intended_target_chr"] == merged["gene_chr"]
    # Distance: closest point of element to closest point of gene
    elem_mid = (merged["intended_target_start"] + merged["intended_target_end"]) / 2
    gene_mid = (merged["gene_start"] + merged["gene_end"]) / 2
    distance = (elem_mid - gene_mid).abs()

    is_cis = same_chr & (distance <= cis_window)
    logger.info(f"  {is_cis.sum():,} cis pairs (within {cis_window/1000:.0f}kb)")

    return is_cis


def annotate_direct_target(trans_df: pd.DataFrame) -> pd.Series:
    """
    Annotate which rows test the element's own intended target.

    Direct target = gene_id matches intended_target_name.
    """
    return trans_df["gene_id"] == trans_df["intended_target_name"]


# ---------------------------------------------------------------------
# Main calibration pipeline
# ---------------------------------------------------------------------
def run_calibration(
    trans_results_path: str,
    mudata_path: str,
    outdir: str,
    prefix: str,
    null_method: Literal["ecdf", "t-fit"] = "t-fit",
    cis_window: int = 100_000,
    non_targeting_label: str = "non_targeting",
    guide_mod_key: str = "guide",
    two_sided: bool = True,
) -> None:
    """
    Run the full calibration pipeline.

    Parameters
    ----------
    trans_results_path : str
        Path to perturbo_trans_per_element_output.tsv.gz.
    mudata_path : str
        Path to inference_mudata.h5mu.
    outdir : str
        Output directory.
    prefix : str
        Prefix for output filenames.
    null_method : str
        Calibration method: "ecdf" or "t-fit".
    cis_window : int
        Distance (bp) for cis annotation (default: 500kb).
    non_targeting_label : str
        Label for non-targeting elements in guide.var["type"].
    guide_mod_key : str
        Modality key for guide data.
    two_sided : bool
        Whether to compute two-sided p-values.
    """
    os.makedirs(outdir, exist_ok=True)

    # --- Load data ---
    trans_df = load_trans_results(trans_results_path)
    guide_var, gene_var, n_cells_per_element = load_guide_metadata(mudata_path, guide_mod_key)

    # --- Identify NTC elements ---
    ntc_elements = get_ntc_elements(guide_var, non_targeting_label)

    # --- Build lookup maps ---
    gene_symbol_map = build_gene_symbol_map(gene_var)
    element_symbol_map = build_element_symbol_map(guide_var, gene_var)
    element_type_map = build_element_type_map(guide_var)

    # --- Split NTC vs targeting results ---
    is_ntc = trans_df["intended_target_name"].isin(ntc_elements)
    ntc_results = trans_df[is_ntc].copy()
    # Keep all results for output (including NTC)
    logger.info(f"  NTC results: {len(ntc_results):,} rows")
    logger.info(f"  Targeting results: {(~is_ntc).sum():,} rows")

    # --- Compute z-values (log2fc / log2fc_std) as test statistic ---
    if "log2_fc_std" not in trans_df.columns:
        logger.error(
            "log2_fc_std column not found in trans results. "
            "Calibration requires z-values (log2fc / log2fc_std). "
            "This may indicate a pipeline version that doesn't output standard errors "
            "(e.g., sceptre_v11 is known to be missing this column)."
        )
        return

    valid_se = trans_df["log2_fc_std"].notna() & (trans_df["log2_fc_std"] > 0)
    trans_df["z_value"] = np.where(
        valid_se,
        trans_df["log2_fc"] / trans_df["log2_fc_std"],
        np.nan,
    )
    stat_col = "z_value"
    logger.info(f"  Using z-values (log2fc / log2fc_std) as test statistic")
    logger.info(f"  {valid_se.sum():,} rows with valid z-values out of {len(trans_df):,}")

    # --- Build null distribution from NTC ---
    null_stats = trans_df.loc[is_ntc, stat_col].dropna().values
    logger.info(f"  Null distribution: {len(null_stats):,} NTC test statistics")

    if len(null_stats) == 0:
        logger.error("No valid NTC test statistics found. Cannot calibrate.")
        return

    # --- Compute empirical p-values for ALL rows ---
    test_stats = trans_df[stat_col].values

    if null_method == "ecdf":
        logger.info("Computing empirical p-values (eCDF method)")
        emp_pvals = empirical_pvals_ecdf(null_stats, test_stats, two_sided=two_sided)
    elif null_method == "t-fit":
        logger.info("Computing empirical p-values (t-distribution fit)")
        emp_pvals = empirical_pvals_tfit(null_stats, test_stats, two_sided=two_sided)
    else:
        raise ValueError(f"Unknown null_method: {null_method}")

    trans_df["empirical_pval"] = emp_pvals

    # Set NaN for rows where test statistic was NaN
    trans_df.loc[trans_df[stat_col].isna(), "empirical_pval"] = np.nan

    # --- BH correction on the discovery set (targeting elements only) ---
    logger.info("Applying BH FDR correction (targeting elements only)")
    trans_df["empirical_pval_adj"] = np.nan

    # Only correct targeting element tests (NTC tests are not part of discovery)
    discovery_mask = (~is_ntc) & trans_df["empirical_pval"].notna()
    if discovery_mask.sum() > 0:
        discovery_pvals = trans_df.loc[discovery_mask, "empirical_pval"].values
        _, pvals_adj, _, _ = multipletests(discovery_pvals, alpha=0.05, method="fdr_bh")
        trans_df.loc[discovery_mask, "empirical_pval_adj"] = pvals_adj

    n_sig_01 = (trans_df["empirical_pval_adj"] < 0.10).sum()
    n_sig_05 = (trans_df["empirical_pval_adj"] < 0.05).sum()
    logger.info(f"  FDR < 0.10: {n_sig_01:,} tests")
    logger.info(f"  FDR < 0.05: {n_sig_05:,} tests")

    # --- Annotate ---
    logger.info("Annotating results")

    # Element metadata
    trans_df["element_symbol"] = trans_df["intended_target_name"].map(element_symbol_map)
    trans_df["element_label"] = trans_df["intended_target_name"].map(element_type_map)

    # Tested gene symbol
    trans_df["tested_gene_symbol"] = trans_df["gene_id"].map(gene_symbol_map)

    # Number of cells per element — match on coordinates if available
    coord_cols = ["intended_target_name", "intended_target_chr",
                  "intended_target_start", "intended_target_end"]
    if isinstance(n_cells_per_element.index, pd.MultiIndex) and all(c in trans_df.columns for c in coord_cols):
        # Merge on coordinate tuple
        n_cells_lookup = n_cells_per_element.reset_index()
        n_cells_lookup.columns = coord_cols + ["n_cells"]
        trans_df = trans_df.merge(n_cells_lookup, on=coord_cols, how="left")
    else:
        trans_df["n_cells"] = trans_df["intended_target_name"].map(n_cells_per_element)

    # Cis / direct target flags
    trans_df["is_direct_target"] = annotate_direct_target(trans_df)
    trans_df["is_cis"] = annotate_cis(trans_df, gene_var, cis_window=cis_window)

    # Direct target is always cis by definition
    trans_df.loc[trans_df["is_direct_target"], "is_cis"] = True

    logger.info(f"  Direct target tests: {trans_df['is_direct_target'].sum():,}")
    logger.info(f"  Cis tests: {trans_df['is_cis'].sum():,}")

    # --- Build output table ---
    # Construct unique element_id: ENSG...|chr:start-end
    # Multiple genomic elements can target the same gene (e.g. different promoters)
    if "intended_target_chr" in trans_df.columns:
        trans_df["element_id"] = (
            trans_df["intended_target_name"] + "|"
            + trans_df["intended_target_chr"] + ":"
            + trans_df["intended_target_start"].astype(str) + "-"
            + trans_df["intended_target_end"].astype(str)
        )
    else:
        # No coordinates available — fall back to just the gene ID
        trans_df["element_id"] = trans_df["intended_target_name"]
        logger.warning("No coordinate columns found; element_id will not be unique for multi-element targets")

    # Rename columns to match plan schema
    output_df = trans_df.rename(columns={
        "gene_id": "tested_gene_id",
        "log2_fc": "log2fc",
        "log2_fc_std": "log2fc_se",
        "p_value": "posterior_pval",
    })

    output_cols = [
        "element_id",
        "element_symbol",
        "element_label",
        "tested_gene_id",
        "tested_gene_symbol",
        "n_cells",
        "log2fc",
        "log2fc_se",
        "is_cis",
        "is_direct_target",
        "posterior_pval",
        "empirical_pval",
        "empirical_pval_adj",
    ]
    # Keep only columns that exist
    output_cols = [c for c in output_cols if c in output_df.columns]
    output_df = output_df[output_cols]

    # --- Write output tables ---
    # (a) All calibrated results
    all_path = os.path.join(outdir, f"{prefix}_calibrated_all_results.tsv")
    output_df.to_csv(all_path, sep="\t", index=False)
    logger.info(f"Written: {all_path} ({len(output_df):,} rows)")

    # (b) Direct target results
    direct_df = output_df[output_df["is_direct_target"]].copy()
    direct_path = os.path.join(outdir, f"{prefix}_calibrated_direct_target_results.tsv")
    direct_df.to_csv(direct_path, sep="\t", index=False)
    logger.info(f"Written: {direct_path} ({len(direct_df):,} rows)")

    # (c) Cis results
    cis_out_df = output_df[output_df["is_cis"]].copy()
    cis_out_path = os.path.join(outdir, f"{prefix}_calibrated_cis_results.tsv")
    cis_out_df.to_csv(cis_out_path, sep="\t", index=False)
    logger.info(f"Written: {cis_out_path} ({len(cis_out_df):,} rows)")

    # (d) Trans results (excluding cis)
    trans_df_out = output_df[~output_df["is_cis"]].copy()
    trans_path = os.path.join(outdir, f"{prefix}_calibrated_trans_results.tsv")
    trans_df_out.to_csv(trans_path, sep="\t", index=False)
    logger.info(f"Written: {trans_path} ({len(trans_df_out):,} rows)")

    # --- Summary ---
    logger.info("=" * 50)
    logger.info("Calibration Summary:")
    logger.info(f"  Method: {null_method}")
    logger.info(f"  Total tests: {len(output_df):,}")
    logger.info(f"  NTC null size: {len(null_stats):,}")
    logger.info(f"  Direct target tests: {len(direct_df):,}")
    logger.info(f"  Cis tests: {len(cis_out_df):,}")
    logger.info(f"  Trans tests: {len(trans_df_out):,}")
    logger.info(f"  FDR < 0.10: {n_sig_01:,}")
    logger.info(f"  FDR < 0.05: {n_sig_05:,}")
    logger.info("=" * 50)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calibrate PerTurbo per-element results using empirical null."
    )
    parser.add_argument(
        "--trans-results", required=True,
        help="Path to perturbo_trans_per_element_output.tsv.gz",
    )
    parser.add_argument(
        "--mudata", required=True,
        help="Path to inference_mudata.h5mu",
    )
    parser.add_argument(
        "--outdir", "-o", required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--prefix", required=True,
        help="Prefix for output filenames (e.g., dataset name)",
    )
    parser.add_argument(
        "--null-method", default="t-fit", choices=["ecdf", "t-fit"],
        help="Calibration method (default: t-fit)",
    )
    parser.add_argument(
        "--cis-window", type=int, default=100_000,
        help="Distance (bp) for cis annotation (default: 100000)",
    )
    parser.add_argument(
        "--non-targeting-label", default="non_targeting",
        help="Label for NTC elements in guide.var['type'] (default: non_targeting)",
    )
    parser.add_argument(
        "--guide-mod-key", default="guide",
        help="Modality key for guide data (default: guide)",
    )
    parser.add_argument(
        "--one-sided", action="store_true",
        help="Compute one-sided p-values (default: two-sided)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_calibration(
        trans_results_path=args.trans_results,
        mudata_path=args.mudata,
        outdir=args.outdir,
        prefix=args.prefix,
        null_method=args.null_method,
        cis_window=args.cis_window,
        non_targeting_label=args.non_targeting_label,
        guide_mod_key=args.guide_mod_key,
        two_sided=not args.one_sided,
    )
