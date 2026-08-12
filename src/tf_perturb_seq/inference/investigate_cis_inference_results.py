#!/usr/bin/env python3
"""
investigate_cis_inference_results.py

Investigates cis inference results from SCEPTRE and Perturbo on a TF perturb-seq screen.
Produces three analyses:
  1. Guide efficiency: how many guides per target perturbed the intended target gene (p < 0.05)
  2. Top downregulated TFs by magnitude of knockdown (cis, significant)
  3. UMAP of cells colored by SERPINA1 and APOE (hepatocyte markers) expression

File format notes (from actual headers):
  mergedresults/per_guide_output.tsv.gz:
      gene_id | guide_id | sceptre_log2_fc | sceptre_p_value | perturbo_log2_fc | perturbo_p_value
  mergedresults/per_element_output.tsv.gz:
      gene_id | intended_target_name | sceptre_log2_fc | sceptre_p_value | perturbo_log2_fc | perturbo_p_value
  inference/sceptre_per_guide_output.tsv.gz:
      gene_id | guide_id | p_value | log2_fc
  inference/sceptre_per_element_output.tsv.gz:
      gene_id | intended_target_name | p_value | log2_fc

  - All files have DUPLICATE ROWS (same result repeated per guide/element); must drop_duplicates.
  - No pre-computed adjusted p-values in files; BH correction is applied here over ALL tests
    (full genome-wide universe) before subsetting to cis pairs.
  - guide_id format: "GENENAME#chr:start-end(strand)"  e.g. "ZNF91#chr19:23395447-23395465(-)"
  - Cis pairs: gene_id == intended_target_name (both Ensembl IDs) in element files.
  - For guide files: cis means the guide's target gene (from guide_id prefix) == gene_id.
    We use the element-level files to map Ensembl ID -> gene symbol for labelling.
  - MuData modalities: "gene" (var_names=Ensembl, var["symbol"]=gene symbol, log-normalized)
    and "guide" (layers["guide_assignment"]=binary); we extract mdata["gene"] for the UMAP.

Usage:
    python investigate_cis_inference_results.py [--downsample N] [--output_dir ./figures]
    python investigate_cis_inference_results.py --use_sceptre_only  # use separate sceptre files
    python investigate_cis_inference_results.py --run_umap           # run UMAP analysis
"""

import argparse
import glob
import os
import re
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

warnings.filterwarnings("ignore")

# --- Default paths (overridden by CLI arguments) ------------------------------
_DEFAULT_RESULTS_DIR          = (
    "/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/"
    "helen_output_poolabcd_prod_v10"
)
_DEFAULT_PIPELINE_OUTPUTS_DIR = _DEFAULT_RESULTS_DIR + "/pipeline_outputs"
_DEFAULT_INFERENCE_DIR        = _DEFAULT_RESULTS_DIR + "/inference"
_DEFAULT_GUIDE_META           = (
    "/hpc/group/gersbachlab/seg95/crispr-pipeline-personal/"
    "example-data/production_guide_metadata_v6.tsv"
)
_DEFAULT_MUDATA_NAME          = "inference_mudata.h5mu"

P_THRESH           = 0.05
DOWNSAMPLE_DEFAULT = 50_000


def _find_calibrated_file(directory: str, suffix: str) -> str | None:
    """
    Auto-detect a calibrated result file in `directory` by suffix convention.
    Looks for any file matching *{suffix} (case-insensitive).
    Returns the full path of the first match, or None if not found.
    """
    pattern = os.path.join(directory, f"*{suffix}")
    matches = glob.glob(pattern)
    if not matches:
        return None
    if len(matches) > 1:
        print(f"  WARNING: multiple files match *{suffix} in {directory}; "
              f"using first: {os.path.basename(matches[0])}")
    return matches[0]


def _build_paths(args):
    """Construct all file paths from the root directory arguments."""
    global PIPELINE_OUTPUTS_DIR, RESULTS_DIR, INFERENCE_DIR
    global MERGED_CIS_GUIDE, MERGED_CIS_ELEM
    global CALIB_DIRECT, CALIB_CIS
    global SCEPTRE_GUIDE, SCEPTRE_ELEM
    global MUDATA_PATH, RNA_H5AD
    global GUIDE_METADATA_PATH

    RESULTS_DIR          = args.results_dir
    PIPELINE_OUTPUTS_DIR = args.pipeline_outputs_dir
    INFERENCE_DIR        = args.inference_dir
    GUIDE_METADATA_PATH  = args.guide_metadata

    MERGED_CIS_GUIDE = os.path.join(RESULTS_DIR,          "cis_per_guide_results.tsv.gz")
    MERGED_CIS_ELEM  = os.path.join(RESULTS_DIR,          "cis_per_element_results.tsv.gz")
    SCEPTRE_GUIDE    = os.path.join(INFERENCE_DIR,        "sceptre_per_guide_output.tsv.gz")
    SCEPTRE_ELEM     = os.path.join(INFERENCE_DIR,        "sceptre_per_element_output.tsv.gz")
    MUDATA_PATH      = os.path.join(PIPELINE_OUTPUTS_DIR, args.mudata_name)
    RNA_H5AD         = os.path.join(RESULTS_DIR,          "mergedresults", "rna.h5ad")

    # Auto-detect calibrated files from the inference directory
    calib_dir = INFERENCE_DIR
    CALIB_DIRECT = _find_calibrated_file(calib_dir, "_calibrated_direct_target_results.tsv")
    CALIB_CIS    = _find_calibrated_file(calib_dir, "_calibrated_cis_results.tsv")

    if CALIB_DIRECT:
        print(f"  Auto-detected calibrated direct target file: {os.path.basename(CALIB_DIRECT)}")
    else:
        print(f"  WARNING: no *_calibrated_direct_target_results.tsv found in {calib_dir}")
    if CALIB_CIS:
        print(f"  Auto-detected calibrated cis file:           {os.path.basename(CALIB_CIS)}")
    else:
        print(f"  WARNING: no *_calibrated_cis_results.tsv found in {calib_dir}")


# --- CLI ----------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--downsample", type=int, default=DOWNSAMPLE_DEFAULT,
        help=f"Max cells for UMAP (default {DOWNSAMPLE_DEFAULT:,}; 0 = all cells)",
    )
    p.add_argument(
        "--output_dir", default="./figures",
        help="Directory to save figures and tables (default: ./figures)",
    )
    p.add_argument(
        "--top_n_tfs", type=int, default=20,
        help="Number of top downregulated TFs to show (default: 20)",
    )
    p.add_argument(
        "--use_sceptre_only", action="store_true", default=False,
        help="Use the standalone SCEPTRE files (inference/) instead of merged files",
    )
    p.add_argument(
        "--run_umap", action="store_true", default=False,
        help="Run UMAP analysis (Analysis 3); skipped by default as it is slow",
    )
    p.add_argument(
        "--umap_precomputed", action="store_true", default=False,
        help="Re-use X_umap already stored in the AnnData object (skips recomputation)",
    )
    p.add_argument(
        "--p_thresh", type=float, default=P_THRESH,
        help=f"P-value threshold for significance (default: {P_THRESH})",
    )
    # --- Path arguments -------------------------------------------------------
    p.add_argument(
        "--results_dir",
        default=_DEFAULT_RESULTS_DIR,
        help="Root results directory containing cis_per_*_results.tsv.gz",
    )
    p.add_argument(
        "--pipeline_outputs_dir",
        default=_DEFAULT_PIPELINE_OUTPUTS_DIR,
        help="Directory containing the MuData .h5mu file",
    )
    p.add_argument(
        "--inference_dir",
        default=_DEFAULT_INFERENCE_DIR,
        help="Directory containing sceptre_per_*_output.tsv.gz and "
             "*_calibrated_*.tsv files (auto-detected by suffix)",
    )
    p.add_argument(
        "--mudata_name",
        default=_DEFAULT_MUDATA_NAME,
        help=f"Filename of the MuData .h5mu file inside pipeline_outputs_dir "
             f"(default: {_DEFAULT_MUDATA_NAME})",
    )
    p.add_argument(
        "--guide_metadata",
        default=_DEFAULT_GUIDE_META,
        help="Path to guide metadata TSV (guide_id, type, intended_target_name, ...)",
    )
    # --calib_prefix kept as a no-op for backwards compatibility but ignored
    p.add_argument(
        "--calib_prefix", default=None,
        help="[DEPRECATED] Calibrated file prefix is now auto-detected; this argument is ignored.",
    )
    return p.parse_args()


# --- Helpers ------------------------------------------------------------------
def load_tsv(path: str, label: str) -> pd.DataFrame:
    print(f"  Loading {label}:\n    {path}")
    df = pd.read_csv(path, sep="\t", compression="gzip")
    n_raw = len(df)
    df = df.drop_duplicates()
    print(f"    {n_raw:,} rows -> {len(df):,} after drop_duplicates  |  cols: {list(df.columns)}")
    return df


def parse_target_from_guide_id(guide_id: str) -> str:
    """Extract gene symbol from 'GENENAME#chr...' guide ID format."""
    return guide_id.split("#")[0]


def bh_correct(df: pd.DataFrame, p_cols: list[str]) -> pd.DataFrame:
    """
    Add Benjamini-Hochberg FDR columns to df for each raw p-value column in p_cols.
    Correction is applied over ALL rows in df (the full multiple-testing universe),
    not just cis pairs -- this is intentional: the FDR reflects the true testing burden.

    New columns are named by appending '_bh' to each input column name.
    Requires statsmodels; falls back to a manual BH implementation if unavailable.
    """
    from scipy.stats import rankdata

    def _bh(pvals: np.ndarray) -> np.ndarray:
        """Vectorised BH correction, handling NaNs."""
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
        print(f"    {col}: {n_sig_raw:,} raw p<0.05  ->  {n_sig_bh:,} BH-adjusted FDR<0.05"
              f"  (over {len(df):,} tests)")
    return df


def build_ensg_to_symbol_mygene(ensg_ids, verbose=True) -> dict:
    """
    Map Ensembl gene IDs to gene symbols using a local two-column TSV file.
    """
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


def build_ensg_to_symbol(elem_df) -> dict:
    """Stub kept for API compatibility."""
    return {}


def build_ensg_to_symbol_from_mudata() -> dict:
    """
    Fetch Ensembl->symbol map for all genes in the guide file.
    """
    ensg_ids = set()
    for path in [MERGED_CIS_GUIDE, SCEPTRE_GUIDE]:
        if os.path.exists(path):
            try:
                df = pd.read_csv(path, sep="\t", compression="gzip",
                                 usecols=["gene_id"]).drop_duplicates()
                ensg_ids.update(df["gene_id"].dropna().tolist())
            except Exception:
                pass
    return build_ensg_to_symbol_mygene(list(ensg_ids))


# --- Analysis 1 - Guide efficiency -------------------------------------------
def guide_efficiency(args, out_dir: str) -> dict:
    """
    For each target, count: (a) total unique guides tested, (b) guides with p < threshold.

    Intermediate TSV output (for R/ggplot):
        guide_efficiency_{method}.csv  -- one row per target with total_guides, sig_guides, frac_sig
    """
    print("\n=== Analysis 1: Guide efficiency ===")

    if args.use_sceptre_only:
        guide_path = SCEPTRE_GUIDE
        method_cols = {"SCEPTRE": "p_value"}
    else:
        guide_path = MERGED_CIS_GUIDE
        method_cols = {
            "SCEPTRE":  "sceptre_p_value",
            "Perturbo": "perturbo_p_value",
        }

    if not os.path.exists(guide_path):
        print(f"  Skipping: per-guide file not found at {guide_path}")
        print("  (Only element-level calibrated outputs are available for this dataset.)")
        return {}

    guide_df = load_tsv(guide_path, "per-guide")

    print("  Building Ensembl->symbol map from inference MuData ...")
    ensg2sym = build_ensg_to_symbol_from_mudata()

    guide_df["target_symbol"] = guide_df["guide_id"].apply(parse_target_from_guide_id)
    guide_df["gene_symbol"]   = guide_df["gene_id"].map(ensg2sym).fillna(guide_df["gene_id"])

    print("  Applying BH correction over all guide-level tests ...")
    p_cols_guide = list(method_cols.values())
    guide_df = bh_correct(guide_df, p_cols_guide)
    method_padj = {m: (pc + "_bh") for m, pc in method_cols.items()}

    cis = guide_df[
        guide_df["target_symbol"].str.upper() == guide_df["gene_symbol"].str.upper()
    ].copy()
    print(f"  Cis guide-gene pairs (guide target == tested gene): {len(cis):,}")

    results = {}
    for method, p_col in method_cols.items():
        padj_col = method_padj[method]
        if padj_col not in cis.columns:
            print(f"  WARNING: column '{padj_col}' not found, skipping {method}")
            continue

        df_m = cis[["target_symbol", "guide_id", p_col, padj_col]].dropna(subset=[padj_col]).copy()
        df_m = (
            df_m.groupby(["target_symbol", "guide_id"], as_index=False)
                .agg({p_col: "min", padj_col: "min"})
        )
        df_m["significant"] = df_m[padj_col] < args.p_thresh

        summary = (
            df_m.groupby("target_symbol")
                .agg(
                    total_guides=("guide_id",    "nunique"),
                    sig_guides  =("significant", "sum"),
                )
                .reset_index()
        )
        summary["frac_sig"] = summary["sig_guides"] / summary["total_guides"]
        summary["method"]   = method
        summary["p_thresh"] = args.p_thresh
        results[method] = summary

        avg_all     = summary["total_guides"].mean()
        avg_sig     = summary["sig_guides"].mean()
        pct_any_sig = (summary["sig_guides"] > 0).mean() * 100

        print(f"\n  [{method}] BH FDR < {args.p_thresh}")
        print(f"    Unique targets:                  {len(summary):,}")
        print(f"    Avg guides per target:           {avg_all:.2f}")
        print(f"    Avg significant guides/target:   {avg_sig:.2f}")
        print(f"    Targets with >=1 sig guide:       {pct_any_sig:.1f}%")
        top10 = summary.sort_values("sig_guides", ascending=False).head(10)
        print(top10.to_string(index=False))

    # -- Plot ------------------------------------------------------------------
    n_methods = len(results)
    fig, axes = plt.subplots(1, n_methods, figsize=(7 * n_methods, 5), squeeze=False)
    fig.suptitle(
        f"Guide efficiency: significant cis knockdown per target  (BH FDR < {args.p_thresh})",
        fontsize=13, y=1.02,
    )

    for ax, (method, summary) in zip(axes[0], results.items()):
        counts = summary["sig_guides"].value_counts().sort_index()
        max_guides = int(summary["total_guides"].max())
        all_counts = pd.Series(0, index=range(0, max_guides + 1))
        all_counts.update(counts)

        ax.bar(all_counts.index, all_counts.values,
               color="#4C72B0", edgecolor="white", linewidth=0.5)

        avg = summary["sig_guides"].mean()
        ax.axvline(avg, color="firebrick", linestyle="--", linewidth=1.5,
                   label=f"Mean = {avg:.2f}")
        ax.set_xlabel("Number of significant guides per target")
        ax.set_ylabel("Number of targets")
        ax.set_title(method, fontsize=11)
        ax.legend(fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    out_path = os.path.join(out_dir, "guide_efficiency.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\n  Saved: {out_path}")

    # -- Intermediate TSV outputs for R/ggplot ---------------------------------
    for method, summary in results.items():
        csv_path = os.path.join(out_dir, f"guide_efficiency_{method.lower()}.csv")
        summary.to_csv(csv_path, index=False)
        print(f"  Saved R-ready table: {csv_path}")

    # Combined long-format TSV (all methods, one file) for easy ggplot faceting
    if results:
        combined = pd.concat(results.values(), ignore_index=True)
        combined_path = os.path.join(out_dir, "guide_efficiency_all_methods.tsv")
        combined.to_csv(combined_path, sep="\t", index=False)
        print(f"  Saved R-ready combined table: {combined_path}")

    return results


# --- Analysis 2 - Top downregulated TFs --------------------------------------
def top_downregulated_tfs(args, out_dir: str):
    """
    Uses the element-level file (one row per elementxgene pair after dedup).
    Cis = gene_id == intended_target_name (both Ensembl IDs).
    Shows top N most downregulated significant cis hits for each method.

    Intermediate TSV outputs (for R/ggplot):
        top_downregulated_tfs_{method}.tsv  -- top N rows with symbol, lfc, padj, neg_log10_fdr
        cis_element_all_{method}.tsv        -- ALL significant cis hits (not just top N)
    """
    print("\n=== Analysis 2: Top downregulated TFs ===")

    if args.use_sceptre_only:
        if not os.path.exists(SCEPTRE_ELEM):
            print(f"  Skipping SCEPTRE element analysis: file not found at {SCEPTRE_ELEM}")
            elem_df = None
            method_specs = []
        else:
            elem_df = load_tsv(SCEPTRE_ELEM, "SCEPTRE per-element")
            method_specs = [("SCEPTRE", "p_value", "log2_fc")]
        calib_direct_df = None
        calib_cis_df    = None
    else:
        if os.path.exists(MERGED_CIS_ELEM):
            elem_df = load_tsv(MERGED_CIS_ELEM, "Merged cis per-element")
            method_specs = [("SCEPTRE", "sceptre_p_value", "sceptre_log2_fc")]
        else:
            print(f"  No merged element file found at {MERGED_CIS_ELEM}; "
                  f"will use calibrated outputs only.")
            elem_df = None
            method_specs = []

        calib_direct_df = None
        calib_cis_df    = None
        if CALIB_DIRECT and os.path.exists(CALIB_DIRECT):
            print("  Loading calibrated Perturbo direct target results ...")
            calib_direct_df = pd.read_csv(CALIB_DIRECT, sep="\t")
            print(f"    {len(calib_direct_df):,} rows")
        else:
            print("  WARNING: calibrated direct target file not found; skipping Perturbo (direct)")

        if CALIB_CIS and os.path.exists(CALIB_CIS):
            print("  Loading calibrated Perturbo cis results ...")
            calib_cis_df = pd.read_csv(CALIB_CIS, sep="\t")
            print(f"    {len(calib_cis_df):,} rows")
        else:
            print("  WARNING: calibrated cis file not found; skipping Perturbo (cis 100kb)")

    if elem_df is not None:
        print("  Applying BH correction over all element-level tests ...")
        p_cols_elem = [spec[1] for spec in method_specs]
        elem_df = bh_correct(elem_df, p_cols_elem)
        method_specs_bh = [(m, pc, pc + "_bh", lfc) for m, pc, lfc in method_specs]
        cis = elem_df[elem_df["gene_id"] == elem_df["intended_target_name"]].copy()
        print(f"  Cis pairs (gene_id == intended_target_name): {len(cis):,}")
    else:
        method_specs_bh = []
        cis = pd.DataFrame()

    print("  Building Ensembl->symbol map from inference MuData ...")
    ensg2sym = build_ensg_to_symbol_from_mudata()
    if not cis.empty:
        cis["symbol"] = cis["gene_id"].map(ensg2sym).fillna(cis["gene_id"])

    calib_dfs = {}
    if not args.use_sceptre_only:
        for calib_label, calib_df, cis_flag, direct_flag in [
            ("Perturbo (direct target)", calib_direct_df, None,   True),
            ("Perturbo (cis 100kb)",     calib_cis_df,    True,   None),
        ]:
            if calib_df is None:
                continue
            calib_df = calib_df.copy()

            # Use pre-computed is_cis / is_direct_target flags where available
            # rather than re-deriving cis from gene_id == intended_target_name.
            if direct_flag is not None and "is_direct_target" in calib_df.columns:
                n_before = len(calib_df)
                calib_df = calib_df[calib_df["is_direct_target"] == True].copy()
                print(f"  [{calib_label}] is_direct_target filter: "
                      f"{n_before:,} -> {len(calib_df):,} rows")
            elif cis_flag is not None and "is_cis" in calib_df.columns:
                n_before = len(calib_df)
                calib_df = calib_df[calib_df["is_cis"] == True].copy()
                print(f"  [{calib_label}] is_cis filter: "
                      f"{n_before:,} -> {len(calib_df):,} rows")

            # Normalise column names to match the rest of the pipeline
            calib_df["symbol"]  = calib_df["element_symbol"]
            calib_df["gene_id"] = calib_df["tested_gene_id"]
            method_specs_bh.append(
                (calib_label, "empirical_pval", "empirical_pval_adj", "log2fc")
            )
            calib_dfs[calib_label] = calib_df

    # -- Plot (one panel per method, side by side) -----------------------------
    n_methods = len(method_specs_bh)
    fig, axes = plt.subplots(1, n_methods, figsize=(9 * n_methods, max(6, args.top_n_tfs * 0.38)))
    if n_methods == 1:
        axes = [axes]

    all_top = {}
    for ax, (method, p_col, padj_col, lfc_col) in zip(axes, method_specs_bh):
        if method in calib_dfs:
            df_method = calib_dfs[method]
        else:
            df_method = cis

        if padj_col not in df_method.columns or lfc_col not in df_method.columns:
            print(f"  WARNING: columns not found for {method}; skipping")
            ax.set_visible(False)
            continue

        sig = df_method[(df_method[padj_col] < args.p_thresh) & (df_method[lfc_col] < 0)].copy()
        print(f"\n  [{method}] Significant downregulated cis TFs (BH FDR < {args.p_thresh}): {len(sig):,}")

        top_n = sig.nsmallest(args.top_n_tfs, padj_col).copy()
        if "symbol" not in top_n.columns:
            top_n["symbol"] = top_n["gene_id"].map(ensg2sym).fillna(top_n["gene_id"])
        top_n["neg_log10_fdr"] = -np.log10(top_n[padj_col].clip(lower=1e-300))
        top_n["method"]        = method
        all_top[method]        = top_n
        print(top_n[["symbol", lfc_col, padj_col]].to_string(index=False))

        lfc_min = top_n[lfc_col].min()
        lfc_max = top_n[lfc_col].max()
        norm  = mcolors.Normalize(vmin=lfc_min, vmax=min(lfc_max, 0))
        cmap  = plt.cm.Blues_r

        top_n_sorted = top_n.sort_values(padj_col, ascending=False)
        labels = top_n_sorted["symbol"].values
        bar_lengths = top_n_sorted["neg_log10_fdr"].values
        lfcs   = top_n_sorted[lfc_col].values
        pvals  = top_n_sorted[padj_col].values
        clrs   = [cmap(norm(v)) for v in lfcs]

        bars = ax.barh(labels, bar_lengths, color=clrs, edgecolor="white", linewidth=0.4)
        ax.set_xlabel("-log10(BH FDR)", fontsize=11)
        ax.set_title(
            f"{method}: top {len(top_n)} downregulated TFs\n"
            f"(cis, BH FDR < {args.p_thresh}, ranked by FDR)",
            fontsize=11,
        )
        ax.spines[["top", "right"]].set_visible(False)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label="log2 fold change", shrink=0.5, pad=0.01)

    fig.suptitle("Top downregulated TFs (cis knockdown)", fontsize=13, y=1.01)
    plt.tight_layout()
    out_path = os.path.join(out_dir, "top_downregulated_tfs.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\n  Saved: {out_path}")

    # -- Intermediate TSV outputs for R/ggplot ---------------------------------
    for method, top_n in all_top.items():
        safe = method.lower().replace(" ", "_").replace("(", "").replace(")", "")
        # Top N table
        tsv_path = os.path.join(out_dir, f"top_downregulated_tfs_{safe}.tsv")
        top_n.to_csv(tsv_path, sep="\t", index=False)
        print(f"  Saved R-ready table: {tsv_path}")

    # Combined long-format TSV for faceting
    if all_top:
        combined = pd.concat(all_top.values(), ignore_index=True)
        combined_path = os.path.join(out_dir, "top_downregulated_tfs_all_methods.tsv")
        combined.to_csv(combined_path, sep="\t", index=False)
        print(f"  Saved R-ready combined table: {combined_path}")

    # Full significant cis hits table (all methods, not just top N)
    all_sig_frames = []
    for method, p_col, padj_col, lfc_col in method_specs_bh:
        df_method = calib_dfs[method] if method in calib_dfs else cis
        if padj_col not in df_method.columns or lfc_col not in df_method.columns:
            continue
        sig_all = df_method[df_method[padj_col] < args.p_thresh].copy()
        if "symbol" not in sig_all.columns:
            sig_all["symbol"] = sig_all["gene_id"].map(ensg2sym).fillna(sig_all["gene_id"])
        sig_all["neg_log10_fdr"] = -np.log10(sig_all[padj_col].clip(lower=1e-300))
        sig_all["method"] = method
        # Standardise column names for easy R use; include extra calibrated columns
        # where available (n_cells, log2fc_se, tested_gene_symbol, posterior_pval)
        base_cols = ["method", "symbol", "gene_id", lfc_col, padj_col, "neg_log10_fdr"]
        extra_cols = [c for c in
                      ["tested_gene_symbol", "n_cells", "log2fc_se", "posterior_pval",
                       "is_cis", "is_direct_target"]
                      if c in sig_all.columns]
        sig_out = sig_all[base_cols + extra_cols].copy()
        sig_out = sig_out.rename(columns={lfc_col: "log2fc", padj_col: "padj"})
        all_sig_frames.append(sig_out)

    if all_sig_frames:
        all_sig = pd.concat(all_sig_frames, ignore_index=True)
        all_sig_path = os.path.join(out_dir, "cis_significant_all_methods.tsv")
        all_sig.to_csv(all_sig_path, sep="\t", index=False)
        print(f"  Saved R-ready full significant hits: {all_sig_path}")

    return all_top


# --- Analysis 3 - UMAP colored by hepatocyte marker expression ---------------
def umap_hepatocyte_markers(args, out_dir: str):
    """
    Loads the gene modality from the inference MuData and plots UMAPs colored
    by hepatocyte marker expression.

    Intermediate TSV outputs (for R/ggplot):
        umap_coordinates.tsv      -- cell barcode, UMAP1, UMAP2
        umap_marker_expression.tsv -- cell barcode + one column per marker gene
    """
    print("\n=== Analysis 3: UMAP colored by hepatocyte markers ===")

    try:
        import mudata as md
        import scanpy as sc
        import scipy.sparse as sp
    except ImportError:
        print("  ERROR: mudata/scanpy not installed.")
        return

    if not os.path.exists(MUDATA_PATH):
        print(f"  ERROR: {MUDATA_PATH} not found")
        return

    markers = {
        "ALB":      "ENSG00000163631",
        "SERPINA1": "ENSG00000197249",
        "APOE":     "ENSG00000130203",
        "AFP":      "ENSG00000081051",
        "TTR":      "ENSG00000118271",
        "FGB":      "ENSG00000171557",
        "ARG1":     "ENSG00000118520",
    }

    print(f"  Loading gene modality (backed) from {MUDATA_PATH} ...")
    import scipy.sparse as _sp

    mdata = md.read_h5mu(MUDATA_PATH, backed="r")
    adata = mdata["gene"]
    n_obs = adata.n_obs
    print(f"  Gene modality: {n_obs:,} cells x {adata.n_vars:,} genes (backed)")

    all_ensg   = list(adata.var_names)
    ensg_index = {e: i for i, e in enumerate(all_ensg)}
    sym_index  = {}
    if "symbol" in adata.var.columns:
        sym_index = {s: i for i, s in enumerate(adata.var["symbol"])}

    marker_col_idx = {}
    marker_keys    = {}
    for name, ensg in markers.items():
        if ensg in ensg_index:
            marker_col_idx[name] = ensg_index[ensg]
            marker_keys[name]    = ensg
            print(f"  Found {name} via Ensembl ID {ensg} (col {ensg_index[ensg]})")
        elif name.upper() in {s.upper() for s in sym_index}:
            sym_match = next(s for s in sym_index if s.upper() == name.upper())
            marker_col_idx[name] = sym_index[sym_match]
            marker_keys[name]    = sym_match
            print(f"  Found {name} via symbol match (col {sym_index[sym_match]})")
        else:
            marker_keys[name] = None
            print(f"  WARNING: {name} ({ensg}) not found in gene modality")

    print(f"  Extracting {len(marker_col_idx)} marker gene columns from X ...")

    exprs_full = {}
    for name, col_idx in marker_col_idx.items():
        X   = adata.X
        col = X[:, col_idx]
        if _sp.issparse(col):
            vec = np.asarray(col.todense()).flatten()
        elif hasattr(col, "A"):
            vec = col.A.flatten()
        else:
            vec = np.asarray(col).flatten()
        vec = vec.astype(float)
        exprs_full[name] = vec
        print(f"    {name}: mean={vec.mean():.3f}, nonzero={(vec > 0).sum():,}")

    exprs_full.update({
        name: np.zeros(n_obs)
        for name in markers if name not in exprs_full
    })

    np.random.seed(42)
    target_total = args.downsample if (args.downsample and args.downsample < n_obs) else n_obs

    if target_total >= n_obs:
        sampled_idx = np.arange(n_obs)
        print(f"  No downsampling needed ({n_obs:,} cells)")
    else:
        print(f"  Balanced downsample: target {target_total:,} cells from {n_obs:,} ...")

        guide_mod    = mdata["guide"]
        guide_layer  = guide_mod.layers["guide_assignment"]
        guide_var    = guide_mod.var
        guide_names  = list(guide_mod.var_names)

        if "targeting" in guide_var.columns:
            raw = guide_var["targeting"]
            is_tgt = (raw.str.upper().str.strip() == "TRUE") \
                     if raw.dtype == object else raw.astype(bool)
        else:
            is_tgt = pd.Series(True, index=guide_var.index)

        tgt_cols = [i for i, g in enumerate(guide_names) if is_tgt.iloc[i]]
        nt_cols  = [i for i, g in enumerate(guide_names) if not is_tgt.iloc[i]]

        if "intended_target_name" in guide_var.columns:
            tgt_col_to_target = {
                i: guide_var["intended_target_name"].iloc[i]
                for i in tgt_cols
            }
        else:
            tgt_col_to_target = {
                i: guide_names[i].split("#")[0]
                for i in tgt_cols
            }

        import scipy.sparse as _sp2
        gl = guide_layer if _sp2.issparse(guide_layer) else _sp2.csr_matrix(guide_layer)
        gl_csc = gl.tocsc()

        cell_target = np.full(n_obs, "non-targeting", dtype=object)
        for col_i in tgt_cols:
            col_data = gl_csc[:, col_i]
            assigned = col_data.nonzero()[0]
            target   = tgt_col_to_target[col_i]
            cell_target[assigned] = target

        target_to_cells = {}
        for cell_i, tgt in enumerate(cell_target):
            target_to_cells.setdefault(tgt, []).append(cell_i)

        rare_thresh   = 50
        rare_targets  = {t: c for t, c in target_to_cells.items() if len(c) <= rare_thresh}
        abund_targets = {t: c for t, c in target_to_cells.items() if len(c) > rare_thresh}

        n_rare   = sum(len(c) for c in rare_targets.values())
        n_budget = target_total - n_rare
        n_abund  = len(abund_targets)

        print(f"    Rare targets (<=50 cells): {len(rare_targets):,}  "
              f"| Abundant targets: {n_abund:,}  | Rare cells kept: {n_rare:,}")

        if n_budget <= 0:
            print("    WARNING: rare cells alone exceed target_total; keeping all rare cells only")
            sampled_idx = np.sort(np.concatenate([c for c in rare_targets.values()]))
        else:
            per_target = n_budget // max(1, n_abund)
            print(f"    Per-abundant-target budget: {per_target:,}")

            selected = list(np.concatenate(list(rare_targets.values())) if rare_targets else [])
            sample_sizes = {}
            for tgt, cells in abund_targets.items():
                n = min(len(cells), per_target)
                chosen = np.random.choice(cells, size=n, replace=False).tolist()
                selected.extend(chosen)
                sample_sizes[tgt] = n

            shortfall = target_total - len(selected)
            if shortfall > 0:
                print(f"    Redistributing {shortfall:,} leftover slots ...")
                spare = {t: len(c) - sample_sizes[t]
                         for t, c in abund_targets.items()
                         if len(c) > sample_sizes[t]}
                selected_set = set(selected)
                abund_list   = list(abund_targets.keys())
                while shortfall > 0 and spare:
                    for tgt in abund_list:
                        if shortfall <= 0 or tgt not in spare:
                            continue
                        available = [c for c in abund_targets[tgt] if c not in selected_set]
                        if not available:
                            del spare[tgt]
                            continue
                        extra = np.random.choice(available, size=1)[0]
                        selected.append(extra)
                        selected_set.add(extra)
                        sample_sizes[tgt] += 1
                        spare[tgt] -= 1
                        if spare[tgt] == 0:
                            del spare[tgt]
                        shortfall -= 1

            sampled_idx = np.sort(np.array(selected, dtype=int))
            print(f"    Final sample: {len(sampled_idx):,} cells "
                  f"across {len(target_to_cells):,} guide targets")

    exprs = {name: vec[sampled_idx] for name, vec in exprs_full.items()}
    del exprs_full

    print(f"  Loading expression matrix for {len(sampled_idx):,} sampled cells ...")
    import anndata as ad
    X_sub  = adata.X[sampled_idx, :]
    cell_barcodes = list(adata.obs_names[sampled_idx])
    adata  = ad.AnnData(X=X_sub, obs=adata.obs.iloc[sampled_idx].copy(),
                        var=adata.var.copy())
    mdata.file.close()
    print(f"  Subsetted AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # -- UMAP ------------------------------------------------------------------
    print("  Computing UMAP (data already log-normalized, skipping renormalization) ...")

    if adata.n_vars > 2000:
        print(f"  Selecting 2,000 HVGs from {adata.n_vars:,} genes ...")
        expressed = np.asarray((adata.X > 0).sum(axis=0)).flatten() > 0
        print(f"  Removing {(~expressed).sum():,} zero-count genes before HVG selection")
        adata_hvg = adata[:, expressed].copy()
        sc.pp.highly_variable_genes(adata_hvg, n_top_genes=min(2000, expressed.sum()),
                                    flavor="seurat_v3")
        hvg_mask = adata_hvg.var["highly_variable"]
        adata_pp  = adata_hvg[:, hvg_mask].copy()
    else:
        adata_pp = adata.copy()

    sc.pp.scale(adata_pp, max_value=10)
    sc.tl.pca(adata_pp, n_comps=50, svd_solver="arpack")
    sc.pp.neighbors(adata_pp, n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata_pp, min_dist=0.3)
    adata.obsm["X_umap"] = adata_pp.obsm["X_umap"]

    umap_xy = adata.obsm["X_umap"]

    # -- Intermediate TSV outputs for R/ggplot ---------------------------------
    # 1. UMAP coordinates
    umap_df = pd.DataFrame(
        umap_xy,
        index=cell_barcodes,
        columns=["UMAP1", "UMAP2"],
    )
    umap_df.index.name = "cell_barcode"

    coords_path = os.path.join(out_dir, "umap_coordinates.tsv")
    umap_df.to_csv(coords_path, sep="\t")
    print(f"  Saved R-ready UMAP coordinates: {coords_path}")

    # 2. Marker expression (long format: cell_barcode, marker, expression)
    expr_df = pd.DataFrame(exprs, index=cell_barcodes)
    expr_df.index.name = "cell_barcode"
    expr_wide_path = os.path.join(out_dir, "umap_marker_expression_wide.tsv")
    expr_df.to_csv(expr_wide_path, sep="\t")
    print(f"  Saved R-ready marker expression (wide): {expr_wide_path}")

    # Long format (easier for ggplot faceting)
    expr_long = expr_df.reset_index().melt(
        id_vars="cell_barcode", var_name="marker", value_name="log_norm_expr"
    )
    # Attach UMAP coords to long format
    expr_long = expr_long.merge(umap_df.reset_index(), on="cell_barcode")
    expr_long_path = os.path.join(out_dir, "umap_marker_expression_long.tsv")
    expr_long.to_csv(expr_long_path, sep="\t", index=False)
    print(f"  Saved R-ready marker expression (long, with UMAP coords): {expr_long_path}")

    # -- Plot UMAPs in a 2x3 grid (one per marker) ----------------------------
    n_markers = len(exprs)
    ncols = 3
    nrows = int(np.ceil(n_markers / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(8 * ncols, 7 * nrows), squeeze=False)
    fig.suptitle(
        f"Hepatocyte marker expression  (n = {adata.n_obs:,} cells)",
        fontsize=13, fontweight="bold",
    )

    for j in range(n_markers, nrows * ncols):
        axes[j // ncols][j % ncols].set_visible(False)

    for i, (name, expr) in enumerate(exprs.items()):
        ax = axes[i // ncols][i % ncols]
        order   = np.argsort(expr)
        scatter = ax.scatter(
            umap_xy[order, 0], umap_xy[order, 1],
            c=expr[order],
            cmap="viridis",
            s=1.5, alpha=0.6, linewidths=0, rasterized=True,
        )
        cbar = plt.colorbar(scatter, ax=ax, pad=0.01, shrink=0.8)
        cbar.set_label("log-normalized expression", fontsize=9)

        found = marker_keys[name] is not None
        ax.set_title(
            f"{name}" + ("" if found else " (not found)"),
            fontsize=12, fontweight="bold",
        )
        ax.set_xlabel("UMAP 1", fontsize=10)
        ax.set_ylabel("UMAP 2", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines[["top", "right", "left", "bottom"]].set_visible(False)

    plt.tight_layout()
    out_path = os.path.join(out_dir, "umap_hepatocyte_markers.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\n  Saved: {out_path}")


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
    Also writes an intermediate TSV for R/ggplot if out_dir is provided.
    """
    if guide_meta.empty:
        print(f"  Skipping label volcano (no guide metadata)")
        return

    keep_cols = [col for col in
                 [p_col, lfc_col, "guide_id", "intended_target_name",
                  "gene_id", "type", "tested_gene_id"]
                 if col in df_results.columns]
    df = df_results[keep_cols].copy().dropna(subset=[p_col, lfc_col])

    if title_suffix == "cis":
        if level == "guide" and "guide_id" in df.columns:
            meta_cols = (
                guide_meta[["guide_id", "intended_target_name", "type"]]
                .drop_duplicates("guide_id")
            )
            df = df.merge(meta_cols, on="guide_id", how="left")

            targeting_mask = df["type"].isin(["targeting"]) & df["intended_target_name"].notna()
            tgt_cis = df[targeting_mask & (df["gene_id"] == df["intended_target_name"])].copy()
            ctrl = df[~targeting_mask].copy()
            ctrl = ctrl.sort_values(p_col).drop_duplicates("guide_id", keep="first")

            n_before = df["guide_id"].nunique()
            df = pd.concat([tgt_cis, ctrl], ignore_index=True)
            print(f"  Cis filter (guide): {n_before:,} guides -> "
                  f"{len(df):,} rows "
                  f"(targeting: cis row; controls: min-p row)")

        elif level == "element" and "intended_target_name" in df.columns:
            n_before = len(df)
            df = df[df["gene_id"] == df["intended_target_name"]].copy()
            print(f"  Cis filter (element): {n_before:,} -> {len(df):,} rows")

    if level == "guide":
        if "type" not in df.columns:
            if "guide_id" in df.columns:
                df = df.merge(guide_meta[["guide_id", "type"]], on="guide_id", how="left")
            else:
                print(f"  WARNING: cannot merge type for guide level -- no guide_id col")
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
                print(f"  WARNING: cannot merge type for element level -- "
                      f"df cols: {list(df.columns)}")
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

    # -- Intermediate TSV for R/ggplot -----------------------------------------
    if out_dir is not None:
        safe_method = method.lower().replace(" ", "_").replace("(", "").replace(")", "")
        tsv_name = f"volcano_{title_suffix}_{level}_{safe_method}.tsv"
        tsv_path = os.path.join(out_dir, tsv_name)
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
               linewidth=0.9, linestyle="--", zorder=5,
               label=f"p = {p_thresh}")
    ax.axvline( fc_thresh, color="#666666", linewidth=0.9, linestyle="--", zorder=5)
    ax.axvline(-fc_thresh, color="#666666", linewidth=0.9, linestyle="--", zorder=5)
    ax.axvline(0,          color="#444444", linewidth=0.6, zorder=5)

    ax.set_xlabel("log2 Fold Change", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)

    present_types = sorted(df["type"].dropna().unique())
    absent_types  = [t for t in LABEL_ORDER if t not in present_types]
    absent_note   = f"  (no {', '.join(absent_types)} in this file)" if absent_types else ""

    ax.set_title(
        f"{method} {title_suffix} volcano  (per {level})\n"
        f"nominal p < {p_thresh} (dashed line), |log2FC| > {fc_thresh}{absent_note}",
        fontsize=10, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper left", fontsize=8, frameon=False,
              bbox_to_anchor=(1.01, 1))

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


def plot_cis_window_volcano(df_raw, p_col, lfc_col, method, out_path,
                             guide_meta, p_thresh=0.01, fc_thresh=0.5,
                             out_dir=None):
    """
    Element-level cis volcano showing ALL cis pairs.
    Also writes an intermediate TSV for R/ggplot if out_dir is provided.
    """
    keep_cols = [col for col in
                 [p_col, lfc_col, "gene_id", "intended_target_name",
                  "type", "tested_gene_id", "is_self_target"]
                 if col in df_raw.columns]
    df = df_raw[keep_cols].copy().dropna(subset=[p_col, lfc_col])

    if "type" not in df.columns:
        elem_labels = (
            guide_meta[["intended_target_name", "type"]]
            .dropna(subset=["intended_target_name"])
            .drop_duplicates("intended_target_name")
        )
        if "intended_target_name" in df.columns:
            df = df.merge(elem_labels, on="intended_target_name", how="left")
        else:
            print(f"  WARNING: cannot merge type for cis window volcano")
            return

    # Respect a pre-set is_self_target column (e.g. from calibrated files where
    # every row is a self-target by definition and the ID formats don't match).
    if "is_self_target" not in df.columns:
        itn_col = "intended_target_name"
        gid_col = "gene_id"
        if itn_col in df.columns and gid_col in df.columns:
            df["is_self_target"] = df[gid_col] == df[itn_col]
        else:
            df["is_self_target"] = False

    df["-log10p"] = -np.log10(df[p_col].clip(lower=1e-350))
    df.loc[~np.isfinite(df["-log10p"]), "-log10p"] = 350
    df["method"] = method

    sig  = df[p_col] < p_thresh
    up   = sig & (df[lfc_col] >  fc_thresh)
    down = sig & (df[lfc_col] < -fc_thresh)
    df["direction"] = np.where(up, "up", np.where(down, "down", "ns"))

    # -- Intermediate TSV for R/ggplot -----------------------------------------
    if out_dir is not None:
        safe_method = method.lower().replace(" ", "_").replace("(", "").replace(")", "")
        tsv_path = os.path.join(out_dir, f"volcano_cis_window_{safe_method}.tsv")
        df.rename(columns={p_col: "p_value", lfc_col: "log2fc"}).to_csv(
            tsv_path, sep="\t", index=False
        )
        print(f"  Saved R-ready cis window volcano data: {tsv_path}")

    fig, ax = plt.subplots(figsize=(8, 7))

    bg = df[~df["is_self_target"]].copy()
    bg_up   = bg[(bg[p_col] < p_thresh) & (bg[lfc_col] >  fc_thresh)]
    bg_down = bg[(bg[p_col] < p_thresh) & (bg[lfc_col] < -fc_thresh)]
    bg_ns   = bg[~((bg[p_col] < p_thresh) & (bg[lfc_col].abs() > fc_thresh))]
    ax.scatter(bg_ns[lfc_col],   bg_ns["-log10p"],
               color="#2ca02c", s=8, alpha=0.25, linewidths=0, rasterized=True, zorder=1,
               label=f"other cis pairs  (up={len(bg_up):,}, down={len(bg_down):,}, ns={len(bg_ns):,})")
    ax.scatter(bg_up[lfc_col],   bg_up["-log10p"],
               color="#2ca02c", s=12, alpha=0.6, linewidths=0, rasterized=True, zorder=2)
    ax.scatter(bg_down[lfc_col], bg_down["-log10p"],
               color="#2ca02c", s=12, alpha=0.6, linewidths=0, rasterized=True, zorder=2)

    self_df = df[df["is_self_target"]]
    for lbl in LABEL_ORDER:
        sub = self_df[self_df["type"] == lbl]
        if len(sub) == 0:
            continue
        color  = LABEL_COLORS.get(lbl, "#aaaaaa")
        n_up   = ((sub[p_col] < p_thresh) & (sub[lfc_col] >  fc_thresh)).sum()
        n_down = ((sub[p_col] < p_thresh) & (sub[lfc_col] < -fc_thresh)).sum()
        n_ns   = len(sub) - n_up - n_down
        ax.scatter(sub[lfc_col], sub["-log10p"],
                   color=color, s=22, alpha=0.85, linewidths=0, rasterized=True, zorder=10,
                   label=f"{lbl} (self-target)  (up={n_up}, down={n_down}, ns={n_ns})")

    ax.axhline(-np.log10(p_thresh), color="#666666",
               linewidth=0.9, linestyle="--", zorder=5, label=f"p = {p_thresh}")
    ax.axvline( fc_thresh, color="#666666", linewidth=0.9, linestyle="--", zorder=5)
    ax.axvline(-fc_thresh, color="#666666", linewidth=0.9, linestyle="--", zorder=5)
    ax.axvline(0,          color="#444444", linewidth=0.6, zorder=5)

    n_self = df["is_self_target"].sum()
    n_all  = len(df)
    ax.set_xlabel("log2 Fold Change", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    ax.set_title(
        f"{method} cis window volcano  (per element, all pairs)\n"
        f"grey = other cis pairs  |  colored = self-target "
        f"(gene == element)  |  n_total={n_all:,}, n_self={n_self:,}\n"
        f"nominal p < {p_thresh} (dashed line), |log2FC| > {fc_thresh}",
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
    _build_paths(args)
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Output directory        : {args.output_dir}")
    print(f"Results dir             : {args.results_dir}")
    print(f"Pipeline outputs dir    : {args.pipeline_outputs_dir}")
    print(f"Inference dir           : {args.inference_dir}")
    print(f"MuData filename         : {args.mudata_name}")
    print(f"P-value threshold       : {args.p_thresh}")
    print(f"Run UMAP                : {args.run_umap}")
    print(f"Mode                    : {'SCEPTRE-only files' if args.use_sceptre_only else 'Merged files (both methods)'}")

    guide_efficiency(args, args.output_dir)
    top_downregulated_tfs(args, args.output_dir)

    if args.run_umap:
        umap_hepatocyte_markers(args, args.output_dir)
    else:
        print("\n=== Analysis 3: UMAP skipped (pass --run_umap to enable) ===")

    # --- Analysis 4: Label-colored volcano plots (cis) ------------------------
    print("\n=== Analysis 4: Label-colored volcano plots (cis) ===")
    guide_meta = load_guide_metadata()
    if args.use_sceptre_only:
        method_file_specs = [
            ("SCEPTRE", SCEPTRE_GUIDE, SCEPTRE_ELEM, "p_value", "log2_fc"),
        ]
    else:
        method_file_specs = [
            ("SCEPTRE",  MERGED_CIS_GUIDE, MERGED_CIS_ELEM, "sceptre_p_value",  "sceptre_log2_fc"),
            ("Perturbo", MERGED_CIS_GUIDE, MERGED_CIS_ELEM, "perturbo_p_value", "perturbo_log2_fc"),
        ]
    for method, guide_path, elem_path, p_col, lfc_col in method_file_specs:
        for level, path in [("guide", guide_path), ("element", elem_path)]:
            if not os.path.exists(path):
                print(f"  Skipping {method} {level} volcano: file not found at {path}")
                continue
            try:
                df_raw = load_tsv(path, f"{method} {level} (volcano)")
                out_path = os.path.join(
                    args.output_dir,
                    f"volcano_cis_{level}_{method.lower()}.png",
                )
                plot_label_volcano(
                    df_results=df_raw,
                    p_col=p_col, lfc_col=lfc_col,
                    level=level, method=method,
                    out_path=out_path,
                    guide_meta=guide_meta,
                    title_suffix="cis",
                    p_thresh=0.01,
                    fc_thresh=0.5,
                    out_dir=args.output_dir,
                )
            except Exception as e:
                print(f"  WARNING: could not plot {method} {level} volcano: {e}")

            if level == "element":
                try:
                    window_out = os.path.join(
                        args.output_dir,
                        f"volcano_cis_window_{method.lower()}.png",
                    )
                    plot_cis_window_volcano(
                        df_raw=df_raw,
                        p_col=p_col, lfc_col=lfc_col,
                        method=method,
                        out_path=window_out,
                        guide_meta=guide_meta,
                        p_thresh=0.01,
                        fc_thresh=0.5,
                        out_dir=args.output_dir,
                    )
                except Exception as e:
                    print(f"  WARNING: could not plot {method} cis window volcano: {e}")

    # --- Volcano plots from calibrated files (element level only) ---------------
    # calib_volcano_specs: (method, path, p_col, lfc_col)
    # For the label volcano we use whichever file is available (direct or cis).
    # For the window volcano we always use CALIB_CIS because it contains ALL
    # cis-window pairs (is_direct_target=True and False), giving the green
    # background of non-self pairs. CALIB_DIRECT only has self-target rows.
    calib_volcano_specs = []
    if not args.use_sceptre_only:
        if CALIB_DIRECT and os.path.exists(CALIB_DIRECT):
            calib_volcano_specs.append(
                ("Perturbo (direct target)", CALIB_DIRECT, "empirical_pval", "log2fc")
            )
        if CALIB_CIS and os.path.exists(CALIB_CIS):
            calib_volcano_specs.append(
                ("Perturbo (cis 100kb)", CALIB_CIS, "empirical_pval", "log2fc")
            )

    def _prep_calib_df(path):
        """Load and normalise a calibrated TSV for the volcano functions."""
        df = pd.read_csv(path, sep="\t")
        df = df.copy()
        df["intended_target_name"] = df["element_id"]
        df["gene_id"]              = df["tested_gene_id"]
        if "element_label" in df.columns and "type" not in df.columns:
            df["type"] = df["element_label"]
        # is_self_target from the pre-computed is_direct_target flag; falls back
        # to True (for direct-target-only files where every row is self-target).
        if "is_direct_target" in df.columns:
            df["is_self_target"] = df["is_direct_target"].astype(bool)
        else:
            df["is_self_target"] = True
        return df

    for method, path, p_col, lfc_col in calib_volcano_specs:
        print(f"  Loading {method} calibrated file for volcano ...")
        df_raw = _prep_calib_df(path)
        print(f"    {len(df_raw):,} rows  "
              f"(self-target: {df_raw['is_self_target'].sum():,}, "
              f"other cis: {(~df_raw['is_self_target']).sum():,})")

        safe = method.lower().replace(" ", "_").replace("(", "").replace(")", "")

        # Label volcano: for the direct-target file every row is a self-target,
        # which is the correct view for a cis knockdown QC plot.
        out_path = os.path.join(args.output_dir, f"volcano_cis_element_{safe}.png")
        try:
            plot_label_volcano(
                df_results=df_raw,
                p_col=p_col, lfc_col=lfc_col,
                level="element", method=method,
                out_path=out_path,
                guide_meta=guide_meta,
                title_suffix="calibrated",   # skip internal cis filter
                p_thresh=0.01,
                fc_thresh=0.5,
                out_dir=args.output_dir,
            )
        except Exception as e:
            print(f"  WARNING: could not plot {method} calibrated volcano: {e}")

        # Window volcano: use CALIB_CIS so non-self cis pairs appear as green
        # background. For the direct-target entry, swap in the cis file.
        if method == "Perturbo (direct target)" and CALIB_CIS and os.path.exists(CALIB_CIS):
            print(f"  Using cis 100kb file for window volcano background ...")
            df_window = _prep_calib_df(CALIB_CIS)
        else:
            df_window = df_raw

        window_out = os.path.join(args.output_dir, f"volcano_cis_window_{safe}.png")
        try:
            plot_cis_window_volcano(
                df_raw=df_window,
                p_col=p_col, lfc_col=lfc_col,
                method=method,
                out_path=window_out,
                guide_meta=guide_meta,
                p_thresh=0.01,
                fc_thresh=0.5,
                out_dir=args.output_dir,
            )
        except Exception as e:
            print(f"  WARNING: could not plot {method} cis window volcano: {e}")

    print("\n=== All analyses complete ===")
    print(f"Results saved to: {args.output_dir}/")


if __name__ == "__main__":
    main()