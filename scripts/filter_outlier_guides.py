#!/usr/bin/env python3
"""
filter_outlier_guides.py

Filter outlier guides — and all cells that received any of those guides — from
a MuData object, based on p-value outlier tables with BH-FDR correction.

Steps
-----
1. Load outlier tables, apply Benjamini-Hochberg FDR correction independently
   to each table.
2. Collect the union of outlier guide IDs at the given FDR threshold.
3. Identify every cell that received any outlier guide (non-zero entry in the
   ``guide_assignment`` layer for any outlier guide column).
4. Remove those cells from ALL modalities.
5. Remove the outlier guide vars from the guide modality.
6. Recompute guide-count obs fields (``num_expressed_guides``,
   ``total_guide_umis``) for the surviving cells.
7. Write the filtered object to disk.

Usage
-----
    srun python filter_outlier_guides.py \\
        --input  /hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/post_pipeline_processing/inference_mudata.h5mu \\
        --output /hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/post_pipeline_processing/inference_mudata_outlier_guide_filt.h5mu \\
        --non_targeting_outliers /hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/energy_distance/data_hep/non_targeting_outlier_table.csv \\
        --targeting_outliers     /hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/energy_distance/data_hep/targeting_outlier_table.csv \\
        [--fdr_threshold 0.05]
"""

import argparse

import numpy as np
import pandas as pd
import scipy.sparse as sp
import mudata as md
from statsmodels.stats.multitest import multipletests


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_outlier_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    if "pval_outlier" not in df.columns:
        raise ValueError(
            f"Expected a 'pval_outlier' column in {path}. "
            f"Found columns: {list(df.columns)}"
        )
    return df


def apply_fdr(df: pd.DataFrame, alpha: float) -> pd.Series:
    """Benjamini-Hochberg FDR correction. Returns bool Series (True = outlier)."""
    pvals = df["pval_outlier"].values.astype(float)
    reject, _, _, _ = multipletests(pvals, alpha=alpha, method="fdr_bh")
    return pd.Series(reject, index=df.index, name="fdr_reject")


def find_outlier_guides(
    non_targeting_path: str,
    targeting_path: str,
    fdr_threshold: float,
) -> set:
    nt_df = load_outlier_table(non_targeting_path)
    t_df  = load_outlier_table(targeting_path)

    nt_reject = apply_fdr(nt_df, alpha=fdr_threshold)
    t_reject  = apply_fdr(t_df,  alpha=fdr_threshold)

    outlier_guides = set(nt_df.index[nt_reject]) | set(t_df.index[t_reject])

    print(
        f"[FDR={fdr_threshold}]  "
        f"Non-targeting outliers: {nt_reject.sum()}/{len(nt_df)}  |  "
        f"Targeting outliers: {t_reject.sum()}/{len(t_df)}  |  "
        f"Total unique outlier guides: {len(outlier_guides)}"
    )
    return outlier_guides


def cells_with_any_outlier_guide(guide_assignment, outlier_col_idx: np.ndarray) -> np.ndarray:
    """
    Return a boolean array (length = n_cells) that is True for every cell
    with a non-zero assignment in any outlier guide column.

    Works for both sparse and dense matrices without materialising the full
    matrix in a new dense array.
    """
    if sp.issparse(guide_assignment):
        # Slice only the outlier columns — cheap for CSC; convert if needed
        mat = guide_assignment.tocsc()[:, outlier_col_idx]
        # .sum(axis=1) on a sparse matrix returns a dense matrix; flatten to 1-D
        return np.asarray(mat.sum(axis=1)).flatten() > 0
    else:
        return guide_assignment[:, outlier_col_idx].sum(axis=1) > 0


def recompute_guide_obs(guide_adata) -> None:
    """
    Recompute num_expressed_guides and total_guide_umis in-place on a
    (already-filtered) guide AnnData from its layers / X matrix.

    Only updates fields that are actually present in obs.
    """
    obs_fields = set(guide_adata.obs.columns)

    if "num_expressed_guides" in obs_fields:
        assignment = guide_adata.layers["guide_assignment"]
        if sp.issparse(assignment):
            counts = np.asarray((assignment > 0).sum(axis=1)).flatten()
        else:
            counts = (assignment > 0).sum(axis=1)
        guide_adata.obs["num_expressed_guides"] = counts

    if "total_guide_umis" in obs_fields:
        X = guide_adata.X
        if sp.issparse(X):
            totals = np.asarray(X.sum(axis=1)).flatten()
        else:
            totals = X.sum(axis=1)
        guide_adata.obs["total_guide_umis"] = totals


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def filter_mudata(
    input_path: str,
    output_path: str,
    non_targeting_path: str,
    targeting_path: str,
    fdr_threshold: float,
) -> None:

    # ---- Load ---------------------------------------------------------------
    print(f"Loading MuData from: {input_path}")
    mdata = md.read_h5mu(input_path)
    print(mdata)
    print()

    if "guide" not in mdata.mod:
        raise KeyError("MuData object does not contain a 'guide' modality.")

    guide_adata = mdata["guide"]
    n_obs_original  = mdata.n_obs
    n_vars_original = guide_adata.n_vars

    # ---- Identify outlier guides --------------------------------------------
    outlier_ids = find_outlier_guides(non_targeting_path, targeting_path, fdr_threshold)

    guide_var_names   = guide_adata.var_names
    outlier_ids_found = outlier_ids & set(guide_var_names)
    missing           = outlier_ids - outlier_ids_found

    if missing:
        print(
            f"WARNING: {len(missing)} outlier guide ID(s) not present in "
            f"var_names and will be skipped."
        )

    # Boolean / index arrays for guide var axis
    outlier_var_mask = guide_var_names.isin(outlier_ids_found)   # True = remove
    outlier_col_idx  = np.where(outlier_var_mask)[0]
    keep_var_mask    = ~outlier_var_mask
    guides_to_keep   = guide_var_names[keep_var_mask]

    n_guides_removed = int(outlier_var_mask.sum())
    print(f"Outlier guides identified in dataset : {n_guides_removed}")

    # ---- Identify cells that received any outlier guide ---------------------
    print("Scanning guide_assignment layer for cells receiving outlier guides...")
    cells_with_outlier = cells_with_any_outlier_guide(
        guide_adata.layers["guide_assignment"], outlier_col_idx
    )
    cells_to_keep = ~cells_with_outlier

    n_cells_removed  = int(cells_with_outlier.sum())
    n_cells_retained = int(cells_to_keep.sum())
    print(f"Cells removed (received ≥1 outlier guide) : {n_cells_removed:,}  "
          f"({100 * n_cells_removed / n_obs_original:.2f}% of {n_obs_original:,})")
    print(f"Cells retained                            : {n_cells_retained:,}")
    print(f"Guides removed                            : {n_guides_removed:,}  "
          f"({100 * n_guides_removed / n_vars_original:.2f}% of {n_vars_original:,})")
    print(f"Guides retained                           : {len(guides_to_keep):,}")
    print()

    # ---- Filter each modality -----------------------------------------------
    # Process modalities one at a time. Extracting filtered copies before
    # building the new MuData avoids holding more than (original + one copy)
    # in memory at once. The gene copy is unavoidable — we are removing rows.

    filtered_mods = {}

    for mod_name, adata in mdata.mod.items():
        if mod_name == "guide":
            # Subset both cells (obs) and vars
            filtered = adata[cells_to_keep, :][:, guides_to_keep].copy()
            recompute_guide_obs(filtered)
        else:
            # Subset cells only
            filtered = adata[cells_to_keep, :].copy()
        filtered_mods[mod_name] = filtered
        print(f"Filtered modality '{mod_name}': {filtered.n_obs:,} × {filtered.n_vars:,}")

    # Carry over top-level obs/uns; free the original as soon as possible
    top_obs = mdata.obs[cells_to_keep].copy()
    top_uns = mdata.uns          # dict of small summary objects — reference is fine
    del mdata                    # release original; filtered modalities are self-contained

    # ---- Reassemble MuData --------------------------------------------------
    new_mdata = md.MuData(filtered_mods)
    new_mdata.obs = top_obs
    new_mdata.uns = top_uns

    print(f"\nFiltered MuData:\n{new_mdata}")

    # ---- Write --------------------------------------------------------------
    print(f"\nWriting filtered MuData to: {output_path}")
    new_mdata.write_h5mu(output_path)
    print("Done.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Remove outlier guides and all cells that received them from a "
            "MuData object, using BH-FDR correction on p-value outlier tables."
        )
    )
    parser.add_argument("--input",  "-i", required=True,
                        help="Path to the input .h5mu MuData file.")
    parser.add_argument("--output", "-o", required=True,
                        help="Path for the filtered .h5mu MuData file.")
    parser.add_argument("--non_targeting_outliers", required=True,
                        help="Path to non_targeting_outlier_table.csv.")
    parser.add_argument("--targeting_outliers", required=True,
                        help="Path to targeting_outlier_table.csv.")
    parser.add_argument("--fdr_threshold", type=float, default=0.05,
                        help="BH-FDR threshold for calling a guide an outlier "
                             "(default: 0.05).")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    filter_mudata(
        input_path=args.input,
        output_path=args.output,
        non_targeting_path=args.non_targeting_outliers,
        targeting_path=args.targeting_outliers,
        fdr_threshold=args.fdr_threshold,
    )


if __name__ == "__main__":
    main()