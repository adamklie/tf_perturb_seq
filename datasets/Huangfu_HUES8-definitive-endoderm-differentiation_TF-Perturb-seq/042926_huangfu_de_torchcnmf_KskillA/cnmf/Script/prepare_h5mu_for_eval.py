"""Pre-process h5mu files before Stage 2 (Evaluation/Calibration).

Two modifications to each cNMF_<K>_2_0.h5mu under Inference/adata/:

1. Add a single-value categorical column obs['sample'] = 'all' on both
   the rna and cNMF modalities. Used as --categorical_key for the
   evaluation/calibration steps so that perturbation association is
   pooled across all 8 batches (not stratified per-batch).

2. Remap uns['guide_targets'] for the rna and cNMF modalities: any guide
   whose name starts with 'non-targeting' has its target relabeled
   from the literal 'nan' to 'non-targeting'. This matches the
   reference_targets convention expected by compute_perturbation_association
   and the U-test calibration's --guide_annotation_key.

Idempotent: re-running is a no-op once the column / labels are in place.
"""

import argparse
import os

import mudata as mu
import numpy as np
import pandas as pd


def update_h5mu(path: str, sample_value: str, control_token: str) -> None:
    print(f"Loading {path}...")
    mdata = mu.read(path)

    n_cells = mdata["rna"].n_obs

    # 1. Add sample column to both modalities (categorical, single value).
    for mod in ("rna", "cNMF"):
        if mod not in mdata.mod:
            continue
        ad = mdata.mod[mod]
        if "sample" in ad.obs.columns and (ad.obs["sample"].astype(str) == sample_value).all():
            print(f"  [{mod}] obs['sample'] already set; skipping.")
        else:
            ad.obs["sample"] = pd.Categorical([sample_value] * ad.n_obs, categories=[sample_value])
            print(f"  [{mod}] added obs['sample'] = '{sample_value}' for {ad.n_obs} cells.")

    # 2. Remap NT guide_targets from 'nan' to control_token.
    for mod in ("rna", "cNMF"):
        if mod not in mdata.mod:
            continue
        ad = mdata.mod[mod]
        if "guide_names" not in ad.uns or "guide_targets" not in ad.uns:
            continue
        names = np.asarray(ad.uns["guide_names"], dtype=str)
        targets = np.asarray(ad.uns["guide_targets"], dtype=str)
        nt_mask = np.array([n.lower().startswith("non-targeting") for n in names])
        n_remap = int(((targets != control_token) & nt_mask).sum())
        if n_remap == 0:
            print(f"  [{mod}] guide_targets already labeled '{control_token}' for NT guides ({nt_mask.sum()} NT guides).")
        else:
            targets = targets.copy()
            targets[nt_mask] = control_token
            ad.uns["guide_targets"] = targets.tolist()
            print(f"  [{mod}] remapped {n_remap} NT guide_targets -> '{control_token}'.")

    print(f"  Writing back to {path}...")
    mdata.write(path)
    del mdata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--adata_dir", required=True,
        help="Directory containing cNMF_<K>_<sel_thresh>.h5mu files",
    )
    parser.add_argument(
        "--K", nargs="+", type=int, required=True,
        help="K values to process",
    )
    parser.add_argument(
        "--sel_thresh_str", default="2_0",
        help="sel_thresh as filename token (e.g. '2_0' for sel_thresh=2.0)",
    )
    parser.add_argument(
        "--sample_value", default="all",
        help="Value for the single-category obs['sample'] column",
    )
    parser.add_argument(
        "--control_token", default="non-targeting",
        help="Target name to use for non-targeting guides",
    )
    args = parser.parse_args()

    for k in args.K:
        path = os.path.join(args.adata_dir, f"cNMF_{k}_{args.sel_thresh_str}.h5mu")
        if not os.path.exists(path):
            print(f"WARNING: {path} not found; skipping.")
            continue
        update_h5mu(path, args.sample_value, args.control_token)

    print("Done.")


if __name__ == "__main__":
    main()
