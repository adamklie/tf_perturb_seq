"""Subsample cells in a cNMF h5mu so Stage 3b plotting fits in memory.

Stage 3b (cNMF_program_analysis.py) OOMs on production-scale datasets
(~270k cells × K=200 × 2k targets) at any tested memory level. The plots
themselves don't need 269k cells for UMAP density; ~27k cells (10%) is
plenty visually. Subsample once here so the upstream Stage 3b script
can run unchanged.

Preserves:
  - obs columns (sample, batch, n_counts, etc.) + UMAP/PCA obsm
  - guide_assignment obsm (subsetted to the kept cells)
  - var (genes/programs) untouched
  - uns intact (guide_names, guide_targets)

Output is `cNMF_<K>_<sel>_thinned.h5mu` next to the input.

Usage:
    python thin_h5mu_for_3b.py --in-h5mu .../cNMF_200_2_0.h5mu --frac 0.1 --seed 0
"""

import argparse
import os

import muon as mu
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--in-h5mu", required=True)
    p.add_argument("--frac", type=float, default=0.1)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out-suffix", default="_thinned")
    args = p.parse_args()

    print(f"Reading {args.in_h5mu}...")
    mdata = mu.read(args.in_h5mu)
    n_cells = mdata["rna"].n_obs
    print(f"  full: {n_cells} cells x {mdata['rna'].n_vars} genes (rna), {mdata['cNMF'].n_vars} programs (cNMF)")

    rng = np.random.default_rng(args.seed)
    n_keep = int(round(n_cells * args.frac))
    keep_idx = np.sort(rng.choice(n_cells, size=n_keep, replace=False))
    print(f"  keeping {n_keep} cells ({100 * args.frac:.0f}%)")

    # Use mudata's slicing: mdata[bool_mask, :] returns a new MuData
    # Need positional index → bool mask aligned to mdata.obs
    mask = np.zeros(n_cells, dtype=bool)
    mask[keep_idx] = True

    thinned = mdata[mask].copy()
    print(f"  thinned: {thinned['rna'].n_obs} cells")

    base, ext = os.path.splitext(args.in_h5mu)
    out_path = base + args.out_suffix + ext
    print(f"  writing {out_path}...")
    thinned.write(out_path)
    print(f"  size: {os.path.getsize(out_path):,} B ({os.path.getsize(out_path) / 1024 / 1024:.0f} MB)")
    print("Done.")


if __name__ == "__main__":
    main()
