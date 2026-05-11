"""Compute a gene-based UMAP and inject it into ONE OR MORE cNMF h5mu files.

Standard scanpy recipe applied to a temporary copy of mdata[rna]:
  normalize_total -> log1p -> highly_variable_genes -> scale -> PCA -> neighbors -> UMAP

The result is written into BOTH modalities of every target h5mu:
  - mdata[rna].obsm['X_pca'], ['X_umap']   (canonical home)
  - mdata[cNMF].obsm['X_pca'], ['X_umap']  (where Stage 3b/3c look)

Across different K h5mu files the rna matrix is identical (cNMF differs only
in the cNMF modality), so UMAP is computed ONCE from the first h5mu's rna and
the same embedding is broadcast into every additional target. Pass multiple
--h5mu_path arguments (or use --h5mu_glob) to inject the same UMAP across a
whole K sweep.

Going forward we should compute UMAP on the input h5ad / inference_mudata
BEFORE Stage 1 inference (e.g. in Convert_file_adata.py, right before
writing the cNMF AnnData input) so the embedding propagates through cNMF
naturally and this post-hoc inject becomes unnecessary. This script is the
remediation for h5mu files produced before that convention.

Because UMAP depends only on the rna matrix (identical across K), for a
given dataset you only need to inject the *selected K* h5mu — defer the
inject until Stage 3a K-selection has picked a k, then run this on the
single h5mu Stage 3 will use.

Idempotent: if a target h5mu already has X_umap+X_pca in both modalities,
it is skipped. Pass --force to overwrite.
"""

import argparse
import glob
import os
import sys

import mudata as mu
import numpy as np
import scanpy as sc


def compute_embedding(rna, *, n_top_hvg, n_comps_pca, n_neighbors, target_sum):
    """Run the canonical scanpy recipe and return (X_pca, X_umap)."""
    print(f"  rna: {rna.n_obs} cells x {rna.n_vars} genes")
    adata = rna.copy()
    print(f"  normalize_total(target_sum={target_sum})...")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    print("  log1p...")
    sc.pp.log1p(adata)
    print(f"  highly_variable_genes(n_top={n_top_hvg}, subset=True)...")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_hvg, subset=True)
    print(f"  -> retained {adata.n_vars} HVGs")
    print("  scale(max_value=10)...")
    sc.pp.scale(adata, max_value=10)
    print(f"  PCA (n_comps={n_comps_pca})...")
    sc.tl.pca(adata, n_comps=n_comps_pca)
    print(f"  neighbors (n_neighbors={n_neighbors})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    print("  UMAP...")
    sc.tl.umap(adata)
    return np.asarray(adata.obsm["X_pca"]), np.asarray(adata.obsm["X_umap"])


def inject_embedding(h5mu_path, X_pca, X_umap, *, rna_key, prog_key, force):
    """Write the precomputed X_pca + X_umap into both modalities of an h5mu."""
    print(f"\n[{h5mu_path}]")
    mdata = mu.read(h5mu_path)
    rna = mdata[rna_key]
    prog = mdata[prog_key]
    if (
        not force
        and "X_umap" in prog.obsm
        and "X_umap" in rna.obsm
        and "X_pca" in prog.obsm
        and "X_pca" in rna.obsm
    ):
        print("  X_umap+X_pca already present in both modalities; skipping.")
        return False
    if rna.n_obs != X_pca.shape[0]:
        sys.exit(
            f"Cell count mismatch: {h5mu_path} rna has {rna.n_obs} cells, "
            f"embedding has {X_pca.shape[0]} rows. UMAP was computed on a "
            "different h5mu — the rna matrix is supposed to be identical across "
            "K h5mu files but this one disagrees."
        )
    rna.obsm["X_pca"] = X_pca
    rna.obsm["X_umap"] = X_umap
    prog.obsm["X_pca"] = X_pca
    prog.obsm["X_umap"] = X_umap
    print(f"  writing back...")
    mdata.write(h5mu_path)
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--h5mu_path",
        nargs="+",
        default=None,
        help="One or more h5mu files to inject. UMAP is computed from the "
        "first listed file's rna and broadcast to the rest.",
    )
    parser.add_argument(
        "--h5mu_glob",
        default=None,
        help="Glob pattern (e.g. '/.../cNMF_*_2_0.h5mu') matched in addition "
        "to --h5mu_path.",
    )
    parser.add_argument("--rna_key", default="rna")
    parser.add_argument("--prog_key", default="cNMF")
    parser.add_argument("--n_top_hvg", type=int, default=2000)
    parser.add_argument("--n_comps_pca", type=int, default=50)
    parser.add_argument("--n_neighbors", type=int, default=30)
    parser.add_argument("--target_sum", type=float, default=1e4)
    parser.add_argument("--force", action="store_true",
                        help="Recompute / overwrite even if X_umap already present.")
    args = parser.parse_args()

    paths = list(args.h5mu_path or [])
    if args.h5mu_glob:
        paths += sorted(glob.glob(args.h5mu_glob))
    # dedup, preserve order
    seen = set()
    paths = [p for p in paths if not (p in seen or seen.add(p))]
    if not paths:
        sys.exit("No h5mu paths given (use --h5mu_path and/or --h5mu_glob).")
    for p in paths:
        if not os.path.exists(p):
            sys.exit(f"Missing h5mu: {p}")

    print(f"Will inject UMAP into {len(paths)} h5mu file(s):")
    for p in paths:
        print(f"  - {p}")

    print(f"\nComputing UMAP from {paths[0]} ...")
    mdata0 = mu.read(paths[0])
    X_pca, X_umap = compute_embedding(
        mdata0[args.rna_key],
        n_top_hvg=args.n_top_hvg,
        n_comps_pca=args.n_comps_pca,
        n_neighbors=args.n_neighbors,
        target_sum=args.target_sum,
    )
    print(f"  computed embedding: X_pca={X_pca.shape}, X_umap={X_umap.shape}")

    # Inject into first h5mu using the in-memory mdata0 (avoid re-read)
    print(f"\n[{paths[0]}] (in-memory)")
    rna0 = mdata0[args.rna_key]
    prog0 = mdata0[args.prog_key]
    already = (
        "X_umap" in prog0.obsm and "X_umap" in rna0.obsm
        and "X_pca" in prog0.obsm and "X_pca" in rna0.obsm
    )
    if already and not args.force:
        print("  X_umap+X_pca already present in both modalities; skipping write.")
    else:
        rna0.obsm["X_pca"] = X_pca
        rna0.obsm["X_umap"] = X_umap
        prog0.obsm["X_pca"] = X_pca
        prog0.obsm["X_umap"] = X_umap
        print(f"  writing back...")
        mdata0.write(paths[0])
    del mdata0

    n_written = 1 if (not already or args.force) else 0
    n_skipped = 1 - n_written
    for p in paths[1:]:
        wrote = inject_embedding(
            p, X_pca, X_umap,
            rna_key=args.rna_key, prog_key=args.prog_key, force=args.force,
        )
        if wrote:
            n_written += 1
        else:
            n_skipped += 1

    print(f"\nDone. wrote={n_written}, skipped={n_skipped}, total={len(paths)}")


if __name__ == "__main__":
    main()
