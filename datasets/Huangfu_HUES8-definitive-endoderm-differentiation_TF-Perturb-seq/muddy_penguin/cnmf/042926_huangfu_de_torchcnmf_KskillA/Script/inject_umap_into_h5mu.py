"""Compute a gene-based UMAP and inject into the cNMF h5mu.

Standard scanpy recipe applied to a temporary copy of mdata[rna]:
  normalize_total -> log1p -> highly_variable_genes -> scale -> PCA -> neighbors -> UMAP

Result is written into BOTH:
  - mdata[rna].obsm['X_pca'], ['X_umap']   (canonical home)
  - mdata[cNMF].obsm['X_pca'], ['X_umap']  (where Stage 3b/3c look)

Going forward we should compute UMAP on the input h5ad before Stage 1
inference so this post-hoc inject becomes unnecessary; this script is the
remediation for h5mu files produced before that convention.

Idempotent: if both modalities already have X_umap, exits without recomputing.
Pass --force to recompute anyway (e.g., if you want to overwrite a prior bad inject).
"""

import argparse
import os

import mudata as mu
import scanpy as sc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5mu_path", required=True)
    parser.add_argument("--rna_key", default="rna")
    parser.add_argument("--prog_key", default="cNMF")
    parser.add_argument("--n_top_hvg", type=int, default=2000)
    parser.add_argument("--n_comps_pca", type=int, default=50)
    parser.add_argument("--n_neighbors", type=int, default=30)
    parser.add_argument("--target_sum", type=float, default=1e4)
    parser.add_argument("--force", action="store_true",
                        help="Recompute even if X_umap already present.")
    args = parser.parse_args()

    print(f"Reading {args.h5mu_path}...")
    mdata = mu.read(args.h5mu_path)

    rna = mdata[args.rna_key]
    prog = mdata[args.prog_key]

    if (
        not args.force
        and "X_umap" in prog.obsm
        and "X_umap" in rna.obsm
        and "X_pca" in prog.obsm
        and "X_pca" in rna.obsm
    ):
        print("  X_umap and X_pca already present in both modalities; skipping.")
        return

    print(f"  Building tmp anndata from rna ({rna.n_obs} cells x {rna.n_vars} genes)...")
    adata = rna.copy()

    print(f"  normalize_total(target_sum={args.target_sum})...")
    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    print("  log1p...")
    sc.pp.log1p(adata)

    print(f"  highly_variable_genes(n_top={args.n_top_hvg}, subset=True)...")
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_hvg, subset=True)
    print(f"  -> retained {adata.n_vars} HVGs")

    print("  scale(max_value=10)...")
    sc.pp.scale(adata, max_value=10)

    print(f"  PCA (n_comps={args.n_comps_pca})...")
    sc.tl.pca(adata, n_comps=args.n_comps_pca)
    print(f"  neighbors (n_neighbors={args.n_neighbors})...")
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors)
    print("  UMAP...")
    sc.tl.umap(adata)

    # Write into both modalities
    print("  Injecting X_pca + X_umap into mdata[rna].obsm and mdata[cNMF].obsm...")
    rna.obsm["X_pca"] = adata.obsm["X_pca"]
    rna.obsm["X_umap"] = adata.obsm["X_umap"]
    prog.obsm["X_pca"] = adata.obsm["X_pca"]
    prog.obsm["X_umap"] = adata.obsm["X_umap"]

    print(f"  Writing back to {args.h5mu_path}...")
    mdata.write(args.h5mu_path)
    print("Done.")


if __name__ == "__main__":
    main()
