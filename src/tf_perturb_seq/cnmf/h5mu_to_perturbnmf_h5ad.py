#!/usr/bin/env python3
"""
Convert a CRISPR pipeline `inference_mudata.h5mu` into a PerturbNMF-compatible
single-modality `.h5ad` file.

The CRISPR pipeline output has separate `gene` and `guide` modalities. PerturbNMF
expects ONE AnnData with the gene matrix as `X` and the guide assignment / names /
targets packed into `obsm` and `uns`. See:
  external/PerturbNMF/.claude/skills/perturbNMF-runner/references/data-format-spec.md

Mapping:

  Source (inference_mudata.h5mu)              ->  Destination (output .h5ad)
  -------------------------------------------     ----------------------------
  mod['gene'].X                                  X (cells x genes counts)
  mod['gene'].var['symbol']  (column)            var.index (gene SYMBOLS, required by spec)
  mod['gene'].var (whole)                        var (preserves mt/ribo flags etc.)
  mod['gene'].obs                                obs  (must contain a categorical key, default 'batch')
  mod['guide'].layers['guide_assignment']        obsm['guide_assignment']  (sparse cells x guides)
  mod['guide'].var[<guide_id_col>]               uns['guide_names']
  mod['guide'].var[<target_col>]                 uns['guide_targets']

Validates with the perturbNMF-runner skill's `validate_data.py` after writing.

Usage:
    python src/tf_perturb_seq/cnmf/h5mu_to_perturbnmf_h5ad.py \
        --in_h5mu  datasets/<dataset>/runs/<run>/pipeline_outputs/inference_mudata.h5mu \
        --out_h5ad datasets/<dataset>/PerturbNMF/Data/<dataset>_<run>_perturbnmf.h5ad
"""

import argparse
import os
import warnings

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.filterwarnings("ignore")


def filter_cells_zero_hvg_counts(adata, num_highvar_genes):
    """Match cnmf.cnmf.get_norm_counts: drop cells with zero counts in the
    overdispersed-gene set that cnmf would select. This is the exact failure
    mode cnmf raises ('%d cells have zero counts of overdispersed genes').

    Replicates cnmf's HVG selection by:
      1. TPM-normalizing the count matrix (each cell -> 1e6)
      2. Calling cnmf's get_highvar_genes_sparse with numgenes=num_highvar_genes
      3. Identifying cells whose RAW counts in HVG columns sum to 0
      4. Returning a boolean mask of cells to KEEP.
    """
    from cnmf.cnmf import get_highvar_genes_sparse

    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    total = np.array(X.sum(axis=1)).flatten()
    # avoid div-by-zero (cells with literally zero counts will be dropped anyway)
    total_safe = np.where(total > 0, total, 1.0)
    tpm = X.multiply(1e6 / total_safe[:, None]).tocsr()
    gene_stats, _ = get_highvar_genes_sparse(tpm, numgenes=num_highvar_genes)
    hvg_mask = gene_stats["high_var"].values
    hvg_counts_per_cell = np.array(X[:, hvg_mask].sum(axis=1)).flatten()
    keep = hvg_counts_per_cell > 0
    n_drop = (~keep).sum()
    print(f"  HVG filter (n_top={num_highvar_genes}): dropping {n_drop} cells with zero counts in the {hvg_mask.sum()} overdispersed genes")
    return keep


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in_h5mu", required=True, help="Path to inference_mudata.h5mu")
    p.add_argument("--out_h5ad", required=True, help="Output .h5ad path")
    p.add_argument("--guide_id_col", default="guide_id",
                   help="Column in guide.var holding guide IDs (default: guide_id)")
    p.add_argument("--target_col", default="gene_name",
                   help="Column in guide.var holding target gene SYMBOLS (default: gene_name)")
    p.add_argument("--symbol_col", default="symbol",
                   help="Column in gene.var holding gene SYMBOLS (default: symbol)")
    p.add_argument("--num_highvar_genes", type=int, default=2000,
                   help="Number of HVGs cnmf will select (must match torch-cNMF --numhvgenes; "
                        "default 2000). Cells with zero counts in the resulting HVG set are dropped, "
                        "matching cnmf's get_norm_counts() pre-flight check.")
    p.add_argument("--no_hvg_filter", action="store_true",
                   help="Skip the HVG-zero-counts cell filter (not recommended — cnmf will fail).")
    args = p.parse_args()

    print(f"Loading {args.in_h5mu}")
    m = md.read_h5mu(args.in_h5mu)
    g, gd = m.mod["gene"], m.mod["guide"]
    print(f"  gene:  {g.shape}")
    print(f"  guide: {gd.shape}, layers: {list(gd.layers.keys())}")

    # var_names = gene SYMBOLS (spec requirement). Drop duplicates by keeping the first row.
    print("Mapping var_names from Ensembl IDs to gene symbols...")
    new_var = g.var.copy()
    new_var.index = new_var[args.symbol_col].astype(str).values
    new_var.index.name = args.symbol_col
    keep = ~new_var.index.duplicated()
    n_dup = (~keep).sum()
    if n_dup:
        print(f"  dropping {n_dup} rows with duplicate symbols (keeping first)")
    X = g.X[:, keep] if hasattr(g.X, "shape") else g.X[:, keep]
    new_var = new_var[keep]

    a = ad.AnnData(X=X, obs=g.obs.copy(), var=new_var)
    print(f"  result: {a.shape}")

    # Pack guide info
    print("Packing guide modality into obsm/uns...")
    ga = gd.layers["guide_assignment"]
    if not sp.issparse(ga):
        ga = sp.csr_matrix(ga)
    print(f"  guide_assignment: shape={ga.shape}, nnz={ga.nnz}, sparsity={1 - ga.nnz / np.prod(ga.shape):.2%}")
    a.obsm["guide_assignment"] = ga
    a.uns["guide_names"] = np.array(gd.var[args.guide_id_col].astype(str).tolist())
    a.uns["guide_targets"] = np.array(gd.var[args.target_col].astype(str).tolist())

    n_unique_targets = pd.Series(a.uns["guide_targets"]).nunique()
    print(f"  n_guides: {len(a.uns['guide_names'])}; unique targets: {n_unique_targets}")

    if "batch" not in a.obs.columns:
        print("  WARNING: no 'batch' column in obs — PerturbNMF expects a categorical key (default 'batch'/'sample')")
    else:
        print(f"  batch: {a.obs['batch'].nunique()} unique values")

    # Pre-flight HVG-zero-counts filter — matches cnmf.cnmf.get_norm_counts behavior.
    # Without this, cnmf raises mid-run after the user's already burned compute on prepare/factorize.
    if not args.no_hvg_filter:
        print(f"\nRunning HVG pre-flight filter (num_highvar_genes={args.num_highvar_genes})...")
        keep = filter_cells_zero_hvg_counts(a, args.num_highvar_genes)
        a = a[keep].copy()  # anndata auto-subsets obsm/uns layers along obs axis
        print(f"  result after HVG filter: {a.shape}")

    print(f"\nWriting {args.out_h5ad}")
    os.makedirs(os.path.dirname(args.out_h5ad), exist_ok=True)
    a.write_h5ad(args.out_h5ad, compression="gzip")
    print(f"Done. {os.path.getsize(args.out_h5ad) / 1e9:.2f} GB on disk.")


if __name__ == "__main__":
    main()
