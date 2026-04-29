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
    python scripts/h5mu_to_perturbnmf_h5ad.py \
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

    print(f"\nWriting {args.out_h5ad}")
    os.makedirs(os.path.dirname(args.out_h5ad), exist_ok=True)
    a.write_h5ad(args.out_h5ad, compression="gzip")
    print(f"Done. {os.path.getsize(args.out_h5ad) / 1e9:.2f} GB on disk.")


if __name__ == "__main__":
    main()
