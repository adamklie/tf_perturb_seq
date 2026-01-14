#!/usr/bin/env python3
"""
Convert a MuData (.h5mu) object into AnnData compatible format (.h5ad) for gene program discovery.
"""

import argparse

import anndata as ad
import mudata as md
import scanpy as sc
from scipy.sparse import issparse


def build_anndata_from_mudata(
    mudata_path: str,
    out_h5ad: str,
    n_pcs: int = 50,
    dense_guide_assignment: bool = True,
) -> None:
    
    # Load MuData
    mdata = md.read_h5mu(mudata_path)
    gene = mdata.mod["gene"]
    guide = mdata.mod["guide"]

    # Create a gene-expression AnnData
    adata = ad.AnnData(
        X=gene.X,
        obs=gene.obs.copy(),
        var=gene.var.copy(),
    )

    # Add guide assignment (optionally densify)
    GA = guide.layers["guide_assignment"]
    if dense_guide_assignment and issparse(GA):
        GA = GA.toarray()
    adata.obsm["guide_assignment"] = GA

    # Pull guide metadata
    guide_names = guide.var["guide_id"].astype(str).to_numpy()
    guide_targets = guide.var["intended_target_name"].astype(str).to_numpy()

    # Merge .uns (guide first, then mudata; mudata overwrites on key collision)
    for k in guide.uns.keys():
        adata.uns[k] = guide.uns[k]
    for k in mdata.uns.keys():
        adata.uns[k] = mdata.uns[k]

    # Store explicit arrays in .uns (as in your snippet)
    adata.uns["guide_names"] = guide_names
    adata.uns["guide_targets"] = guide_targets

    # Compute embeddings
    sc.tl.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Save
    adata.write_h5ad(out_h5ad)
    print(f"Wrote AnnData to: {out_h5ad}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Create an AnnData (.h5ad) from MuData (.h5mu) and add guide assignment metadata."
    )
    p.add_argument("--mudata", required=True, help="Path to input .h5mu file (MuData).")
    p.add_argument("--out_h5ad", required=True, help="Path to output .h5ad file.")
    p.add_argument("--n_pcs", type=int, default=50, help="Number of PCs for PCA (default: 50).")
    p.add_argument(
        "--dense_guide_assignment",
        action="store_true",
        help="If set, convert guide_assignment to a dense ndarray before storing in .obsm.",
    )
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build_anndata_from_mudata(
        mudata_path=args.mudata,
        out_h5ad=args.out_h5ad,
        n_pcs=args.n_pcs,
        dense_guide_assignment=args.dense_guide_assignment,
    )
