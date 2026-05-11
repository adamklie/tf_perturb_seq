"""Convert an inference MuData (gene + guide modalities) to AnnData for cNMF.

Mirrors the pattern in Hon's PerturbNMF/Script/Convert_file_adata.py — copies
the gene mod as the output adata, attaches `guide_assignment` from
`guide.layers`, and stages `guide_names` / `guide_targets` in `adata.uns`.

The .X of the output adata must stay as raw counts (cNMF input). UMAP is
computed on a temporary copy and only the embeddings are written back to
`adata.obsm`.

Optional: `--compute_umap` runs the canonical scanpy recipe
(normalize_total → log1p → HVG → scale → PCA → neighbors → UMAP) on a temp
copy of the gene mod and stashes `obsm['X_pca']` + `obsm['X_umap']` in the
output AnnData. Stage 1 carries them through into every per-K h5mu's rna
modality, which removes the need to run `inject_umap_into_h5mu.py` post-hoc
on Stage 1 outputs. Recommended for any new dataset run.
"""

import argparse

import muon as mu
import scanpy as sc

DEFAULT_INPUT = "/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/PerturbNMF/Data/inference_mudata.h5mu"
DEFAULT_OUTPUT = "/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/PerturbNMF/Data/inference_mudata_cleaned.h5ad"


def compute_umap_into_obsm(adata, *, n_top_hvg, n_comps_pca, n_neighbors, target_sum):
    """Run the canonical recipe on a temp copy; write X_pca + X_umap into adata.obsm in place.

    adata.X is NOT modified — only obsm is updated. cNMF still gets raw counts in .X.
    """
    tmp = adata.copy()
    print(f"  rna temp: {tmp.n_obs} cells x {tmp.n_vars} genes")
    print(f"  normalize_total(target_sum={target_sum}) -> log1p")
    sc.pp.normalize_total(tmp, target_sum=target_sum)
    sc.pp.log1p(tmp)
    print(f"  highly_variable_genes(n_top={n_top_hvg}, subset=True)")
    sc.pp.highly_variable_genes(tmp, n_top_genes=n_top_hvg, subset=True)
    print(f"  -> retained {tmp.n_vars} HVGs")
    print("  scale(max_value=10)")
    sc.pp.scale(tmp, max_value=10)
    print(f"  PCA (n_comps={n_comps_pca})")
    sc.tl.pca(tmp, n_comps=n_comps_pca)
    print(f"  neighbors (n_neighbors={n_neighbors})")
    sc.pp.neighbors(tmp, n_neighbors=n_neighbors)
    print("  UMAP")
    sc.tl.umap(tmp)
    adata.obsm["X_pca"] = tmp.obsm["X_pca"]
    adata.obsm["X_umap"] = tmp.obsm["X_umap"]
    print(f"  injected: X_pca={adata.obsm['X_pca'].shape}, X_umap={adata.obsm['X_umap'].shape}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default=DEFAULT_INPUT, help="Input inference MuData (.h5mu)")
    ap.add_argument("--output", default=DEFAULT_OUTPUT, help="Output AnnData (.h5ad)")
    ap.add_argument("--gene_mod", default="gene", help="MuData modality name for the gene matrix (default: gene)")
    ap.add_argument("--guide_mod", default="guide", help="MuData modality name for guides (default: guide)")
    ap.add_argument(
        "--compute_umap",
        action="store_true",
        help="Compute PCA + UMAP on the gene matrix and stash in obsm. Recommended "
        "going forward so Stage 1 carries embeddings into every output h5mu and "
        "post-hoc inject_umap_into_h5mu is not needed.",
    )
    ap.add_argument("--n_top_hvg", type=int, default=2000)
    ap.add_argument("--n_comps_pca", type=int, default=50)
    ap.add_argument("--n_neighbors", type=int, default=30)
    ap.add_argument("--target_sum", type=float, default=1e4)
    args = ap.parse_args()

    print(f"Reading {args.input}...")
    mudata = mu.read(args.input)
    adata = mudata[args.gene_mod].copy()

    adata.obsm["guide_assignment"] = mudata[args.guide_mod].layers["guide_assignment"].copy()
    adata.uns["guide_names"] = list(mudata[args.guide_mod].var["guide_id"])
    adata.uns["guide_targets"] = list(mudata[args.guide_mod].var["gene_name"])

    if args.compute_umap:
        print("Computing UMAP from gene matrix (cNMF will inherit it through Stage 1)...")
        compute_umap_into_obsm(
            adata,
            n_top_hvg=args.n_top_hvg,
            n_comps_pca=args.n_comps_pca,
            n_neighbors=args.n_neighbors,
            target_sum=args.target_sum,
        )

    print(f"Writing {args.output}...")
    adata.write(args.output)
    print(f"  adata shape: {adata.shape}")
    print(f"  guide_assignment shape: {adata.obsm['guide_assignment'].shape}")
    print(f"  guide_names={len(adata.uns['guide_names'])}, guide_targets={len(adata.uns['guide_targets'])}")
    print(f"  obsm keys: {list(adata.obsm.keys())}")


if __name__ == "__main__":
    main()
