#!/usr/bin/env python3
"""
Integrate gene expression across benchmark datasets.

Reads gene modality from each dataset's MuData, merges, normalizes,
finds HVGs, runs PCA, Harmony batch correction, Leiden clustering, and UMAP.

Usage:
    python 1_integrate_gene.py --run-label cleanser_unified
"""

import argparse
import logging
from pathlib import Path

import numpy as np
import scanpy as sc
import mudata as mu
import harmonypy as hm

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

PROJECT_ROOT = Path("/cellar/users/aklie/projects/tf_perturb_seq")
BENCHMARK_DIR = PROJECT_ROOT / "datasets" / "technology-benchmark_WTC11_TF-Perturb-seq"

DATASETS = [
    "Hon_WTC11-benchmark_TF-Perturb-seq",
    "Huangfu_WTC11-benchmark_TF-Perturb-seq",
    "Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3",
    "Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2",
    "Engreitz_WTC11-benchmark_TF-Perturb-seq",
]

SUFFIX = "_WTC11-benchmark_TF-Perturb-seq"


def short_name(ds: str) -> str:
    return ds.replace(SUFFIX, "").replace("_", " ")


def load_gene_adata(dataset: str, run_label: str) -> sc.AnnData:
    """Load gene modality from a dataset's MuData."""
    run_dir = PROJECT_ROOT / "datasets" / dataset / "runs" / run_label
    mudata_path = run_dir / "pipeline_dashboard" / "inference_mudata.h5mu"
    if not mudata_path.exists():
        mudata_path = run_dir / "pipeline_outputs" / "inference_mudata.h5mu"

    logger.info(f"Loading {dataset} from {mudata_path}")
    mdata = mu.read_h5mu(str(mudata_path))
    adata = mdata.mod["gene"].copy()

    # Add sample metadata
    adata.obs["sample"] = dataset
    adata.obs["sample_short"] = short_name(dataset)

    # Propagate batch info
    if "batch" in adata.obs.columns:
        adata.obs["lane"] = adata.obs["batch"].astype(str)

    # Make barcodes unique across samples
    adata.obs_names = [f"{dataset}__{bc}" for bc in adata.obs_names]
    adata.obs_names_make_unique()

    logger.info(f"  {adata.n_obs} cells, {adata.n_vars} genes")
    return adata


def main():
    parser = argparse.ArgumentParser(description="Integrate gene expression across benchmark datasets")
    parser.add_argument("--run-label", default="cleanser_unified", help="Pipeline run label")
    parser.add_argument("--n-hvgs", type=int, default=3000, help="Number of highly variable genes")
    parser.add_argument("--n-pcs", type=int, default=50, help="Number of PCs")
    parser.add_argument("--resolutions", type=float, nargs="+", default=[0.2, 0.5, 0.8, 1.0, 1.5],
                        help="Leiden clustering resolutions")
    parser.add_argument("--subsample", type=int, default=None,
                        help="Subsample N cells per dataset for testing")
    args = parser.parse_args()

    outdir = BENCHMARK_DIR / "results" / "integration" / args.run_label
    outdir.mkdir(parents=True, exist_ok=True)

    # --- Step 1: Load and merge ---
    logger.info("Step 1: Loading gene modalities")
    adatas = []
    for ds in DATASETS:
        try:
            adata = load_gene_adata(ds, args.run_label)
            adatas.append(adata)
        except Exception as e:
            logger.warning(f"Failed to load {ds}: {e}")

    # Subsample if requested
    if args.subsample:
        logger.info(f"Subsampling to {args.subsample} cells per dataset")
        subsampled = []
        for a in adatas:
            n = min(args.subsample, a.n_obs)
            sc.pp.subsample(a, n_obs=n, random_state=42)
            subsampled.append(a)
        adatas = subsampled

    logger.info(f"Merging {len(adatas)} datasets")
    adata = sc.concat(adatas, join="inner", label="sample_concat", keys=[short_name(ds) for ds in DATASETS[:len(adatas)]])
    logger.info(f"Merged: {adata.n_obs} cells, {adata.n_vars} genes")

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    # --- Step 2: QC metrics ---
    logger.info("Step 2: QC metrics")
    adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # --- Step 3: Normalize ---
    logger.info("Step 3: Normalize")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["log1p_norm"] = adata.X.copy()

    # --- Step 4: HVG selection ---
    logger.info(f"Step 4: Finding {args.n_hvgs} HVGs")
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvgs, batch_key="sample")
    logger.info(f"  {adata.var['highly_variable'].sum()} HVGs selected")

    # --- Step 5: Scale and PCA ---
    logger.info(f"Step 5: PCA ({args.n_pcs} components)")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=args.n_pcs, svd_solver="arpack")

    # Transfer PCA to full adata
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.uns["pca"] = adata_hvg.uns["pca"]
    adata.varm["PCs"] = np.zeros((adata.n_vars, args.n_pcs))
    hvg_idx = np.where(adata.var["highly_variable"])[0]
    adata.varm["PCs"][hvg_idx] = adata_hvg.varm["PCs"]

    del adata_hvg

    # --- Step 6: Harmony batch correction ---
    logger.info("Step 6: Harmony batch correction")
    harmony_out = hm.run_harmony(adata.obsm["X_pca"], adata.obs, "sample", max_iter_harmony=20)
    Z = harmony_out.result()
    logger.info(f"  Harmony output shape: {Z.shape}")
    n_nan = np.isnan(Z).sum()
    if n_nan > 0:
        logger.warning(f"  {n_nan} NaN values in Harmony output — replacing with 0")
        Z = np.nan_to_num(Z, nan=0.0)
    assert Z.shape == (adata.n_obs, args.n_pcs), f"Shape mismatch: {Z.shape} vs ({adata.n_obs}, {args.n_pcs})"
    adata.obsm["X_harmony"] = Z

    # --- Step 7: Neighbors, clustering, UMAP ---
    logger.info("Step 7: Neighbors + Leiden + UMAP")

    # Uncorrected
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=30, key_added="pca_neighbors")
    sc.tl.umap(adata, neighbors_key="pca_neighbors")
    adata.obsm["X_umap_uncorrected"] = adata.obsm["X_umap"].copy()

    # Harmony-corrected
    sc.pp.neighbors(adata, use_rep="X_harmony", n_neighbors=30)
    sc.tl.umap(adata, min_dist=0.3, spread=1.0)

    for res in args.resolutions:
        key = f"leiden_{res}"
        logger.info(f"  Clustering at resolution {res}")
        sc.tl.leiden(adata, resolution=res, key_added=key)

    # --- Step 8: Save ---
    out_path = outdir / "integrated_gene.h5ad"
    logger.info(f"Step 8: Saving to {out_path}")
    # Restore counts as X for saving
    adata.X = adata.layers["counts"]
    adata.write_h5ad(str(out_path))
    logger.info(f"Done: {adata.n_obs} cells, {adata.n_vars} genes")

    # Also save a slim version
    slim = adata.copy()
    del slim.layers
    for key in list(slim.obsm.keys()):
        if key not in ["X_umap", "X_harmony", "X_pca"]:
            del slim.obsm[key]
    slim_path = outdir / "integrated_gene_slim.h5ad"
    slim.write_h5ad(str(slim_path))
    logger.info(f"Slim version: {slim_path}")


if __name__ == "__main__":
    main()
