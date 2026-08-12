import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import mudata as md
import anndata as ad

KEY_PROGRAMS = [32, 59, 25, 29, 5, 37]

def main(mudata_path, usage_path, downsample=200000, out_prefix="program_umap"):

    print("Opening MuData (backed)...")
    mdata = md.read_h5mu(mudata_path, backed="r")
    adata_gene = mdata["gene"]

    n_obs = adata_gene.n_obs
    print(f"Dataset: {n_obs:,} cells")

    # Downsample cells
    if downsample and downsample < n_obs:
        np.random.seed(0)
        sampled_idx = np.random.choice(n_obs, downsample, replace=False)
    else:
        sampled_idx = np.arange(n_obs)

    print(f"Using {len(sampled_idx):,} cells")

    # Load expression matrix for sampled cells
    print("Loading expression subset...")
    X_sub = adata_gene.X[sampled_idx, :]

    adata = ad.AnnData(
        X=X_sub,
        obs=adata_gene.obs.iloc[sampled_idx].copy(),
        var=adata_gene.var.copy()
    )

    mdata.file.close()

    print(f"Subset AnnData: {adata.n_obs:,} x {adata.n_vars:,}")

    # Compute UMAP (same procedure as other script)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")
    adata = adata[:, adata.var.highly_variable].copy()

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata, min_dist=0.3)

    print("Loading deltaNMF usages...")
    usages = pd.read_csv(usage_path, sep="\t", index_col=0)

    common = adata.obs_names.intersection(usages.index)
    adata = adata[common].copy()
    usages = usages.loc[common]

    for p in KEY_PROGRAMS:
        col = f"case_{p}"
        adata.obs[col] = usages[col]

    print("Plotting...")
    sc.pl.umap(
        adata,
        color=[f"case_{p}" for p in KEY_PROGRAMS],
        cmap="viridis",
        ncols=3,
        save=f"_{out_prefix}.png"
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mudata")
    parser.add_argument("usages")
    parser.add_argument("--downsample", type=int, default=200000)
    parser.add_argument("--out_prefix", default="program_umap")

    args = parser.parse_args()

    main(args.mudata, args.usages, args.downsample, args.out_prefix)

# Sample run: srun -c 4 --account=singhlab --mem=100G --time=02:00:00 python visualize_programs_umap.py \
"/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/helen_output_poolabcd_prod_v9/mergedresults/inference_mudata.h5mu" \
"/hpc/group/gersbachlab/agk21/hep_perturbseq/deltanmf_fm/H_stage2_usages_barcoded.tsv"