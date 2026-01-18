import argparse
import time
import scanpy as sc
import logging
import os
import rapids_singlecell as rsc
import rmm
import cupy as cp
from rmm.allocators.cupy import rmm_cupy_allocator
import matplotlib.pyplot as plt

def setup_logging(log_file):
    open(log_file, 'w').close()  # Overwrite existing log file
    
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def setup_gpu():
    rmm.reinitialize(
        managed_memory=True,
        pool_allocator=False,
        devices=0,
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)

def cluster(input_path, output_path, plot_dir=None, n_pcs=40, integrate=None, n_neighbors=15, resolution=0.6):
    start_time = time.time()
    logging.info("Starting clustering on AnnData file.")
    
    setup_gpu()
    
    adata = sc.read_h5ad(input_path)
    
    # Perform PCA
    rsc.tl.pca(adata, n_comps=n_pcs)

    # Optionally integrate data
    if integrate:
        logging.info(f"Integrating data with Harmony using key {integrate}.")
        rsc.pp.harmony_integrate(adata, key=integrate)
        use_rep = "X_pca_harmony"
    else:
        use_rep = "X_pca"
    
    # Compute neighborhood graph
    rsc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep)
    
    # Compute UMAP
    rsc.tl.umap(adata)
    
    # Perform clustering
    rsc.tl.louvain(adata, resolution=resolution)
    
    # Move AnnData back to CPU
    rsc.get.anndata_to_CPU(adata)
    
    if plot_dir:
        os.makedirs(plot_dir, exist_ok=True)
        with plt.rc_context({'axes.facecolor':'white', 'figure.facecolor':'white'}):
            sc.pl.umap(adata, color=["louvain"], legend_loc="on data", show=False)
            plt.savefig(os.path.join(plot_dir, "umap_clustering.png"), bbox_inches="tight")
            plt.close()
    
    adata.write(output_path)
    
    logging.info(f"Saved clustered AnnData file to {output_path}.")
    logging.info(f"Clustering completed in {time.time() - start_time:.2f} seconds.")

def main():
    parser = argparse.ArgumentParser(description="Perform clustering on AnnData file.")
    parser.add_argument("--input", type=str, required=True, help="Path to input AnnData file.")
    parser.add_argument("--output", type=str, required=True, help="Path to save output AnnData file.")
    parser.add_argument("--plot_dir", type=str, default=None, help="Directory to save clustering plots.")
    parser.add_argument("--n_pcs", type=int, default=40, help="Number of principal components to use.")
    parser.add_argument("--integrate", type=str, default=None, help="Key for Harmony integration.")
    parser.add_argument("--n_neighbors", type=int, default=15, help="Number of neighbors for graph computation.")
    parser.add_argument("--resolution", type=float, default=0.6, help="Resolution parameter for Louvain clustering.")
    parser.add_argument("--log_file", type=str, default="clustering.log", help="Path to log file.")
    
    args = parser.parse_args()
    setup_logging(args.log_file)
    cluster(args.input, args.output, args.plot_dir, args.n_pcs, args.integrate, args.n_neighbors, args.resolution)

if __name__ == "__main__":
    main()
