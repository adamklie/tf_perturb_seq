import argparse
import time
import scanpy as sc
import logging
import os
import rapids_singlecell as rsc
import matplotlib.pyplot as plt
import rmm
import cupy as cp
from rmm.allocators.cupy import rmm_cupy_allocator

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

def preprocess(input_path, output_path, plot_dir=None, target_sum=1e4, n_top_genes=5000, max_value=10):
    start_time = time.time()
    logging.info("Starting preprocessing on AnnData file.")
    
    setup_gpu()
    
    adata = sc.read_h5ad(input_path)
    
    # Move AnnData to GPU
    rsc.get.anndata_to_GPU(adata)
    
    # Normalize
    rsc.pp.normalize_total(adata, target_sum=target_sum)
    rsc.pp.log1p(adata)
    
    # Select highly variable genes
    rsc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="cell_ranger")
    if plot_dir:
        os.makedirs(plot_dir, exist_ok=True)
        with plt.rc_context({'axes.facecolor':'white', 'figure.facecolor':'white'}):
            sc.pl.highly_variable_genes(adata, show=False)
            plt.savefig(os.path.join(plot_dir, "highly_variable_genes.png"), bbox_inches="tight")
            plt.close()
    adata.raw = adata.copy()
    rsc.pp.filter_highly_variable(adata)
    
    # Regress out confounders
    rsc.pp.regress_out(adata, keys=["total_gene_umis", "percent_mito"])
    
    # Scale data
    rsc.pp.scale(adata, max_value=max_value)
    
    # Move AnnData back to CPU
    rsc.get.anndata_to_CPU(adata)
    
    adata.write(output_path)
    
    logging.info(f"Saved preprocessed AnnData file to {output_path}.")
    logging.info(f"Preprocessing completed in {time.time() - start_time:.2f} seconds.")

def main():
    parser = argparse.ArgumentParser(description="Preprocess AnnData file.")
    parser.add_argument("--input", type=str, required=True, help="Path to input AnnData file.")
    parser.add_argument("--output", type=str, required=True, help="Path to save output AnnData file.")
    parser.add_argument("--plot_dir", type=str, default=None, help="Directory to save preprocessing plots.")
    parser.add_argument("--target_sum", type=float, default=1e4, help="Target sum for normalization.")
    parser.add_argument("--n_top_genes", type=int, default=5000, help="Number of top variable genes to select.")
    parser.add_argument("--max_value", type=float, default=10, help="Maximum value for scaling.")
    parser.add_argument("--log_file", type=str, default="preprocess.log", help="Path to log file.")
    
    args = parser.parse_args()
    setup_logging(args.log_file)
    preprocess(args.input, args.output, args.plot_dir, args.target_sum, args.n_top_genes, args.max_value)

if __name__ == "__main__":
    main()
