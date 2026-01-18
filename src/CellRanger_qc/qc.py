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
        managed_memory=True,  # Allows oversubscription
        pool_allocator=False,  # default is False
        devices=0,  # GPU device IDs to register. By default registers only GPU 0.
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)

def perform_qc(input_path, output_path, plot_dir=None, min_genes=500, max_mt_pct=20, min_gene_count=3):
    start_time = time.time()
    logging.info("Starting quality control on AnnData file.")
    
    setup_gpu()
    
    adata = sc.read_h5ad(input_path)
    
    # Move AnnData to GPU
    rsc.get.anndata_to_GPU(adata)
    
    # Flag mitochondrial and ribosomal genes
    rsc.pp.flag_gene_family(adata, gene_family_name="MT", gene_family_prefix="MT")
    rsc.pp.flag_gene_family(adata, gene_family_name="RIBO", gene_family_prefix="RPS")
    
    # Compute QC metrics
    rsc.pp.calculate_qc_metrics(adata, qc_vars=["MT", "RIBO"])
    
    # Filtering
    adata = adata[adata.obs["n_genes_by_counts"] >= min_genes]
    adata = adata[adata.obs["pct_counts_MT"] <= max_mt_pct]
    rsc.pp.filter_genes(adata, min_count=min_gene_count)
    logging.info(f"Filtered cells with >={min_genes} genes and <={max_mt_pct}% mitochondrial content.")
    logging.info(f"Filtered genes with <{min_gene_count} counts.")
    
    if plot_dir:
        os.makedirs(plot_dir, exist_ok=True)
        with plt.rc_context({'axes.facecolor':'white', 'figure.facecolor':'white'}):
            sc.pl.scatter(adata, x="total_counts", y="pct_counts_MT", show=False)
            plt.savefig(os.path.join(plot_dir, "total_counts_vs_pct_MT.png"), bbox_inches="tight")
            plt.close()
            
            sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", show=False)
            plt.savefig(os.path.join(plot_dir, "total_counts_vs_n_genes.png"), bbox_inches="tight")
            plt.close()
            
            sc.pl.violin(adata, "n_genes_by_counts", jitter=0.4, show=False)
            plt.savefig(os.path.join(plot_dir, "violin_n_genes.png"), bbox_inches="tight")
            plt.close()
            
            sc.pl.violin(adata, "total_counts", jitter=0.4, show=False)
            plt.savefig(os.path.join(plot_dir, "violin_total_counts.png"), bbox_inches="tight")
            plt.close()
            
            sc.pl.violin(adata, "pct_counts_MT", jitter=0.4, show=False)
            plt.savefig(os.path.join(plot_dir, "violin_pct_MT.png"), bbox_inches="tight")
            plt.close()
    
    adata.write(output_path)
    
    logging.info(f"Saved QC-processed AnnData file to {output_path}.")
    logging.info(f"Quality control completed in {time.time() - start_time:.2f} seconds.")

def main():
    parser = argparse.ArgumentParser(description="Perform quality control on AnnData file.")
    parser.add_argument("--input", type=str, required=True, help="Path to input AnnData file.")
    parser.add_argument("--output", type=str, required=True, help="Path to save output AnnData file.")
    parser.add_argument("--plot_dir", type=str, default=None, help="Directory to save QC plots.")
    parser.add_argument("--min_genes", type=int, default=500, help="Minimum number of genes detected per cell.")
    parser.add_argument("--max_mt_pct", type=float, default=20, help="Maximum mitochondrial gene percentage.")
    parser.add_argument("--min_gene_count", type=int, default=3, help="Minimum count threshold for genes.")
    parser.add_argument("--log_file", type=str, default="qc.log", help="Path to log file.")
    
    args = parser.parse_args()
    setup_logging(args.log_file)
    perform_qc(args.input, args.output, args.plot_dir, args.min_genes, args.max_mt_pct, args.min_gene_count)

if __name__ == "__main__":
    main()
