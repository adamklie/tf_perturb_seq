import argparse
import time
import scanpy as sc
import logging
import os

def setup_logging(log_file):
    open(log_file, 'w').close()  # Overwrite existing log file
    
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def merge_anndata(input_files, output_file, names=None, names_key="sample", make_bcs_unique=True):
    start_time = time.time()
    logging.info("Starting merging of AnnData files.")
    
    if names is None:
        names = [str(i) for i in range(len(input_files))]
        logging.info("No sample names provided. Using indices as names.")
    
    adata_list = []
    samples = []
    
    for i, input_file in enumerate(input_files):
        sample = names[i]
        logging.info(f"Reading in AnnData for sample {sample}")
        adata = sc.read_h5ad(input_file)
        
        if make_bcs_unique:
            logging.info(f"Making barcodes unique by prepending sample name {sample}.")
            adata.obs.index = sample + "#" + adata.obs.index
        
        adata_list.append(adata)
        samples.append(sample)
    
    logging.info("Concatenating AnnData objects.")
    adata_merged = sc.concat(adata_list, label=names_key, keys=samples)
    
    adata_merged.write(output_file)
    logging.info(f"Saved merged AnnData file to {output_file}.")
    logging.info(f"Merging completed in {time.time() - start_time:.2f} seconds.")

def main():
    parser = argparse.ArgumentParser(description="Merge multiple AnnData files.")
    parser.add_argument("--inputs", type=str, nargs='+', required=True, help="Paths to input AnnData files.")
    parser.add_argument("--output", type=str, required=True, help="Path to save merged AnnData file.")
    parser.add_argument("--names", type=str, nargs='+', default=None, help="Sample names to use when merging.")
    parser.add_argument("--names_key", type=str, default="sample", help="Key to use for sample labels in metadata.")
    parser.add_argument("--make_bcs_unique", action='store_true', help="Make barcodes unique by prepending sample names.")
    parser.add_argument("--log_file", type=str, default="merge.log", help="Path to log file.")
    
    args = parser.parse_args()
    setup_logging(args.log_file)
    merge_anndata(args.inputs, args.output, args.names, args.names_key, args.make_bcs_unique)

if __name__ == "__main__":
    main()
