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

def convert_10x_to_anndata(input_path, output_path, subset=None, random_seed=1234):
    start_time = time.time()
    logging.info("Starting conversion of 10X h5 file to AnnData.")
    
    adata = sc.read_10x_h5(input_path)
    adata.var_names_make_unique()
    
    if subset:
        logging.info(f"Subsampling {subset} cells with random seed {random_seed}.")
        sc.pp.subsample(adata, n_obs=subset, random_state=random_seed)
    
    adata.write(output_path)
    
    logging.info(f"Saved AnnData file to {output_path}.")
    logging.info(f"Conversion completed in {time.time() - start_time:.2f} seconds.")

def main():
    parser = argparse.ArgumentParser(description="Convert 10X h5 file to AnnData format.")
    parser.add_argument("--input", type=str, required=True, help="Path to 10X h5 file.")
    parser.add_argument("--output", type=str, required=True, help="Path to save output AnnData file.")
    parser.add_argument("--subset", type=int, default=None, help="Number of cells to randomly select.")
    parser.add_argument("--random_seed", type=int, default=1234, help="Random seed for subsampling.")
    parser.add_argument("--log_file", type=str, default="anndata.log", help="Path to log file.")
    
    args = parser.parse_args()
    setup_logging(args.log_file)
    convert_10x_to_anndata(args.input, args.output, args.subset, args.random_seed)

if __name__ == "__main__":
    main()
