import argparse
import os
import shutil
import pickle
from collections import defaultdict

import scanpy as sc
import muon
import pandas as pd
import numpy as np
from tqdm import tqdm
import synapseclient

import json

def get_promoter_name(row_info):
    # Determine the promoter name based on type
    if row_info["type"] == "non-targeting":
        return "non-targeting"
    else:
        # Fixed nested quotation marks inside the f-string
        promoter_name = f"{row_info['intended_target_name']}|"
        promoter_name += f"{row_info['intended_target_chr']}:{row_info['intended_target_start']:.0f}-{row_info['intended_target_end']:.0f}"
        return promoter_name

def main():
    # Set up argument parser for command line execution
    parser = argparse.ArgumentParser(description="Process single-cell RNA and sgRNA data from Synapse.")
    parser.add_argument("--synapse_id", required=True, help="Synapse ID to download (e.g., syn70753570)")
    parser.add_argument("--auth_token", required=True, help="Synapse personal access token for login")
    parser.add_argument("--out_dir", default="./data", help="Output directory for processed files")
    parser.add_argument("--config1_2_path", default="./data/config1_2.json",
                        help="path for the config file of step 1&2.")
    parser.add_argument("--config3_path", default="./data/config3.json",
                        help="path for the config file of step 3.")
    args = parser.parse_args()

    synapse_id = args.synapse_id
    out_dir = args.out_dir

    config1_2_path = args.config1_2_path
    config3_path = args.config3_path
    
    # Define output file paths
    final_path_name = os.path.join(out_dir, f"{synapse_id}.h5mu")
    dict_file = os.path.join(out_dir, "gRNA_dict.pickle")
    pca_file = os.path.join(out_dir, "pca_dataframe.pickle")
    annotation_file_path = os.path.join(out_dir, "annotation_table.csv")

    # Prepare data folder
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    # Login to Synapse using the provided auth token
    syn = synapseclient.Synapse()
    syn.login(authToken=args.auth_token)

    # Download mudata from Synapse
    print(f"Downloading {synapse_id} from Synapse...")
    entity = syn.get(synapse_id)
    print("Moving downloaded file to target directory...")
    shutil.move(entity.path, final_path_name)

    # Load mudata
    print("Loading MuData...")
    adata_mu = muon.read_h5mu(final_path_name)
    print(adata_mu)

    ### Process gene expression
    print("Processing gene expression...")
    adata_exp = adata_mu["gene"].copy()
    sc.pp.filter_genes(adata_exp, min_counts=1) # only consider genes with more than 1 count
    sc.pp.normalize_total(adata_exp)

    # Normalize and scale the matrix
    sc.pp.log1p(adata_exp)
    sc.pp.scale(adata_exp)

    # PCA
    print("Running PCA...")
    sc.tl.pca(adata_exp, random_state=0, n_comps=50)

    X = pd.DataFrame(adata_exp.obsm["X_pca"].copy(), index=adata_exp.obs.index)
    print(f"Saving PCA data to '{pca_file}'...")
    X.to_pickle(pca_file)

    ### Process gRNA
    print("Processing gRNA...")
    adata_guide = adata_mu["guide"].copy()

    gRNA_info_df = adata_guide.var.copy()
    promoter_inteded_names = gRNA_info_df.apply(lambda x: get_promoter_name(x), axis=1)
    gRNA_info_df["intended_target_promoter"] = promoter_inteded_names

    print(f"Saving annotation table to '{annotation_file_path}'...")
    gRNA_info_df.to_csv(annotation_file_path)

    # Modify gRNA dict
    print("Extracting non-zero counts for gRNA dictionary...")
    non_zero_rows, non_zero_cols = np.where((adata_guide.X != 0).todense())

    if len(non_zero_rows) == 0:
        print("Warning: No non-zero counts found in the sgRNA data. Returning an empty dictionary.")
        gRNA_dict = {}
    else:
        gRNA_name_list = adata_guide.var.index.to_list()
        cell_name_list = adata_guide.obs.index.to_list()

        print("Creating gRNA dictionary (improved method)...")
        # Build the dictionary efficiently using defaultdict
        gRNA_dict_default = defaultdict(list)
        num_pairs = len(non_zero_rows)
        try:
            for i in tqdm(range(num_pairs), desc="Processing gRNA pairs"):
                row_idx = non_zero_rows[i]
                col_idx = non_zero_cols[i]
                cell_name = cell_name_list[row_idx]
                gRNA_name = gRNA_name_list[col_idx]
                gRNA_dict_default[gRNA_name].append(cell_name)

            # Convert defaultdict back to a regular dict
            gRNA_dict = dict(gRNA_dict_default)
            print(f"gRNA dictionary creation complete. Found {len(gRNA_dict)} types of gRNAs.")
        except Exception as e:
            print(f"An error occurred during gRNA dictionary creation loop: {e}")

    # Save dictionary
    print(f"Saving gRNA dictionary to '{dict_file}'...")
    try:
        with open(dict_file, mode='wb') as fo:
            pickle.dump(gRNA_dict, fo)
        print("Processing complete!")
    except IOError as e:
        print(f"Error: Failed to write gRNA dictionary file '{dict_file}': {e}")
    except Exception as e:
        print(f"An unexpected error occurred while saving the gRNA dictionary: {e}")

    # Prepare config file for pipeline
    with open("./energy_dist_pipeline/config.json", 'r') as f:
        step12_config = json.load(f)
    
    with open("./energy_dist_pipeline/config_clustering.json", 'r') as f:
        step3_config = json.load(f)
    
    #Modify config file for this pipeline
    step12_config["output_file_name_list"]['OUTPUT_FOLDER'] = "./data"
    step12_config["output_file_name_list"]['pca_table'] = "pca_dataframe.pickle"
    step12_config["output_file_name_list"]['gRNA_dict'] = "gRNA_dict.pickle"
    step12_config["output_file_name_list"]['OVERWRITE_PCA_DICT'] = False
    
    step12_config["input_data"]["annotation_file"]["file_path"] = "./data/annotation_table.csv"
    step12_config["input_data"]["annotation_file"]["concatenate_key"] = "intended_target_promoter"
    
    step12_config["gRNA_filtering"]["perform_targeting_filtering"] = False
    step12_config["gRNA_filtering"]["perform_nontargeting_filtering"] = False
    
    #As this pipeline doesn't use h5ad or sgRNA dataframe fill with dummy
    step12_config["input_data"]["h5ad_file"]["file_path"] = "dummy"
    step12_config["input_data"]["sgRNA_file"]["file_path"] = "dummy"
    
    with open(config1_2_path, 'w') as f:
        json.dump(step12_config,f)
    with open(config3_path, 'w') as f:
        json.dump(step3_config,f)
    
if __name__ == "__main__":
    main()