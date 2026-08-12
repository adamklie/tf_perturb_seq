import sys
import os

import pandas as pd
import numpy as np
import scanpy as sc

from scipy.sparse import coo_matrix

from tqdm import tqdm

import gc
import warnings
import pickle
import json
from collections import defaultdict

### Helper Functions
def load_sgRNA_sparse(sgRNA_file, cell_name_list=None):
    """
    Loads sgRNA data and processes it into a dictionary with maximum optimization.
    This version avoids transposing the data and filters cells first.
    """
    print("\n--- Processing gRNA dictionary (Fully Optimized Version) ---")
    print(f"Loading sgRNA file '{sgRNA_file}'...")
    
    # Load the initial data. Assume rows are gRNAs and columns are cells.
    sgRNA_data = pd.read_pickle(sgRNA_file)
    
    # --- Optimization 1: Filter columns (cells) first, before any heavy operations ---
    if cell_name_list is not None:
        print("Filtering by cells (columns) present in input data...")
        # Find which of the requested cells are actually in the columns
        available_cells = [cell for cell in cell_name_list if cell in sgRNA_data.columns]
        if not available_cells:
            print("Warning: None of the provided cells were found. Returning an empty dictionary.")
            return {}
        # Reduce the DataFrame to only the necessary columns
        sgRNA_data = sgRNA_data[available_cells]

    if sgRNA_data.empty:
        print("Warning: DataFrame is empty. Returning an empty dictionary.")
        return {}
        
    print("Converting filtered data to sparse matrix (COO format)...")
    try:
        # --- Optimization 2: No transpose needed. Directly use the data structure ---
        # Get gRNA names from the index and cell names from the columns
        gRNA_names = sgRNA_data.index.to_list()
        cell_names = sgRNA_data.columns.to_list()
        
        # Convert the filtered DataFrame values to a COO sparse matrix
        sparse_sgRNA = coo_matrix(sgRNA_data.values)
        
        # Free up memory by deleting the DataFrame
        print("Releasing dense sgRNA data from memory.")
        del sgRNA_data
        gc.collect()

    except Exception as e:
        print(f"An error occurred during sparse matrix conversion: {e}")
        return None

    if sparse_sgRNA.nnz == 0:
        print("Warning: No non-zero counts found in the sgRNA data. Returning an empty dictionary.")
        return {}

    print("Creating gRNA dictionary from sparse matrix...")
    gRNA_dict_default = defaultdict(list)
    
    try:
        # Iterate over non-zero elements.
        # Note the mapping change: row index maps to gRNA, col index maps to cell.
        for row_idx, col_idx in tqdm(zip(sparse_sgRNA.row, sparse_sgRNA.col), total=sparse_sgRNA.nnz, desc="Processing gRNA pairs"):
            gRNA_name = gRNA_names[row_idx]
            cell_name = cell_names[col_idx]
            gRNA_dict_default[gRNA_name].append(cell_name)
            
        gRNA_dict = dict(gRNA_dict_default)
        print(f"gRNA dictionary creation complete. Found {len(gRNA_dict)} types of gRNAs.")

    except Exception as e:
        print(f"An error occurred during gRNA dictionary creation: {e}")
        return None
    finally:
        # Clean up
        del sparse_sgRNA, gRNA_dict_default, gRNA_names, cell_names
        gc.collect()

    print("\n--- Processing finished ---")
    return gRNA_dict

def gini(x):
    mad = np.abs(np.subtract.outer(x, x)).astype(float).mean()
    # Relative mean absolute difference
    rmad = mad/np.mean(x)
    # Gini coefficient
    g = 0.5 * rmad
    return g

### Load configs
config_filepath = sys.argv[1]

try:
    with open(config_filepath, 'r', encoding='utf-8') as f:
        # json.load() reads from a file object
        config_data = json.load(f)

    # --- Store values into variables ---

    # Extract data from the 'input_data' section
    input_config = config_data.get('input_data', {}) # Use .get() for safety
    input_file = input_config.get('input_file')
    sgRNA_file = input_config.get('sgRNA_file')
    
    # Extract data from the 'output_file' section
    output_config = config_data.get('output_file', {})
    
    output_folder = output_config.get('output_folder',"./pipeline_output")
    top50_clonal_info = output_config.get('top50_clonal_info',True)
    top50_clonal_info_filename = output_config.get('top50_clonal_info_filename',"top50_clones.csv")
    gini_index_table = output_config.get('gini_index_table',True)
    gini_index_table_filename = output_config.get('gini_index_table_filename',"gini_index_analysis.csv")
    output_transcriptome_h5ad = output_config.get('output_transcriptome_h5ad',True)
    output_transcriptome_h5ad_filename = output_config.get('output_transcriptome_h5ad_filename',"clonal_cell_removal_transcriptome.h5ad")
    
    os.makedirs(output_folder, exist_ok=True)
    
    # --- Print the variables to verify ---
    print("--- Input Data ---")
    print(f"Input H5AD file: {input_file}")
    print(f"sgRNA file: {sgRNA_file}")
    print("\n--- Output Settings ---")
    print(f"Output folder: {output_folder}")
    print(f"Generate top 50 clonal info: {top50_clonal_info}")
    print(f"Top 50 clonal info filename: {top50_clonal_info_filename}")
    print(f"Generate Gini index table: {gini_index_table}")
    print(f"Gini index table filename: {gini_index_table_filename}")
    print(f"Output transcriptome H5AD: {output_transcriptome_h5ad}")
    print(f"Output transcriptome filename: {output_transcriptome_h5ad_filename}")


except FileNotFoundError:
    print(f"Error: The file '{config_filepath}' was not found.")
except json.JSONDecodeError:
    print(f"Error: Could not decode JSON from the file '{config_filepath}'.")

### Main functions
print("Loading h5ad")
adata = sc.read_h5ad(input_file)
valid_cell_name = adata.obs.index.to_list()

print("Loading gRNA df")
gRNA_dict = load_sgRNA_sparse(sgRNA_file,valid_cell_name)

gRNA_per_cell_dict = {}
for key in tqdm(gRNA_dict):
    for cell_name in gRNA_dict[key]:
        if not cell_name in gRNA_per_cell_dict.keys():
            gRNA_per_cell_dict[cell_name] = [key]
        else:
            gRNA_per_cell_dict[cell_name] += [key]
            
print(f"Number of cells before filtering: {len(gRNA_per_cell_dict.keys())}")

cell_per_clone_dict = {}

for key in tqdm(gRNA_per_cell_dict.keys()):
    gRNA_list = gRNA_per_cell_dict[key]
    gRNA_list.sort()
    
    gRNA_name_total = ""
    for gRNA_name in gRNA_list:
        gRNA_name_total += gRNA_name+"$$"
        
    if gRNA_name_total in cell_per_clone_dict.keys():
        cell_per_clone_dict[gRNA_name_total] += [key]
    else:
        cell_per_clone_dict[gRNA_name_total] = [key]
    
clonal_list = np.array(list(cell_per_clone_dict.keys()))
print(f"Number of clones: {len(cell_per_clone_dict.keys())}")

if top50_clonal_info:
    unique_clone_list = cell_per_clone_dict.keys()
    unique_clone_count = [len(cell_per_clone_dict[key]) for key in cell_per_clone_dict.keys()]

    gRNA_df = pd.DataFrame([unique_clone_list,unique_clone_count]).T
    gRNA_df.columns=["clone_name","count"]

    gRNA_df = gRNA_df.sort_values("count",ascending=False)
    gRNA_df.head(50).to_csv(os.path.join(output_folder,top50_clonal_info_filename))
    
    
### Statistics for each clones
res_output = {}
for target in tqdm(gRNA_dict.keys()):
    res_output[target] = {}
    target_cell = gRNA_dict[target]
    clonal_list = []
        
    for key in target_cell:
        gRNA_list = gRNA_per_cell_dict[key]
        gRNA_list.sort()

        gRNA_name_total = ""
        for gRNA_name in gRNA_list:
            gRNA_name_total += str(gRNA_name+"$$")
        
        if gRNA_name_total!="":
            clonal_list += [gRNA_name_total]
    clonal_list = np.array(clonal_list)

    gRNA_name_list,gRNA_cell_count_list = np.unique(clonal_list,return_counts=True)
    gRNA_num_count_list = np.array([len(label.split("$$"))-1 for label in gRNA_name_list])
    
    res_output[target]["num_total_cells"] = len(target_cell)
    res_output[target]["num_clones"] = len(gRNA_name_list)
    res_output[target]["num_single_gRNA"] = np.sum(gRNA_num_count_list==1)
    
    gini_coef = gini(gRNA_cell_count_list[gRNA_num_count_list>1])
    res_output[target]["gini_coef"] = gini_coef

gRNA_analysis_df = pd.DataFrame(res_output).T
gRNA_analysis_df = gRNA_analysis_df.sort_values("gini_coef",ascending=False)

if gini_index_table:
    gRNA_analysis_df.to_csv(os.path.join(output_folder,gini_index_table_filename))
    
### Filterin 1 cell per clone
normalized_cell_list = []

for clone_name in cell_per_clone_dict:
    gRNA_num_count = len(clone_name.split("$$"))-1
    if gRNA_num_count == 1:
        normalized_cell_list += cell_per_clone_dict[clone_name]
    else:
        normalized_cell_list += [cell_per_clone_dict[clone_name][0]]
        
print(f"Number of cells after filtering: (1 cell per clone): {len(normalized_cell_list)}")
###output data
if output_transcriptome_h5ad:
    adata_non_clone = adata[normalized_cell_list].copy()
    adata_non_clone.write_h5ad(os.path.join(output_folder,output_transcriptome_h5ad_filename))