"""
Rename gene spectra score file columns from Ensembl IDs to gene symbols.

Reads the Ensembl ID -> gene symbol mapping from cNMF_5_0_2.h5mu,
saves the mapping as a CSV, then renames columns in all 60
gene_spectra_score TSV files in the result directory.
"""

import glob
import os

import muon as mu
import pandas as pd

# --- Paths ---
PROJECT = "/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Hon_WTC11"
H5MU_PATH = os.path.join(
    PROJECT,
    "Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/adata/cNMF_5_0_2.h5mu",
)
RESULT_DIR = os.path.join(
    PROJECT, "Result/030726_20iter_5KHVG_torch_halsvar_batch_e7"
)
MAPPING_CSV = os.path.join(PROJECT, "Data/ensembl_to_symbol.csv")

# --- Step 1: Load mapping from h5mu ---
print("Loading mapping from h5mu...")
mdata = mu.read_h5mu(H5MU_PATH, backed="r")
rna_var = mdata.mod["rna"].var
mapping = dict(zip(rna_var.index, rna_var["symbol"]))
mdata.file.close()
print(f"  Loaded {len(mapping)} gene mappings")

# --- Step 2: Save mapping as CSV ---
mapping_df = pd.DataFrame(
    list(mapping.items()), columns=["ensembl_id", "symbol"]
)
mapping_df.to_csv(MAPPING_CSV, index=False)
print(f"  Saved mapping to {MAPPING_CSV}")

# --- Step 3: Find all gene_spectra_score files ---
pattern = os.path.join(RESULT_DIR, "*.gene_spectra_score.*.txt")
files = sorted(glob.glob(pattern))
print(f"  Found {len(files)} gene_spectra_score files")

# --- Step 4: Rename columns in each file ---
for i, fpath in enumerate(files, 1):
    df = pd.read_csv(fpath, sep="\t", index_col=0)
    df.columns = [mapping.get(col, col) for col in df.columns]
    df.to_csv(fpath, sep="\t")
    if i % 10 == 0 or i == len(files):
        print(f"  Processed {i}/{len(files)} files")

print("Done.")
