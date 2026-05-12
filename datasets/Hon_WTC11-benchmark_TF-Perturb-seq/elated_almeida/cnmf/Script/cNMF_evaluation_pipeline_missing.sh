#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=030726_eval_missing
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Hon_WTC11/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/Eval/logs/%j.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Hon_WTC11/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/Eval/logs/%j.err
#SBATCH --partition=engreitz
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G

# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ymo@stanford.edu

# Define the cNMF case
OUT_DIR="/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Hon_WTC11/Result"
RUN_NAME="030726_20iter_5KHVG_torch_halsvar_batch_e7"
LOG_DIR="$OUT_DIR/$RUN_NAME"

# Store start time
START_TIME=$(date +%s)

# Print some job information
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Working directory: $(pwd)"
echo "Number of CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Partition: $SLURM_JOB_PARTITION"
echo "Log directory: $LOG_DIR"


# Create logs directory if it doesn't exist
mkdir -p "$LOG_DIR/Eval/logs"

# Activate conda environment
echo "Activating conda environment..."
eval "$(conda shell.bash hook)"
conda activate NMF_Benchmarking


echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"


# Add 'WTC' batch column to h5mu files for perturbation association
echo "Adding 'WTC' batch column to h5mu files..."
python3 -c "
import mudata as mu

adata_dir = '$OUT_DIR/$RUN_NAME/Inference/adata'

K_values = [5,6,7,8,9,10,11,12,13,14,15,17,19,21,23,25,27,30,35,40,45,50,55,60,70,80,90,100,150,200]
thresh_strs = ['0_2', '2_0']
files = [f'{adata_dir}/cNMF_{k}_{t}.h5mu' for k in K_values for t in thresh_strs]

for f in files:
    print(f'Processing {f}...')
    mdata = mu.read(f)
    mdata['rna'].obs['WTC'] = 'WTC'
    mdata.write(f)
    print(f'  Added WTC column and saved.')
print('Done.')
"

# Run the Python script — only missing K values at thresh 2.0
echo "Running Python script for missing K=25,27,30,35,40 at thresh 2.0..."
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation/Slurm_Version/cNMF_evaluation_pipeline.py \
        --out_dir "$OUT_DIR" \
        --run_name "$RUN_NAME" \
        --X_normalized_path "$LOG_DIR/cnmf_tmp/$RUN_NAME.norm_counts.h5ad" \
        --Perform_perturbation \
        --data_key 'rna' \
        --prog_key 'cNMF' \
        --categorical_key 'WTC' \
        --organism 'human' \
        --data_guide_path "/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Hon_WTC11/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/adata/cNMF_5_0_2.h5mu" \
        --gwas_data_path '/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation/Resources/OpenTargets_L2G_Filtered.csv.gz' \
        --sel_thresh 0.2 2.0 \
        --K 5 6 7 8 9 10 11 12 13 14 15 17 19 21 23 25 27 30 35 40 45 50 55 60 70 80 90 100 150 200 \
        --FDR_method "StoreyQ"


# Calculate and print elapsed time at the end
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
HOURS=$((ELAPSED_TIME / 3600))
MINUTES=$(((ELAPSED_TIME % 3600) / 60))
SECONDS=$((ELAPSED_TIME % 60))

echo "Job completed at: $(date)"
echo "Total elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED_TIME} seconds)"
