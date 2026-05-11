#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=Program_50_2_0
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/Plot/Program_50_2_0/logs/%j.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/Plot/Program_50_2_0/logs/%j.err
#SBATCH --partition=engreitz,owners,bigmem
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=256G

# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ymo@stanford.edu


# Define paths
BASE_DIR="/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7"
LOG_DIR="${BASE_DIR}/Plot/Program_50_2_0"

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
mkdir -p "$LOG_DIR/logs"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate NMF_Benchmarking

echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"


# Run the Python script
echo "Running Python script..."
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/PerturbNMF/src/Stage3_Interpretation/A_Plotting/Slurm_Version/cNMF_program_analysis.py \
        --mdata_path "${BASE_DIR}/adata/cNMF_50_2_0.h5mu" \
        --perturb_path_base "${BASE_DIR}/Eval/50_2_0/50_perturbation_association_results" \
        --GO_path "${BASE_DIR}/Eval/50_2_0/50_GO_term_enrichment.txt" \
        --top_program 5 \
        --p_value 0.05 \
        --pdf_save_path "$LOG_DIR" \
        --PDF \
        --sample WTC \
        --square_plots \
        --figsize 35 20\
        --categorical_key "batch" \
        --subsample_frac 1.0 \
        --gene_name_key "symbol"


# Calculate and print elapsed time at the end
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
HOURS=$((ELAPSED_TIME / 3600))
MINUTES=$(((ELAPSED_TIME % 3600) / 60))
SECONDS=$((ELAPSED_TIME % 60))

echo "Job completed at: $(date)"
echo "Total elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED_TIME} seconds)"
