#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=030726_k_sel_2_0_Huangfu
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/Plot/k_selection_2_0/logs/%j.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/Plot/k_selection_2_0/logs/%j.err
#SBATCH --partition=engreitz
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G

# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ymo@stanford.edu

# Define paths
OUT_DIR="/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Result"
RUN_NAME="030726_20iter_5KHVG_torch_halsvar_batch_e7"
RUN_DIR="$OUT_DIR/$RUN_NAME"

# Store start time
START_TIME=$(date +%s)

# Print some job information
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Working directory: $(pwd)"
echo "Number of CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Partition: $SLURM_JOB_PARTITION"
echo "Run directory: $RUN_DIR"

# Create logs directory if it doesn't exist
mkdir -p "$RUN_DIR/Plot/k_selection_2_0/logs"

# Activate conda environment
echo "Activating conda base environment..."
source activate torch-nmf-dl

echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"

# Run the Python script
echo "Running Python script..."
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/PerturbNMF/src/Stage3_Interpretation/A_Plotting/Slurm_Version/cNMF_k_selection.py \
        --output_directory "$OUT_DIR" \
        --run_name "$RUN_NAME" \
        --save_folder_name "$RUN_DIR/Plot/k_selection_2_0" \
        --eval_folder_name "$RUN_DIR/Eval" \
        --stability_file "$RUN_DIR/$RUN_NAME.k_selection_stats.df.npz" \
        --groupby "batch" \
        --K 5 6 7 8 9 10 11 12 13 14 15 17 19 21 23 25 27 30 35 40 45 50 55 60 70 80 90 100 150 200 \
        --sel_threshs 2.0 \
        --samples IGVFDS0471AYHF IGVFDS1260YCMC IGVFDS1889TBEY IGVFDS5642SPLX

# Calculate and print elapsed time at the end
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
HOURS=$((ELAPSED_TIME / 3600))
MINUTES=$(((ELAPSED_TIME % 3600) / 60))
SECONDS=$((ELAPSED_TIME % 60))

echo "Job completed at: $(date)"
echo "Total elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED_TIME} seconds)"
