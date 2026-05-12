#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=Perturb_gene_Huangfu
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/Plot/Perturb_gene_2_0/logs/%j.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/Plot/Perturb_gene_2_0/logs/%j.err
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
LOG_DIR="/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/Plot/Perturb_gene_2_0"

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

# Activate conda base environment
echo "Activating conda base environment..."
eval "$(conda shell.bash hook)" && conda activate NMF_Benchmarking

echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"

# Run the Python script
echo "Running Python script..."
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/PerturbNMF/src/Stage3_Interpretation/A_Plotting/Slurm_Version/cNMF_perturbed_gene_analysis.py \
        --mdata_path "/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/adata/cNMF_50_2_0.h5mu" \
        --perturb_path_base "/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/Eval/50_2_0/50_perturbation_association_results" \
        --top_n_programs 10 \
        --perturb_target_col "target_name" \
        --perturb_program_col "program_name" \
        --perturb_log2fc_col "log2FC" \
        --top_corr_genes 5 \
        --significance_threshold 0.05 \
        --volcano_log2fc_min -0.00 \
        --volcano_log2fc_max 0.00 \
        --save_path "$LOG_DIR" \
        --square_plots \
        --figsize 35 20 \
        --sample WTC \
        --PDF \
        --n_processes -1 \
        --umap_dot_size 10 \
        --data_key 'rna' \
        --prog_key 'cNMF' \
        --categorical_key 'batch' \
        --gene_name_key 'symbol' \
        --control_target_name 'non-targeting' \
        --subsample_frac 1.0 \
        --parallel

# Calculate and print elapsed time at the end
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
HOURS=$((ELAPSED_TIME / 3600))
MINUTES=$(((ELAPSED_TIME % 3600) / 60))
SECONDS=$((ELAPSED_TIME % 60))

echo "Job completed at: $(date)"
echo "Total elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED_TIME} seconds)"
