#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=030726_20iter_5KHVG_torch_halsvar_batch_e7           # Job name
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/Eval/logs/%j.out      # Output file (%j = job ID)
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/Eval/logs/%j.err       # Error file
#SBATCH --partition=engreitz            # partition name
#SBATCH --time=48:00:00                 # Time limit
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=20              # CPUs per task
#SBATCH --mem=128G                       # Memory per node

# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL      # Send email at start, end, and on failure
#SBATCH --mail-user=ymo@stanford.edu    # Email address

# Define the cNMF case
OUT_DIR="/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Result"
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

# Activate conda base environment
echo "Activating conda environment..."
eval "$(conda shell.bash hook)"
conda activate NMF_Benchmarking


echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"


# Run the Python script
echo "Running Python script..."
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation/Slurm_Version/cNMF_evaluation_pipeline.py \
        --out_dir "$OUT_DIR" \
        --run_name "$RUN_NAME" \
        --X_normalized_path "$LOG_DIR/cnmf_tmp/$RUN_NAME.norm_counts.h5ad" \
        --Perform_explained_variance \
        --Perform_categorical \
        --Perform_perturbation \
        --Perform_geneset \
        --Perform_trait \
        --data_key 'rna' \
        --prog_key 'cNMF' \
        --categorical_key 'batch' \
        --organism 'human' \
        --data_guide_path "$LOG_DIR/adata/cNMF_5_0_2.h5mu" \
        --gwas_data_path '/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation/Resources/OpenTargets_L2G_Filtered.csv.gz' \
        --reference_gtf_path "/oak/stanford/groups/engreitz/Users/opushkar/genome/IGVFFI9573KOZR.gtf.gz" \
        --sel_thresh 0.2 2.0 \
        --K 5 6 7 8 9 10 11 12 50 55 60 70 80 90 100 150 200 \
        --FDR_method "StoreyQ"
        #--Perform_motif \
        #--check_format \
        #--guide_annotation_path "/oak/stanford/groups/engreitz/Users/ymo/cc-perturb-seq/Data/guide/guide_metadata_v43.tsv" \





# Calculate and print elapsed time at the end
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
HOURS=$((ELAPSED_TIME / 3600))
MINUTES=$(((ELAPSED_TIME % 3600) / 60))
SECONDS=$((ELAPSED_TIME % 60))

echo "Job completed at: $(date)"
echo "Total elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED_TIME} seconds)"
