#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=CM_ED_SC_overlap
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/Testing_Scripts/Meta_program/Result/CM_ED_SC_analysis/logs/%j.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/Testing_Scripts/Meta_program/Result/CM_ED_SC_analysis/logs/%j.err
#SBATCH --partition=engreitz,owners,bigmem
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ymo@stanford.edu


LOG_DIR="/oak/stanford/groups/engreitz/Users/ymo/Testing_Scripts/Meta_program/Result/CM_ED_SC_analysis"
mkdir -p "$LOG_DIR/logs"

START_TIME=$(date +%s)

echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Partition: $SLURM_JOB_PARTITION"
echo "Log directory: $LOG_DIR"

# Activate conda env
eval "$(conda shell.bash hook)" && conda activate NMF_Benchmarking
echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python: $(which python)"

# Run analysis (unbuffered so progress prints stream into the .out file)
python -u /oak/stanford/groups/engreitz/Users/ymo/Testing_Scripts/Meta_program/scripts/CM_ED_SC_analysis/Script/CM_ED_SC_overlap.py

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo "Job completed at: $(date)"
echo "Total elapsed time: $((ELAPSED / 60))m $((ELAPSED % 60))s"
