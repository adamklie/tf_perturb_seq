#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# ENVIRONMENT
# =============================================================================
# Nextflow lives in a conda env (not on default PATH).
# Activate before running. Same env Adam used for DE/ESC.
source /cellar/users/aklie/opt/miniconda3/etc/profile.d/conda.sh
conda activate nextflow

# =============================================================================
# CONFIGURATION
# =============================================================================

# Dataset name (used for log file naming)
DATASET_NAME=Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq

# Base directory for the dataset (local path)
BASE_DIR=/carter/users/aklie/projects/tf_perturb_seq/datasets/${DATASET_NAME}

# Sample metadata with GCS paths (patched, all 311 GCS paths verified 2026-04-29)
SAMPLE_METADATA=$BASE_DIR/sample_metadata_gcp_2026_04_15_patched.csv

# CRISPR Pipeline path
PIPELINE_PATH=/cellar/users/aklie/opt/CRISPR_Pipeline

# Dataset-specific config (adapted from Hon benchmark; HTO multiplexed production)
CONFIG=$BASE_DIR/${DATASET_NAME}_2026_04_15.config

# Output directory on GCS
# Name following the DE/ESC convention (sceptre_v1). Spacer-tagged Hon run #1.
OUTDIR=gs://igvf-pertub-seq-pipeline-data/${DATASET_NAME}/2026_04_15/outs/sceptre_v1

# Log file with dataset name and timestamp
LOG_FILE=$BASE_DIR/logs/${DATASET_NAME}_crispr_pipeline_$(date +%Y%m%d_%H%M%S).log

# Run in background? Set to true to run with nohup
RUN_IN_BACKGROUND=${RUN_IN_BACKGROUND:-false}

# =============================================================================
# RUN PIPELINE
# =============================================================================

# Create logs directory if it doesn't exist
mkdir -p $BASE_DIR/logs

cd $PIPELINE_PATH

# Build the nextflow command. Using absolute path to nextflow so background nohup
# subshell (which doesn't inherit conda activation) can still run it.
NEXTFLOW=/cellar/users/aklie/opt/miniconda3/envs/nextflow/bin/nextflow

NF_CMD="$NEXTFLOW run main.nf \
    -profile google \
    -c $CONFIG \
    --input $SAMPLE_METADATA \
    --outdir $OUTDIR \
    -with-tower \
    -resume"

echo "============================================="
echo " Hon cardio CRISPR pipeline launch"
echo "============================================="
echo " Dataset:           $DATASET_NAME"
echo " Sample metadata:   $SAMPLE_METADATA  (112 measurement sets)"
echo " Config:            $CONFIG"
echo " OUTDIR (GCS):      $OUTDIR"
echo " Pipeline path:     $PIPELINE_PATH"
echo " Nextflow:          $NEXTFLOW"
echo "============================================="

if [ "$RUN_IN_BACKGROUND" = true ]; then
    echo "Running CRISPR pipeline in background..."
    echo "Log file: $LOG_FILE"
    echo "Monitor with: tail -f $LOG_FILE"
    nohup bash -c "$NF_CMD" > $LOG_FILE 2>&1 &
    echo "Pipeline started with PID: $!"
else
    echo "Running CRISPR pipeline in foreground..."
    echo "To run in background, use: RUN_IN_BACKGROUND=true $0"
    $NF_CMD
fi
