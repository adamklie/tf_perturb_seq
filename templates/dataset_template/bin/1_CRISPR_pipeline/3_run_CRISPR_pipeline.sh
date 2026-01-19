#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

# Dataset name (used for log file naming)
DATASET_NAME=TEMPLATE_DATASET  # <-- CHANGE THIS

# Base directory for the dataset (local path)
BASE_DIR=/Users/adamklie/Desktop/projects/tf_perturb_seq/datasets/${DATASET_NAME}

# Sample metadata with GCS paths (update date as needed)
SAMPLE_METADATA=$BASE_DIR/sample_metadata_gcp_YYYY_MM_DD.csv  # <-- CHANGE DATE

# CRISPR Pipeline path
PIPELINE_PATH=/Users/adamklie/Desktop/projects/CRISPR_Pipeline

# Output directory on GCS (update date as needed)
OUTDIR=gs://igvf-pertub-seq-pipeline-data/${DATASET_NAME}/YYYY_MM_DD/outs  # <-- CHANGE DATE

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

# Build the nextflow command
# -profile google: Use Google Cloud Batch for execution
# --input: Sample metadata CSV with GCS paths
# --outdir: Output directory on GCS
# -with-tower: Enable Nextflow Tower monitoring (requires TOWER_ACCESS_TOKEN)
NF_CMD="nextflow run main.nf \
    -profile google \
    --input $SAMPLE_METADATA \
    --outdir $OUTDIR \
    -with-tower"

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
