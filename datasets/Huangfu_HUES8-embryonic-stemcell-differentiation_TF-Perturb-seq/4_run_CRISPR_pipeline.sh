#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

# Dataset name (used for log file naming)
DATASET_NAME=Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq

# Base directory for the dataset (local path)
BASE_DIR=/data4/yyang117/tf_perturb_seq/datasets/${DATASET_NAME}

# Sample metadata with GCS paths
SAMPLE_METADATA=$BASE_DIR/sample_metadata_gcp_2026_04_13_patched.csv

# CRISPR Pipeline path
# TODO: Update to the local clone of the CRISPR pipeline
PIPELINE_PATH=/data4/yyang117/CRISPR_Pipeline

# Dataset-specific config (adapted from Huangfu WTC11 benchmark)
CONFIG=$BASE_DIR/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq_2026_04_13.config

# Output directory on GCS
OUTDIR=gs://igvf-pertub-seq-pipeline-data/${DATASET_NAME}/2026_04_13/outs/sceptre_v1

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
# Uses the Huangfu benchmark config as a starting point
# TODO: Create a dataset-specific .config if parameters differ
#       (e.g. different guide library, MOI, capture method)
#       The Huangfu benchmark config is at:
#       datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/Huangfu_WTC11-benchmark_TF-Perturb-seq_2026_03_11.config
NF_CMD="nextflow run main.nf \
    -profile google \
    -c $CONFIG \
    --input $SAMPLE_METADATA \
    --outdir $OUTDIR \
    -with-tower \
    -resume"

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
