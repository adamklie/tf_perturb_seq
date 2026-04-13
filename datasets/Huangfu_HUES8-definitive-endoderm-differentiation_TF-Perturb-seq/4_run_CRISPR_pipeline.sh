#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

# Dataset name (used for log file naming)
DATASET_NAME=Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq

# Base directory for the dataset (local path)
BASE_DIR=/carter/users/aklie/projects/tf_perturb_seq/datasets/${DATASET_NAME}

# Sample metadata with GCS paths (patched = decompressed barcode_onlist, guide_design, seqspec)
SAMPLE_METADATA=$BASE_DIR/sample_metadata_gcp_2026_04_09_patched.csv

# CRISPR Pipeline path
PIPELINE_PATH=/cellar/users/aklie/opt/CRISPR_Pipeline

# Dataset-specific config (adapted from Huangfu WTC11 benchmark)
CONFIG=$BASE_DIR/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq_2026_04_09.config

# Output directory on GCS
OUTDIR=gs://igvf-pertub-seq-pipeline-data/${DATASET_NAME}/2026_04_09/outs/sceptre_v1

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
