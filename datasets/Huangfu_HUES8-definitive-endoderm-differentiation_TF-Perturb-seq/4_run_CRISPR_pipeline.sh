#!/usr/bin/env bash
set -euo pipefail

# Pin Nextflow version (Nextflow downloads this at startup)
export NXF_VER=26.03.2-edge

# GCS auth for the pipeline (service account key for igvf-pertub-seq-pipeline)
export GOOGLE_APPLICATION_CREDENTIALS=/Users/adamklie/Desktop/tfp3/igvf-pertub-seq-pipeline-aa37ac78d6f1.json

# =============================================================================
# CONFIGURATION
# =============================================================================

# Dataset name (used for log file naming)
DATASET_NAME=Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq

# Base directory for the dataset (local path)
BASE_DIR=/Users/adamklie/Desktop/tfp3/tf_perturb_seq/datasets/${DATASET_NAME}

# Data date (GCS subfolder under the dataset — matches sample_metadata_gcp_<DATA_DATE>.csv)
DATA_DATE=2026_04_09

# Sample metadata with GCS paths (patched = decompressed barcode_onlist, guide_design, seqspec)
SAMPLE_METADATA=$BASE_DIR/sample_metadata_gcp_${DATA_DATE}_patched.csv

# CRISPR Pipeline path
PIPELINE_PATH=/Users/adamklie/Desktop/tfp3/CRISPR_Pipeline

# Run label (Nextflow run name) — bump this per run
RUN_LABEL=muddy_penguin

# Dataset-specific config (adapted from Huangfu WTC11 benchmark)
CONFIG=$BASE_DIR/${DATASET_NAME}_${RUN_LABEL}.config

# Output directory on GCS
OUTDIR=gs://igvf-pertub-seq-pipeline-data/${DATASET_NAME}/${DATA_DATE}/outs/${RUN_LABEL}

# Log file with dataset name, run label, and timestamp
LOG_FILE=$BASE_DIR/logs/${DATASET_NAME}_${RUN_LABEL}_$(date +%Y%m%d_%H%M%S).log

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
