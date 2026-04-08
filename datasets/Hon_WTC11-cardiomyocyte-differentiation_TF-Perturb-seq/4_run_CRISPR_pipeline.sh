#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

# Dataset name (used for log file naming)
DATASET_NAME=Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq

# Base directory for the dataset (local path)
BASE_DIR=/carter/users/aklie/projects/tf_perturb_seq/datasets/${DATASET_NAME}

# Sample metadata with GCS paths
# TODO: Update to the correct patched metadata CSV after steps 2-3
SAMPLE_METADATA=$BASE_DIR/sample_metadata_gcp_YYYY_MM_DD_patched.csv

# CRISPR Pipeline path
# TODO: Update to the local clone of the CRISPR pipeline
PIPELINE_PATH=/path/to/CRISPR_Pipeline

# Output directory on GCS
# TODO: Update date and run name
OUTDIR=gs://igvf-pertub-seq-pipeline-data/${DATASET_NAME}/YYYY_MM_DD/outs/run_name

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
# NOTE: This is a production dataset (~26 measurement sets, full TF library)
#       Config needs to handle HTO multiplexing and larger scale
# TODO: Create dataset-specific .config
#       Start from Hon benchmark config:
#       datasets/Hon_WTC11-benchmark_TF-Perturb-seq/Hon_WTC11-benchmark_TF-Perturb-seq_2026_03_11.config
#       Key differences:
#         - ENABLE_DATA_HASHING = true (HTO multiplexed)
#         - Potentially higher resource limits for the larger dataset
#         - guide_design file is IGVFFI8270UPKB (production library, same as Huangfu production)
NF_CMD="nextflow run main.nf \
    -profile google \
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
