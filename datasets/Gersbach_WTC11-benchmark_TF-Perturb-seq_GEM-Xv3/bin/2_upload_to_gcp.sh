#!/usr/bin/env bash
#
# Upload files from IGVF portal and local paths to GCP
#
# This script wraps upload_to_gcp.py with dataset-specific configuration.
# It handles both IGVF accessions (transferred directly from S3 to GCS)
# and local files (uploaded via gsutil).
#
# Prerequisites:
#   - gcloud CLI authenticated and configured
#   - gsutil available
#   - Python 3.8+ with requests library
#   - Access to IGVF portal API
#
# Usage:
#   ./2_upload_to_gcp.sh              # Normal run
#   DRY_RUN=true ./2_upload_to_gcp.sh # Dry run
#
set -euo pipefail

# =============================================================================
# CONFIGURATION - Only change BASE_DIR for your environment
# =============================================================================

# Base directory (change this for local vs remote)
BASE_DIR=/Users/adamklie/Desktop/projects/tf_perturb_seq

# Dataset name
DATASET_NAME=Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3

# GCP project ID
PROJECT=igvf-pertub-seq-pipeline

# GCS bucket (without gs:// prefix)
GCS_BUCKET=igvf-pertub-seq-pipeline-data

# =============================================================================
# DERIVED PATHS (no need to change)
# =============================================================================

DATASET_DIR=${BASE_DIR}/datasets/${DATASET_NAME}
UPLOAD_SCRIPT=${BASE_DIR}/scripts/upload_to_gcp.py

# Input sample metadata file (CSV format)
INPUT_FILE=${DATASET_DIR}/sample_metadata.csv

# Output file with updated GCS paths
OUTPUT_FILE=${DATASET_DIR}/sample_metadata_gcp_$(date +%Y_%m_%d).csv

# GCS prefix/folder path within the bucket
GCS_PREFIX=${DATASET_NAME}/$(date +%Y_%m_%d)/

# Columns to process for file references
# Default: R1_path,R2_path,seqspec,barcode_onlist,guide_design,barcode_hashtag_map
# Uncomment to customize:
# COLUMNS="R1_path R2_path"

# =============================================================================
# EXECUTION
# =============================================================================

echo "=========================================="
echo "IGVF to GCP Upload Script"
echo "=========================================="
echo ""
echo "Configuration:"
echo "  Input:      ${INPUT_FILE}"
echo "  Output:     ${OUTPUT_FILE}"
echo "  GCS Dest:   gs://${GCS_BUCKET}/${GCS_PREFIX}"
echo "  Project:    ${PROJECT}"
echo ""

# Check prerequisites
if ! command -v gcloud &> /dev/null; then
    echo "Error: gcloud CLI not found. Please install Google Cloud SDK."
    exit 1
fi

if ! command -v gsutil &> /dev/null; then
    echo "Error: gsutil not found. Please install Google Cloud SDK."
    exit 1
fi

if ! python3 -c "import requests" &> /dev/null; then
    echo "Error: Python requests library not found."
    echo "Install with: pip install requests"
    exit 1
fi

# Authenticate if needed
echo "Checking GCP authentication..."
if ! gcloud auth print-access-token &> /dev/null; then
    echo "Not authenticated. Running gcloud auth login..."
    gcloud auth login
fi

# Set project
gcloud config set project "${PROJECT}"

# Build command
CMD=(
    python3 "${UPLOAD_SCRIPT}"
    --input "${INPUT_FILE}"
    --output "${OUTPUT_FILE}"
    --gcs-bucket "${GCS_BUCKET}"
    --gcs-prefix "${GCS_PREFIX}"
    --project "${PROJECT}"
)

# Add columns if specified
if [[ -n "${COLUMNS:-}" ]]; then
    CMD+=(--columns ${COLUMNS})
fi

# Add dry-run flag if DRY_RUN is set
if [[ "${DRY_RUN:-false}" == "true" ]]; then
    CMD+=(--dry-run)
    echo "DRY RUN MODE - No files will be uploaded"
fi

# Run the upload script
echo ""
echo "Running upload script..."
echo ""
"${CMD[@]}"

echo ""
echo "=========================================="
echo "Upload complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "  1. Review the output file: ${OUTPUT_FILE}"
echo "  2. Verify files in GCS: gsutil ls gs://${GCS_BUCKET}/${GCS_PREFIX}"
echo "  3. Use the updated metadata for the CRISPR pipeline"
