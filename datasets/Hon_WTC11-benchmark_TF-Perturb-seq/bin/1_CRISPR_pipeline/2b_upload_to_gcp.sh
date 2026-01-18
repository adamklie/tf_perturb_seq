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
#   ./2b_upload_to_gcp.sh
#
set -euo pipefail

# =============================================================================
# CONFIGURATION - Modify these values as needed
# =============================================================================

# Base directory for the dataset
BASE_DIR=/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq

# Input sample metadata file
INPUT_FILE=${BASE_DIR}/sample_metadata.tsv

# Output file with updated GCS paths
# Tip: Use a date suffix to track versions
OUTPUT_FILE=${BASE_DIR}/sample_metadata_gcp_$(date +%Y_%m_%d).tsv

# GCP project ID
PROJECT=igvf-pertub-seq-pipeline

# GCS bucket (without gs:// prefix)
GCS_BUCKET=igvf-pertub-seq-pipeline-data

# GCS prefix/folder path within the bucket
# Tip: Include a date to organize uploads
GCS_PREFIX=Hon_WTC11-benchmark_TF-Perturb-seq/$(date +%Y_%m_%d)/

# Columns to process for file references
# Default: R1_path,R2_path,seqspec,barcode_onlist,guide_design,barcode_hashtag_map
# Uncomment to customize:
# COLUMNS="R1_path R2_path"

# Path to the upload script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UPLOAD_SCRIPT="${SCRIPT_DIR}/../upload_to_gcp.py"

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
