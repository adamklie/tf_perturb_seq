#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# CONFIGURATION - Only change BASE_DIR for your environment
# =============================================================================

# Base directory (change this for local vs remote)
BASE_DIR=/Users/adamklie/Desktop/projects/tf_perturb_seq

# Dataset name
DATASET_NAME=Engreitz_WTC11-benchmark_TF-Perturb-seq

# =============================================================================
# DERIVED PATHS (no need to change)
# =============================================================================

DATASET_DIR=${BASE_DIR}/datasets/${DATASET_NAME}
PATCH_SCRIPT=${BASE_DIR}/scripts/patch_gcp_files.py

# Input: the GCP samplesheet with .gz paths
INPUT_FILE=${DATASET_DIR}/sample_metadata_gcp_2026_02_26.csv

# Output: patched samplesheet with decompressed paths
OUTPUT_FILE=${DATASET_DIR}/sample_metadata_gcp_2026_02_26_patched.csv

# =============================================================================
# RUN
# =============================================================================

echo "=========================================="
echo "Patch GCP Files (decompress .gz)"
echo "=========================================="
echo "Dataset:    ${DATASET_NAME}"
echo "Input:      ${INPUT_FILE}"
echo "Output:     ${OUTPUT_FILE}"
echo ""

CMD=(
    python3 "${PATCH_SCRIPT}"
    --input "${INPUT_FILE}"
    --output "${OUTPUT_FILE}"
)

# Add dry-run flag if DRY_RUN is set
if [[ "${DRY_RUN:-false}" == "true" ]]; then
    CMD+=(--dry-run)
    echo "DRY RUN MODE - No files will be decompressed"
    echo ""
fi

"${CMD[@]}"

echo ""
echo "Done! Patched samplesheet: ${OUTPUT_FILE}"
