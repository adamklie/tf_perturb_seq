#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# CONFIGURATION - Update these for your dataset
# =============================================================================

# Base directory (change this for local vs remote)
BASE_DIR=/Users/adamklie/Desktop/projects/tf_perturb_seq

# Dataset name (should match folder name)
DATASET_NAME=TEMPLATE_DATASET  # <-- CHANGE THIS

# Analysis set ID from IGVF portal
ACCESSION=IGVFDSXXXXXXXX  # <-- CHANGE THIS

# =============================================================================
# DERIVED PATHS (no need to change)
# =============================================================================

DATASET_DIR=${BASE_DIR}/datasets/${DATASET_NAME}
SCRIPT=${BASE_DIR}/scripts/generate_per_sample.py

# =============================================================================
# RUN
# =============================================================================

echo "=========================================="
echo "Generate Per-Sample Metadata"
echo "=========================================="
echo "Dataset:    ${DATASET_NAME}"
echo "Accession:  ${ACCESSION}"
echo "Output:     ${DATASET_DIR}/sample_metadata.csv"
echo ""

python3 ${SCRIPT} \
  --accession ${ACCESSION} \
  --output ${DATASET_DIR}/sample_metadata.csv

echo ""
echo "Done! Output: ${DATASET_DIR}/sample_metadata.csv"
