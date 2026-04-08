#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# CONFIGURATION - Only change BASE_DIR for your environment
# =============================================================================

# Base directory (change this for local vs remote)
BASE_DIR=/carter/users/aklie/projects/tf_perturb_seq

# Dataset name
DATASET_NAME=Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq

# Analysis set ID in IGVF portal
# 78 inputs: 26 scRNA measurement sets + 52 gRNA auxiliary sets
ACCESSION=IGVFDS6332VCTO

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
