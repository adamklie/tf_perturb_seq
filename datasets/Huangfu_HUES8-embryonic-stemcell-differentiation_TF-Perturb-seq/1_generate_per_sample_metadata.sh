#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# CONFIGURATION - Only change BASE_DIR for your environment
# =============================================================================

# Base directory (change this for local vs remote)
BASE_DIR=/data4/yyang117/tf_perturb_seq

# Dataset name
DATASET_NAME=Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq

# Analysis set ID in IGVF portal
# TODO: Find the correct accession ID on https://data.igvf.org/
# The Huangfu benchmark accession is IGVFDS8556LYRW for reference.
# Search for HUES8 embryonic stem cell analysis sets on the portal.
ACCESSION=IGVFDS1216AEWT  # <-- CHANGE THIS

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
