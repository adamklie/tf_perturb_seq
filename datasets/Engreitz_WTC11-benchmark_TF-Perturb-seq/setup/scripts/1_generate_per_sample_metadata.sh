#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# CONFIGURATION - Only change BASE_DIR for your environment
# =============================================================================

# Base directory (change this for local vs remote)
BASE_DIR=/Users/adamklie/Desktop/projects/tf_perturb_seq

# Dataset name
DATASET_NAME=Engreitz_WTC11-benchmark_TF-Perturb-seq

# Analysis set ID in IGVF portal
ACCESSION=IGVFDS5057HJKP

# =============================================================================
# DERIVED PATHS (no need to change)
# =============================================================================

DATASET_DIR=${BASE_DIR}/datasets/${DATASET_NAME}
SCRIPT=${BASE_DIR}/src/tf_perturb_seq/portal/generate_per_sample.py

# =============================================================================
# RUN
# =============================================================================

echo "=========================================="
echo "Generate Per-Sample Metadata"
echo "=========================================="
echo "Dataset:    ${DATASET_NAME}"
echo "Accession:  ${ACCESSION}"
echo "Output:     ${DATASET_DIR}/cleanser_800_mito_15pc/crispr_pipeline/inputs/sample_metadata.csv"
echo ""

python3 ${SCRIPT} \
  --accession ${ACCESSION} \
  --output ${DATASET_DIR}/cleanser_800_mito_15pc/crispr_pipeline/inputs/sample_metadata.csv \
  --rna_seqspec /Users/adamklie/Downloads/tf_50genes/all_GEX_sublibraries/IGVFDS2550EQNC.yaml

echo ""
echo "Done! Output: ${DATASET_DIR}/cleanser_800_mito_15pc/crispr_pipeline/inputs/sample_metadata.csv"
