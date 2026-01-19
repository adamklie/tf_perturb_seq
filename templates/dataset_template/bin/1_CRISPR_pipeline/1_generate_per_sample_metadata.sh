#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# CONFIGURATION - Only change BASE_DIR for your environment
# =============================================================================

# Base directory (change this for local vs remote)
BASE_DIR=/Users/adamklie/Desktop/projects/tf_perturb_seq

# Dataset name
DATASET_NAME=TEMPLATE_DATASET  # <-- CHANGE THIS

# Analysis set ID in IGVF portal
ACCESSION=IGVFDSXXXXXXXX  # <-- CHANGE THIS

# =============================================================================
# DERIVED PATHS (no need to change)
# =============================================================================

DATASET_DIR=${BASE_DIR}/datasets/${DATASET_NAME}
SCRIPT=${BASE_DIR}/scripts/generate_per_sample.py
SEQSPEC_DIR=${DATASET_DIR}/bin/1_CRISPR_pipeline/seqspec/yaml_files

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
  --output ${DATASET_DIR}/sample_metadata.csv \
  --hash_seqspec ${SEQSPEC_DIR}/hash_seqspec.yml \
  --rna_seqspec ${SEQSPEC_DIR}/rna_seqspec.yml \
  --sgrna_seqspec ${SEQSPEC_DIR}/guide_seqspec.yml

echo ""
echo "Done! Output: ${DATASET_DIR}/sample_metadata.csv"
