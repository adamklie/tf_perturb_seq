#!/bin/bash
# Run gene expression mapping QC
#
# Usage:
#   bash 1_mapping_gene.sh
#
# This script runs mapping_gene.py to compute gene expression QC metrics
# and generate distribution plots.

set -euo pipefail

###############################################
# PATHS
###############################################

# Project root (adjust if needed)
PROJECT_ROOT="/Users/adamklie/Desktop/projects/tf_perturb_seq"

###############################################
# ENVIRONMENT
###############################################

# Activate project venv (managed by uv)
source "${PROJECT_ROOT}/.venv/bin/activate"

###############################################
# INPUT/OUTPUT
###############################################

INPUT="${PROJECT_ROOT}/datasets/test/results/inference_mudata.h5mu"
OUTDIR="${PROJECT_ROOT}/datasets/test/results/qc/mapping_gene"
SCRIPT="${PROJECT_ROOT}/src/tf_perturb_seq/qc/mapping_gene.py"

###############################################
# RUN
###############################################

echo "============================================="
echo " Gene Expression Mapping QC"
echo " Input:  ${INPUT}"
echo " Output: ${OUTDIR}"
echo "============================================="

mkdir -p "${OUTDIR}"

python "${SCRIPT}" \
    --input "${INPUT}" \
    --outdir "${OUTDIR}" \
    --prefix "gene"

echo "============================================="
echo " Done! Results in: ${OUTDIR}"
echo "============================================="
