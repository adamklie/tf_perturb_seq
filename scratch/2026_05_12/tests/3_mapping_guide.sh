#!/bin/bash
# Run guide mapping QC
#
# Usage:
#   bash 2_mapping_guide.sh
#
# This script runs mapping_guide.py to compute guide assignment QC metrics
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
OUTDIR="${PROJECT_ROOT}/datasets/test/results/qc/mapping_guide"
SCRIPT="${PROJECT_ROOT}/src/tf_perturb_seq/qc/mapping_guide.py"

###############################################
# RUN
###############################################

echo "============================================="
echo " Guide Mapping QC"
echo " Input:  ${INPUT}"
echo " Output: ${OUTDIR}"
echo "============================================="

mkdir -p "${OUTDIR}"

python "${SCRIPT}" \
    --input "${INPUT}" \
    --outdir "${OUTDIR}" \
    --prefix "guide"

echo "============================================="
echo " Done! Results in: ${OUTDIR}"
echo "============================================="
