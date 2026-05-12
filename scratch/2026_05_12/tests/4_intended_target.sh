#!/bin/bash
# Run intended target inference QC
#
# Usage:
#   bash 3_intended_target.sh
#
# This script runs intended_target.py to evaluate guide knockdown efficiency
# on their intended target genes using trans-regulatory test results.
# Trans results are used (instead of cis) because they include all guide-gene
# pairs, including non-targeting guides needed for AUROC/AUPRC evaluation.

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
OUTDIR="${PROJECT_ROOT}/datasets/test/results/qc/intended_target"
SCRIPT="${PROJECT_ROOT}/src/tf_perturb_seq/qc/intended_target.py"

###############################################
# RUN
###############################################

echo "============================================="
echo " Intended Target Inference QC"
echo " Input:  ${INPUT}"
echo " Output: ${OUTDIR}"
echo "============================================="

mkdir -p "${OUTDIR}"

python "${SCRIPT}" \
    --input "${INPUT}" \
    --outdir "${OUTDIR}" \
    --prefix "intended_target"

echo "============================================="
echo " Done! Results in: ${OUTDIR}"
echo "============================================="
