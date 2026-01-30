#!/bin/bash
# Run trans-regulatory inference QC
#
# Usage:
#   bash 4_trans.sh
#
# This script runs trans.py to evaluate trans-regulatory effects across the
# genome, measuring how many genes each perturbation significantly affects
# beyond its direct target.
#
# Optionally evaluates recovery of validated trans regulatory relationships
# using AUROC/AUPRC if a validated links file is provided.

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
OUTDIR="${PROJECT_ROOT}/datasets/test/results/qc/trans"
SCRIPT="${PROJECT_ROOT}/src/tf_perturb_seq/qc/trans.py"

# Optional: Path to validated trans links TSV (uncomment to use)
# VALIDATED_LINKS="${PROJECT_ROOT}/datasets/test/validated_trans_links.tsv"

###############################################
# RUN
###############################################

echo "============================================="
echo " Trans-Regulatory Inference QC"
echo " Input:  ${INPUT}"
echo " Output: ${OUTDIR}"
echo "============================================="

mkdir -p "${OUTDIR}"

# Run without validated links (general trans QC only)
python "${SCRIPT}" \
    --input "${INPUT}" \
    --outdir "${OUTDIR}" \
    --prefix "trans"

# Uncomment below to run with validated trans links for AUROC/AUPRC evaluation
# python "${SCRIPT}" \
#     --input "${INPUT}" \
#     --outdir "${OUTDIR}" \
#     --validated-links "${VALIDATED_LINKS}" \
#     --prefix "trans"

echo "============================================="
echo " Done! Results in: ${OUTDIR}"
echo "============================================="
