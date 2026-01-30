#!/bin/bash
# Generate QC summary PDF report
#
# Usage:
#   bash 4_generate_report.sh
#
# This script combines QC metrics and plots into a single PDF report.

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

QC_BASE="${PROJECT_ROOT}/datasets/test/results/qc"
OUTPUT="${QC_BASE}/qc_summary_report.pdf"
SCRIPT="${PROJECT_ROOT}/src/tf_perturb_seq/qc/generate_report.py"

###############################################
# RUN
###############################################

echo "============================================="
echo " Generate QC Summary Report"
echo " Output: ${OUTPUT}"
echo "============================================="

python "${SCRIPT}" \
    --output "${OUTPUT}" \
    --sample-name "Test Dataset" \
    --mapping-gene-dir "${QC_BASE}/mapping_gene" \
    --mapping-guide-dir "${QC_BASE}/mapping_guide" \
    --intended-target-dir "${QC_BASE}/intended_target"

echo "============================================="
echo " Done! Report: ${OUTPUT}"
echo "============================================="
