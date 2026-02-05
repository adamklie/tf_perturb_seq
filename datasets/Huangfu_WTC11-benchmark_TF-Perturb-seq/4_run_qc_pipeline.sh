#!/usr/bin/env bash
set -euo pipefail

###############################################
# Huangfu Dataset: Run QC Pipeline
#
# USAGE (local):
#   ./4_run_qc_pipeline.sh
#
# USAGE (dry-run):
#   ./4_run_qc_pipeline.sh --dry-run
###############################################

# Parse args
DRY_RUN_FLAG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY_RUN_FLAG="--dry-run"; shift ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

###############################################
# Paths - edit these as needed
###############################################

# Project root (where .venv lives)
PROJECT_ROOT="/Users/adamklie/Desktop/projects/tf_perturb_seq"

# Dataset root
DATASET_ROOT="${PROJECT_ROOT}/datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq"

# Input: inference MuData from CRISPR pipeline
INPUT="${DATASET_ROOT}/runs/triangular_squid/pipeline_dashboard/inference_mudata.h5mu"

# Output directory
OUTDIR="${DATASET_ROOT}/runs/triangular_squid/pipeline_qc"

# Run name (used as prefix for output files)
RUN_NAME="Huangfu_WTC11-benchmark_TF-Perturb-seq"

# Optional: validated trans links file (uncomment if available)
# VALIDATED_TRANS_LINKS="${DATASET_ROOT}/resources/validated_trans_links.tsv"

###############################################
# Run QC pipeline
###############################################

echo "============================================="
echo " Huangfu Dataset: QC Pipeline"
echo " Input:  ${INPUT}"
echo " Output: ${OUTDIR}"
echo "============================================="

# Build command
CMD=(
  "${PROJECT_ROOT}/scripts/run_qc_pipeline.sh"
  --project-root "${PROJECT_ROOT}"
  --input "${INPUT}"
  --outdir "${OUTDIR}"
  --run-name "${RUN_NAME}"
)

# Add validated trans links if defined
if [[ -n "${VALIDATED_TRANS_LINKS:-}" ]]; then
  CMD+=(--validated-trans-links "${VALIDATED_TRANS_LINKS}")
fi

# Add dry-run flag if requested
if [[ -n "${DRY_RUN_FLAG}" ]]; then
  CMD+=(${DRY_RUN_FLAG})
fi

# Execute
"${CMD[@]}"
