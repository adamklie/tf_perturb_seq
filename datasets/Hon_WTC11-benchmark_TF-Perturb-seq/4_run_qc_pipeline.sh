#!/usr/bin/env bash
set -euo pipefail

###############################################
# Hon Dataset: Run QC Pipeline
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
DATASET_ROOT="${PROJECT_ROOT}/datasets/Hon_WTC11-benchmark_TF-Perturb-seq"

# Input: inference MuData from CRISPR pipeline
INPUT="${DATASET_ROOT}/results/1_CRISPR_pipeline/2026_01_22/inference_mudata.h5mu"

# Output directory
OUTDIR="${DATASET_ROOT}/results/2_qc_pipeline/2026_01_22"

# Run name (used as prefix for output files)
RUN_NAME="Hon_WTC11-benchmark"

# Optional: validated trans links file (uncomment if available)
# VALIDATED_TRANS_LINKS="${DATASET_ROOT}/resources/validated_trans_links.tsv"

###############################################
# Run QC pipeline
###############################################

echo "============================================="
echo " Hon Dataset: QC Pipeline"
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
