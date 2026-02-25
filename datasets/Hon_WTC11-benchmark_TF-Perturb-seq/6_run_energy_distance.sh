#!/usr/bin/env bash
set -euo pipefail

###############################################
# Hon Dataset: Run Energy Distance Pipeline
#
# USAGE (local):
#   ./5_run_energy_distance.sh
#
# USAGE (dry-run):
#   ./5_run_energy_distance.sh --dry-run
#
# USAGE (with Harmony batch correction):
#   ./5_run_energy_distance.sh --use-harmony
#
# USAGE (force re-run all steps):
#   ./5_run_energy_distance.sh --force
###############################################

# Parse args
DRY_RUN_FLAG=""
HARMONY_FLAG=""
FORCE_FLAG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY_RUN_FLAG="--dry-run"; shift ;;
    --use-harmony) HARMONY_FLAG="--use-harmony"; shift ;;
    --force) FORCE_FLAG="--force"; shift ;;
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
INPUT="${DATASET_ROOT}/runs/elated_almeida/inference_mudata_cleaned.h5mu"

# Output directory
OUTDIR="${DATASET_ROOT}/runs/elated_almeida/energy_distance"

# Run name (used as prefix for output files)
RUN_NAME="Hon_WTC11-benchmark"

# Batch key for Harmony integration (optional)
BATCH_KEY="batch"

###############################################
# Pipeline parameters
###############################################

# Number of highly variable genes
N_TOP_GENES=5000

# Number of PCA components
N_PCS=50

###############################################
# Run Energy Distance Pipeline
###############################################

echo "============================================="
echo " Hon Dataset: Energy Distance Pipeline"
echo " Input:  ${INPUT}"
echo " Output: ${OUTDIR}"
echo "============================================="

# Build command
CMD=(
  "${PROJECT_ROOT}/scripts/run_energy_distance_pipeline.sh"
  --project-root "${PROJECT_ROOT}"
  --input "${INPUT}"
  --outdir "${OUTDIR}"
  --run-name "${RUN_NAME}"
  --n-top-genes "${N_TOP_GENES}"
  --n-pcs "${N_PCS}"
  --batch-key "${BATCH_KEY}"
)

# Add Harmony flag if requested
if [[ -n "${HARMONY_FLAG}" ]]; then
  CMD+=(${HARMONY_FLAG})
fi

# Add force flag if requested
if [[ -n "${FORCE_FLAG}" ]]; then
  CMD+=(${FORCE_FLAG})
fi

# Add dry-run flag if requested
if [[ -n "${DRY_RUN_FLAG}" ]]; then
  CMD+=(${DRY_RUN_FLAG})
fi

# Execute
"${CMD[@]}"
