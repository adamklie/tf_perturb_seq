#!/usr/bin/env bash
set -euo pipefail

###############################################
# Huangfu Dataset: Run cNMF Gene Program Discovery
#
# USAGE (local):
#   ./6_run_cnmf.sh
#
# USAGE (dry-run):
#   ./6_run_cnmf.sh --dry-run
#
# USAGE (generate SLURM script):
#   ./6_run_cnmf.sh --slurm
#
# USAGE (run specific phase only):
#   ./6_run_cnmf.sh --phase 2
#
# All extra flags are forwarded to the
# orchestrator (run_cnmf_pipeline.sh).
###############################################

###############################################
# Paths - edit these as needed
###############################################

# Project root (where .venv lives)
PROJECT_ROOT="/Users/adamklie/Desktop/projects/tf_perturb_seq"

# Dataset root
DATASET_ROOT="${PROJECT_ROOT}/datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq"

# Input: inference MuData from CRISPR pipeline
# TODO: Update path once CRISPR pipeline run completes
INPUT="${DATASET_ROOT}/runs/triangular_squid/inference_mudata.h5mu"

# Output directory
OUTDIR="${DATASET_ROOT}/runs/triangular_squid/cnmf"

# Run name (used as prefix for output files)
RUN_NAME="Huangfu_WTC11-benchmark"

###############################################
# Run cNMF Pipeline
###############################################

echo "============================================="
echo " Huangfu Dataset: cNMF Gene Program Discovery"
echo " Input:  ${INPUT}"
echo " Output: ${OUTDIR}"
echo "============================================="

# Build command
CMD=(
  "${PROJECT_ROOT}/scripts/run_cnmf_pipeline.sh"
  --project-root "${PROJECT_ROOT}"
  --input "${INPUT}"
  --outdir "${OUTDIR}"
  --run-name "${RUN_NAME}"
  --use-gpu
)

# Forward all extra arguments (--dry-run, --slurm, --phase, etc.)
"${CMD[@]}" "$@"
