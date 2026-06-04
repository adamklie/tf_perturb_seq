#!/usr/bin/env bash
set -euo pipefail

###############################################
# Gersbach HTv2 Dataset: Run cNMF Gene Program Discovery
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
DATASET_ROOT="${PROJECT_ROOT}/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2"

# Input: inference MuData from CRISPR pipeline
INPUT="${DATASET_ROOT}/runs/red_lettuce/inference_mudata.h5mu"

# Output directory
OUTDIR="${DATASET_ROOT}/runs/red_lettuce/cnmf"

# Run name (used as prefix for output files)
RUN_NAME="Gersbach_HTv2"

###############################################
# Run cNMF Pipeline
###############################################

echo "============================================="
echo " Gersbach HTv2: cNMF Gene Program Discovery"
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
