#!/usr/bin/env bash
#SBATCH -J inference_array
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --partition=carter-compute
#SBATCH -t 01:00:00
#SBATCH -o /cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/technology-benchmark_WTC11_TF-Perturb-seq/bin/slurm_logs/inference_array.%A_%a.out
#SBATCH -e /cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/technology-benchmark_WTC11_TF-Perturb-seq/bin/slurm_logs/inference_array.%A_%a.err

set -euo pipefail

# Runs calibration + pathway enrichment for one dataset per SLURM array task.
#
# Usage:
#   sbatch --array=1-N scripts/inference_array.sh /path/to/manifest.tsv /path/to/project_root [--dry-run]
#
# Manifest TSV format (tab-separated, with header):
#   dataset  run_label  trans_results  cis_results  mudata  gene_sets  outdir

TSV="${1:?ERROR: provide manifest TSV as arg1}"
PROJECT_ROOT="${2:?ERROR: provide project root as arg2}"

DRY_RUN_FLAG=""
if [[ "${3:-}" == "--dry-run" ]]; then
  DRY_RUN_FLAG="--dry-run"
elif [[ -n "${3:-}" ]]; then
  echo "ERROR: Unknown optional arg: ${3}" >&2
  exit 1
fi

SLURM_LOG_DIR="/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/technology-benchmark_WTC11_TF-Perturb-seq/bin/slurm_logs"
mkdir -p "${SLURM_LOG_DIR}"

[[ -f "${TSV}" ]] || { echo "ERROR: TSV not found: ${TSV}" >&2; exit 1; }

LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))   # skip header
LINE="$(sed -n "${LINE_NUM}p" "${TSV}" || true)"
[[ -n "${LINE}" ]] || { echo "ERROR: No TSV line ${LINE_NUM}" >&2; exit 1; }

DATASET="$(echo "${LINE}" | cut -f1)"
RUN_LABEL="$(echo "${LINE}" | cut -f2)"
TRANS_RESULTS="$(echo "${LINE}" | cut -f3)"
CIS_RESULTS="$(echo "${LINE}" | cut -f4)"
MUDATA="$(echo "${LINE}" | cut -f5)"
GENE_SETS="$(echo "${LINE}" | cut -f6)"
OUTDIR="$(echo "${LINE}" | cut -f7)"

CALIBRATE_SCRIPT="${PROJECT_ROOT}/scripts/run_calibration.sh"
PATHWAYS_SCRIPT="${PROJECT_ROOT}/scripts/run_pathways.sh"

echo "============================================="
echo " SLURM job:   ${SLURM_JOB_ID}"
echo " Task:        ${SLURM_ARRAY_TASK_ID}"
echo " Dataset:     ${DATASET}"
echo " Run label:   ${RUN_LABEL}"
echo " Outdir:      ${OUTDIR}"
echo " Dry-run:     ${DRY_RUN_FLAG:-false}"
echo "============================================="

# Step 1: Calibration
echo ""
echo "--- Step 1: Calibration ---"
bash "${CALIBRATE_SCRIPT}" \
  --project-root "${PROJECT_ROOT}" \
  --trans-results "${TRANS_RESULTS}" \
  --cis-results "${CIS_RESULTS}" \
  --mudata "${MUDATA}" \
  --outdir "${OUTDIR}" \
  --prefix "${DATASET}" \
  ${DRY_RUN_FLAG}

# Step 2: Pathway enrichment (if gene_sets provided and calibration produced output)
ALL_RESULTS="${OUTDIR}/${DATASET}_calibrated_all_results.tsv"
if [[ -f "${ALL_RESULTS}" && -n "${GENE_SETS}" && -f "${GENE_SETS}" ]]; then
  echo ""
  echo "--- Step 2: Pathway enrichment ---"
  bash "${PATHWAYS_SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --calibrated-results "${ALL_RESULTS}" \
    --gene-sets "${GENE_SETS}" \
    --outdir "${OUTDIR}" \
    --prefix "${DATASET}" \
    ${DRY_RUN_FLAG}
elif [[ -n "${DRY_RUN_FLAG}" ]]; then
  echo ""
  echo "[DRY-RUN] Step 2: Pathway enrichment would run after calibration"
else
  echo ""
  echo "[SKIP] Pathway enrichment — calibrated results or gene sets not available"
fi

echo ""
echo "Done."
