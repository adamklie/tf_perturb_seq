#!/usr/bin/env bash
#SBATCH -J qc_array
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --partition=carter-compute
#SBATCH -t 02:00:00
#SBATCH -o slurm_logs/qc_array.%A_%a.out
#SBATCH -e slurm_logs/qc_array.%A_%a.err

set -euo pipefail

# Usage:
#   sbatch --array=1-N qc_array.sbatch /path/to/samples.tsv /path/to/project_root [--dry-run]

TSV="${1:?ERROR: provide TSV path as arg1}"
PROJECT_ROOT="${2:?ERROR: provide project root as arg2}"

DRY_RUN_FLAG=""
if [[ "${3:-}" == "--dry-run" ]]; then
  DRY_RUN_FLAG="--dry-run"
elif [[ -n "${3:-}" ]]; then
  echo "ERROR: Unknown optional arg: ${3}" >&2
  exit 1
fi

mkdir -p slurm_logs

[[ -f "${TSV}" ]] || { echo "ERROR: TSV not found: ${TSV}" >&2; exit 1; }

LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))   # skip header
LINE="$(sed -n "${LINE_NUM}p" "${TSV}" || true)"
[[ -n "${LINE}" ]] || { echo "ERROR: No TSV line ${LINE_NUM}" >&2; exit 1; }

INPUT="$(echo "${LINE}" | cut -f1)"
OUTDIR="$(echo "${LINE}" | cut -f2)"
RUN_NAME="$(echo "${LINE}" | cut -f3)"

[[ -n "${INPUT}" && -n "${OUTDIR}" && -n "${RUN_NAME}" ]] || {
  echo "ERROR: Failed to parse TSV line ${LINE_NUM}: ${LINE}" >&2
  exit 1
}

PIPELINE_SCRIPT="${PROJECT_ROOT}/scripts/run_qc_pipeline.sh"
[[ -x "${PIPELINE_SCRIPT}" ]] || {
  echo "ERROR: Pipeline script not executable: ${PIPELINE_SCRIPT}" >&2
  exit 1
}

echo "============================================="
echo " SLURM job:  ${SLURM_JOB_ID}"
echo " Task:       ${SLURM_ARRAY_TASK_ID}"
echo " Input:      ${INPUT}"
echo " Outdir:     ${OUTDIR}"
echo " Run name:   ${RUN_NAME}"
echo " Dry-run:    ${DRY_RUN_FLAG:-false}"
echo "============================================="

bash "${PIPELINE_SCRIPT}" \
  --project-root "${PROJECT_ROOT}" \
  --input "${INPUT}" \
  --outdir "${OUTDIR}" \
  --run-name "${RUN_NAME}" \
  ${DRY_RUN_FLAG}
