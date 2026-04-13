#!/usr/bin/env bash
#SBATCH -J pathway_array
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH --partition=carter-compute
#SBATCH -t 01:00:00

set -euo pipefail

TSV="${1:?ERROR: provide TSV as arg1}"
PROJECT_ROOT="${2:?ERROR: provide project root as arg2}"

LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))
LINE="$(sed -n "${LINE_NUM}p" "${TSV}" || true)"
[[ -n "${LINE}" ]] || { echo "ERROR: No TSV line ${LINE_NUM}" >&2; exit 1; }

DATASET="$(echo "${LINE}" | cut -f1)"
RUN_LABEL="$(echo "${LINE}" | cut -f2)"
CAL_RESULTS="$(echo "${LINE}" | cut -f3)"
GENE_SETS="$(echo "${LINE}" | cut -f4)"
OUTDIR="$(echo "${LINE}" | cut -f5)"

echo "Dataset: ${DATASET}, Run: ${RUN_LABEL}"

bash "${PROJECT_ROOT}/scripts/run_pathways.sh" \
  --project-root "${PROJECT_ROOT}" \
  --calibrated-results "${CAL_RESULTS}" \
  --gene-sets "${GENE_SETS}" \
  --outdir "${OUTDIR}" \
  --prefix "${DATASET}"
