#!/usr/bin/env bash
#SBATCH -J seqspec_val
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH --partition=carter-compute
#SBATCH -t 00:30:00
#SBATCH -o /carter/users/aklie/projects/tf_perturb_seq/scratch/seqspec_validation_logs/seqspec_val.%A_%a.out
#SBATCH -e /carter/users/aklie/projects/tf_perturb_seq/scratch/seqspec_validation_logs/seqspec_val.%A_%a.err

# DRAFT — NOT yet submitted. One array task per dataset; each runs the crispr_validator
# (Mode 1 samplesheet) against that dataset's tiny representative-lane CSV.
#
# NRNB note: `sbatch`/`sinfo` are not on PATH on the login node by default — run
#   `module load slurm/nrnb/23.02.7`
# before submitting. carter-compute is the project partition (confirmed via `sinfo -s`).
#
# Usage:
#   N=$(($(wc -l < manifest.tsv) - 1))
#   sbatch --array=1-${N} submit_seqspec_validation.sh manifest.tsv /carter/users/aklie/projects/tf_perturb_seq [--dry-run]
#
# manifest.tsv columns (tab, header row): dataset  validation_csv  analysis_root  downloads_dir
# (produced by build_validation_samplesheets.py --manifest ...)

set -euo pipefail

MANIFEST="${1:?ERROR: provide manifest TSV as arg1}"
REPO_ROOT="${2:?ERROR: provide repo root as arg2}"
DRY_RUN=0
if [[ "${3:-}" == "--dry-run" ]]; then
  DRY_RUN=1
elif [[ -n "${3:-}" ]]; then
  echo "ERROR: unknown optional arg: ${3}" >&2; exit 1
fi

# Validator + tool paths (verified present on NRNB).
VALIDATOR_DIR="${REPO_ROOT}/external/crispr_validator"
KEYPAIR="${IGVF_KEYPAIR:-${REPO_ROOT}/config/igvf_keypair.json}"   # adjust to taste
EXTRA_PATH="/cellar/users/aklie/.local/bin:/carter/users/aklie/opt/google-cloud-sdk/bin"
BARCODE_SAMPLE_READS="${BARCODE_SAMPLE_READS:-3000}"
FEATURE_SAMPLE_READS="${FEATURE_SAMPLE_READS:-5000}"
CHUNK_BYTES="${CHUNK_BYTES:-8000000}"

LOG_DIR="${REPO_ROOT}/scratch/seqspec_validation_logs"
mkdir -p "${LOG_DIR}"

[[ -f "${MANIFEST}" ]] || { echo "ERROR: manifest not found: ${MANIFEST}" >&2; exit 1; }
[[ -d "${VALIDATOR_DIR}" ]] || { echo "ERROR: validator dir not found: ${VALIDATOR_DIR}" >&2; exit 1; }

LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))   # skip header
LINE="$(sed -n "${LINE_NUM}p" "${MANIFEST}" || true)"
[[ -n "${LINE}" ]] || { echo "ERROR: no manifest line ${LINE_NUM}" >&2; exit 1; }

DATASET="$(echo "${LINE}" | cut -f1)"
VALIDATION_CSV="$(echo "${LINE}" | cut -f2)"
ANALYSIS_ROOT="$(echo "${LINE}" | cut -f3)"
DOWNLOADS_DIR="$(echo "${LINE}" | cut -f4)"

[[ -n "${DATASET}" && -n "${VALIDATION_CSV}" && -n "${ANALYSIS_ROOT}" && -n "${DOWNLOADS_DIR}" ]] || {
  echo "ERROR: failed to parse manifest line ${LINE_NUM}: ${LINE}" >&2; exit 1; }

echo "============================================="
echo " SLURM job:      ${SLURM_JOB_ID:-NA}"
echo " Task:           ${SLURM_ARRAY_TASK_ID:-NA}"
echo " Dataset:        ${DATASET}"
echo " Validation CSV: ${VALIDATION_CSV}"
echo " Analysis root:  ${ANALYSIS_ROOT}"
echo " Downloads dir:  ${DOWNLOADS_DIR}"
echo " Dry-run:        ${DRY_RUN}"
echo "============================================="

mkdir -p "${ANALYSIS_ROOT}" "${DOWNLOADS_DIR}"

# Build the verbatim working invocation.
run_cmd() {
  cd "${VALIDATOR_DIR}"
  PATH="${EXTRA_PATH}:${PATH}" \
  uv run --no-project --with pyyaml --with certifi --with seqspec \
    python -u seqspec_parser.py samplesheet \
    --samplesheet "${VALIDATION_CSV}" \
    --analysis-root "${ANALYSIS_ROOT}" \
    --downloads-dir "${DOWNLOADS_DIR}" \
    --group-by measurement_sets \
    --igvf-keypair "${KEYPAIR}" \
    --barcode-sample-reads "${BARCODE_SAMPLE_READS}" \
    --feature-sample-reads "${FEATURE_SAMPLE_READS}" \
    --chunk-bytes "${CHUNK_BYTES}"
}

if [[ "${DRY_RUN}" -eq 1 ]]; then
  echo "[dry-run] would run validator for ${DATASET}:"
  echo "  cd ${VALIDATOR_DIR} && PATH=${EXTRA_PATH}:\$PATH uv run --no-project \\"
  echo "    --with pyyaml --with certifi --with seqspec python -u seqspec_parser.py samplesheet \\"
  echo "    --samplesheet ${VALIDATION_CSV} --analysis-root ${ANALYSIS_ROOT} \\"
  echo "    --downloads-dir ${DOWNLOADS_DIR} --group-by measurement_sets \\"
  echo "    --igvf-keypair ${KEYPAIR} --barcode-sample-reads ${BARCODE_SAMPLE_READS} \\"
  echo "    --feature-sample-reads ${FEATURE_SAMPLE_READS} --chunk-bytes ${CHUNK_BYTES}"
  exit 0
fi

run_cmd
echo "DONE ${DATASET}: see ${ANALYSIS_ROOT}/<group>/analysis_summary.json"
