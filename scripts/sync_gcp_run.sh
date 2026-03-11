#!/usr/bin/env bash
set -euo pipefail

###############################################
# sync_gcp_run.sh
#
# Download a single dataset's outputs from a GCP CRISPR pipeline run to a
# local runs/<label>/ directory.
#
# By default downloads the full pipeline_dashboard/ (preferred — includes
# pre-computed QC, plots, h5mu, and dashboard.html).  Falls back to
# pipeline_outputs/ when no dashboard is present (e.g. a partially-complete
# run).
#
# USAGE:
#   sync_gcp_run.sh \
#     --gcp-dataset-uri gs://bucket/scratch/.../RunName/DatasetName \
#     --local-dir /path/to/datasets/DatasetName/runs/run-label \
#     [--mode dashboard|outputs|auto]  # default: auto
#     [--force] [--dry-run]
#
# Downloaded layout in local-dir:
#
#   pipeline_dashboard/             (mode=dashboard or auto with dashboard)
#     dashboard.html
#     inference_mudata.h5mu
#     additional_qc/
#       gene/gene_metrics.tsv
#       guide/guide_metrics.tsv
#       intended_target/intended_target_{results,metrics}.tsv
#       trans/trans_{results,metrics,per_guide_summary,significant_results}.tsv
#     benchmark_output/...
#     evaluation_output/...
#     figures/...
#     ...
#
#   inference_mudata.h5mu           (mode=outputs or auto fallback)
#   perturbo_*.tsv.gz               (mode=outputs or auto fallback)
###############################################

GCP_DATASET_URI=""
LOCAL_DIR=""
MODE="auto"     # auto | dashboard | outputs
FORCE=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gcp-dataset-uri) GCP_DATASET_URI="$2"; shift 2 ;;
    --local-dir)       LOCAL_DIR="$2";       shift 2 ;;
    --mode)            MODE="$2";            shift 2 ;;
    --force)           FORCE=1;              shift 1 ;;
    --dry-run)         DRY_RUN=1;            shift 1 ;;
    -h|--help)
      sed -n '1,45p' "$0"
      exit 0
      ;;
    *)
      echo "ERROR: Unknown arg: $1" >&2
      exit 1
      ;;
  esac
done

[[ -n "${GCP_DATASET_URI}" && -n "${LOCAL_DIR}" ]] || {
  echo "ERROR: --gcp-dataset-uri and --local-dir are required." >&2
  exit 1
}

[[ "${MODE}" == "auto" || "${MODE}" == "dashboard" || "${MODE}" == "outputs" ]] || {
  echo "ERROR: --mode must be auto, dashboard, or outputs (got: ${MODE})" >&2
  exit 1
}

GCP_DATASET_URI="${GCP_DATASET_URI%/}"
LOCAL_DIR="${LOCAL_DIR%/}"

GCP_DASHBOARD="${GCP_DATASET_URI}/pipeline_dashboard"
GCP_OUTPUTS="${GCP_DATASET_URI}/pipeline_outputs"

###############################################
# Resolve mode
###############################################
RESOLVED_MODE="${MODE}"

if [[ "${MODE}" == "auto" ]]; then
  if gsutil ls "${GCP_DASHBOARD}/" &>/dev/null; then
    RESOLVED_MODE="dashboard"
  elif gsutil ls "${GCP_OUTPUTS}/" &>/dev/null; then
    RESOLVED_MODE="outputs"
  else
    echo "ERROR: Neither pipeline_dashboard/ nor pipeline_outputs/ found under ${GCP_DATASET_URI}" >&2
    exit 1
  fi
fi

echo "============================================="
echo " sync_gcp_run.sh"
echo " GCP dataset:  ${GCP_DATASET_URI}"
echo " Local dest:   ${LOCAL_DIR}"
echo " Mode:         ${RESOLVED_MODE}"
echo " Force:        ${FORCE}"
echo " Dry-run:      ${DRY_RUN}"
echo "============================================="

if [[ "${DRY_RUN}" -eq 0 ]]; then
  mkdir -p "${LOCAL_DIR}"
fi

###############################################
# Mode: dashboard — rsync the full directory
###############################################
if [[ "${RESOLVED_MODE}" == "dashboard" ]]; then
  LOCAL_DASHBOARD="${LOCAL_DIR}/pipeline_dashboard"

  if [[ "${FORCE}" -eq 0 && -f "${LOCAL_DASHBOARD}/inference_mudata.h5mu" ]]; then
    echo "[SKIP] pipeline_dashboard/ already present (use --force to re-sync)"
    exit 0
  fi

  if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "[DRY-RUN] gsutil -m rsync -r ${GCP_DASHBOARD} ${LOCAL_DASHBOARD}"
    exit 0
  fi

  echo "[SYNC] pipeline_dashboard/ → ${LOCAL_DASHBOARD}"
  mkdir -p "${LOCAL_DASHBOARD}"
  gsutil -m rsync -r "${GCP_DASHBOARD}" "${LOCAL_DASHBOARD}"

  echo "============================================="
  echo " Sync complete: ${LOCAL_DASHBOARD}"
  echo "============================================="
  exit 0
fi

###############################################
# Mode: outputs — download individual files
###############################################
REQUIRED_FILES=(
  "inference_mudata.h5mu"
)
OPTIONAL_FILES=(
  "perturbo_cis_per_element_output.tsv.gz"
  "perturbo_cis_per_guide_output.tsv.gz"
  "perturbo_trans_per_element_output.tsv.gz"
  "perturbo_trans_per_guide_output.tsv.gz"
)

DOWNLOADED=0
SKIPPED=0
WARNINGS=0

for f in "${REQUIRED_FILES[@]}"; do
  src="${GCP_OUTPUTS}/${f}"
  dst="${LOCAL_DIR}/${f}"

  if [[ "${FORCE}" -eq 0 && -f "${dst}" ]]; then
    echo "[SKIP]     ${f}  (already exists)"
    SKIPPED=$((SKIPPED + 1))
    continue
  fi

  if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "[DRY-RUN]  gsutil cp ${src} ${dst}"
    continue
  fi

  echo "[DOWNLOAD] ${f}"
  gsutil cp "${src}" "${dst}"
  DOWNLOADED=$((DOWNLOADED + 1))
done

for f in "${OPTIONAL_FILES[@]}"; do
  src="${GCP_OUTPUTS}/${f}"
  dst="${LOCAL_DIR}/${f}"

  if [[ "${FORCE}" -eq 0 && -f "${dst}" ]]; then
    echo "[SKIP]     ${f}  (already exists)"
    SKIPPED=$((SKIPPED + 1))
    continue
  fi

  if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "[DRY-RUN]  gsutil cp ${src} ${dst}  (optional)"
    continue
  fi

  if gsutil ls "${src}" &>/dev/null; then
    echo "[DOWNLOAD] ${f}  (optional)"
    gsutil cp "${src}" "${dst}"
    DOWNLOADED=$((DOWNLOADED + 1))
  else
    echo "[WARN]     ${f} not on GCP — skipping (optional)"
    WARNINGS=$((WARNINGS + 1))
  fi
done

echo "============================================="
echo " Downloaded: ${DOWNLOADED}"
echo " Skipped:    ${SKIPPED}"
echo " Warnings:   ${WARNINGS}  (optional files missing on GCP)"
echo "============================================="
