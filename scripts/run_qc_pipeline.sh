#!/usr/bin/env bash
set -euo pipefail

###############################################
# Args
###############################################
PROJECT_ROOT=""
INPUT=""
OUTDIR=""
RUN_NAME=""
DRY_RUN=0

# NEW (optional)
VALIDATED_TRANS_LINKS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root)            PROJECT_ROOT="$2"; shift 2 ;;
    --input)                   INPUT="$2"; shift 2 ;;
    --outdir)                  OUTDIR="$2"; shift 2 ;;
    --run-name)                RUN_NAME="$2"; shift 2 ;;
    --validated-trans-links)   VALIDATED_TRANS_LINKS="$2"; shift 2 ;;  # optional
    --dry-run)                 DRY_RUN=1; shift 1 ;;
    -h|--help)
      sed -n '1,220p' "$0"
      exit 0
      ;;
    *)
      echo "ERROR: Unknown arg: $1" >&2
      exit 1
      ;;
  esac
done

[[ -n "${PROJECT_ROOT}" && -n "${INPUT}" && -n "${OUTDIR}" && -n "${RUN_NAME}" ]] || {
  echo "ERROR: Missing required args." >&2
  exit 1
}

[[ -f "${INPUT}" ]] || { echo "ERROR: Input does not exist: ${INPUT}" >&2; exit 1; }

if [[ -n "${VALIDATED_TRANS_LINKS}" && ! -f "${VALIDATED_TRANS_LINKS}" ]]; then
  echo "ERROR: --validated-trans-links provided but file not found: ${VALIDATED_TRANS_LINKS}" >&2
  exit 1
fi

###############################################
# Environment
###############################################
VENV="${PROJECT_ROOT}/.venv/bin/activate"
[[ -f "${VENV}" ]] || { echo "ERROR: venv not found: ${VENV}" >&2; exit 1; }
# shellcheck disable=SC1090
source "${VENV}"

###############################################
# Helpers
###############################################
run_cmd () {
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    printf '[DRY-RUN] %q' "$1"
    shift
    for arg in "$@"; do printf ' %q' "${arg}"; done
    printf '\n'
  else
    "$@"
  fi
}

###############################################
# Scripts
###############################################
SCRIPT_GENE="${PROJECT_ROOT}/src/tf_perturb_seq/qc/mapping_gene.py"
SCRIPT_GUIDE="${PROJECT_ROOT}/src/tf_perturb_seq/qc/mapping_guide.py"
SCRIPT_INTENDED="${PROJECT_ROOT}/src/tf_perturb_seq/qc/intended_target.py"
SCRIPT_TRANS="${PROJECT_ROOT}/src/tf_perturb_seq/qc/trans.py"               # NEW
SCRIPT_REPORT="${PROJECT_ROOT}/src/tf_perturb_seq/qc/generate_report.py"

for s in "${SCRIPT_GENE}" "${SCRIPT_GUIDE}" "${SCRIPT_INTENDED}" "${SCRIPT_TRANS}" "${SCRIPT_REPORT}"; do
  [[ -f "${s}" ]] || { echo "ERROR: Missing script: ${s}" >&2; exit 1; }
done

###############################################
# Output layout
###############################################
DIR_GENE="${OUTDIR}/mapping_gene"
DIR_GUIDE="${OUTDIR}/mapping_guide"
DIR_INTENDED="${OUTDIR}/intended_target"
DIR_TRANS="${OUTDIR}/trans"                              # NEW
REPORT_PDF="${OUTDIR}/qc_summary_report.pdf"

###############################################
# Prefixes
###############################################
PREFIX_GENE="${RUN_NAME}_gene"
PREFIX_GUIDE="${RUN_NAME}_guide"
PREFIX_INTENDED="${RUN_NAME}_intended_target"
PREFIX_TRANS="${RUN_NAME}_trans"                          # NEW
PREFIX_REPORT="${RUN_NAME}"

echo "============================================="
echo " QC pipeline"
echo " Project:    ${PROJECT_ROOT}"
echo " Input:      ${INPUT}"
echo " Outdir:     ${OUTDIR}"
echo " Run name:   ${RUN_NAME}"
echo " Dry-run:    ${DRY_RUN}"
if [[ -n "${VALIDATED_TRANS_LINKS}" ]]; then
  echo " Trans links:${VALIDATED_TRANS_LINKS}"
else
  echo " Trans links:<none>"
fi
echo "============================================="

run_cmd mkdir -p "${DIR_GENE}" "${DIR_GUIDE}" "${DIR_INTENDED}" "${DIR_TRANS}"

###############################################
# Step 1: gene mapping QC
###############################################
echo "[1/5] Gene mapping QC"
run_cmd python "${SCRIPT_GENE}" \
  --input "${INPUT}" \
  --outdir "${DIR_GENE}" \
  --prefix "${PREFIX_GENE}"

###############################################
# Step 2: guide mapping QC
###############################################
echo "[2/5] Guide mapping QC"
run_cmd python "${SCRIPT_GUIDE}" \
  --input "${INPUT}" \
  --outdir "${DIR_GUIDE}" \
  --prefix "${PREFIX_GUIDE}"

###############################################
# Step 3: intended target QC
###############################################
echo "[3/5] Intended target QC"
run_cmd python "${SCRIPT_INTENDED}" \
  --input "${INPUT}" \
  --outdir "${DIR_INTENDED}" \
  --prefix "${PREFIX_INTENDED}"

###############################################
# Step 4: trans-regulatory QC (NEW)
###############################################
echo "[4/5] Trans-regulatory inference QC"
if [[ -n "${VALIDATED_TRANS_LINKS}" ]]; then
  run_cmd python "${SCRIPT_TRANS}" \
    --input "${INPUT}" \
    --outdir "${DIR_TRANS}" \
    --validated-links "${VALIDATED_TRANS_LINKS}" \
    --prefix "${PREFIX_TRANS}"
else
  run_cmd python "${SCRIPT_TRANS}" \
    --input "${INPUT}" \
    --outdir "${DIR_TRANS}" \
    --prefix "${PREFIX_TRANS}"
fi

###############################################
# Step 5: report
###############################################
echo "[5/5] Generate report"
run_cmd python "${SCRIPT_REPORT}" \
  --output "${REPORT_PDF}" \
  --sample-name "${RUN_NAME}" \
  --prefix "${PREFIX_REPORT}" \
  --mapping-gene-dir "${DIR_GENE}" \
  --mapping-guide-dir "${DIR_GUIDE}" \
  --intended-target-dir "${DIR_INTENDED}" \
  --trans-dir "${DIR_TRANS}"

echo "============================================="
[[ "${DRY_RUN}" -eq 1 ]] && echo " Dry-run complete (no files written)." || echo " Done."
echo "============================================="
