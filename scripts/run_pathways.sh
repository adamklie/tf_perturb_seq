#!/usr/bin/env bash
set -euo pipefail

###############################################
# run_pathways.sh
#
# Run per-element pathway enrichment on a single dataset's calibrated results.
#
# USAGE:
#   run_pathways.sh \
#     --project-root /path/to/tf_perturb_seq \
#     --calibrated-results /path/to/calibrated_all_results.tsv \
#     --gene-sets /path/to/gene_sets.gmt \
#     --outdir /path/to/output \
#     --prefix dataset_name \
#     [--fdr-threshold 0.10] \
#     [--min-genes 3] \
#     [--max-pathway-size 500] \
#     [--min-pathway-size 10] \
#     [--dry-run]
###############################################

PROJECT_ROOT=""
CALIBRATED_RESULTS=""
GENE_SETS=""
OUTDIR=""
PREFIX=""
FDR_THRESHOLD="0.10"
MIN_GENES="3"
MAX_PATHWAY_SIZE="500"
MIN_PATHWAY_SIZE="10"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root)         PROJECT_ROOT="$2";         shift 2 ;;
    --calibrated-results)   CALIBRATED_RESULTS="$2";   shift 2 ;;
    --gene-sets)            GENE_SETS="$2";            shift 2 ;;
    --outdir)               OUTDIR="$2";               shift 2 ;;
    --prefix)               PREFIX="$2";               shift 2 ;;
    --fdr-threshold)        FDR_THRESHOLD="$2";        shift 2 ;;
    --min-genes)            MIN_GENES="$2";            shift 2 ;;
    --max-pathway-size)     MAX_PATHWAY_SIZE="$2";     shift 2 ;;
    --min-pathway-size)     MIN_PATHWAY_SIZE="$2";     shift 2 ;;
    --dry-run)              DRY_RUN=1;                 shift 1 ;;
    -h|--help)
      sed -n '1,20p' "$0"
      exit 0
      ;;
    *)
      echo "ERROR: Unknown arg: $1" >&2
      exit 1
      ;;
  esac
done

# Validate required args
[[ -n "${PROJECT_ROOT}" ]] || { echo "ERROR: --project-root required" >&2; exit 1; }
[[ -n "${CALIBRATED_RESULTS}" ]] || { echo "ERROR: --calibrated-results required" >&2; exit 1; }
[[ -n "${GENE_SETS}" ]] || { echo "ERROR: --gene-sets required" >&2; exit 1; }
[[ -n "${OUTDIR}" ]] || { echo "ERROR: --outdir required" >&2; exit 1; }
[[ -n "${PREFIX}" ]] || { echo "ERROR: --prefix required" >&2; exit 1; }

[[ -f "${CALIBRATED_RESULTS}" ]] || { echo "ERROR: Calibrated results not found: ${CALIBRATED_RESULTS}" >&2; exit 1; }
[[ -f "${GENE_SETS}" ]] || { echo "ERROR: Gene sets not found: ${GENE_SETS}" >&2; exit 1; }

PATHWAYS_SCRIPT="${PROJECT_ROOT}/src/tf_perturb_seq/inference/pathways.py"
[[ -f "${PATHWAYS_SCRIPT}" ]] || { echo "ERROR: pathways.py not found: ${PATHWAYS_SCRIPT}" >&2; exit 1; }

# Activate venv
source "${PROJECT_ROOT}/.venv/bin/activate"

echo "============================================="
echo " run_pathways.sh"
echo " Results:     ${CALIBRATED_RESULTS}"
echo " Gene sets:   ${GENE_SETS}"
echo " Outdir:      ${OUTDIR}"
echo " Prefix:      ${PREFIX}"
echo " FDR:         ${FDR_THRESHOLD}"
echo "============================================="

CMD="python ${PATHWAYS_SCRIPT}"
CMD+=" --calibrated-results ${CALIBRATED_RESULTS}"
CMD+=" --gene-sets ${GENE_SETS}"
CMD+=" --outdir ${OUTDIR}"
CMD+=" --prefix ${PREFIX}"
CMD+=" --fdr-threshold ${FDR_THRESHOLD}"
CMD+=" --min-genes ${MIN_GENES}"
CMD+=" --max-pathway-size ${MAX_PATHWAY_SIZE}"
CMD+=" --min-pathway-size ${MIN_PATHWAY_SIZE}"

if [[ "${DRY_RUN}" -eq 1 ]]; then
  echo "[DRY-RUN] ${CMD}"
  exit 0
fi

mkdir -p "${OUTDIR}"
eval "${CMD}"
