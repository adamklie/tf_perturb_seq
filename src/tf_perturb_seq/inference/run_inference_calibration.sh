#!/usr/bin/env bash
set -euo pipefail

###############################################
# run_calibration.sh
#
# Run empirical null calibration on a single dataset's PerTurbo results.
#
# USAGE:
#   run_calibration.sh \
#     --project-root /path/to/tf_perturb_seq \
#     --trans-results /path/to/perturbo_trans_per_element_output.tsv.gz \
#     --mudata /path/to/inference_mudata.h5mu \
#     --outdir /path/to/output \
#     --prefix dataset_name \
#     [--cis-results /path/to/perturbo_cis_per_element_output.tsv.gz] \
#     [--null-method t-fit|ecdf] \
#     [--dry-run]
###############################################

PROJECT_ROOT=""
TRANS_RESULTS=""
MUDATA=""
OUTDIR=""
PREFIX=""
NULL_METHOD="t-fit"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root)   PROJECT_ROOT="$2";   shift 2 ;;
    --trans-results)  TRANS_RESULTS="$2";   shift 2 ;;
    --mudata)         MUDATA="$2";         shift 2 ;;
    --outdir)         OUTDIR="$2";         shift 2 ;;
    --prefix)         PREFIX="$2";         shift 2 ;;
    --null-method)    NULL_METHOD="$2";    shift 2 ;;
    --dry-run)        DRY_RUN=1;           shift 1 ;;
    -h|--help)
      sed -n '1,25p' "$0"
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
[[ -n "${TRANS_RESULTS}" ]] || { echo "ERROR: --trans-results required" >&2; exit 1; }
[[ -n "${MUDATA}" ]] || { echo "ERROR: --mudata required" >&2; exit 1; }
[[ -n "${OUTDIR}" ]] || { echo "ERROR: --outdir required" >&2; exit 1; }
[[ -n "${PREFIX}" ]] || { echo "ERROR: --prefix required" >&2; exit 1; }

[[ -f "${TRANS_RESULTS}" ]] || { echo "ERROR: Trans results not found: ${TRANS_RESULTS}" >&2; exit 1; }
[[ -f "${MUDATA}" ]] || { echo "ERROR: MuData not found: ${MUDATA}" >&2; exit 1; }

CALIBRATE_SCRIPT="${PROJECT_ROOT}/calibrate.py"
[[ -f "${CALIBRATE_SCRIPT}" ]] || { echo "ERROR: calibrate.py not found: ${CALIBRATE_SCRIPT}" >&2; exit 1; }

# Build command
CMD="python ${CALIBRATE_SCRIPT}"
CMD+=" --trans-results ${TRANS_RESULTS}"
CMD+=" --mudata ${MUDATA}"
CMD+=" --outdir ${OUTDIR}"
CMD+=" --prefix ${PREFIX}"
CMD+=" --null-method ${NULL_METHOD}"

echo "============================================="
echo " run_calibration.sh"
echo " Trans:       ${TRANS_RESULTS}"
echo " MuData:      ${MUDATA}"
echo " Outdir:      ${OUTDIR}"
echo " Prefix:      ${PREFIX}"
echo " Method:      ${NULL_METHOD}"
echo "============================================="

if [[ "${DRY_RUN}" -eq 1 ]]; then
  echo "[DRY-RUN] ${CMD}"
  exit 0
fi

mkdir -p "${OUTDIR}"
eval "${CMD}"


# Sample run:
#srun --mem=200G --account=crawfordlab --pty bash -c "
#  source $(conda info --base)/etc/profile.d/conda.sh && \
#  conda activate scanpy_env && \
#  bash /hpc/group/gersbachlab/seg95/crispr-pipeline-personal/run_inference_calibration.sh \
#    --project-root /hpc/group/gersbachlab/seg95/crispr-pipeline-personal \
#    --trans-results /hpc/group/gersbachlab/seg95/CRISPR_Pipeline/helen_output_poolabcd_prod_v10/pipeline_outputs/perturbo_trans_per_element_output.tsv.gz \
#    --mudata /hpc/group/gersbachlab/seg95/CRISPR_Pipeline/helen_output_poolabcd_prod_v10/pipeline_outputs/inference_mudata.h5mu \
#    --outdir /hpc/group/gersbachlab/seg95/CRISPR_Pipeline/helen_output_poolabcd_prod_v10/pipeline_outputs \
#    --prefix gersbach_ihep \
#    --null-method t-fit 
#"

#srun --mem=200G --account=crawfordlab --pty bash -c "
#  source $(conda info --base)/etc/profile.d/conda.sh && \
#  conda activate scanpy_env && \
#  bash /hpc/group/gersbachlab/seg95/crispr-pipeline-personal/run_inference_calibration.sh \
#    --project-root /hpc/group/gersbachlab/seg95/crispr-pipeline-personal \
#    --trans-results /hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/perturbo_trans_per_element_output.tsv.gz \
#    --mudata /hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/inference_mudata.h5mu \
#    --outdir /hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq \
#    --prefix hon_cardio \
#    --null-method t-fit 
#"
