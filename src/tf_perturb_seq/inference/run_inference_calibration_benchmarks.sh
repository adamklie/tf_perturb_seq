#!/bin/bash
#SBATCH --job-name=inference_calibration_benchmarks
#SBATCH --output=logs/inference_calibration_benchmarks_%j.out
#SBATCH --error=logs/inference_calibration_benchmarks_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=common

###############################################
# run_inference_calibration_benchmarks.sh
#
# Runs calibrate.py in --from-uns mode on all 5 TF Perturb-seq benchmark
# datasets. Each dataset's inference_mudata.h5mu already carries
# cis_per_element_results / trans_per_element_results directly in .uns
# (unlike calibrate.py's original standalone-TSV workflow), so this calls
# calibrate.py --from-uns for each one and writes results to a Calibration/
# folder under that dataset's own directory.
#
# Trans is calibrated via NTC empirical-null (z-value based, ecdf/t-fit).
# Cis is calibrated via plain BH-FDR on its existing sceptre_p_value /
# perturbo_p_value columns (see calibrate.py's calibrate_cis_bh() for why:
# no std-error column and too few rows for an empirical null).
###############################################

set -euo pipefail

# - Paths -
PROJECT_ROOT="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/inference"
CALIBRATE_SCRIPT="${PROJECT_ROOT}/calibrate.py"
BASE_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets"
NULL_METHOD="t-fit"
LOG_DIR="${PROJECT_ROOT}/logs"

mkdir -p "${LOG_DIR}"

# - Environment -
# conda's activation scripts (e.g. activate.d/env_vars.sh) reference variables
# like LD_LIBRARY_PATH without a default, which trips `set -u` if they're not
# already set in the job's environment. Disable nounset just for activation.
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scanpy_env
set -u

echo "============================================="
echo " run_inference_calibration_benchmarks.sh"
echo " Job:         ${SLURM_JOB_ID:-local}"
echo " Node:        $(hostname)"
echo " Started:     $(date)"
echo " Python:      $(which python)"
echo " Script:      ${CALIBRATE_SCRIPT}"
echo " Null method: ${NULL_METHOD}"
echo "============================================="

[[ -f "${CALIBRATE_SCRIPT}" ]] || { echo "ERROR: calibrate.py not found: ${CALIBRATE_SCRIPT}" >&2; exit 1; }

# One entry per benchmark dataset: "<dir under BASE_DIR>:<output prefix>"
DATASETS=(
  "Engreitz_WTC11-benchmark_TF-Perturb-seq:engreitz"
  "Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3:gersbach_gemx"
  "Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2:gersbach_htv2"
  "Hon_WTC11-benchmark_TF-Perturb-seq:hon"
  "Huangfu_WTC11-benchmark_TF-Perturb-seq:huangfu"
)

N_FAILED=0

for entry in "${DATASETS[@]}"; do
  DIR_NAME="${entry%%:*}"
  PREFIX="${entry##*:}"

  MUDATA="${BASE_DIR}/${DIR_NAME}/sceptre/inference_mudata.h5mu"
  OUTDIR="${BASE_DIR}/${DIR_NAME}/Calibration"

  echo ""
  echo "--- ${PREFIX} ---"
  echo "MuData: ${MUDATA}"
  echo "Outdir: ${OUTDIR}"

  if [[ ! -f "${MUDATA}" ]]; then
    echo "WARNING: MuData not found, skipping ${PREFIX}: ${MUDATA}" >&2
    N_FAILED=$((N_FAILED + 1))
    continue
  fi

  mkdir -p "${OUTDIR}"

  if python "${CALIBRATE_SCRIPT}" \
      --mudata "${MUDATA}" \
      --outdir "${OUTDIR}" \
      --prefix "${PREFIX}" \
      --from-uns \
      --null-method "${NULL_METHOD}"; then
    echo "Done: ${PREFIX}"
  else
    echo "WARNING: calibration failed for ${PREFIX} (see traceback above); continuing with remaining datasets" >&2
    N_FAILED=$((N_FAILED + 1))
  fi
done

echo ""
echo "============================================="
echo " Finished: $(date)"
if [[ "${N_FAILED}" -gt 0 ]]; then
  echo " ${N_FAILED} dataset(s) failed or were skipped -- see WARNINGs above"
else
  echo " All datasets calibrated successfully"
fi
echo "============================================="
