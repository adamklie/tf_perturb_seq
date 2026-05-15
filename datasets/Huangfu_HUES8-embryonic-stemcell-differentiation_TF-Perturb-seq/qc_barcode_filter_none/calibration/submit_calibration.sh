#!/usr/bin/env bash
#SBATCH -J calib_ESC_nobc
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH --partition=carter-compute
#SBATCH -t 02:00:00
#SBATCH -o /cellar/users/aklie/projects/tf_perturb_seq/scratch/calibration_logs/calib_ESC_nobc.%j.out
#SBATCH -e /cellar/users/aklie/projects/tf_perturb_seq/scratch/calibration_logs/calib_ESC_nobc.%j.err

set -euo pipefail

PROJECT_ROOT="/cellar/users/aklie/projects/tf_perturb_seq"
RUN_DIR="${PROJECT_ROOT}/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/qc_barcode_filter_none"

bash "${PROJECT_ROOT}/scripts/run_calibration.sh" \
  --project-root "${PROJECT_ROOT}" \
  --trans-results "${RUN_DIR}/crispr_pipeline/pipeline_outputs/perturbo_trans_per_element_output.tsv.gz" \
  --mudata "${RUN_DIR}/crispr_pipeline/pipeline_dashboard/inference_mudata.h5mu" \
  --outdir "${RUN_DIR}/calibration" \
  --prefix "Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq_qc_barcode_filter_none"
