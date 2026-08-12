#!/bin/bash
#SBATCH --job-name=trans_inference
#SBATCH --partition=common
#SBATCH --account=singhlab
#SBATCH --mem=200G
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --output=logs/trans_inference_%j.out
#SBATCH --error=logs/trans_inference_%j.err
# =============================================================================
# run_trans_inference.sh
#
# Usage:
#   sbatch run_trans_inference.sh               # default
#   sbatch run_trans_inference.sh --sceptre_only
#   sbatch run_trans_inference.sh --keep_nontargeting
#
# Calibrated trans file is auto-detected from INFERENCE_DIR by suffix
# (*_calibrated_trans_results.tsv). --calib_prefix is no longer needed.
# Guide-level analysis is skipped automatically if no per-guide file exists.
# =============================================================================
set -eo pipefail
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# --- Edit these for a new dataset --------------------------------------------
CONDA_ENV="scanpy_env"
SCRIPT_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/inference"
LOG_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/inference/logs"

DATASET_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq"
PIPELINE_OUTPUTS_DIR="${DATASET_DIR}/post_pipeline_processing"
#INFERENCE_DIR="${DATASET_DIR}/calibrated_outs"
INFERENCE_DIR="/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/helen_output_poolabcd_prod_v10/pipeline_outputs/"

# MuData filename inside PIPELINE_OUTPUTS_DIR
MUDATA_NAME="inference_mudata_outlier_guide_filt.h5mu"

GUIDE_METADATA="/hpc/group/gersbachlab/seg95/crispr-pipeline-personal/example-data/production_guide_metadata_v6.tsv"
OUTPUT_DIR="${DATASET_DIR}/calibrated_outs_manual/trans_inference_visualizations"

# Analysis parameters
FDR_THRESH=0.05
TOP_N=20
# -----------------------------------------------------------------------------

mkdir -p "${LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

echo "============================================="
echo " run_trans_inference.sh"
echo " Job ID             : ${SLURM_JOB_ID:-local}"
echo " Node               : $(hostname)"
echo " Started            : $(date)"
echo " Dataset dir        : ${DATASET_DIR}"
echo " Pipeline outputs   : ${PIPELINE_OUTPUTS_DIR}"
echo " Inference/calib dir: ${INFERENCE_DIR}"
echo " MuData filename    : ${MUDATA_NAME}"
echo " Output dir         : ${OUTPUT_DIR}"
echo " FDR threshold      : ${FDR_THRESH}"
echo "============================================="

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
echo "Python: $(which python) ($(python --version))"
echo ""

python "${SCRIPT_DIR}/investigate_trans_inference_results.py" \
    --pipeline_outputs_dir  "${PIPELINE_OUTPUTS_DIR}" \
    --inference_dir         "${INFERENCE_DIR}"         \
    --mudata_name           "${MUDATA_NAME}"           \
    --guide_metadata        "${GUIDE_METADATA}"        \
    --output_dir            "${OUTPUT_DIR}"            \
    --fdr_thresh            "${FDR_THRESH}"            \
    --top_n                 "${TOP_N}"                 \
    "$@"

echo ""
echo "============================================="
echo " Finished: $(date)"
echo "============================================="