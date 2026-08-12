#!/bin/bash
#SBATCH --job-name=cis_inference
#SBATCH --partition=common
#SBATCH --account=singhlab
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --output=logs/cis_inference_%j.out
#SBATCH --error=logs/cis_inference_%j.err
# =============================================================================
# run_cis_inference.sh
#
# Runs investigate_cis_inference_results.py for the Gersbach dataset.
#
# Usage:
#   sbatch run_cis_inference.sh                 # default (no UMAP)
#   sbatch run_cis_inference.sh --run_umap      # include UMAP analysis
#   sbatch run_cis_inference.sh --use_sceptre_only
#
# Any arguments after the script name are passed through to the Python script.
# Calibrated result files are auto-detected from INFERENCE_DIR by suffix
# (*_calibrated_cis_results.tsv, *_calibrated_direct_target_results.tsv).
# --calib_prefix is no longer needed.
# =============================================================================
set -eo pipefail
# LD_LIBRARY_PATH may be unset in the SLURM environment; default it to empty
# so that conda's activate.d scripts don't trip the -u (nounset) flag.
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# --- Edit these for a new dataset --------------------------------------------
CONDA_ENV="scanpy_env"
SCRIPT_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/inference"
LOG_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/inference/logs"

# Gersbach dataset
DATASET_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq"
PIPELINE_OUTPUTS_DIR="${DATASET_DIR}/post_pipeline_processing"
#INFERENCE_DIR="${DATASET_DIR}/calibrated_outs"
INFERENCE_DIR="/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/helen_output_poolabcd_prod_v10/pipeline_outputs/"

# MuData filename inside PIPELINE_OUTPUTS_DIR
MUDATA_NAME="inference_mudata_outlier_guide_filt.h5mu"

# RESULTS_DIR: where cis_per_guide_results.tsv.gz / cis_per_element_results.tsv.gz live.
# If these files are also in calibrated_outs, set RESULTS_DIR=INFERENCE_DIR.
#RESULTS_DIR="${DATASET_DIR}/calibrated_outs"
RESULTS_DIR=INFERENCE_DIR

GUIDE_METADATA="/hpc/group/gersbachlab/seg95/crispr-pipeline-personal/example-data/production_guide_metadata_v6.tsv"
OUTPUT_DIR="${DATASET_DIR}/calibrated_outs_manual/cis_inference_visualizations"

# Analysis parameters
TOP_N_TFS=20
DOWNSAMPLE=50000
P_THRESH=0.05
# -----------------------------------------------------------------------------

mkdir -p "${LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

echo "============================================="
echo " run_cis_inference.sh"
echo " Job ID             : ${SLURM_JOB_ID:-local}"
echo " Node               : $(hostname)"
echo " Started            : $(date)"
echo " Dataset dir        : ${DATASET_DIR}"
echo " Pipeline outputs   : ${PIPELINE_OUTPUTS_DIR}"
echo " Inference/calib dir: ${INFERENCE_DIR}"
echo " MuData filename    : ${MUDATA_NAME}"
echo " Results dir        : ${RESULTS_DIR}"
echo " Output dir         : ${OUTPUT_DIR}"
echo " Downsample (UMAP)  : ${DOWNSAMPLE}"
echo "============================================="

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
echo "Python: $(which python) ($(python --version))"
echo ""

python "${SCRIPT_DIR}/investigate_cis_inference_results.py" \
    --results_dir           "${RESULTS_DIR}"          \
    --pipeline_outputs_dir  "${PIPELINE_OUTPUTS_DIR}" \
    --inference_dir         "${INFERENCE_DIR}"         \
    --mudata_name           "${MUDATA_NAME}"           \
    --guide_metadata        "${GUIDE_METADATA}"        \
    --output_dir            "${OUTPUT_DIR}"            \
    --top_n_tfs             "${TOP_N_TFS}"             \
    --downsample            "${DOWNSAMPLE}"            \
    --p_thresh              "${P_THRESH}"              \
    --run_umap 					       \
    "$@"

#

echo ""
echo "============================================="
echo " Finished: $(date)"
echo "============================================="