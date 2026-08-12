#!/bin/bash
#SBATCH --job-name=edist_viz
#SBATCH --partition=common
#SBATCH --account=crawfordlab
#SBATCH --mem=100G
#SBATCH --cpus-per-task=4
#SBATCH --time=0:30:00
#SBATCH --output=logs/edist_viz_%j.out
#SBATCH --error=logs/edist_viz_%j.err
# =============================================================================
# run_energy_dist_visualize.sh
#
# Runs visualize_energy_dist_results.py.
#
# Usage:
#   sbatch run_energy_dist_visualize.sh
#   sbatch run_energy_dist_visualize.sh --fdr_thresh 1e-3 --leiden_resolution 1.0
#
# Any arguments after the script name are forwarded to the Python script.
# =============================================================================
set -eo pipefail
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# --- Edit these for a new dataset --------------------------------------------
CONDA_ENV="scanpy_env"
SCRIPT_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/energy_distance"
LOG_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/energy_distance/logs"

ENERGY_DIST_DIR="/hpc/group/gersbachlab/agk21/hep_perturbseq/energy_distance"
OUTPUT_DIR="${SCRIPT_DIR}/visualizations"

FDR_THRESH=1e-5
LEIDEN_RESOLUTION=2.0
# -----------------------------------------------------------------------------

mkdir -p "${LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

echo "============================================="
echo " run_energy_dist_visualize.sh"
echo " Job ID             : ${SLURM_JOB_ID:-local}"
echo " Node               : $(hostname)"
echo " Started            : $(date)"
echo " Energy dist dir    : ${ENERGY_DIST_DIR}"
echo " Output dir         : ${OUTPUT_DIR}"
echo " FDR threshold      : ${FDR_THRESH}"
echo " Leiden resolution  : ${LEIDEN_RESOLUTION}"
echo "============================================="

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
echo "Python: $(which python) ($(python --version))"
echo ""

python "${SCRIPT_DIR}/visualize_energy_dist_results.py" \
    "${ENERGY_DIST_DIR}"                    \
    --outdir            "${OUTPUT_DIR}"     \
    --fdr_thresh        "${FDR_THRESH}"     \
    --leiden_resolution "${LEIDEN_RESOLUTION}" \
    "$@"

echo ""
echo "============================================="
echo " Finished: $(date)"
echo "============================================="