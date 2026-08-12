#!/bin/bash
#SBATCH --job-name=cluster_enrichment
#SBATCH --partition=common
#SBATCH --account=singhlab
#SBATCH --mem=400G
#SBATCH --cpus-per-task=8
#SBATCH --time=20:00:00
#SBATCH --output=cluster_enrichment_%j.out
#SBATCH --error=cluster_enrichment_%j.err

# =============================================================================
# cluster_enrichment.sh
#
# Leiden-clusters all ~1M cells from the TF perturb-seq screen and tests
# whether cells assigned to each TF perturbation are enriched or depleted
# in particular clusters.
#
# Usage:
#   sbatch cluster_enrichment.sh
#   sbatch cluster_enrichment.sh --background all
#   sbatch cluster_enrichment.sh --skip_clustering   # reuse existing cluster_labels.csv
#
# All arguments after the script name are passed through to the Python script.
# =============================================================================

set -euo pipefail

# --- Paths -------------------------------------------------------------------
CONDA_ENV="scanpy_env"
SCRIPT_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/inference"
OUTPUT_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/clustering"
LOG_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/logs"

# --- Parameters (edit here or override via sbatch --export) ------------------
BACKGROUND="nt"          # "nt" or "all"
RESOLUTION="0.5"         # Leiden resolution (higher = more clusters)
N_HVGS="2000"            # HVGs for PCA
N_PCS="50"               # PCA components
N_NEIGHBORS="15"         # k-NN neighbors
MIN_CELLS="50"           # Skip TFs with fewer cells
FDR_THRESH="0.05"        # BH FDR threshold
TOP_N_TFS="40"           # TFs shown in heatmap
RANDOM_SEED="42"

# -----------------------------------------------------------------------------

echo "============================================="
echo " cluster_enrichment.sh"
echo " Job ID      : ${SLURM_JOB_ID}"
echo " Node        : $(hostname)"
echo " Started     : $(date)"
echo " Output dir  : ${OUTPUT_DIR}"
echo " Background  : ${BACKGROUND}"
echo " Resolution  : ${RESOLUTION}"
echo " HVGs        : ${N_HVGS}"
echo " PCs         : ${N_PCS}"
echo " k-NN        : ${N_NEIGHBORS}"
echo " Min cells   : ${MIN_CELLS}"
echo " FDR thresh  : ${FDR_THRESH}"
echo " Top N TFs   : ${TOP_N_TFS}"
echo "============================================="

# Create output and log directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

echo ""
echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo ""

# Check for pynndescent (approximate k-NN -- much faster at 1M cells)
python -c "import pynndescent; print('pynndescent available: ' + pynndescent.__version__)"     2>/dev/null || echo "WARNING: pynndescent not installed -- will fall back to exact k-NN (slower)"
echo ""

# Run the Python script, passing through any extra arguments from the command line
python "${SCRIPT_DIR}/cluster_enrichment.py" \
    --output_dir  "${OUTPUT_DIR}"  \
    --background  "${BACKGROUND}"  \
    --resolution  "${RESOLUTION}"  \
    --n_hvgs      "${N_HVGS}"      \
    --n_pcs       "${N_PCS}"       \
    --n_neighbors "${N_NEIGHBORS}" \
    --min_cells   "${MIN_CELLS}"   \
    --fdr_thresh  "${FDR_THRESH}"  \
    --top_n_tfs   "${TOP_N_TFS}"   \
    --random_seed "${RANDOM_SEED}" \
    "$@"

echo ""
echo "============================================="
echo " Finished: $(date)"
echo "============================================="