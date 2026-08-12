#!/bin/bash
#SBATCH --job-name=cluster_vis
#SBATCH --partition=common
#SBATCH --account=singhlab
#SBATCH --mem=100G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --output=logs/cluster_vis_%j.out
#SBATCH --error=logs/cluster_vis_%j.err

# =============================================================================
# run_cluster_enrichment_visualization.sh
#
# Generates UMAP plots from the outputs of cluster_enrichment.py:
#   - umap_leiden_clusters.png          (all cells, colored by cluster)
#   - umap_per_cluster_highlights.png   (one panel per cluster)
#   - umap_per_tf_page*.png             (paginated, one panel per TF)
#
# Requires cluster_labels.csv and pca_coords.npz from cluster_enrichment.py.
# UMAP coordinates are computed from PCA and saved to umap_coords.npz so
# subsequent runs can use --skip_umap to replot without recomputing.
#
# Usage:
#   sbatch run_cluster_enrichment_visualization.sh
#   sbatch run_cluster_enrichment_visualization.sh --skip_umap   # replot only
# =============================================================================

set -euo pipefail

# --- Edit these for a new dataset --------------------------------------------
CONDA_ENV="scanpy_env"
SCRIPT_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/inference"
LOG_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/logs"

CLUSTERING_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/clustering"
OUTPUT_DIR="${CLUSTERING_DIR}"   # save plots alongside cluster outputs

MUDATA_PATH="/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/helen_output_poolabcd_prod_v9/pipeline_outputs/inference_mudata.h5mu"
GUIDE_METADATA="/hpc/group/gersbachlab/seg95/crispr-pipeline-personal/example-data/production_guide_metadata_v6.tsv"

# UMAP / plot parameters
DOWNSAMPLE=200000    # cells to plot (balanced across clusters)
N_NEIGHBORS=15
MIN_DIST=0.3
TFS_PER_PAGE=12
MIN_CELLS_TF=50
RANDOM_SEED=42
# -----------------------------------------------------------------------------

mkdir -p "${LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

echo "============================================="
echo " run_cluster_enrichment_visualization.sh"
echo " Job ID          : ${SLURM_JOB_ID}"
echo " Node            : $(hostname)"
echo " Started         : $(date)"
echo " Clustering dir  : ${CLUSTERING_DIR}"
echo " Output dir      : ${OUTPUT_DIR}"
echo " Downsample      : ${DOWNSAMPLE}"
echo " n_neighbors     : ${N_NEIGHBORS}"
echo " min_dist        : ${MIN_DIST}"
echo "============================================="

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

echo "Python: $(which python) ($(python --version))"

# Check for umap-learn (preferred) vs scanpy fallback
python -c "import umap; print('umap-learn available: ' + umap.__version__)" \
    2>/dev/null || echo "WARNING: umap-learn not installed -- will use scanpy UMAP fallback"
echo ""

python "${SCRIPT_DIR}/cluster_enrichment_visualization.py" \
    --clustering_dir  "${CLUSTERING_DIR}"  \
    --output_dir      "${OUTPUT_DIR}"      \
    --mudata_path     "${MUDATA_PATH}"     \
    --guide_metadata  "${GUIDE_METADATA}"  \
    --downsample      "${DOWNSAMPLE}"      \
    --n_neighbors     "${N_NEIGHBORS}"     \
    --min_dist        "${MIN_DIST}"        \
    --tfs_per_page    "${TFS_PER_PAGE}"    \
    --min_cells_tf    "${MIN_CELLS_TF}"    \
    --random_seed     "${RANDOM_SEED}"     \
    "$@"

echo ""
echo "============================================="
echo " Finished: $(date)"
echo "============================================="