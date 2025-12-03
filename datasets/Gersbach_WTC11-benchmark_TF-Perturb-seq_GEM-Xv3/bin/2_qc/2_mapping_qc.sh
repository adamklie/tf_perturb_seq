#!/bin/bash
#SBATCH --job-name=Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3_mapping_qc
#SBATCH --partition=carter-gpu
#SBATCH --account=carter-gpu
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=14-00:00:00
#SBATCH --output=/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/bin/2_qc/%x.%A.out

#####
# USAGE:
#   sbatch mapping_qc.sh
#####

set -euo pipefail

###############################################
# USER-DEFINED PARAMETERS (EDIT THESE ONLY)
###############################################

# Dataset name
SAMPLE="Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3"

# Input MuData from CRISPR pipeline
H5MU="/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/results/1_CRISPR_pipeline/2025_11_26/inference_mudata.h5mu"

# Output directory
OUT="/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/results/2_qc/2025_12_01"

# Location of scripts
SCRIPTS="/cellar/users/aklie/data/datasets/tf_perturb_seq/scripts/qc"

# Conda envs
ENV_QC="scverse-lite-py39"
ENV_GPU="rapids_singlecell"

###############################################
# DERIVED PATHS (NO NEED TO EDIT)
###############################################
MAPPINGQC_DIR="${OUT}/mapping_qc"
PRE_DIR="${OUT}/preprocess"
CLUST_DIR="${OUT}/clustering"

GENE_FILTERED="${MAPPINGQC_DIR}/gene.filtered.h5ad"
GENE_PREPRO="${PRE_DIR}/gene.preprocessed.h5ad"
GENE_CLUSTERED="${CLUST_DIR}/gene.clustered.h5ad"

PRE_LOG="${PRE_DIR}/preprocess.log"
CLUST_LOG="${CLUST_DIR}/clustering.log"

###############################################
# START PIPELINE
###############################################
echo "============================================="
echo " Running TF Perturb-seq QC Pipeline"
echo " SAMPLE: ${SAMPLE}"
echo " OUTDIR: ${OUT}"
echo "============================================="

mkdir -p "${MAPPINGQC_DIR}"

###############################################
# 1) QC + FILTERING
###############################################
echo "---- [1/3] Running mapping_qc.py ----"

source activate /cellar/users/aklie/opt/miniconda3/envs/${ENV_QC}
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

python "${SCRIPTS}/mapping_qc.py" \
    --input-h5mu "${H5MU}" \
    --outdir "${MAPPINGQC_DIR}" \
    --sample-name "${SAMPLE}" \
    --plot-gene-by-batch \
    --write-filtered-h5ad

echo "QC complete → ${GENE_FILTERED}"
echo

###############################################
# 2) PREPROCESSING
###############################################
echo "---- [2/3] Running preprocess.py ----"

conda activate /cellar/users/aklie/opt/miniconda3/envs/rapids_singlecell

mkdir -p "${PRE_DIR}"

python "${SCRIPTS}/preprocess.py" \
    --input "${GENE_FILTERED}" \
    --output "${GENE_PREPRO}" \
    --plot_dir "${PRE_DIR}" \
    --target_sum 10000 \
    --n_top_genes 5000 \
    --max_value 10 \
    --log_file "${PRE_LOG}"

echo "Preprocessing complete → ${GENE_PREPRO}"
echo

###############################################
# 3) CLUSTERING
###############################################
echo "---- [3/3] Running cluster.py ----"

mkdir -p "${CLUST_DIR}"

python "${SCRIPTS}/cluster.py" \
    --input "${GENE_PREPRO}" \
    --output "${GENE_CLUSTERED}" \
    --plot_dir "${CLUST_DIR}" \
    --n_pcs 40 \
    --n_neighbors 15 \
    --resolution 0.6 \
    --log_file "${CLUST_LOG}"

echo "Clustering complete → ${GENE_CLUSTERED}"
echo

###############################################
# DONE
###############################################
echo "============================================="
echo "  Pipeline finished for sample: ${SAMPLE}"
echo "  Outputs written to: ${OUT}"
echo "============================================="
