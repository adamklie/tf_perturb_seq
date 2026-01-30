#!/bin/bash
# Format MuData to AnnData for gene program discovery
#
# Usage:
#   bash 1_format_anndata.sh
#
# This script converts a MuData (.h5mu) to AnnData (.h5ad) format,
# adding guide assignment to obsm and computing PCA/UMAP embeddings.

set -euo pipefail

###############################################
# PATHS
###############################################

# Project root (adjust if needed)
PROJECT_ROOT="/Users/adamklie/Desktop/projects/tf_perturb_seq"

###############################################
# ENVIRONMENT
###############################################

# Activate project venv (managed by uv)
source "${PROJECT_ROOT}/.venv/bin/activate"

###############################################
# INPUT/OUTPUT
###############################################

INPUT="${PROJECT_ROOT}/datasets/test/results/inference_mudata.selected_guides_1pc_5nt_5tf_2kgenes.h5mu"
OUTDIR="${PROJECT_ROOT}/datasets/test/results/gene_program_discovery"
OUTPUT="${OUTDIR}/adata.h5ad"
SCRIPT="${PROJECT_ROOT}/src/tf_perturb_seq/gene_program_discovery/format_anndata.py"

###############################################
# RUN
###############################################

echo "============================================="
echo " Format AnnData for Gene Program Discovery"
echo " Input:  ${INPUT}"
echo " Output: ${OUTPUT}"
echo "============================================="

mkdir -p "${OUTDIR}"

python "${SCRIPT}" \
    --mudata "${INPUT}" \
    --out_h5ad "${OUTPUT}" \
    --n_pcs 50 \
    --dense_guide_assignment

echo "============================================="
echo " Done! Output: ${OUTPUT}"
echo "============================================="
