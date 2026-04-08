#!/bin/bash
#
# Decompress .gz files on GCS to a patch directory
#
# The pipeline requires uncompressed files for:
# - barcode_onlist (.tsv)
# - guide_design (.tsv)
# - seqspec (.yaml)
#
# This script decompresses them via streaming (no local temp files).
#
# Usage:
#   ./3_patch_gcp_files.sh
#
# After running, create a patched sample metadata CSV with updated paths.
#
set -e

# =============================================================================
# CONFIGURATION
# =============================================================================

BUCKET="gs://igvf-pertub-seq-pipeline-data"
DATASET="Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq"
# TODO: Update DATE to match the upload date from step 2
DATE="YYYY_MM_DD"
PATCH_DIR="${BUCKET}/${DATASET}/${DATE}/patch"

# =============================================================================
# TODO: Fill in the actual file paths after step 1 generates the metadata
# Use the Huangfu benchmark (3_patch_gcp_files.sh) as a reference for structure
# =============================================================================

# BARCODE_ONLIST_SRC="${BUCKET}/${DATASET}/${DATE}/path/to/barcode_onlist.tsv.gz"
# GUIDE_DESIGN_SRC="${BUCKET}/${DATASET}/${DATE}/path/to/guide_design.tsv.gz"

# echo "=== Decompressing barcode_onlist ==="
# gsutil cat "${BARCODE_ONLIST_SRC}" | gunzip | gsutil cp - "${PATCH_DIR}/barcode_onlist.tsv"

# echo "=== Decompressing guide_design ==="
# gsutil cat "${GUIDE_DESIGN_SRC}" | gunzip | gsutil cp - "${PATCH_DIR}/guide_design.tsv"

# =============================================================================
# SEQSPEC YAML FILES
# TODO: Fill in after step 1
# =============================================================================

# SEQSPEC_FILES=(
#     "path/to/seqspec1.yaml.gz"
#     "path/to/seqspec2.yaml.gz"
# )

# echo ""
# echo "=== Decompressing seqspec YAML files ==="
# for file in "${SEQSPEC_FILES[@]}"; do
#     filename=$(basename "${file}" .gz)
#     echo "  ${filename}"
#     gsutil cat "${BUCKET}/${DATASET}/${DATE}/${file}" | gunzip | gsutil cp - "${PATCH_DIR}/${filename}"
# done

echo ""
echo "TODO: Fill in actual file paths from step 1 metadata output"
echo "See datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/3_patch_gcp_files.sh for reference"
