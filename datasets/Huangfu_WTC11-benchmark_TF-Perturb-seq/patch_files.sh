#!/bin/bash
# Script to create patch directory and copy/unzip barcode_onlist and guide_design files

set -e

# Define paths
BUCKET="gs://igvf-pertub-seq-pipeline-data"
DATASET="Huangfu_WTC11-benchmark_TF-Perturb-seq"
DATE="2026_01_30"
PATCH_DIR="${BUCKET}/${DATASET}/${DATE}/patch"

# Source files (gzipped)
BARCODE_ONLIST_SRC="${BUCKET}/${DATASET}/${DATE}/2024/12/05/996cf65e-5b5c-494d-b8d3-f1f1e8b0181c/IGVFFI9487JPEN.tsv.gz"
GUIDE_DESIGN_SRC="${BUCKET}/${DATASET}/${DATE}/2025/08/04/cc3d00eb-118c-4230-945f-3cdb847cd849/IGVFFI5765HMZH.tsv.gz"

# Destination files (unzipped)
BARCODE_ONLIST_DST="${PATCH_DIR}/IGVFFI9487JPEN.tsv"
GUIDE_DESIGN_DST="${PATCH_DIR}/IGVFFI5765HMZH.tsv"

echo "Decompressing and uploading barcode_onlist file..."
gsutil cat "${BARCODE_ONLIST_SRC}" | gunzip | gsutil cp - "${BARCODE_ONLIST_DST}"

echo "Decompressing and uploading guide_design file..."
gsutil cat "${GUIDE_DESIGN_SRC}" | gunzip | gsutil cp - "${GUIDE_DESIGN_DST}"

echo "Done! Files uploaded to:"
echo "  ${BARCODE_ONLIST_DST}"
echo "  ${GUIDE_DESIGN_DST}"
