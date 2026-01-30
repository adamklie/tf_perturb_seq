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
#   ./patch_files.sh
#
# After running, create a patched sample metadata CSV with updated paths.
#
set -e

# =============================================================================
# CONFIGURATION
# =============================================================================

BUCKET="gs://igvf-pertub-seq-pipeline-data"
DATASET="Huangfu_WTC11-benchmark_TF-Perturb-seq"
DATE="2026_01_30"
PATCH_DIR="${BUCKET}/${DATASET}/${DATE}/patch"

# =============================================================================
# BARCODE ONLIST & GUIDE DESIGN
# =============================================================================

BARCODE_ONLIST_SRC="${BUCKET}/${DATASET}/${DATE}/2024/12/05/996cf65e-5b5c-494d-b8d3-f1f1e8b0181c/IGVFFI9487JPEN.tsv.gz"
GUIDE_DESIGN_SRC="${BUCKET}/${DATASET}/${DATE}/2025/08/04/cc3d00eb-118c-4230-945f-3cdb847cd849/IGVFFI5765HMZH.tsv.gz"

echo "=== Decompressing barcode_onlist ==="
gsutil cat "${BARCODE_ONLIST_SRC}" | gunzip | gsutil cp - "${PATCH_DIR}/IGVFFI9487JPEN.tsv"

echo "=== Decompressing guide_design ==="
gsutil cat "${GUIDE_DESIGN_SRC}" | gunzip | gsutil cp - "${PATCH_DIR}/IGVFFI5765HMZH.tsv"

# =============================================================================
# SEQSPEC YAML FILES
# =============================================================================

SEQSPEC_FILES=(
    "2025/12/11/0b5f464c-5b89-4f63-ad9d-565758d4b22b/IGVFFI5456UAIA.yaml.gz"
    "2025/12/11/1b3b1159-7b86-490b-a159-d61abeaf2d85/IGVFFI6493YTYC.yaml.gz"
    "2025/12/11/396706d4-d8a7-4eb9-831e-44047334d098/IGVFFI6307WAIF.yaml.gz"
    "2025/12/11/417cf9db-e348-4a2a-b4f2-0e6648810dbb/IGVFFI4722FIVQ.yaml.gz"
    "2025/12/11/c051a41d-9976-4c91-90ae-053ae19acc1b/IGVFFI9032ETVM.yaml.gz"
    "2025/12/11/d164081f-0b9a-47f7-8a9f-300bdd30932f/IGVFFI2881ZSXN.yaml.gz"
    "2025/12/11/d403091b-a398-456d-84cf-4bda7ab15849/IGVFFI8819SESP.yaml.gz"
    "2025/12/11/f79ad7b4-28a3-4b03-bd14-e8e46178983b/IGVFFI8798JBXN.yaml.gz"
)

echo ""
echo "=== Decompressing seqspec YAML files ==="
for file in "${SEQSPEC_FILES[@]}"; do
    filename=$(basename "${file}" .gz)
    echo "  ${filename}"
    gsutil cat "${BUCKET}/${DATASET}/${DATE}/${file}" | gunzip | gsutil cp - "${PATCH_DIR}/${filename}"
done

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "=== Done! Files uploaded to ${PATCH_DIR}/ ==="
echo ""
echo "Patched files:"
gsutil ls "${PATCH_DIR}/"
echo ""
echo "Next: Update sample_metadata_gcp_*.csv with paths pointing to patch/ directory"
