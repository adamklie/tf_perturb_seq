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
DATE="2026_04_13"
PATCH_DIR="${BUCKET}/${DATASET}/${DATE}/patch"

# =============================================================================
# GUIDE DESIGN, BARCODE
# =============================================================================

BARCODE_ONLIST_SRC="${BUCKET}/${DATASET}/${DATE}/2025/03/20/3a502f60-e1a1-4c82-9b41-c1201f9d92ed/IGVFFI4695IKAL.tsv.gz"
GUIDE_DESIGN_SRC="${BUCKET}/${DATASET}/${DATE}/2024/08/21/2dca8d2e-e59f-425d-b0d4-835ec3fb61f9/IGVFFI8270UPKB.csv.gz"

echo "=== Decompressing barcode_onlist ==="
gsutil cat "${BARCODE_ONLIST_SRC}" | gunzip | gsutil cp - "${PATCH_DIR}/IGVFFI4695IKAL.tsv"

echo "=== Decompressing guide_design ==="
gsutil cat "${GUIDE_DESIGN_SRC}" | gunzip | gsutil cp - "${PATCH_DIR}/IGVFFI8270UPKB.tsv"

# =============================================================================
# SEQSPEC YAML FILES
# =============================================================================

SEQSPEC_FILES=(
    "michael-beer_TF-perturb-seq-ES-rep1.yaml.gz"
    "michael-beer_TF-perturb-seq-ES-rep2.yaml.gz"
    "michael-beer_TF-perturb-seq-ES-rep3.yaml.gz"
    "michael-beer_TF-perturb-seq-ES-rep4.yaml.gz"
    "michael-beer_TF-perturb-seq-ES-rep5.yaml.gz"
    "michael-beer_TF-perturb-seq-ES-rep6.yaml.gz"
    "michael-beer_TF-perturb-seq-ES-rep7.yaml.gz"
    "michael-beer_TF-perturb-seq-ES-rep8.yaml.gz"
    "michael-beer_TF-perturb-seq-ES-guide-sequencing-rep1.yaml.gz"
    "michael-beer_TF-perturb-seq-ES-guide-sequencing-rep2.yaml.gz"
    "michael-beer_TF-perturb-seq-ES-guide-sequencing-rep3.yaml.gz"
    "michael-beer_TF-perturb-seq-ES-guide-sequencing-rep4.yaml.gz"
    "michael-beer_TF-perturb-seq-ES-guide-sequencing-rep5.yaml.gz"
    "michael-beer_TF-perturb-seq-ES-guide-sequencing-rep6.yaml.gz"
    "michael-beer_TF-perturb-seq-ES-guide-sequencing-rep7.yaml.gz"
    "michael-beer_TF-perturb-seq-ES-guide-sequencing-rep8.yaml.gz"
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
