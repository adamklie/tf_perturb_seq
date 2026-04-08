#!/bin/bash
#
# Decompress .gz files on GCS to a patch directory
#
# The pipeline requires uncompressed files for:
# - barcode_onlist (.tsv)
# - guide_design (.tsv/.csv)
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
DATASET="Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq"
# TODO: Update DATE to match the upload date from step 2
DATE="YYYY_MM_DD"
PATCH_DIR="${BUCKET}/${DATASET}/${DATE}/patch"

# =============================================================================
# TODO: Fill in the actual file paths after step 1 generates the metadata
# =============================================================================
# Reference: The old GCP upload had these files:
#   seqspec:       .../2025/07/21/.../IGVFFI5565FKGG.yaml.gz
#   barcode_onlist: .../2024/12/05/.../IGVFFI9487JPEN.tsv.gz
#   guide_design:   .../2024/08/21/.../IGVFFI8270UPKB.csv.gz
#
# These may have changed -- use the fresh metadata from step 1.
# See datasets/Hon_WTC11-benchmark_TF-Perturb-seq/3_patch_gcp_files.sh for reference.

echo ""
echo "TODO: Fill in actual file paths from step 1 metadata output"
echo "See datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/3_patch_gcp_files.sh for reference"
echo ""
echo "Files that likely need decompression:"
echo "  - seqspec YAML files"
echo "  - barcode_onlist TSV"
echo "  - guide_design CSV/TSV"
