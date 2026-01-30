#!/bin/bash
#
# Decompress .tsv.gz files on GCS to a patch directory
#
# Some pipeline components require uncompressed TSV files for barcode_onlist
# and guide_design. This script decompresses them and creates a patched
# sample metadata file with updated paths.
#
# Usage:
#   ./patch_files.sh
#
# After running, use sample_metadata_gcp_*_patched.csv for the pipeline.
#
set -e

# =============================================================================
# CONFIGURATION - Update these for your dataset
# =============================================================================

BUCKET="gs://igvf-pertub-seq-pipeline-data"
DATASET="TEMPLATE_DATASET"  # <-- CHANGE THIS
DATE="YYYY_MM_DD"           # <-- CHANGE THIS (match upload date)

# Source files (gzipped) - update these paths based on your sample_metadata_gcp file
BARCODE_ONLIST_SRC="${BUCKET}/${DATASET}/${DATE}/path/to/IGVFFIXXXXXX.tsv.gz"  # <-- CHANGE THIS
GUIDE_DESIGN_SRC="${BUCKET}/${DATASET}/${DATE}/path/to/IGVFFIXXXXXX.tsv.gz"    # <-- CHANGE THIS

# =============================================================================
# DERIVED PATHS (no need to change)
# =============================================================================

PATCH_DIR="${BUCKET}/${DATASET}/${DATE}/patch"

# Extract filenames without .gz extension
BARCODE_FILENAME=$(basename "${BARCODE_ONLIST_SRC}" .gz)
GUIDE_FILENAME=$(basename "${GUIDE_DESIGN_SRC}" .gz)

BARCODE_ONLIST_DST="${PATCH_DIR}/${BARCODE_FILENAME}"
GUIDE_DESIGN_DST="${PATCH_DIR}/${GUIDE_FILENAME}"

# =============================================================================
# EXECUTION
# =============================================================================

echo "Decompressing and uploading barcode_onlist file..."
gsutil cat "${BARCODE_ONLIST_SRC}" | gunzip | gsutil cp - "${BARCODE_ONLIST_DST}"

echo "Decompressing and uploading guide_design file..."
gsutil cat "${GUIDE_DESIGN_SRC}" | gunzip | gsutil cp - "${GUIDE_DESIGN_DST}"

echo ""
echo "Done! Files uploaded to:"
echo "  ${BARCODE_ONLIST_DST}"
echo "  ${GUIDE_DESIGN_DST}"
echo ""
echo "Next steps:"
echo "  1. Create patched sample metadata with updated paths"
echo "  2. Replace .tsv.gz paths with .tsv paths in the patch/ directory"
echo "  3. Use the patched CSV for 3_run_CRISPR_pipeline.sh"
