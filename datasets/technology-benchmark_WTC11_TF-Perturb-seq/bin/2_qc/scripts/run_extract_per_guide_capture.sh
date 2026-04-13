#!/bin/bash
# Run per-guide capture extraction for all datasets in a manifest.
#
# Usage:
#   bash run_extract_per_guide_capture.sh <qc_paths.tsv>
#
# Example:
#   bash run_extract_per_guide_capture.sh ../../manifests/cleanser_unified_qc_paths.tsv
#   bash run_extract_per_guide_capture.sh ../../manifests/sceptre_v11_qc_paths.tsv
#
# The script derives the MuData path from each row's qc_dir:
#   qc_dir = .../pipeline_dashboard/additional_qc
#   mudata  = .../pipeline_dashboard/inference_mudata.h5mu
#
# Output is written to: <qc_dir>/guide/per_guide_capture.tsv

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QC_PATHS="${1:?Usage: $0 <qc_paths.tsv>}"

if [[ ! -f "$QC_PATHS" ]]; then
    echo "ERROR: File not found: $QC_PATHS"
    exit 1
fi

echo "Reading manifest: $QC_PATHS"
echo "---"

# Activate the project venv
PROJECT_ROOT="/cellar/users/aklie/projects/tf_perturb_seq"
source "$PROJECT_ROOT/.venv/bin/activate"

# Read the TSV (skip header), columns: dataset qc_dir ...
tail -n +2 "$QC_PATHS" | while IFS=$'\t' read -r dataset qc_dir _rest; do
    # Derive MuData path from qc_dir
    # qc_dir ends in .../pipeline_dashboard/additional_qc
    mudata_path="$(dirname "$qc_dir")/inference_mudata.h5mu"
    output_file="${qc_dir}/guide/per_guide_capture.tsv"

    if [[ ! -f "$mudata_path" ]]; then
        echo "SKIP: MuData not found for $dataset: $mudata_path"
        continue
    fi

    if [[ -f "$output_file" ]]; then
        echo "SKIP: Already exists for $dataset: $output_file"
        continue
    fi

    echo "Processing: $dataset"
    echo "  MuData: $mudata_path"
    echo "  Output: $output_file"

    python "$SCRIPT_DIR/extract_per_guide_capture.py" \
        "$mudata_path" \
        "$output_file" \
        "$dataset"

    echo "  Done!"
    echo "---"
done

echo "All datasets processed!"
