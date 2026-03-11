#!/bin/bash
# Run per-guide capture extraction for all datasets
# Execute this on the cluster where MuData files are accessible

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MUDATA_PATHS="/Users/adamklie/Desktop/projects/tf_perturb_seq/datasets/technology-benchmark_WTC11_TF-Perturb-seq/latest_mudata_paths.tsv"

# Read the TSV and process each dataset
tail -n +2 "$MUDATA_PATHS" | while IFS=$'\t' read -r input outdir run_name; do
    echo "Processing: $run_name"

    output_file="${outdir}/per_guide_capture.tsv"

    python "$SCRIPT_DIR/extract_per_guide_capture.py" \
        "$input" \
        "$output_file" \
        "$run_name"

    echo "Done: $output_file"
    echo "---"
done

echo "All datasets processed!"
