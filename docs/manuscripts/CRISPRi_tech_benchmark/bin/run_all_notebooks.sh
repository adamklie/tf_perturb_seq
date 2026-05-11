#!/usr/bin/env bash
set -uo pipefail

# Run all QC and inference notebooks for a given RUN_LABEL
# Usage: bash run_all_notebooks.sh <run_label>
#   e.g.: bash run_all_notebooks.sh cleanser_unified
#         bash run_all_notebooks.sh sceptre_v11

PROJECT_ROOT="/cellar/users/aklie/projects/tf_perturb_seq"
BENCHMARK_DIR="${PROJECT_ROOT}/datasets/technology-benchmark_WTC11_TF-Perturb-seq"
VENV="${PROJECT_ROOT}/.venv/bin/activate"

RUN_LABEL="${1:?Usage: $0 <run_label>}"

source "${VENV}"

OUTPUT_DIR="${BENCHMARK_DIR}/results/cross_tech_comparison/${RUN_LABEL}/executed_notebooks"
mkdir -p "${OUTPUT_DIR}"

run_notebook() {
    local nb_path="$1"
    local nb_name
    nb_name="$(basename "${nb_path}" .ipynb)"
    local out_path="${OUTPUT_DIR}/${nb_name}_${RUN_LABEL}.ipynb"

    echo "=========================================="
    echo "Running: ${nb_name} (${RUN_LABEL})"
    echo "=========================================="

    # Inject RUN_LABEL via papermill parameter
    if python -m papermill "${nb_path}" "${out_path}" \
        -p RUN_LABEL "${RUN_LABEL}" \
        --no-progress-bar \
        2>&1 | tail -20; then
        echo "Done: ${out_path}"
    else
        echo "FAILED: ${nb_name} (${RUN_LABEL}) — continuing with remaining notebooks"
    fi
    echo ""
}

echo "============================================="
echo " Running all notebooks for: ${RUN_LABEL}"
echo "============================================="

# 2_qc notebooks
for nb in \
    "${BENCHMARK_DIR}/bin/2_qc/1_metrics_comparison.ipynb" \
    "${BENCHMARK_DIR}/bin/2_qc/2_repression_efficiency_comparison.ipynb" \
    "${BENCHMARK_DIR}/bin/2_qc/3_guide_capture_comparison.ipynb" \
; do
    run_notebook "${nb}"
done

# 3_inference notebooks (only if calibration data exists)
CAL_CHECK="${PROJECT_ROOT}/datasets/Hon_WTC11-benchmark_TF-Perturb-seq/runs/${RUN_LABEL}/calibration"
if [[ -d "${CAL_CHECK}" ]] && ls "${CAL_CHECK}"/*.tsv &>/dev/null; then
    echo "Calibration data found for ${RUN_LABEL} — running inference notebooks"
    for nb in \
        "${BENCHMARK_DIR}/bin/3_inference/3_calibration_summary.ipynb" \
        "${BENCHMARK_DIR}/bin/3_inference/4_trans_comparison.ipynb" \
        "${BENCHMARK_DIR}/bin/3_inference/5_tf_enrichment_comparison.ipynb" \
    ; do
        run_notebook "${nb}"
    done
else
    echo "SKIP: No calibration data for ${RUN_LABEL} — skipping 3_inference notebooks"
fi

echo "============================================="
echo " All notebooks complete for: ${RUN_LABEL}"
echo " Output: ${OUTPUT_DIR}"
echo "============================================="
