PROJECT_ROOT="/cellar/users/aklie/projects/tf_perturb_seq"
DATASETS_ROOT="${PROJECT_ROOT}/datasets"
BENCHMARK_DIR="${DATASETS_ROOT}/technology-benchmark_WTC11_TF-Perturb-seq"
SCRIPT="${PROJECT_ROOT}/scripts/run_calibration.sh"

DATASETS=(
  Hon_WTC11-benchmark_TF-Perturb-seq
  Huangfu_WTC11-benchmark_TF-Perturb-seq
  Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3
  Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2
  Engreitz_WTC11-benchmark_TF-Perturb-seq
)

RUN_LABELS=(
  cleanser_unified
  sceptre_v11
)

### Run all
for run_label in "${RUN_LABELS[@]}"; do
  for dataset in "${DATASETS[@]}"; do
    run_dir="${DATASETS_ROOT}/${dataset}/runs/${run_label}"
    mudata="${run_dir}/pipeline_dashboard/inference_mudata.h5mu"
    [[ -f "${mudata}" ]] || mudata="${run_dir}/pipeline_outputs/inference_mudata.h5mu"

    # bash "${SCRIPT}" \
    #   --project-root "${PROJECT_ROOT}" \
    #   --trans-results "${run_dir}/pipeline_outputs/perturbo_trans_per_element_output.tsv.gz" \
    #   --mudata "${mudata}" \
    #   --outdir "${run_dir}/calibration" \
    #   --prefix "${dataset}"
  done
done
