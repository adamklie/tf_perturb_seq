PROJECT_ROOT="/cellar/users/aklie/data/datasets/tf_perturb_seq"
DATASETS_ROOT="${PROJECT_ROOT}/datasets"
SCRIPT="${PROJECT_ROOT}/scripts/run_pathways.sh"
GMT="${PROJECT_ROOT}/ref/gene_sets/c5.go.bp.v2024.1.Hs.symbols.gmt"

DATASETS=(
  Hon_WTC11-benchmark_TF-Perturb-seq
  Huangfu_WTC11-benchmark_TF-Perturb-seq
  Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3
  Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2
  Engreitz_WTC11-benchmark_TF-Perturb-seq
)

RUN_LABELS=(
  cleanser_unified
)

### Run all
for run_label in "${RUN_LABELS[@]}"; do
  for dataset in "${DATASETS[@]}"; do
    cal_dir="${DATASETS_ROOT}/${dataset}/runs/${run_label}/calibration"

    # bash "${SCRIPT}" \
    #   --project-root "${PROJECT_ROOT}" \
    #   --calibrated-results "${cal_dir}/${dataset}_calibrated_all_results.tsv" \
    #   --gene-sets "${GMT}" \
    #   --outdir "${cal_dir}" \
    #   --prefix "${dataset}"
  done
done
