#!/usr/bin/env bash
set -euo pipefail

###############################################
# sync_benchmark.sh
#
# Sync GCP pipeline outputs for one benchmark run and generate manifests:
#
#   1. Check GCP for available pipeline_dashboard/ per dataset
#   2. Download pipeline_dashboard/ (preferred) or pipeline_outputs/ (fallback)
#   3. Generate manifests/
#        <run-label>_qc_input.tsv   — for 2_qc/scripts/qc_array.sh (fallback datasets only)
#        <run-label>_qc_paths.tsv   — for 2_qc/ comparison notebooks
#
# USAGE:
#   cd datasets/technology-benchmark_WTC1TF-Perturb-seq/bin/get_data
#
#   sync_benchmark.sh \
#     --gcp-prefix gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/FirstBenchmarkRun_sceptre_v11 \
#     --run-label sceptre_v11
#
#   sync_benchmark.sh \
#     --gcp-prefix gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/Benchmark_cleanser_unified \
#     --run-label cleanser_unified
#
#   # Dry-run (see what would happen)
#   sync_benchmark.sh --gcp-prefix gs://... --run-label sceptre_v11 --dry-run
#
#   # Skip download (already synced), regenerate manifests only
#   sync_benchmark.sh --gcp-prefix gs://... --run-label sceptre_v11 --skip-download
###############################################

###############################################
# Dataset mapping: GCP name → local dir name
# Format: "gcp_name:local_name:dataset_label"
# For Gersbach HTv2: GCP name varies by run (Gersbach_HTv2 vs Gersbach_HTV2).
# The alt name is tried if the primary isn't found.
###############################################
declare -a DATASET_MAP=(
  "GaryHonDataset:Hon_WTC11-benchmark_TF-Perturb-seq:Hon_WTC11-benchmark_TF-Perturb-seq"
  "HuangFuDataset:Huangfu_WTC11-benchmark_TF-Perturb-seq:Huangfu_WTC11-benchmark_TF-Perturb-seq"
  "Gersbach_gemX_v3:Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3:Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3"
  "Gersbach_HTv2:Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2:Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2"
  "Engreitz_ccPerturb:Engreitz_WTC11-benchmark_TF-Perturb-seq:Engreitz_WTC11-benchmark_TF-Perturb-seq"
)

# Alternate GCP names to try if primary is not found
declare -A DATASET_ALT_NAMES=(
  ["Gersbach_HTv2"]="Gersbach_HTV2"
  ["Gersbach_HTV2"]="Gersbach_HTv2"
)

###############################################
# Paths
###############################################
PROJECT_ROOT="/cellar/users/aklie/projects/tf_perturb_seq"
DATASETS_ROOT="${PROJECT_ROOT}/datasets"
BENCHMARK_DIR="${DATASETS_ROOT}/technology-benchmark_WTC11_TF-Perturb-seq"

SYNC_SCRIPT="${PROJECT_ROOT}/scripts/sync_gcp_run.sh"
MANIFESTS_DIR="${BENCHMARK_DIR}/manifests"

###############################################
# Args
###############################################
GCP_PREFIX=""
RUN_LABEL=""
SKIP_DOWNLOAD=0
DRY_RUN=0
FORCE_DOWNLOAD=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gcp-prefix)     GCP_PREFIX="$2";   shift 2 ;;
    --run-label)      RUN_LABEL="$2";    shift 2 ;;
    --skip-download)  SKIP_DOWNLOAD=1;   shift 1 ;;
    --force-download) FORCE_DOWNLOAD=1;  shift 1 ;;
    --dry-run)        DRY_RUN=1;         shift 1 ;;
    -h|--help)
      sed -n '1,50p' "$0"
      exit 0
      ;;
    *)
      echo "ERROR: Unknown arg: $1" >&2
      exit 1
      ;;
  esac
done

[[ -n "${GCP_PREFIX}" && -n "${RUN_LABEL}" ]] || {
  echo "ERROR: --gcp-prefix and --run-label are required." >&2
  exit 1
}

GCP_PREFIX="${GCP_PREFIX%/}"

###############################################
# Validate dependencies
###############################################
[[ -f "${SYNC_SCRIPT}" ]] || { echo "ERROR: sync script not found: ${SYNC_SCRIPT}" >&2; exit 1; }

echo "============================================="
echo " sync_benchmark.sh"
echo " GCP prefix:    ${GCP_PREFIX}"
echo " Run label:     ${RUN_LABEL}"
echo " Skip download: ${SKIP_DOWNLOAD}"
echo " Dry-run:       ${DRY_RUN}"
echo "============================================="

if [[ "${DRY_RUN}" -eq 0 ]]; then
  mkdir -p "${MANIFESTS_DIR}"
fi

###############################################
# Step 1: Check GCP availability
#
# Preferred: pipeline_dashboard/  (has pre-computed additional_qc/ + h5mu + plots)
# Fallback:  pipeline_outputs/    (h5mu only; run 2_qc/scripts/qc_array.sh separately)
###############################################
echo ""
echo "--- Step 1/3: Checking GCP availability ---"

declare -a AVAILABLE_DASHBOARD=()   # entries: "dataset_label:gcp_name:local_name"
declare -a AVAILABLE_OUTPUTS=()
declare -a UNAVAILABLE=()

for entry in "${DATASET_MAP[@]}"; do
  IFS=':' read -r gcp_name local_name dataset_label <<< "${entry}"

  # Try primary GCP name, then alternate if it exists
  resolved_gcp_name="${gcp_name}"
  gcp_dashboard="${GCP_PREFIX}/${resolved_gcp_name}/pipeline_dashboard"
  gcp_mudata="${GCP_PREFIX}/${resolved_gcp_name}/pipeline_outputs/inference_mudata.h5mu"

  if ! gsutil ls "${gcp_dashboard}/" &>/dev/null && ! gsutil ls "${gcp_mudata}" &>/dev/null; then
    # Try alternate name if available
    alt_name="${DATASET_ALT_NAMES[${gcp_name}]:-}"
    if [[ -n "${alt_name}" ]]; then
      gcp_dashboard="${GCP_PREFIX}/${alt_name}/pipeline_dashboard"
      gcp_mudata="${GCP_PREFIX}/${alt_name}/pipeline_outputs/inference_mudata.h5mu"
      if gsutil ls "${gcp_dashboard}/" &>/dev/null || gsutil ls "${gcp_mudata}" &>/dev/null; then
        resolved_gcp_name="${alt_name}"
        echo "[NOTE]      Using alternate GCP name: ${alt_name} (instead of ${gcp_name})"
      fi
    fi
  fi

  gcp_dashboard="${GCP_PREFIX}/${resolved_gcp_name}/pipeline_dashboard"
  gcp_mudata="${GCP_PREFIX}/${resolved_gcp_name}/pipeline_outputs/inference_mudata.h5mu"

  if gsutil ls "${gcp_dashboard}/" &>/dev/null; then
    AVAILABLE_DASHBOARD+=("${dataset_label}:${resolved_gcp_name}:${local_name}")
    echo "[DASHBOARD] ${resolved_gcp_name} → ${dataset_label}"
  elif gsutil ls "${gcp_mudata}" &>/dev/null; then
    AVAILABLE_OUTPUTS+=("${dataset_label}:${resolved_gcp_name}:${local_name}")
    echo "[OUTPUTS]   ${resolved_gcp_name} → ${dataset_label}  (no dashboard; run qc_array.sh after sync)"
  else
    UNAVAILABLE+=("${dataset_label}")
    echo "[MISSING]   ${gcp_name} — not ready on GCP (skipping)"
  fi
done

N_DASHBOARD="${#AVAILABLE_DASHBOARD[@]}"
N_OUTPUTS="${#AVAILABLE_OUTPUTS[@]}"
N_UNAVAILABLE="${#UNAVAILABLE[@]}"

echo ""
echo "  Dashboard (${N_DASHBOARD}): $(printf '%s ' "${AVAILABLE_DASHBOARD[@]+"${AVAILABLE_DASHBOARD[@]}"}" | sed 's/:[^:]*:[^ ]*//g')"
echo "  Outputs   (${N_OUTPUTS}):   $(printf '%s ' "${AVAILABLE_OUTPUTS[@]+"${AVAILABLE_OUTPUTS[@]}"}" | sed 's/:[^:]*:[^ ]*//g')"
[[ "${N_UNAVAILABLE}" -gt 0 ]] && \
  echo "  Missing   (${N_UNAVAILABLE}):  ${UNAVAILABLE[*]}"

###############################################
# Step 2: Download
###############################################
if [[ "${SKIP_DOWNLOAD}" -eq 0 ]]; then
  echo ""
  echo "--- Step 2/3: Downloading ---"

  force_flag=""; [[ "${FORCE_DOWNLOAD}" -eq 1 ]] && force_flag="--force"
  dry_flag="";   [[ "${DRY_RUN}" -eq 1 ]]        && dry_flag="--dry-run"

  for entry in "${AVAILABLE_DASHBOARD[@]+"${AVAILABLE_DASHBOARD[@]}"}"; do
    IFS=':' read -r dataset_label gcp_name local_name <<< "${entry}"
    local_run_dir="${DATASETS_ROOT}/${local_name}/runs/${RUN_LABEL}"
    echo ""
    echo "Syncing ${gcp_name} (dashboard) → ${local_run_dir}"
    bash "${SYNC_SCRIPT}" \
      --gcp-dataset-uri "${GCP_PREFIX}/${gcp_name}" \
      --local-dir "${local_run_dir}" \
      --mode dashboard \
      ${force_flag} ${dry_flag}
  done

  for entry in "${AVAILABLE_OUTPUTS[@]+"${AVAILABLE_OUTPUTS[@]}"}"; do
    IFS=':' read -r dataset_label gcp_name local_name <<< "${entry}"
    local_run_dir="${DATASETS_ROOT}/${local_name}/runs/${RUN_LABEL}"
    echo ""
    echo "Syncing ${gcp_name} (outputs fallback) → ${local_run_dir}"
    bash "${SYNC_SCRIPT}" \
      --gcp-dataset-uri "${GCP_PREFIX}/${gcp_name}" \
      --local-dir "${local_run_dir}" \
      --mode outputs \
      ${force_flag} ${dry_flag}
  done
else
  echo ""
  echo "--- Step 2/3: Download skipped (--skip-download) ---"
fi

###############################################
# Step 3: Generate manifests
#
# qc_input.tsv — only fallback datasets (those without pipeline_dashboard/)
#                used by 2_qc/scripts/qc_array.sh
# qc_paths.tsv — ALL available datasets
#   dashboard datasets → pipeline_dashboard/additional_qc/
#   fallback datasets  → pipeline_qc/  (populated after running qc_array.sh)
###############################################
echo ""
echo "--- Step 3/3: Generating manifests ---"

QC_INPUT_TSV="${MANIFESTS_DIR}/${RUN_LABEL}_qc_input.tsv"
QC_PATHS_TSV="${MANIFESTS_DIR}/${RUN_LABEL}_qc_paths.tsv"

if [[ "${DRY_RUN}" -eq 1 ]]; then
  echo "[DRY-RUN] Would write: ${QC_INPUT_TSV}"
  echo "[DRY-RUN] Would write: ${QC_PATHS_TSV}"
else
  # ---- qc_input.tsv ----
  printf "input\toutdir\trun_name\n" > "${QC_INPUT_TSV}"
  for entry in "${AVAILABLE_OUTPUTS[@]+"${AVAILABLE_OUTPUTS[@]}"}"; do
    IFS=':' read -r dataset_label gcp_name local_name <<< "${entry}"
    local_run_dir="${DATASETS_ROOT}/${local_name}/runs/${RUN_LABEL}"
    printf "%s\t%s\t%s\n" \
      "${local_run_dir}/pipeline_outputs/inference_mudata.h5mu" \
      "${local_run_dir}/pipeline_qc" \
      "${dataset_label}" \
      >> "${QC_INPUT_TSV}"
  done
  echo "Written: ${QC_INPUT_TSV}"

  # ---- qc_paths.tsv ----
  printf "dataset\tqc_dir\tgene_metrics\tguide_metrics\tintended_target_results\tintended_target_metrics\ttrans_results\ttrans_metrics\n" \
    > "${QC_PATHS_TSV}"

  # Dashboard datasets → pipeline_dashboard/additional_qc/
  for entry in "${AVAILABLE_DASHBOARD[@]+"${AVAILABLE_DASHBOARD[@]}"}"; do
    IFS=':' read -r dataset_label gcp_name local_name <<< "${entry}"
    qc_dir="${DATASETS_ROOT}/${local_name}/runs/${RUN_LABEL}/pipeline_dashboard/additional_qc"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${dataset_label}" \
      "${qc_dir}" \
      "gene/gene_metrics.tsv" \
      "guide/guide_metrics.tsv" \
      "intended_target/intended_target_results.tsv" \
      "intended_target/intended_target_metrics.tsv" \
      "trans/trans_results.tsv" \
      "trans/trans_metrics.tsv" \
      >> "${QC_PATHS_TSV}"
  done

  # Fallback datasets → pipeline_qc/ (populated after running qc_array.sh)
  for entry in "${AVAILABLE_OUTPUTS[@]+"${AVAILABLE_OUTPUTS[@]}"}"; do
    IFS=':' read -r dataset_label gcp_name local_name <<< "${entry}"
    qc_dir="${DATASETS_ROOT}/${local_name}/runs/${RUN_LABEL}/pipeline_qc"
    p="${dataset_label}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${dataset_label}" \
      "${qc_dir}" \
      "mapping_gene/${p}_gene_metrics.tsv" \
      "mapping_guide/${p}_guide_metrics.tsv" \
      "intended_target/${p}_intended_target_results.tsv" \
      "intended_target/${p}_intended_target_metrics.tsv" \
      "trans/${p}_trans_results.tsv" \
      "trans/${p}_trans_metrics.tsv" \
      >> "${QC_PATHS_TSV}"
  done
  echo "Written: ${QC_PATHS_TSV}"

  if [[ "${N_OUTPUTS}" -gt 0 ]]; then
    echo ""
    echo "NOTE: ${N_OUTPUTS} dataset(s) lack pipeline_dashboard/ and need local QC."
    echo "      Run: cd .2_qc && sbatch --array=1-${N_OUTPUTS} scripts/qc_array.sh ${QC_INPUT_TSV} ${PROJECT_ROOT}"
  fi
fi

echo ""
echo "============================================="
echo " Done."
echo " Run label:  ${RUN_LABEL}"
echo " Manifests:  ${MANIFESTS_DIR}/"
echo "============================================="
