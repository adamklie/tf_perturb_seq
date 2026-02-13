#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Energy Distance Pipeline Runner
#
# This script runs the full energy distance analysis pipeline:
#   1. Convert MuData to pipeline inputs (adapter)
#   2. Preprocess annotations (0_preprocess.py)
#   3. Filter gRNAs (1_filtereing_gRNA.py)
#   4. Compute energy distances (2_e_distance_nontargeting.py)
#   5. Generate plots (2_1_Plot_figure.py)
#   6. Cluster phenotypes (3_e_distance_among_regions.py)
#
# USAGE:
#   run_energy_distance_pipeline.sh \
#     --project-root /path/to/tf_perturb_seq \
#     --input /path/to/inference_mudata.h5mu \
#     --outdir /path/to/output \
#     --run-name dataset_name \
#     [--config /path/to/config.json] \
#     [--batch-key batch] \
#     [--use-harmony] \
#     [--n-top-genes 5000] \
#     [--n-pcs 50] \
#     [--dry-run]
#
###############################################################################

# Default values
PROJECT_ROOT=""
INPUT=""
OUTDIR=""
RUN_NAME=""
CONFIG=""
BATCH_KEY=""
USE_HARMONY=""
N_TOP_GENES=5000
N_PCS=50
DRY_RUN=0

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root)     PROJECT_ROOT="$2"; shift 2 ;;
    --input)            INPUT="$2"; shift 2 ;;
    --outdir)           OUTDIR="$2"; shift 2 ;;
    --run-name)         RUN_NAME="$2"; shift 2 ;;
    --config)           CONFIG="$2"; shift 2 ;;
    --batch-key)        BATCH_KEY="$2"; shift 2 ;;
    --use-harmony)      USE_HARMONY="--use-harmony"; shift 1 ;;
    --n-top-genes)      N_TOP_GENES="$2"; shift 2 ;;
    --n-pcs)            N_PCS="$2"; shift 2 ;;
    --dry-run)          DRY_RUN=1; shift 1 ;;
    -h|--help)
      sed -n '1,40p' "$0"
      exit 0
      ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

# Validate required arguments
if [[ -z "${PROJECT_ROOT}" || -z "${INPUT}" || -z "${OUTDIR}" || -z "${RUN_NAME}" ]]; then
  echo "ERROR: Missing required arguments." >&2
  echo "Required: --project-root, --input, --outdir, --run-name" >&2
  exit 1
fi

[[ -f "${INPUT}" ]] || { echo "ERROR: Input file not found: ${INPUT}" >&2; exit 1; }

###############################################################################
# Environment Setup
###############################################################################

VENV="${PROJECT_ROOT}/.venv/bin/activate"
[[ -f "${VENV}" ]] || { echo "ERROR: venv not found: ${VENV}" >&2; exit 1; }
# shellcheck disable=SC1090
source "${VENV}"

# Paths to scripts
ADAPTER_SCRIPT="${PROJECT_ROOT}/src/tf_perturb_seq/energy_distance/mudata_adapter.py"
EDIST_PIPELINE="${PROJECT_ROOT}/external/energy_dist_pipeline/bin"

[[ -f "${ADAPTER_SCRIPT}" ]] || { echo "ERROR: Adapter script not found: ${ADAPTER_SCRIPT}" >&2; exit 1; }
[[ -d "${EDIST_PIPELINE}" ]] || { echo "ERROR: Energy distance pipeline not found: ${EDIST_PIPELINE}" >&2; echo "Run: git submodule add https://github.com/Chikara-Takeuchi/energy_dist_pipeline.git external/energy_dist_pipeline" >&2; exit 1; }

###############################################################################
# Helper Functions
###############################################################################

run_cmd() {
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    printf '[DRY-RUN] %q' "$1"
    shift
    for arg in "$@"; do printf ' %q' "${arg}"; done
    printf '\n'
  else
    "$@"
  fi
}

###############################################################################
# Output Layout
###############################################################################

DIR_INPUTS="${OUTDIR}/inputs"
DIR_FILTER="${OUTDIR}/filtering"
DIR_EDIST="${OUTDIR}/energy_distance"
DIR_CLUSTER="${OUTDIR}/clustering"
DIR_PLOTS="${OUTDIR}/plots"

echo "============================================="
echo " Energy Distance Pipeline"
echo " Project:    ${PROJECT_ROOT}"
echo " Input:      ${INPUT}"
echo " Output:     ${OUTDIR}"
echo " Run name:   ${RUN_NAME}"
echo " Dry-run:    ${DRY_RUN}"
echo "============================================="

run_cmd mkdir -p "${DIR_INPUTS}" "${DIR_FILTER}" "${DIR_EDIST}" "${DIR_CLUSTER}" "${DIR_PLOTS}"

###############################################################################
# Step 1: Convert MuData to Pipeline Inputs
###############################################################################

echo ""
echo "[1/6] Converting MuData to pipeline inputs..."

ADAPTER_ARGS=(
  python "${ADAPTER_SCRIPT}"
  --input "${INPUT}"
  --outdir "${DIR_INPUTS}"
  --n-top-genes "${N_TOP_GENES}"
  --n-pcs "${N_PCS}"
)

if [[ -n "${BATCH_KEY}" ]]; then
  ADAPTER_ARGS+=(--batch-key "${BATCH_KEY}")
fi

if [[ -n "${USE_HARMONY}" ]]; then
  ADAPTER_ARGS+=(${USE_HARMONY})
fi

run_cmd "${ADAPTER_ARGS[@]}"

# Paths to generated inputs
ADATA="${DIR_INPUTS}/adata.h5ad"
GRNA_DICT="${DIR_INPUTS}/grna_dict.pkl"
ANNOTATION="${DIR_INPUTS}/annotation.csv"

###############################################################################
# Step 2: Create Pipeline Config
###############################################################################

echo ""
echo "[2/6] Creating pipeline configuration..."

# Use provided config or create default
if [[ -n "${CONFIG}" && -f "${CONFIG}" ]]; then
  CONFIG_FILE="${CONFIG}"
  echo "Using provided config: ${CONFIG_FILE}"
else
  CONFIG_FILE="${OUTDIR}/config.json"
  run_cmd bash -c "cat > '${CONFIG_FILE}' << 'EOF'
{
  \"output_file_name_list\": {
    \"OUTPUT_FOLDER\": \"${DIR_EDIST}\",
    \"targeting_outlier_table\": \"targeting_outlier_table.csv\",
    \"non_targeting_outlier_table\": \"non_targeting_outlier_table.csv\",
    \"edist_pvalue_table\": \"pval_edist_full.csv\",
    \"edist_target_by_target_matrix\": \"target_by_target_matrix.csv\",
    \"edist_embedding_info\": \"edist_embedding_info.csv\",
    \"pca_table\": \"pca_dataframe.pickle\",
    \"gRNA_dict\": \"gRNA_dictionary.pickle\",
    \"discordance_gRNA_table\": \"discordance_gRNA.csv\",
    \"OVERWRITE_PCA_DICT\": true
  },
  \"input_data\": {
    \"annotation_file\": {
      \"file_path\": \"${ANNOTATION}\",
      \"concatenate_key\": \"intended_target_name\"
    },
    \"h5ad_file\": {
      \"file_path\": \"${ADATA}\",
      \"obsm_key\": \"X_pca\"
    },
    \"sgRNA_file\": {
      \"file_path\": \"${GRNA_DICT}\"
    }
  },
  \"gRNA_filtering\": {
    \"perform_targeting_filtering\": true,
    \"perform_nontargeting_filtering\": true,
    \"threshold_gRNA_num\": 6,
    \"combi_count\": 4,
    \"total_permute_disco\": 1000,
    \"combi_cell_num_max\": 1000,
    \"batch_num_basic\": 120
  },
  \"permutation_test\": {
    \"permute_per_bg\": 1000,
    \"num_of_bg\": 20,
    \"non_target_pick\": 2000,
    \"target_cell_num_max\": 2000,
    \"batch_num_basic\": 200,
    \"use_matched_bg\": false
  },
  \"aggregate\": {
    \"downsampling_maximum\": 10000
  }
}
EOF"
  echo "Created default config: ${CONFIG_FILE}"
fi

# Clustering config
CLUSTER_CONFIG="${OUTDIR}/config_clustering.json"
run_cmd bash -c "cat > '${CLUSTER_CONFIG}' << 'EOF'
{
  \"cutoff\": {
    \"pval_cutoff\": 0.05,
    \"distance_cutoff\": 0
  },
  \"clustering\": {
    \"method\": \"Affinity\"
  }
}
EOF"

###############################################################################
# Step 3: Run gRNA Filtering
###############################################################################

echo ""
echo "[3/6] Running gRNA filtering..."

if [[ "${DRY_RUN}" -eq 0 ]]; then
  python "${EDIST_PIPELINE}/1_filtereing_gRNA.py" "${CONFIG_FILE}" 2>&1 | tee "${DIR_FILTER}/filtering.log" || {
    echo "WARNING: gRNA filtering step had issues. Check ${DIR_FILTER}/filtering.log"
  }
else
  echo "[DRY-RUN] python ${EDIST_PIPELINE}/1_filtereing_gRNA.py ${CONFIG_FILE}"
fi

###############################################################################
# Step 4: Run Energy Distance Computation
###############################################################################

echo ""
echo "[4/6] Computing energy distances..."

if [[ "${DRY_RUN}" -eq 0 ]]; then
  python "${EDIST_PIPELINE}/2_e_distance_nontargeting.py" "${CONFIG_FILE}" 2>&1 | tee "${DIR_EDIST}/edist.log" || {
    echo "WARNING: Energy distance step had issues. Check ${DIR_EDIST}/edist.log"
  }
else
  echo "[DRY-RUN] python ${EDIST_PIPELINE}/2_e_distance_nontargeting.py ${CONFIG_FILE}"
fi

###############################################################################
# Step 5: Generate Plots
###############################################################################

echo ""
echo "[5/6] Generating plots..."

if [[ "${DRY_RUN}" -eq 0 ]]; then
  python "${EDIST_PIPELINE}/2_1_Plot_figure.py" "${CONFIG_FILE}" 2>&1 | tee "${DIR_PLOTS}/plots.log" || {
    echo "WARNING: Plot generation had issues. Check ${DIR_PLOTS}/plots.log"
  }
else
  echo "[DRY-RUN] python ${EDIST_PIPELINE}/2_1_Plot_figure.py ${CONFIG_FILE}"
fi

###############################################################################
# Step 6: Run Phenotype Clustering
###############################################################################

echo ""
echo "[6/6] Running phenotype clustering..."

if [[ "${DRY_RUN}" -eq 0 ]]; then
  python "${EDIST_PIPELINE}/3_e_distance_among_regions.py" "${CONFIG_FILE}" "${CLUSTER_CONFIG}" 2>&1 | tee "${DIR_CLUSTER}/clustering.log" || {
    echo "WARNING: Clustering step had issues. Check ${DIR_CLUSTER}/clustering.log"
  }
else
  echo "[DRY-RUN] python ${EDIST_PIPELINE}/3_e_distance_among_regions.py ${CONFIG_FILE} ${CLUSTER_CONFIG}"
fi

###############################################################################
# Done
###############################################################################

echo ""
echo "============================================="
if [[ "${DRY_RUN}" -eq 1 ]]; then
  echo " Dry-run complete (no files written)."
else
  echo " Energy distance pipeline complete!"
  echo ""
  echo " Outputs:"
  echo "   Inputs:     ${DIR_INPUTS}"
  echo "   Filtering:  ${DIR_FILTER}"
  echo "   E-distance: ${DIR_EDIST}"
  echo "   Clustering: ${DIR_CLUSTER}"
  echo "   Plots:      ${DIR_PLOTS}"
fi
echo "============================================="
