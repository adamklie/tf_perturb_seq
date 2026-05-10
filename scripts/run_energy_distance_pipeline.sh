#!/usr/bin/env bash
###############################################################################
# Run the TF Perturb-seq energy distance pipeline (steps 1, 2, 2.1) on a
# single dataset's inference_mudata.h5mu — UCSD nrnb edition.
#
# Adapted from Sara's Duke HPC template, which itself is based on the canonical
# upstream wrapper:
#   https://github.com/Chikara-Takeuchi/energy_dist_TFperturb
#     -> tracked locally as: external/energy_dist_TFperturb (submodule)
#
# The pipeline bin scripts come from:
#   https://github.com/Chikara-Takeuchi/energy_dist_pipeline
#     -> tracked locally as: external/energy_dist_pipeline (submodule, pinned)
#
# Differences vs the upstream wrapper:
#   - Reads MuData from local disk / GCS / Synapse (upstream is Synapse-only)
#   - Uses our local scripts/preprocess_mudata_local.py (the Synapse-coupled
#     upstream preprocess_mudata.py is not used here; see external/.../preprocess_mudata.py
#     if you want the canonical Synapse-source version)
#   - Uses pinned submodule paths instead of cloning bin/ scripts at runtime
#   - Idempotent: skips download / preprocess / clone / container pull when
#     outputs already exist
#   - SLURM-friendly: assumes /cellar bind mount on UCSD nrnb
#
# Step 3 (target-by-target matrix + clustering) is intentionally SKIPPED here —
# run separately after picking cutoffs in config3.json.
#
# Inputs (one of):
#   --mudata-path      /path/to/inference_mudata.h5mu   (already-downloaded file)
#   --gcs-mudata-path  gs://.../inference_mudata.h5mu   (downloads to OUTPUT_DIR)
#   --synapse-id       synXXXXX                          (downloads via synapseclient)
# And:
#   --output-dir       /path/to/where/results/should/go
#
# Optional:
#   --container-path   Path to .sif (default: shared sif under /cellar/users/aklie/opt/containers)
#   --pipeline-bin     Override path to energy_dist_pipeline/bin
#                      (default: external/energy_dist_pipeline/bin from this repo)
#   --skip-step3 / --run-step3   Step 3 default = SKIP
#
# Designed to be called from a SLURM job (datasets/<dataset>/5_run_energy_distance.sh).
###############################################################################

set -euo pipefail

# Resolve repo root from this script's location (scripts/ → repo)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

CONTAINER_PATH="${CONTAINER_PATH:-/cellar/users/aklie/opt/containers/edist_pipeline.sif}"
PIPELINE_BIN="${PIPELINE_BIN:-${REPO_ROOT}/external/energy_dist_pipeline/bin}"
PREPROCESS_LOCAL="${REPO_ROOT}/scripts/preprocess_mudata_local.py"
MUDATA_PATH=""
GCS_MUDATA_PATH=""
SYNAPSE_ID=""
OUTPUT_DIR=""
SKIP_STEP3=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mudata-path)     MUDATA_PATH="$2"; shift 2 ;;
    --gcs-mudata-path) GCS_MUDATA_PATH="$2"; shift 2 ;;
    --synapse-id)      SYNAPSE_ID="$2"; shift 2 ;;
    --output-dir)      OUTPUT_DIR="$2"; shift 2 ;;
    --container-path)  CONTAINER_PATH="$2"; shift 2 ;;
    --pipeline-bin)    PIPELINE_BIN="$2"; shift 2 ;;
    --skip-step3)      SKIP_STEP3=1; shift ;;
    --run-step3)       SKIP_STEP3=0; shift ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

[[ -z "$OUTPUT_DIR" ]] && { echo "ERROR: --output-dir is required" >&2; exit 2; }
if [[ -z "$MUDATA_PATH" && -z "$GCS_MUDATA_PATH" && -z "$SYNAPSE_ID" ]]; then
  echo "ERROR: one of --mudata-path, --gcs-mudata-path, or --synapse-id is required" >&2; exit 2
fi
[[ ! -f "$PREPROCESS_LOCAL" ]] && { echo "ERROR: preprocess script missing: $PREPROCESS_LOCAL" >&2; exit 2; }
[[ ! -d "$PIPELINE_BIN" ]] && {
  echo "ERROR: pipeline bin dir missing: $PIPELINE_BIN" >&2
  echo "Did you 'git submodule update --init external/energy_dist_pipeline'?" >&2
  exit 2
}

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
OUTPUT_DIR_ABS="$(pwd)"

###############################################################################
# 0. Resolve MuData path (download from GCS / Synapse if needed)
###############################################################################

if [[ -z "$MUDATA_PATH" ]]; then
  MUDATA_LOCAL="${OUTPUT_DIR_ABS}/inference_mudata.h5mu"
  if [[ -f "$MUDATA_LOCAL" ]]; then
    echo "[setup] reusing already-downloaded MuData at $MUDATA_LOCAL"
  elif [[ -n "$GCS_MUDATA_PATH" ]]; then
    echo "[setup] downloading MuData from $GCS_MUDATA_PATH -> $MUDATA_LOCAL"
    gcloud storage cp "$GCS_MUDATA_PATH" "$MUDATA_LOCAL"
  elif [[ -n "$SYNAPSE_ID" ]]; then
    echo "[setup] downloading MuData from Synapse $SYNAPSE_ID -> $MUDATA_LOCAL"
    [[ -z "${SYNAPSE_AUTH_TOKEN:-}" ]] && { echo "ERROR: SYNAPSE_AUTH_TOKEN not set" >&2; exit 2; }
    # Use the project venv explicitly — bare `python` on a SLURM compute node
    # may resolve to base conda without synapseclient.
    /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/python -c "
import os, shutil, synapseclient
syn = synapseclient.Synapse(silent=True)
syn.login(authToken=os.environ['SYNAPSE_AUTH_TOKEN'])
e = syn.get('${SYNAPSE_ID}', downloadLocation='${OUTPUT_DIR_ABS}')
src, dst = e.path, '${MUDATA_LOCAL}'
if src != dst: shutil.move(src, dst)
print(f'fetched {e.id} -> {dst}')
"
  fi
  MUDATA_PATH="$MUDATA_LOCAL"
fi

[[ ! -f "$MUDATA_PATH" ]] && { echo "ERROR: MuData not found: $MUDATA_PATH" >&2; exit 2; }
echo "[setup] MuData: $MUDATA_PATH ($(du -h "$MUDATA_PATH" | awk '{print $1}'))"

###############################################################################
# 1. Container
###############################################################################

if [[ ! -f "$CONTAINER_PATH" ]]; then
  echo "[setup] container not at $CONTAINER_PATH; pulling..."
  mkdir -p "$(dirname "$CONTAINER_PATH")"
  apptainer pull "$CONTAINER_PATH" docker://docker.io/takechikara/energy_distance_env:latest
fi
echo "[setup] container: $CONTAINER_PATH"
echo "[setup] pipeline bin: $PIPELINE_BIN"
echo "[setup] preprocess script: $PREPROCESS_LOCAL"

###############################################################################
# 2. Generate run-specific configs (config1_2.json drives steps 0/1/2;
#    config3.json drives step 3 cutoffs)
###############################################################################

CONFIG12="${OUTPUT_DIR_ABS}/config1_2.json"
CONFIG3="${OUTPUT_DIR_ABS}/config3.json"

cat > "$CONFIG12" <<JSON
{
  "output_file_name_list": {
    "OUTPUT_FOLDER": "${OUTPUT_DIR_ABS}",
    "targeting_outlier_table": "targeting_outlier_table.csv",
    "non_targeting_outlier_table": "non_targeting_outlier_table.csv",
    "edist_pvalue_table": "pval_edist_full.csv",
    "edist_target_by_target_matrix": "target_by_target_matrix.csv",
    "edist_embedding_info": "edist_embedding_info.csv",
    "pca_table": "pca_dataframe.pickle",
    "gRNA_dict": "gRNA_dict.pickle",
    "discordance_gRNA_table": "discordance_gRNA.csv",
    "OVERWRITE_PCA_DICT": false
  },
  "input_data": {
    "annotation_file": {
      "file_path": "${OUTPUT_DIR_ABS}/annotation_table.csv",
      "concatenate_key": "intended_target_name"
    },
    "h5ad_file": {
      "file_path": "${OUTPUT_DIR_ABS}/preprocessed.h5ad",
      "obsm_key": "X_pca"
    },
    "sgRNA_file": {
      "file_path": "${OUTPUT_DIR_ABS}/gRNA_dict.pickle"
    }
  },
  "gRNA_filtering": {
    "perform_targeting_filtering": true,
    "perform_nontargeting_filtering": true,
    "threshold_gRNA_num": 6,
    "combi_count": 4,
    "total_permute_disco": 1000,
    "combi_cell_num_max": 1000,
    "batch_num_basic": 120
  },
  "permutation_test": {
    "permute_per_bg": 1000,
    "num_of_bg": 20,
    "non_target_pick": 2000,
    "target_cell_num_max": 2000,
    "batch_num_basic": 200,
    "use_matched_bg": false
  },
  "aggregate": {
    "downsampling_maximum": 10000
  }
}
JSON

cat > "$CONFIG3" <<JSON
{
  "cutoff": {
    "pvalue_cutoff": 0.05,
    "edist_cutoff": 0.5
  }
}
JSON

echo "[setup] configs written to $OUTPUT_DIR_ABS"

###############################################################################
# 3. Preprocess MuData (local file → preprocessed.h5ad + gRNA_dict + pca + annotation)
###############################################################################

# Bind mounts: data dir, output dir, repo root (for preprocess + bin scripts).
# muon is installed into /tmp/muon_deps inside the container (matches Sara's template).
BIND_ARGS=(
  --bind "$(dirname "$MUDATA_PATH"):$(dirname "$MUDATA_PATH")"
  --bind "${OUTPUT_DIR_ABS}:${OUTPUT_DIR_ABS}"
  --bind "${REPO_ROOT}:${REPO_ROOT}"
)

if [[ ! -f "${OUTPUT_DIR_ABS}/pca_dataframe.pickle" ]] || [[ ! -f "${OUTPUT_DIR_ABS}/gRNA_dict.pickle" ]]; then
  echo "[step0] preprocessing MuData ..."
  # Use container's pre-installed muon 0.1.7 (in /app/.venv). DO NOT pip install muon
  # to /tmp/muon_deps — that pulls a newer numpy (>=2.0) whose pickled artifacts
  # cannot be deserialized by the container's numpy 1.26.4 in steps 1/2/2.1.
  apptainer exec --nv "${BIND_ARGS[@]}" "$CONTAINER_PATH" \
    python ${PREPROCESS_LOCAL} \
      --mudata-path ${MUDATA_PATH} \
      --output-dir ${OUTPUT_DIR_ABS}
else
  echo "[step0] preprocessed files exist — skipping"
fi

###############################################################################
# 4. Pipeline steps 1, 2, 2.1 (and optionally 3)
#
# NOTE: at the currently pinned submodule commit, the step-1 script filename has
# the typo "1_filtereing_gRNA.py" (carried over from upstream). When/if we bump
# external/energy_dist_pipeline past the rename commit, change this to
# "1_filtering_gRNA.py".
###############################################################################

PYTHONPATH_RUN="${PIPELINE_BIN}"
# Prefer the renamed file (upstream main); fall back to the legacy typo'd name only if it's all that exists.
STEP1_SCRIPT="${PIPELINE_BIN}/1_filtering_gRNA.py"
[[ ! -f "$STEP1_SCRIPT" ]] && STEP1_SCRIPT="${PIPELINE_BIN}/1_filtereing_gRNA.py"

echo "[step1] filter outlier gRNAs"
apptainer exec --nv "${BIND_ARGS[@]}" "$CONTAINER_PATH" \
  bash -c "PYTHONPATH=${PYTHONPATH_RUN} python ${STEP1_SCRIPT} ${CONFIG12}"

echo "[step2] energy distance vs non-targeting"
apptainer exec --nv "${BIND_ARGS[@]}" "$CONTAINER_PATH" \
  bash -c "PYTHONPATH=${PYTHONPATH_RUN} python ${PIPELINE_BIN}/2_e_distance_nontargeting.py ${CONFIG12}"

echo "[step2.1] diagnostic plots"
apptainer exec --nv "${BIND_ARGS[@]}" "$CONTAINER_PATH" \
  bash -c "PYTHONPATH=${PYTHONPATH_RUN} python ${PIPELINE_BIN}/2_1_Plot_figure.py ${CONFIG12}"

if [[ "$SKIP_STEP3" -eq 0 ]]; then
  echo "[step3] target-by-target distance + clustering"
  apptainer exec --nv "${BIND_ARGS[@]}" "$CONTAINER_PATH" \
    bash -c "PYTHONPATH=${PYTHONPATH_RUN} python ${PIPELINE_BIN}/3_e_distance_among_regions.py ${CONFIG12} ${CONFIG3}"
else
  echo "[step3] SKIPPED (run separately after picking cutoffs in config3.json)"
fi

echo "[done] energy distance pipeline complete: ${OUTPUT_DIR_ABS}"
