#!/usr/bin/env bash
###############################################################################
# Run the TF Perturb-seq energy distance pipeline (Steps 1 + 2 + 2.1) on a
# single dataset's inference_mudata.h5mu. Built around the wrapper repo
# (https://github.com/Chikara-Takeuchi/energy_dist_TFperturb) but uses a local
# MuData path instead of fetching from Synapse, so we don't redownload the file
# every job.
#
# Required inputs (one of):
#   --mudata-path      /path/to/inference_mudata.h5mu  (use already-downloaded file)
#   --gcs-mudata-path  gs://.../inference_mudata.h5mu  (downloads to working dir first)
#   --synapse-id       synXXXXX                         (downloads via synapseclient)
# And:
#   --output-dir       /path/to/where/results/should/go
#
# Optional:
#   --container-path   Path to .sif (default: shared sif under /cellar/users/aklie/opt/containers)
#   --pipeline-dir     Path to clone of energy_dist_pipeline (default: <output-dir>/energy_dist_pipeline)
#   --skip-step3 / --run-step3   Step 3 default = SKIP (run separately after picking cutoffs)
#
# Designed to be called from a SLURM job (see datasets/<dataset>/5_run_energy_distance.sh).
###############################################################################

set -euo pipefail

CONTAINER_PATH="${CONTAINER_PATH:-/cellar/users/aklie/opt/containers/edist_pipeline.sif}"
PIPELINE_DIR=""
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
    --pipeline-dir)    PIPELINE_DIR="$2"; shift 2 ;;
    --skip-step3)      SKIP_STEP3=1; shift ;;
    --run-step3)       SKIP_STEP3=0; shift ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

[[ -z "$OUTPUT_DIR" ]] && { echo "ERROR: --output-dir is required" >&2; exit 2; }
if [[ -z "$MUDATA_PATH" && -z "$GCS_MUDATA_PATH" && -z "$SYNAPSE_ID" ]]; then
  echo "ERROR: one of --mudata-path, --gcs-mudata-path, or --synapse-id is required" >&2; exit 2
fi

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Absolute paths once cd'd in
OUTPUT_DIR_ABS="$(pwd)"
PIPELINE_DIR="${PIPELINE_DIR:-${OUTPUT_DIR_ABS}/energy_dist_pipeline}"

###############################################################################
# 0. Resolve MuData path (download from GCS if needed)
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
    python -c "
import os, shutil, sys, synapseclient
syn = synapseclient.Synapse(silent=True)
syn.login(authToken=os.environ['SYNAPSE_AUTH_TOKEN'])
e = syn.get('${SYNAPSE_ID}', downloadLocation='${OUTPUT_DIR_ABS}')
src = e.path
dst = '${MUDATA_LOCAL}'
if src != dst:
    shutil.move(src, dst)
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

###############################################################################
# 2. Pipeline source clone
###############################################################################

if [[ ! -d "$PIPELINE_DIR" ]]; then
  echo "[setup] cloning energy_dist_pipeline -> $PIPELINE_DIR"
  git clone https://github.com/Chikara-Takeuchi/energy_dist_pipeline.git "$PIPELINE_DIR"
fi
BIN_PATH="${PIPELINE_DIR}/bin"
echo "[setup] pipeline bin: $BIN_PATH"

###############################################################################
# 3. Generate run-specific configs
#    config1_2.json controls preprocess + step 1 + step 2
#    config3.json   controls step 3 (clustering / aggregation)
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
# 4. Preprocess MuData (local; no Synapse fetch)
#    Writes: preprocessed.h5ad (with obsm["X_pca"]), gRNA_dict.pickle, annotation_table.csv
###############################################################################

PREPROCESS_PY="${OUTPUT_DIR_ABS}/preprocess_mudata_local.py"
cat > "$PREPROCESS_PY" <<'PYEOF'
"""Local-input variant of `preprocess_mudata.py` from Chikara-Takeuchi/energy_dist_TFperturb.

Reads a local inference_mudata.h5mu (no Synapse fetch), runs scanpy preprocessing
+ PCA, builds the gRNA->cells dict, and writes the three artifacts the pipeline
expects in OUTPUT_FOLDER:
  - preprocessed.h5ad  (anndata with obsm["X_pca"])
  - gRNA_dict.pickle   (dict of gRNA_name -> list of cell barcodes)
  - annotation_table.csv (columns: guide_id, intended_target_name, type, spacer)
"""
import argparse, json, os, pickle
from collections import defaultdict
import scanpy as sc, muon, pandas as pd, numpy as np


def get_promoter_name(row):
    if row["type"] == "non-targeting":
        return "non-targeting"
    return f"{row['intended_target_name']}|{row['intended_target_chr']}:{int(row['intended_target_start'])}-{int(row['intended_target_end'])}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mudata-path", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    out = args.output_dir
    os.makedirs(out, exist_ok=True)

    print(f"[preprocess] reading {args.mudata_path}")
    mdata = muon.read_h5mu(args.mudata_path)
    print(mdata)

    # The CRISPR pipeline's MuData uses 'gene' for RNA and 'gRNA' (or 'guide')
    rna_key = "gene" if "gene" in mdata.mod else "rna"
    gRNA_key = "gRNA" if "gRNA" in mdata.mod else "guide"
    print(f"[preprocess] modality keys: rna={rna_key}, guide={gRNA_key}")

    adata_exp = mdata[rna_key].copy()
    sc.pp.filter_genes(adata_exp, min_counts=1)
    sc.pp.normalize_total(adata_exp)
    sc.pp.log1p(adata_exp)
    sc.pp.scale(adata_exp)
    print("[preprocess] running PCA (n_comps=50)")
    sc.tl.pca(adata_exp, random_state=0, n_comps=50)

    # Build gRNA -> cells dict from the guide modality's X (assumed binary or count matrix)
    print("[preprocess] building gRNA -> cells dict")
    adata_g = mdata[gRNA_key]
    g_x = adata_g.X
    try:
        g_x = g_x.tocsc()
    except Exception:
        pass
    gRNA_dict = {}
    for j, gname in enumerate(adata_g.var_names):
        col = g_x[:, j]
        if hasattr(col, "toarray"):
            col = col.toarray().ravel()
        else:
            col = np.asarray(col).ravel()
        cells = list(adata_g.obs_names[col > 0])
        gRNA_dict[gname] = cells

    # Annotation table — pull from var of the guide modality
    print("[preprocess] building annotation table")
    g_var = adata_g.var.reset_index().rename(columns={"index": "guide_id"})
    keep_cols = ["guide_id"]
    for col in ("intended_target_name", "type", "spacer", "intended_target_chr",
                "intended_target_start", "intended_target_end"):
        if col in g_var.columns:
            keep_cols.append(col)
    annotation = g_var[keep_cols].copy()
    if {"intended_target_chr", "intended_target_start", "intended_target_end"} <= set(annotation.columns):
        annotation["intended_target_name"] = annotation.apply(
            lambda r: get_promoter_name(r) if r.get("type") not in ("non-targeting",) else "non-targeting", axis=1
        )
    annotation.to_csv(os.path.join(out, "annotation_table.csv"), index=False)

    # Save preprocessed h5ad
    print("[preprocess] writing preprocessed.h5ad + gRNA_dict.pickle")
    adata_exp.write_h5ad(os.path.join(out, "preprocessed.h5ad"))
    with open(os.path.join(out, "gRNA_dict.pickle"), "wb") as fh:
        pickle.dump(gRNA_dict, fh)

    # PCA dataframe (pickle) — what the pipeline reads via util_functions.load_files
    pca_df = pd.DataFrame(adata_exp.obsm["X_pca"], index=adata_exp.obs_names)
    pca_df.to_pickle(os.path.join(out, "pca_dataframe.pickle"))

    print("[preprocess] done")


if __name__ == "__main__":
    main()
PYEOF

if [[ ! -f "${OUTPUT_DIR_ABS}/pca_dataframe.pickle" ]] || [[ ! -f "${OUTPUT_DIR_ABS}/gRNA_dict.pickle" ]]; then
  echo "[step0] preprocessing MuData ..."
  apptainer exec --nv \
    --bind "$(dirname "$MUDATA_PATH"):$(dirname "$MUDATA_PATH")" \
    --bind "${OUTPUT_DIR_ABS}:${OUTPUT_DIR_ABS}" \
    "$CONTAINER_PATH" \
    bash -c "pip install --quiet --target=/tmp/muon_deps muon && \
             PYTHONPATH=/tmp/muon_deps python ${PREPROCESS_PY} \
               --mudata-path ${MUDATA_PATH} \
               --output-dir ${OUTPUT_DIR_ABS}"
else
  echo "[step0] preprocessed files exist — skipping"
fi

###############################################################################
# 5. Pipeline steps 1, 2, 2.1
###############################################################################

PYTHONPATH_RUN="/tmp/muon_deps:${PIPELINE_DIR}"

echo "[step1] filter outlier gRNAs"
apptainer exec --nv \
  --bind "${OUTPUT_DIR_ABS}:${OUTPUT_DIR_ABS}" \
  --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
  "$CONTAINER_PATH" \
  bash -c "PYTHONPATH=${PYTHONPATH_RUN} python ${BIN_PATH}/1_filtereing_gRNA.py ${CONFIG12}"

echo "[step2] energy distance vs non-targeting"
apptainer exec --nv \
  --bind "${OUTPUT_DIR_ABS}:${OUTPUT_DIR_ABS}" \
  --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
  "$CONTAINER_PATH" \
  bash -c "PYTHONPATH=${PYTHONPATH_RUN} python ${BIN_PATH}/2_e_distance_nontargeting.py ${CONFIG12}"

echo "[step2.1] diagnostic plots"
apptainer exec --nv \
  --bind "${OUTPUT_DIR_ABS}:${OUTPUT_DIR_ABS}" \
  --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
  "$CONTAINER_PATH" \
  bash -c "PYTHONPATH=${PYTHONPATH_RUN} python ${BIN_PATH}/2_1_Plot_figure.py ${CONFIG12}"

if [[ "$SKIP_STEP3" -eq 0 ]]; then
  echo "[step3] target-by-target distance + clustering"
  apptainer exec --nv \
    --bind "${OUTPUT_DIR_ABS}:${OUTPUT_DIR_ABS}" \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    "$CONTAINER_PATH" \
    bash -c "PYTHONPATH=${PYTHONPATH_RUN} python ${BIN_PATH}/3_e_distance_among_regions.py ${CONFIG12} ${CONFIG3}"
else
  echo "[step3] SKIPPED (run separately after picking cutoffs in config3.json)"
fi

echo "[done] energy distance pipeline complete: ${OUTPUT_DIR_ABS}"
