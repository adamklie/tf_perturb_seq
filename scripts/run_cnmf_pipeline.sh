#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# cNMF Gene Program Discovery Pipeline Runner
#
# This script runs the full cNMF analysis pipeline:
#   Phase 1: Convert MuData to AnnData (format_anndata.py)
#   Phase 2: Run torch-cNMF inference
#   Phase 3: Run cNMF evaluation pipeline
#   Phase 4: Generate k-selection plots
#
# Parameters follow docs/analysis/cNMF.md guidelines.
#
# USAGE:
#   run_cnmf_pipeline.sh \
#     --project-root /path/to/tf_perturb_seq \
#     --input /path/to/inference_mudata.h5mu \
#     --outdir /path/to/output \
#     --run-name dataset_name \
#     [--cnmf-repo-path /path/to/cNMF_benchmarking] \
#     [--K "5 6 7 8 9 10 ..."] \
#     [--numiter 100] \
#     [--numhvgenes 1000] \
#     [--seed 123] \
#     [--species human] \
#     [--sel-thresh "0.01 0.05 2.0"] \
#     [--guide-annotation-path /path/to/guide_metadata.tsv] \
#     [--use-gpu] \
#     [--skip-format] \
#     [--skip-eval] \
#     [--phase 1|2|3|4] \
#     [--slurm] \
#     [--slurm-partition GPUp40] \
#     [--slurm-time 99:00:00] \
#     [--slurm-gres gpu:1] \
#     [--slurm-email ""] \
#     [--dry-run]
#
###############################################################################

# Default values
PROJECT_ROOT=""
INPUT=""
OUTDIR=""
RUN_NAME=""
CNMF_REPO=""
K="5 6 7 8 9 10 11 12 13 14 15 17 19 21 23 25 27 30 35 40 45 50 55 60 70 80 90 100 150 200"
NUMITER=100
NUMHVGENES=1000
SEED=123
SPECIES="human"
SEL_THRESH="0.01 0.05 2.0"
GUIDE_ANNOTATION=""
USE_GPU=0
SKIP_FORMAT=0
SKIP_EVAL=0
PHASE=""
SLURM=0
SLURM_PARTITION="GPUp40"
SLURM_TIME="99:00:00"
SLURM_GRES="gpu:1"
SLURM_EMAIL=""
DRY_RUN=0

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root)           PROJECT_ROOT="$2"; shift 2 ;;
    --input)                  INPUT="$2"; shift 2 ;;
    --outdir)                 OUTDIR="$2"; shift 2 ;;
    --run-name)               RUN_NAME="$2"; shift 2 ;;
    --cnmf-repo-path)         CNMF_REPO="$2"; shift 2 ;;
    --K)                      K="$2"; shift 2 ;;
    --numiter)                NUMITER="$2"; shift 2 ;;
    --numhvgenes)             NUMHVGENES="$2"; shift 2 ;;
    --seed)                   SEED="$2"; shift 2 ;;
    --species)                SPECIES="$2"; shift 2 ;;
    --sel-thresh)             SEL_THRESH="$2"; shift 2 ;;
    --guide-annotation-path)  GUIDE_ANNOTATION="$2"; shift 2 ;;
    --use-gpu)                USE_GPU=1; shift 1 ;;
    --skip-format)            SKIP_FORMAT=1; shift 1 ;;
    --skip-eval)              SKIP_EVAL=1; shift 1 ;;
    --phase)                  PHASE="$2"; shift 2 ;;
    --slurm)                  SLURM=1; shift 1 ;;
    --slurm-partition)        SLURM_PARTITION="$2"; shift 2 ;;
    --slurm-time)             SLURM_TIME="$2"; shift 2 ;;
    --slurm-gres)             SLURM_GRES="$2"; shift 2 ;;
    --slurm-email)            SLURM_EMAIL="$2"; shift 2 ;;
    --dry-run)                DRY_RUN=1; shift 1 ;;
    -h|--help)
      sed -n '1,45p' "$0"
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

###############################################################################
# Environment Setup
###############################################################################

VENV="${PROJECT_ROOT}/.venv/bin/activate"
if [[ -f "${VENV}" ]]; then
  # shellcheck disable=SC1090
  source "${VENV}"
fi

# Default to submodule path if not overridden
if [[ -z "${CNMF_REPO}" ]]; then
  CNMF_REPO="${PROJECT_ROOT}/external/cNMF_benchmarking"
fi

# Validate cNMF repo
CNMF_PIPELINE="${CNMF_REPO}/cNMF_benchmarking_pipeline"
INFERENCE_SCRIPT="${CNMF_PIPELINE}/Inference/torch-cNMF/Slurm_Version/torch-cNMF_inference_pipeline.py"
EVAL_SCRIPT="${CNMF_PIPELINE}/Evaluation/Slurm_Version/cNMF_evaluation_pipeline.py"
KSELECT_SCRIPT="${CNMF_PIPELINE}/Plotting/Slurm_Version/cNMF_k_selection.py"

if [[ ! -d "${CNMF_PIPELINE}" ]]; then
  echo "ERROR: cNMF_benchmarking pipeline not found at: ${CNMF_PIPELINE}" >&2
  echo "Run: git submodule update --init --recursive" >&2
  exit 1
fi

# Add to PYTHONPATH so 'import cnmf' and 'from Inference.src import ...' work
export PYTHONPATH="${CNMF_PIPELINE}${PYTHONPATH:+:$PYTHONPATH}"

# Adapter script
FORMAT_SCRIPT="${PROJECT_ROOT}/src/tf_perturb_seq/gene_program_discovery/format_anndata.py"

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

should_run_phase() {
  local phase_num="$1"
  if [[ -n "${PHASE}" ]]; then
    [[ "${PHASE}" == "${phase_num}" ]]
  else
    return 0
  fi
}

###############################################################################
# Output Layout
###############################################################################

DIR_INPUTS="${OUTDIR}/inputs"

echo "============================================="
echo " cNMF Gene Program Discovery Pipeline"
echo " Project:    ${PROJECT_ROOT}"
echo " cNMF repo:  ${CNMF_REPO}"
echo " Input:      ${INPUT}"
echo " Output:     ${OUTDIR}"
echo " Run name:   ${RUN_NAME}"
echo " K values:   ${K}"
echo " Iterations: ${NUMITER}"
echo " HVGs:       ${NUMHVGENES}"
echo " Seed:       ${SEED}"
echo " Species:    ${SPECIES}"
echo " Sel thresh: ${SEL_THRESH}"
echo " GPU:        ${USE_GPU}"
echo " Phase:      ${PHASE:-all}"
echo " SLURM:      ${SLURM}"
echo " Dry-run:    ${DRY_RUN}"
echo "============================================="

run_cmd mkdir -p "${DIR_INPUTS}"

###############################################################################
# SLURM Mode: Generate batch script and exit
###############################################################################

if [[ "${SLURM}" -eq 1 ]]; then
  SLURM_SCRIPT="${OUTDIR}/slurm_cnmf.sh"
  run_cmd mkdir -p "${OUTDIR}"

  echo "Generating SLURM batch script: ${SLURM_SCRIPT}"

  if [[ "${DRY_RUN}" -eq 0 ]]; then
    cat > "${SLURM_SCRIPT}" << SLURM_EOF
#!/bin/bash
#SBATCH --job-name=cnmf_${RUN_NAME}
#SBATCH --output=${OUTDIR}/logs/%j.out
#SBATCH --error=${OUTDIR}/logs/%j.err
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --time=${SLURM_TIME}
#SBATCH --gres=${SLURM_GRES}
SLURM_EOF

    if [[ -n "${SLURM_EMAIL}" ]]; then
      cat >> "${SLURM_SCRIPT}" << SLURM_EMAIL_EOF
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${SLURM_EMAIL}
SLURM_EMAIL_EOF
    fi

    cat >> "${SLURM_SCRIPT}" << 'SLURM_BODY_EOF'

START_TIME=$(date +%s)
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Partition: $SLURM_JOB_PARTITION"

SLURM_BODY_EOF

    # Build the pipeline command (re-invoke this script without --slurm)
    PIPELINE_CMD="${PROJECT_ROOT}/scripts/run_cnmf_pipeline.sh"
    PIPELINE_CMD+=" --project-root ${PROJECT_ROOT}"
    PIPELINE_CMD+=" --input ${INPUT}"
    PIPELINE_CMD+=" --outdir ${OUTDIR}"
    PIPELINE_CMD+=" --run-name ${RUN_NAME}"
    PIPELINE_CMD+=" --cnmf-repo-path ${CNMF_REPO}"
    PIPELINE_CMD+=" --K \"${K}\""
    PIPELINE_CMD+=" --numiter ${NUMITER}"
    PIPELINE_CMD+=" --numhvgenes ${NUMHVGENES}"
    PIPELINE_CMD+=" --seed ${SEED}"
    PIPELINE_CMD+=" --species ${SPECIES}"
    PIPELINE_CMD+=" --sel-thresh \"${SEL_THRESH}\""
    if [[ -n "${GUIDE_ANNOTATION}" ]]; then
      PIPELINE_CMD+=" --guide-annotation-path ${GUIDE_ANNOTATION}"
    fi
    if [[ "${USE_GPU}" -eq 1 ]]; then
      PIPELINE_CMD+=" --use-gpu"
    fi
    if [[ "${SKIP_FORMAT}" -eq 1 ]]; then
      PIPELINE_CMD+=" --skip-format"
    fi
    if [[ "${SKIP_EVAL}" -eq 1 ]]; then
      PIPELINE_CMD+=" --skip-eval"
    fi
    if [[ -n "${PHASE}" ]]; then
      PIPELINE_CMD+=" --phase ${PHASE}"
    fi

    cat >> "${SLURM_SCRIPT}" << SLURM_CMD_EOF
# Create logs directory
mkdir -p ${OUTDIR}/logs

# Run the pipeline
${PIPELINE_CMD}

# Summary
END_TIME=\$(date +%s)
DURATION=\$((END_TIME - START_TIME))
echo "========================================="
echo "EXECUTION SUMMARY"
echo "========================================="
echo "Total execution time: \${DURATION} seconds (\$((\$DURATION / 3600))h \$((\$DURATION % 3600 / 60))m \$((\$DURATION % 60))s)"
echo "Job completed at: \$(date)"
SLURM_CMD_EOF

    chmod +x "${SLURM_SCRIPT}"
    echo ""
    echo "SLURM script generated: ${SLURM_SCRIPT}"
    echo "Submit with: sbatch ${SLURM_SCRIPT}"
  else
    echo "[DRY-RUN] Would generate SLURM script at ${SLURM_SCRIPT}"
  fi
  exit 0
fi

###############################################################################
# Phase 1: Convert MuData to AnnData
###############################################################################

if should_run_phase 1; then
  if [[ "${SKIP_FORMAT}" -eq 0 ]]; then
    echo ""
    echo "[Phase 1/4] Converting MuData to AnnData..."

    [[ -f "${INPUT}" ]] || { echo "ERROR: Input file not found: ${INPUT}" >&2; exit 1; }
    [[ -f "${FORMAT_SCRIPT}" ]] || { echo "ERROR: Format script not found: ${FORMAT_SCRIPT}" >&2; exit 1; }

    COUNTS_H5AD="${DIR_INPUTS}/counts.h5ad"

    run_cmd python "${FORMAT_SCRIPT}" \
      --mudata "${INPUT}" \
      --out_h5ad "${COUNTS_H5AD}" \
      --dense_guide_assignment \
      --skip_embeddings

    echo "Phase 1 complete: ${COUNTS_H5AD}"
  else
    echo ""
    echo "[Phase 1/4] Skipped (--skip-format)"
  fi
fi

COUNTS_H5AD="${DIR_INPUTS}/counts.h5ad"

###############################################################################
# Phase 2: Run torch-cNMF Inference
###############################################################################

if should_run_phase 2; then
  echo ""
  echo "[Phase 2/4] Running torch-cNMF inference..."

  [[ -f "${INFERENCE_SCRIPT}" ]] || { echo "ERROR: Inference script not found: ${INFERENCE_SCRIPT}" >&2; exit 1; }

  # Build inference command
  # shellcheck disable=SC2206
  INFERENCE_CMD=(
    python "${INFERENCE_SCRIPT}"
    --counts_fn "${COUNTS_H5AD}"
    --output_directory "${OUTDIR}"
    --run_name "${RUN_NAME}"
    --K ${K}
    --numiter "${NUMITER}"
    --numhvgenes "${NUMHVGENES}"
    --seed "${SEED}"
    --species "${SPECIES}"
    --algo "halsvar"
    --mode "batch"
    --tol 1e-7
    --sel_thresh ${SEL_THRESH}
    --densify
    --sk_cd_refit
    --shuffle_cells
    --batch_max_iter 1000
    --batch_hals_max_iter 1000
    --batch_hals_tol 0.005
  )

  if [[ "${USE_GPU}" -eq 1 ]]; then
    INFERENCE_CMD+=(--use_gpu)
  fi

  run_cmd "${INFERENCE_CMD[@]}"

  echo "Phase 2 complete: inference results in ${OUTDIR}/${RUN_NAME}/"
fi

###############################################################################
# Phase 3: Run cNMF Evaluation
###############################################################################

if should_run_phase 3; then
  if [[ "${SKIP_EVAL}" -eq 0 ]]; then
    echo ""
    echo "[Phase 3/4] Running cNMF evaluation..."

    [[ -f "${EVAL_SCRIPT}" ]] || { echo "ERROR: Evaluation script not found: ${EVAL_SCRIPT}" >&2; exit 1; }

    # Build evaluation command
    # shellcheck disable=SC2206
    EVAL_CMD=(
      python "${EVAL_SCRIPT}"
      --out_dir "${OUTDIR}"
      --run_name "${RUN_NAME}"
      --K ${K}
      --sel_thresh ${SEL_THRESH}
      --X_normalized_path "${COUNTS_H5AD}"
      --mdata_guide_path "${INPUT}"
      --organism "${SPECIES}"
      --Perform_explained_variance
      --Perform_categorical
      --Perform_perturbation
      --Perform_geneset
      --Perform_trait
    )

    if [[ -n "${GUIDE_ANNOTATION}" ]]; then
      EVAL_CMD+=(--guide_annotation_path "${GUIDE_ANNOTATION}")
    fi

    run_cmd "${EVAL_CMD[@]}"

    echo "Phase 3 complete: evaluation results in ${OUTDIR}/${RUN_NAME}/Eval/"
  else
    echo ""
    echo "[Phase 3/4] Skipped (--skip-eval)"
  fi
fi

###############################################################################
# Phase 4: Generate K-Selection Plots
###############################################################################

if should_run_phase 4; then
  if [[ "${SKIP_EVAL}" -eq 0 ]]; then
    echo ""
    echo "[Phase 4/4] Generating k-selection plots..."

    [[ -f "${KSELECT_SCRIPT}" ]] || { echo "ERROR: K-selection script not found: ${KSELECT_SCRIPT}" >&2; exit 1; }

    DIR_PLOTS="${OUTDIR}/plots"
    run_cmd mkdir -p "${DIR_PLOTS}"

    # Build k-selection command
    # shellcheck disable=SC2206
    KSELECT_CMD=(
      python "${KSELECT_SCRIPT}"
      --output_directory "${OUTDIR}"
      --run_name "${RUN_NAME}"
      --save_folder_name "${DIR_PLOTS}"
      --eval_folder_name "${OUTDIR}/${RUN_NAME}/Eval"
      --K ${K}
      --sel_threshs ${SEL_THRESH}
    )

    run_cmd "${KSELECT_CMD[@]}"

    echo "Phase 4 complete: plots in ${DIR_PLOTS}/"
  else
    echo ""
    echo "[Phase 4/4] Skipped (--skip-eval)"
  fi
fi

###############################################################################
# Done
###############################################################################

echo ""
echo "============================================="
if [[ "${DRY_RUN}" -eq 1 ]]; then
  echo " Dry-run complete (no files written)."
else
  echo " cNMF pipeline complete!"
  echo ""
  echo " Outputs:"
  echo "   Inputs:     ${DIR_INPUTS}"
  echo "   Inference:  ${OUTDIR}/${RUN_NAME}/"
  if [[ "${SKIP_EVAL}" -eq 0 ]]; then
    echo "   Evaluation: ${OUTDIR}/${RUN_NAME}/Eval/"
    echo "   Plots:      ${OUTDIR}/plots/"
  fi
fi
echo "============================================="
