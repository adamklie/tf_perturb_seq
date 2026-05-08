#!/usr/bin/env bash
#SBATCH --job-name=edist_hon_cm
#SBATCH --partition=carter-gpu
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --gres=gpu:a30:1
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/results/energy_distance/2026_04_19_no_spacer/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/results/energy_distance/2026_04_19_no_spacer/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

set -euo pipefail

# As of 2026-05-08, the only completed Hon CM CRISPR pipeline output is on
# Synapse (the GCS prefix /2026_04_15/outs/initial_run/ is missing
# pipeline_outputs/). Pull the inference_mudata directly from Synapse for now.
# The Hon team is uploading a full GCS bundle; once that lands switch to
# --gcs-mudata-path and update OUTPUT_DIR's run-label.
SYNAPSE_ID="syn74522725"   # Hon_CM 2026_04_19_no_spacer pipeline_outputs/inference_mudata.h5mu
OUTPUT_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/results/energy_distance/2026_04_19_no_spacer"

mkdir -p "${OUTPUT_DIR}/logs"
echo "Job: ${SLURM_JOB_ID:-?} @ $(hostname) $(date)"
nvidia-smi || true

# SYNAPSE_AUTH_TOKEN must already be in the environment (from ~/.bashrc).
[[ -z "${SYNAPSE_AUTH_TOKEN:-}" ]] && { echo "ERROR: SYNAPSE_AUTH_TOKEN not set" >&2; exit 2; }

bash /cellar/users/aklie/projects/tf_perturb_seq/scripts/run_energy_distance_pipeline.sh \
  --synapse-id "$SYNAPSE_ID" \
  --output-dir "$OUTPUT_DIR"
