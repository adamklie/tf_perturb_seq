#!/usr/bin/env bash
#SBATCH --job-name=edist_huangfu_de
#SBATCH --partition=carter-gpu
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --gres=gpu:a30:1
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/results/energy_distance/muddy_penguin/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/results/energy_distance/muddy_penguin/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

set -euo pipefail

# Load apptainer (not in default PATH on nrnb compute nodes)
module load apptainer

# MuData lives on GCS; the runner downloads it to OUTPUT_DIR/inference_mudata.h5mu
# (skipped if already present, e.g. on resubmissions).
GCS_MUDATA="gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/2026_04_09/outs/muddy_penguin/inference_mudata.h5mu"
OUTPUT_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/results/energy_distance/muddy_penguin"

mkdir -p "${OUTPUT_DIR}/logs"
echo "Job: ${SLURM_JOB_ID:-?} @ $(hostname) $(date)"
nvidia-smi || true

# gcloud account must be the IGVF service account to read the bucket
gcloud config set account adamklie@igvf-pertub-seq-pipeline.iam.gserviceaccount.com >/dev/null

bash /cellar/users/aklie/projects/tf_perturb_seq/scripts/run_energy_distance_pipeline.sh \
  --gcs-mudata-path "$GCS_MUDATA" \
  --output-dir "$OUTPUT_DIR"
