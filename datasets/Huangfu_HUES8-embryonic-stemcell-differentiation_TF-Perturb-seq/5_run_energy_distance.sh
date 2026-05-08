#!/usr/bin/env bash
#SBATCH --job-name=edist_huangfu_es
#SBATCH --partition=carter-gpu
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --gres=gpu:a30:1
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/results/energy_distance/sceptre_v1/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/results/energy_distance/sceptre_v1/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

set -euo pipefail

# Load apptainer (not in default PATH on nrnb compute nodes)
module load apptainer

GCS_MUDATA="gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/2026_04_13/outs/sceptre_v1/inference_mudata.h5mu"
OUTPUT_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/results/energy_distance/sceptre_v1"

mkdir -p "${OUTPUT_DIR}/logs"
echo "Job: ${SLURM_JOB_ID:-?} @ $(hostname) $(date)"
nvidia-smi || true

gcloud config set account adamklie@igvf-pertub-seq-pipeline.iam.gserviceaccount.com >/dev/null

bash /cellar/users/aklie/projects/tf_perturb_seq/scripts/run_energy_distance_pipeline.sh \
  --gcs-mudata-path "$GCS_MUDATA" \
  --output-dir "$OUTPUT_DIR"
