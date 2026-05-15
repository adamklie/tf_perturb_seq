#!/usr/bin/env bash
#SBATCH --job-name=edist_huangfu_es_qc_bc_none
#SBATCH --partition=carter-gpu
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --gres=gpu:a30:1
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/qc_barcode_filter_none/energy_distance/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/qc_barcode_filter_none/energy_distance/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

set -euo pipefail

module load apptainer

MUDATA="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/qc_barcode_filter_none/crispr_pipeline/pipeline_dashboard/inference_mudata.h5mu"
OUTPUT_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/qc_barcode_filter_none/energy_distance"

mkdir -p "${OUTPUT_DIR}/logs"
echo "Job: ${SLURM_JOB_ID:-?} @ $(hostname) $(date)"
nvidia-smi || true

bash /cellar/users/aklie/projects/tf_perturb_seq/scripts/run_energy_distance_pipeline.sh \
  --mudata-path "$MUDATA" \
  --output-dir "$OUTPUT_DIR"
