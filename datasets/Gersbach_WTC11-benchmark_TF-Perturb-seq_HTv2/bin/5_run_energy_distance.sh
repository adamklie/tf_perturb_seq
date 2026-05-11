#!/usr/bin/env bash
#SBATCH --job-name=edist_htv2
#SBATCH --partition=carter-gpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --gres=gpu:a30:1
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/results/energy_distance/cleanser_800_mito_15pc/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/results/energy_distance/cleanser_800_mito_15pc/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

set -euo pipefail

# Load apptainer (not in default PATH on nrnb compute nodes)
module load apptainer

# HTv2 testbed: same `cleanser_800_mito_15pc` CRISPR pipeline run that's mirrored
# to Synapse (syn74885574). Source MuData lives locally on HPC under that run's
# pipeline_dashboard/. Smaller dataset (~37k cells) — should finish in a few hours,
# unlike the production 270k-cell Huangfu runs.
MUDATA_LOCAL="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/runs/cleanser_800_mito_15pc/pipeline_dashboard/inference_mudata.h5mu"
OUTPUT_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/results/energy_distance/cleanser_800_mito_15pc"

mkdir -p "${OUTPUT_DIR}/logs"
echo "Job: ${SLURM_JOB_ID:-?} @ $(hostname) $(date)"
nvidia-smi || true

bash /cellar/users/aklie/projects/tf_perturb_seq/scripts/run_energy_distance_pipeline.sh \
  --mudata-path "$MUDATA_LOCAL" \
  --output-dir "$OUTPUT_DIR"
