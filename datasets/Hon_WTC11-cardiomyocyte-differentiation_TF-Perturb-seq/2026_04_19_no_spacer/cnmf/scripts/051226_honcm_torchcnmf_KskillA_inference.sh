#!/bin/bash
#SBATCH --job-name=051226_honcm_torchcnmf_KskillA
#SBATCH --partition=carter-gpu
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --gres=gpu:a30:1
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/scratch/cnmf_logs/honcm_infer.%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/scratch/cnmf_logs/honcm_infer.%j.err

set -euo pipefail

DS=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/2026_04_19_no_spacer
RUN_NAME=051226_honcm_torchcnmf_KskillA
OUT_DIR=$DS/cnmf/Result
LOG_DIR=$OUT_DIR/$RUN_NAME/Inference/logs
DATA_H5AD=$DS/cnmf/Data/HonCM_2026_04_19_no_spacer_perturbnmf.h5ad
PIPELINE_PY=/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src/Stage1_Inference/torch-cNMF/Slurm_Version/torch_cnmf_inference_pipeline.py
VENV_PY=/cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/python

mkdir -p "$LOG_DIR"
echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"
nvidia-smi || true

export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

$VENV_PY -u "$PIPELINE_PY" \
    --counts_fn "$DATA_H5AD" \
    --output_directory "$OUT_DIR" \
    --run_name "$RUN_NAME" \
    --species human \
    --K 30 50 60 80 100 200 250 300 \
    --numiter 20 \
    --numhvgenes 5000 \
    --algo halsvar \
    --mode batch \
    --init random \
    --tol 1e-7 \
    --use_gpu \
    --batch_max_epoch 1000 \
    --batch_hals_max_iter 1000 \
    --batch_hals_tol 0.005 \
    --categorical_key batch \
    --gene_names_key symbol \
    --sel_thresh 2.0 \
    --seed 14 \
    --run_factorize --run_refit --run_compile_annotation

echo "Done @ $(date)"
