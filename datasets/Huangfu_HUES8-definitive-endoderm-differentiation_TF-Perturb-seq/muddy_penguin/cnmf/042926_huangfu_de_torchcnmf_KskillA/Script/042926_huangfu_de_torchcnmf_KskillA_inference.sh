#!/bin/bash
#SBATCH --job-name=042926_huangfu_de_torchcnmf_KskillA
#SBATCH --partition=carter-gpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --gres=gpu:a30:1
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_de_torchcnmf_KskillA/Inference/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_de_torchcnmf_KskillA/Inference/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

set -euo pipefail
mkdir -p "/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_de_torchcnmf_KskillA/Inference/logs"
echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"
nvidia-smi || true

source "/cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate"
export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

python -u /cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src/Inference/torch-cNMF/Slurm_Version/torch_cnmf_inference_pipeline.py \
    --counts_fn "/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Data/DE_muddy_penguin_perturbnmf.h5ad" \
    --output_directory "/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result" \
    --run_name "042926_huangfu_de_torchcnmf_KskillA" \
    --species human \
    --K 30 50 60 80 100 200 250 300 \
    --numiter 10 \
    --numhvgenes 2000 \
    --sel_thresh 2.0 \
    --seed 14 \
    --algo halsvar \
    --mode batch \
    --tol 1e-4 \
    --use_gpu \
    --categorical_key batch \
    --gene_names_key symbol \
    --run_factorize --run_refit --run_compile_annotation --run_diagnostic_plots

echo "Done @ $(date)"
