#!/bin/bash
#SBATCH --job-name=051126_honcm_torchcnmf_KskillA
#SBATCH --partition=carter-gpu
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --gres=gpu:a30:1
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/PerturbNMF/Result/051126_honcm_torchcnmf_KskillA/Inference/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/PerturbNMF/Result/051126_honcm_torchcnmf_KskillA/Inference/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

# Hon CM production cNMF Stage 1 (KskillA pattern matching DE/ESC).
# Input: PerturbNMF/Data/inference_mudata_cleaned.h5ad (built by convert_h5mu_to_h5ad.sh
# with --compute_umap, so X_pca/X_umap are pre-injected and propagate through cNMF
# into every output h5mu's rna modality — no post-hoc inject needed for Stage 3).
# ~1M cells x 7,263 genes x 8 K x 10 iter; expect 15-25h on A30 (DE was 3-5h at
# 269k cells, scaling ~linearly with cells).

set -euo pipefail
mkdir -p "/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/PerturbNMF/Result/051126_honcm_torchcnmf_KskillA/Inference/logs"
echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"
nvidia-smi || true

source "/cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate"
export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

BASE=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq

python -u /cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src/Stage1_Inference/torch-cNMF/Slurm_Version/torch_cnmf_inference_pipeline.py \
    --counts_fn "$BASE/PerturbNMF/Data/inference_mudata_cleaned.h5ad" \
    --output_directory "$BASE/PerturbNMF/Result" \
    --run_name "051126_honcm_torchcnmf_KskillA" \
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
    --batch_max_epoch 1000 \
    --batch_hals_max_iter 1000 \
    --batch_hals_tol 0.005 \
    --categorical_key batch \
    --gene_names_key symbol \
    --run_factorize --run_refit --run_compile_annotation --run_diagnostic_plots

echo "Done @ $(date)"
