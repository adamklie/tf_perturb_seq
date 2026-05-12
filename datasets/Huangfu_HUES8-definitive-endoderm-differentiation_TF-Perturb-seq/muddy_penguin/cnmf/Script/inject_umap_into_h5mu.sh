#!/bin/bash
#SBATCH --job-name=huangfu_de_inject_umap_K200
#SBATCH --partition=carter-compute
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_de_torchcnmf_KskillA/Inference/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_de_torchcnmf_KskillA/Inference/logs/%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

set -euo pipefail
echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate

H5MU="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_de_torchcnmf_KskillA/Inference/adata/cNMF_200_2_0.h5mu"

python -u "/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Script/inject_umap_into_h5mu.py" \
    --h5mu_path "$H5MU" \
    --prog_key cNMF \
    --n_comps_pca 50 \
    --n_neighbors 30 \
    --force

echo "Done @ $(date)"
