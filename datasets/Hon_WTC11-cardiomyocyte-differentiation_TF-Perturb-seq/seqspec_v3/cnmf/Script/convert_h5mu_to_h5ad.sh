#!/bin/bash
#SBATCH --job-name=honcm_convert_h5ad
#SBATCH --partition=carter-compute
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/PerturbNMF/Data/convert_%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/PerturbNMF/Data/convert_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

# Convert Hon CM production inference_mudata.h5mu (16.85 GB, ~1M cells, 7,263 genes)
# into an AnnData with PCA + UMAP baked into obsm (per the new convention so
# Stage 1 carries them into every output h5mu without a post-hoc inject).

set -euo pipefail
echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate

BASE=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq
python -u $BASE/PerturbNMF/Script/Convert_file_adata.py \
    --input $BASE/PerturbNMF/Data/inference_mudata.h5mu \
    --output $BASE/PerturbNMF/Data/inference_mudata_cleaned.h5ad \
    --compute_umap \
    --n_top_hvg 2000 \
    --n_comps_pca 50 \
    --n_neighbors 30

echo "Done @ $(date)"
