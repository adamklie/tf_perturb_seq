#!/bin/bash
#SBATCH --job-name=htv2_inject_umap
#SBATCH --partition=carter-compute
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/Inference/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/Inference/logs/%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

# Inject a gene-based UMAP into the SELECTED K h5mu (Stage 3 input).
#
# Run AFTER Stage 3a K-selection picks a k. UMAP only depends on the rna
# matrix (identical across K), so we inject into just the one h5mu Stage 3
# will use — don't fan out to all 30.
#
# Override the K via the SELECTED_K env var:
#   SELECTED_K=50 sbatch inject_umap_into_h5mu.sh
#
# Going forward: compute UMAP in Convert_file_adata.py before Stage 1 so
# this post-hoc step isn't needed at all.

set -euo pipefail
echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate

SELECTED_K="${SELECTED_K:-50}"
ADATA_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/Inference/adata"
H5MU="$ADATA_DIR/cNMF_${SELECTED_K}_2_0.h5mu"

echo "Selected K: $SELECTED_K"
echo "Target h5mu: $H5MU"

python -u "/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Script/inject_umap_into_h5mu.py" \
    --h5mu_path "$H5MU" \
    --prog_key cNMF \
    --n_top_hvg 2000 \
    --n_comps_pca 50 \
    --n_neighbors 30

echo "Done @ $(date)"
