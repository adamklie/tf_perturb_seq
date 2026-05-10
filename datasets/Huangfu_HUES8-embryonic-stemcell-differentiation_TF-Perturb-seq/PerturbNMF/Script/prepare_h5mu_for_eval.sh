#!/bin/bash
#SBATCH --job-name=huangfu_esc_prepare_h5mu
#SBATCH --partition=carter-compute
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_esc_torchcnmf_KskillA/Inference/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_esc_torchcnmf_KskillA/Inference/logs/%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

# Pre-processes the 8 cNMF_<K>_2_0.h5mu files produced by inference:
#   - adds obs['sample'] = 'all' (single-value categorical) so Stage 2
#     perturbation association is NOT stratified by batch
#   - remaps guide_targets for non-targeting guides ('nan' -> 'non-targeting')
# Idempotent. Re-running is safe.

set -euo pipefail
echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate
export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

BASE_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF"
RUN_NAME="042926_huangfu_esc_torchcnmf_KskillA"
ADATA_DIR="$BASE_DIR/Result/$RUN_NAME/Inference/adata"

python -u "$BASE_DIR/Script/prepare_h5mu_for_eval.py" \
    --adata_dir "$ADATA_DIR" \
    --K 30 50 60 80 100 200 250 300 \
    --sel_thresh_str 2_0 \
    --sample_value all \
    --control_token non-targeting

echo "Done @ $(date)"
