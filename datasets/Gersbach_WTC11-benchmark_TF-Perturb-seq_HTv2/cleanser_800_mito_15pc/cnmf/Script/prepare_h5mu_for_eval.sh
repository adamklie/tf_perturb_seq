#!/bin/bash
#SBATCH --job-name=htv2_prepare_h5mu
#SBATCH --partition=carter-compute
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/Inference/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/Inference/logs/%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

set -euo pipefail
echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate
export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

BASE_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF"
RUN_NAME="050926_HTv2_20iter_5KHVG_torch_halsvar_batch"
ADATA_DIR="$BASE_DIR/Result/$RUN_NAME/Inference/adata"

python -u "$BASE_DIR/Script/prepare_h5mu_for_eval.py" \
    --adata_dir "$ADATA_DIR" \
    --K 5 6 7 8 9 10 11 12 13 14 15 17 19 21 23 25 27 30 35 40 45 50 55 60 70 80 90 100 150 200 \
    --sel_thresh_str 2_0 \
    --sample_value all \
    --control_token non-targeting

echo "Done @ $(date)"
