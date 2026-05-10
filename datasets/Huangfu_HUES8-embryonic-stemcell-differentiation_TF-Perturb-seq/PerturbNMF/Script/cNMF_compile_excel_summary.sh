#!/bin/bash
#SBATCH --job-name=huangfu_esc_excel_K200
#SBATCH --partition=carter-compute
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_esc_torchcnmf_KskillA/Interpretation/Summary_table/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_esc_torchcnmf_KskillA/Interpretation/Summary_table/logs/%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

set -euo pipefail
START_TIME=$(date +%s)

OUT_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result"
RUN_NAME="042926_huangfu_esc_torchcnmf_KskillA"
RUN_DIR="$OUT_DIR/$RUN_NAME"
SAVE_PATH="$RUN_DIR/Interpretation/Summary_table/200_2_0/work"

echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"
mkdir -p "$RUN_DIR/Interpretation/Summary_table/logs"
mkdir -p "$SAVE_PATH"

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate
export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

python -u /cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Script/cNMF_compile_excel_summary.py \
    --out_dir "$OUT_DIR" \
    --run_name "$RUN_NAME" \
    --save_path "$SAVE_PATH" \
    --K 200 \
    --sel_thresh 2.0 \
    --samples all \
    --perturbation_file_name perturbation_association_results \
    --non_targeting_key non-targeting \
    --categorical_key sample \
    --prog_key cNMF \
    --data_key rna \
    --guide_targets_key guide_targets \
    --effect_size log2FC

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
printf 'Done @ %s | Elapsed: %dh %dm %ds\n' "$(date)" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
