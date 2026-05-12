#!/bin/bash
#SBATCH --job-name=huangfu_de_kselection
#SBATCH --partition=carter-compute
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_de_torchcnmf_KskillA/Plot/k_selection/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_de_torchcnmf_KskillA/Plot/k_selection/logs/%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

# Stage 3a: K-selection plot
# Combines stability/error from inference with evaluation metrics
# (perturbation, gene-set, GO, explained variance) into per-K diagnostic plots.
# Groups stability/error by --groupby (batch). Depends on Stage 2a outputs in
# Evaluation/<K>_2_0/.

set -euo pipefail
START_TIME=$(date +%s)

OUT_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result"
RUN_NAME="042926_huangfu_de_torchcnmf_KskillA"
RUN_DIR="$OUT_DIR/$RUN_NAME"

echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"
echo "Run dir: $RUN_DIR"

mkdir -p "$RUN_DIR/Plot/k_selection/logs"

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate
export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

python -u /cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src/Stage3_Interpretation/A_Plotting/Slurm_Version/cNMF_k_selection.py \
    --output_directory "$OUT_DIR" \
    --run_name "$RUN_NAME" \
    --save_folder_name "$RUN_DIR/Plot/k_selection" \
    --eval_folder_name "$RUN_DIR/Evaluation" \
    --stability_file "$RUN_DIR/Inference/Inference.k_selection_stats.df.npz" \
    --groupby "batch" \
    --K 30 50 60 80 100 200 250 300 \
    --sel_threshs 2.0 \
    --samples all

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
printf 'Done @ %s | Elapsed: %dh %dm %ds\n' "$(date)" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
