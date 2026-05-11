#!/bin/bash
#SBATCH --job-name=htv2_kselection
#SBATCH --partition=carter-compute
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/Plot/k_selection/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/Plot/k_selection/logs/%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

# Stage 3a: K-selection plot. Depends on Stage 2a Evaluation outputs.
# Uses the same K sweep as Stage 1/2a.

set -euo pipefail
START_TIME=$(date +%s)

OUT_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result"
RUN_NAME="050926_HTv2_20iter_5KHVG_torch_halsvar_batch"
RUN_DIR="$OUT_DIR/$RUN_NAME"

echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"
echo "Run dir: $RUN_DIR"

mkdir -p "$RUN_DIR/Plot/k_selection/logs"

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate
export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

# Required workaround per docs/analysis/PerturbNMF.md issue #6: upstream hardcodes
# <run>/adata/ but our outputs live at <run>/Inference/adata/. Create symlink if absent.
RUN_ADATA="$RUN_DIR/adata"
INF_ADATA="$RUN_DIR/Inference/adata"
[ -L "$RUN_ADATA" ] || [ -d "$RUN_ADATA" ] || ln -sv "$INF_ADATA" "$RUN_ADATA"

python -u /cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src/Stage3_Interpretation/A_Plotting/Slurm_Version/cNMF_k_selection.py \
    --output_directory "$OUT_DIR" \
    --run_name "$RUN_NAME" \
    --save_folder_name "$RUN_DIR/Plot/k_selection" \
    --eval_folder_name "$RUN_DIR/Evaluation" \
    --stability_file "$RUN_DIR/Inference/Inference.k_selection_stats.df.npz" \
    --groupby "batch" \
    --K 5 6 7 8 9 10 11 12 13 14 15 17 19 21 23 25 27 30 35 40 45 50 55 60 70 80 90 100 150 200 \
    --sel_threshs 2.0 \
    --samples all

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
printf 'Done @ %s | Elapsed: %dh %dm %ds\n' "$(date)" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
