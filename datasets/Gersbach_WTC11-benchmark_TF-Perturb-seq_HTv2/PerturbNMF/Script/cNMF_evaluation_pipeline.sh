#!/bin/bash
#SBATCH --job-name=htv2_eval
#SBATCH --partition=carter-compute
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/Evaluation/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/Evaluation/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

# Stage 2a: Evaluation on the full benchmark K sweep at sel_thresh=2.0.
# Mirrors Hon's guiding-light run params. perturbation_association, geneset,
# GO, and explained_variance (trait skipped pending OpenTargets file).

set -euo pipefail
START_TIME=$(date +%s)

OUT_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result"
RUN_NAME="050926_HTv2_20iter_5KHVG_torch_halsvar_batch"
RUN_DIR="$OUT_DIR/$RUN_NAME"

echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"
echo "Run dir: $RUN_DIR"

mkdir -p "$RUN_DIR/Evaluation/logs"

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate
export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

python -u /cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src/Stage2_Evaluation/A_Metrics/Slurm_Version/cNMF_evaluation_pipeline.py \
    --out_dir "$OUT_DIR" \
    --run_name "$RUN_NAME" \
    --X_normalized_path "$RUN_DIR/Inference/cnmf_tmp/Inference.norm_counts.h5ad" \
    --Perform_perturbation \
    --Perform_geneset \
    --Perform_explained_variance \
    --data_key "rna" \
    --prog_key "cNMF" \
    --categorical_key "sample" \
    --gene_names_key "symbol" \
    --guide_names_key "guide_names" \
    --guide_targets_key "guide_targets" \
    --guide_assignment_key "guide_assignment" \
    --guide_annotation_key "non-targeting" \
    --organism "human" \
    --K 5 6 7 8 9 10 11 12 13 14 15 17 19 21 23 25 27 30 35 40 45 50 55 60 70 80 90 100 150 200 \
    --sel_thresh 2.0 \
    --FDR_method "StoreyQ"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
printf 'Done @ %s | Elapsed: %dh %dm %ds\n' "$(date)" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
