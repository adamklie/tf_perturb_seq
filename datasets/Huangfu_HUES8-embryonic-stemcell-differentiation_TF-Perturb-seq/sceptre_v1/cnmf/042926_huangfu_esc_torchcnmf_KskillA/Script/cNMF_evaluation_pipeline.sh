#!/bin/bash
#SBATCH --job-name=huangfu_esc_eval
#SBATCH --partition=carter-compute
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_esc_torchcnmf_KskillA/Evaluation/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_esc_torchcnmf_KskillA/Evaluation/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

# Stage 2a: Evaluation
# Runs perturbation association, gene-set/GO enrichment, and explained variance
# at K = 30,50,60,80,100,200,250,300 with sel_thresh = 2.0.
#
# NOTES:
# - --categorical_key sample (single value 'all'): perturbation association is
#   NOT stratified by batch. Output = single <K>_perturbation_association_results_all.txt
#   per K. Run prepare_h5mu_for_eval.sh first.
# - --Perform_categorical SKIPPED (degenerate against single-value column).
# - --Perform_trait SKIPPED for now: requires OpenTargets_L2G_Filtered.csv.gz,
#   which is only on Sherlock. TODO: obtain file from Yusen Mo and re-run with
#   --Perform_trait + --gwas_data_path.

set -euo pipefail
START_TIME=$(date +%s)

OUT_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result"
RUN_NAME="042926_huangfu_esc_torchcnmf_KskillA"
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
    --K 30 50 60 80 100 200 250 300 \
    --sel_thresh 2.0 \
    --FDR_method "StoreyQ"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
printf 'Done @ %s | Elapsed: %dh %dm %ds\n' "$(date)" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
