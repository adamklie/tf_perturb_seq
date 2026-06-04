#!/bin/bash
#SBATCH --job-name=htv2_eval_trait
#SBATCH --partition=carter-compute
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=96G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/Evaluation/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/Evaluation/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

# Stage 2a — trait-enrichment-only follow-up for HTv2 testbed.
# Original eval (10724837) ran perturbation/geneset/GO/explained_variance.
# This pass only adds <K>_trait_enrichment.txt so Stage 3a kselection can run.

set -euo pipefail
START_TIME=$(date +%s)

OUT_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result"
RUN_NAME="050926_HTv2_20iter_5KHVG_torch_halsvar_batch"
RUN_DIR="$OUT_DIR/$RUN_NAME"
GWAS_DATA="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src/Stage2_Evaluation/Resources/OpenTargets_L2G_Filtered.csv.gz"

echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"

mkdir -p "$RUN_DIR/Evaluation/logs"

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate
export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

python -u /cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src/Stage2_Evaluation/A_Metrics/Slurm_Version/cNMF_evaluation_pipeline.py \
    --out_dir "$OUT_DIR" \
    --run_name "$RUN_NAME" \
    --Perform_trait \
    --gwas_data_path "$GWAS_DATA" \
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
