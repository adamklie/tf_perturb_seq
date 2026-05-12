#!/bin/bash
#SBATCH --job-name=huangfu_de_eval_trait
#SBATCH --partition=carter-compute
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=96G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_de_torchcnmf_KskillA/Evaluation/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_de_torchcnmf_KskillA/Evaluation/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

# Stage 2a — trait-enrichment-only follow-up
# The original eval (10542096 / 10542010) ran perturbation/geneset/GO/explained_variance
# but skipped --Perform_trait because OpenTargets_L2G_Filtered.csv.gz wasn't local.
# That file is now available at external/PerturbNMF/src/Stage2_Evaluation/Resources/
# (downloaded from EngreitzLab/gene_network_evaluation/smk/resources/).
#
# This rerun ONLY runs --Perform_trait. Outputs: <K>_trait_enrichment.txt added
# alongside the existing per-K eval files. After this completes, the K-selection
# plot's trait panel will work.

set -euo pipefail
START_TIME=$(date +%s)

OUT_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result"
RUN_NAME="042926_huangfu_de_torchcnmf_KskillA"
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
    --K 30 50 60 80 100 200 250 300 \
    --sel_thresh 2.0 \
    --FDR_method "StoreyQ"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
printf 'Done @ %s | Elapsed: %dh %dm %ds\n' "$(date)" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
