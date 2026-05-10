#!/bin/bash
#SBATCH --job-name=huangfu_esc_perturbgene_200
#SBATCH --partition=carter-compute
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=256G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_esc_torchcnmf_KskillA/Plot/Perturb_gene_200_2_0/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_esc_torchcnmf_KskillA/Plot/Perturb_gene_200_2_0/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

# Stage 3c: Perturbed gene analysis plots for K=200

set -euo pipefail
START_TIME=$(date +%s)

BASE_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_esc_torchcnmf_KskillA"
LOG_DIR="${BASE_DIR}/Plot/Perturb_gene_200_2_0"
K=200

echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"

mkdir -p "$LOG_DIR/logs"

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate
export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

python -u /cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src/Stage3_Interpretation/A_Plotting/Slurm_Version/cNMF_perturbed_gene_analysis.py \
    --mdata_path "${BASE_DIR}/Inference/adata/cNMF_${K}_2_0.h5mu" \
    --perturb_path_base "${BASE_DIR}/Evaluation/${K}_2_0/${K}_perturbation_association_results" \
    --top_n_programs 10 \
    --perturb_target_col "target_name" \
    --perturb_program_col "program_name" \
    --perturb_log2fc_col "log2FC" \
    --top_corr_genes 5 \
    --significance_threshold 0.05 \
    --volcano_log2fc_min -0.00 \
    --volcano_log2fc_max 0.00 \
    --save_path "$LOG_DIR" \
    --square_plots \
    --figsize 35 20 \
    --sample all \
    --PDF \
    --n_processes -1 \
    --umap_dot_size 10 \
    --data_key "rna" \
    --prog_key "cNMF" \
    --categorical_key "sample" \
    --gene_name_key "symbol" \
    --control_target_name "non-targeting" \
    --subsample_frac 0.1 \
    --parallel

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
printf 'Done @ %s | Elapsed: %dh %dm %ds\n' "$(date)" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
