#!/bin/bash
#SBATCH --job-name=huangfu_de_program_200
#SBATCH --partition=carter-compute
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=700G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_de_torchcnmf_KskillA/Plot/Program_200_2_0/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_de_torchcnmf_KskillA/Plot/Program_200_2_0/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

# Stage 3b: Program analysis plots for K=200 (initial K choice — can revise later)
# Reads h5mu + Stage 2a perturbation_association_results_<sample>.txt + GO enrichment.
# Categorical key 'sample' single-value 'all' so perturbation file is read as
# <K>_perturbation_association_results_all.txt.

set -euo pipefail
START_TIME=$(date +%s)

BASE_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_de_torchcnmf_KskillA"
LOG_DIR="${BASE_DIR}/Plot/Program_200_2_0"
K=200

echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"

mkdir -p "$LOG_DIR/logs"

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate
export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

python -u /cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src/Stage3_Interpretation/A_Plotting/Slurm_Version/cNMF_program_analysis.py \
    --mdata_path "${BASE_DIR}/Inference/adata/cNMF_${K}_2_0.h5mu" \
    --perturb_path_base "${BASE_DIR}/Evaluation/${K}_2_0/${K}_perturbation_association_results" \
    --GO_path "${BASE_DIR}/Evaluation/${K}_2_0/${K}_GO_term_enrichment.txt" \
    --top_program 5 \
    --p_value 0.05 \
    --pdf_save_path "$LOG_DIR" \
    --PDF \
    --sample all \
    --square_plots \
    --figsize 35 20 \
    --categorical_key "sample" \
    --subsample_frac 0.1 \
    --gene_name_key "symbol"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
printf 'Done @ %s | Elapsed: %dh %dm %ds\n' "$(date)" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
