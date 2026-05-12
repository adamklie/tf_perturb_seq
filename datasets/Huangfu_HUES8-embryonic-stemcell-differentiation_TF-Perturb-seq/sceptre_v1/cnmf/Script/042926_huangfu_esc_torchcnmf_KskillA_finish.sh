#!/bin/bash
#SBATCH --job-name=042926_huangfu_esc_torchcnmf_KskillA_finish
#SBATCH --partition=carter-compute
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_esc_torchcnmf_KskillA/Inference/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_esc_torchcnmf_KskillA/Inference/logs/%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

# Re-runs ONLY the compile_annotation + diagnostic_plots stages after the original
# job (9890793) crashed at annotate_genes_to_excel due to missing openpyxl.
# Factorize + refit outputs from the original run are reused as-is.

set -euo pipefail
echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"

source "/cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate"
export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

python -u /cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src/Inference/torch-cNMF/Slurm_Version/torch_cnmf_inference_pipeline.py \
    --counts_fn "/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Data/ESC_sceptre_v1_perturbnmf.h5ad" \
    --output_directory "/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result" \
    --run_name "042926_huangfu_esc_torchcnmf_KskillA" \
    --species human \
    --K 30 50 60 80 100 200 250 300 \
    --numiter 10 \
    --numhvgenes 2000 \
    --sel_thresh 2.0 \
    --seed 14 \
    --algo halsvar \
    --mode batch \
    --tol 1e-4 \
    --categorical_key batch \
    --gene_names_key symbol \
    --run_compile_annotation --run_diagnostic_plots

echo "Done @ $(date)"
