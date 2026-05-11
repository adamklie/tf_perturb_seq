#!/bin/bash
#SBATCH --job-name=convert_esc_perturbnmf
#SBATCH --partition=carter-compute
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Script/convert_esc_perturbnmf.%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

set -euo pipefail
echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate

python -u /cellar/users/aklie/projects/tf_perturb_seq/src/tf_perturb_seq/cnmf/h5mu_to_perturbnmf_h5ad.py \
    --in_h5mu  /cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/runs/sceptre_v1/pipeline_outputs/inference_mudata.h5mu \
    --out_h5ad /cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Data/ESC_sceptre_v1_perturbnmf.h5ad

echo "Done @ $(date)"
