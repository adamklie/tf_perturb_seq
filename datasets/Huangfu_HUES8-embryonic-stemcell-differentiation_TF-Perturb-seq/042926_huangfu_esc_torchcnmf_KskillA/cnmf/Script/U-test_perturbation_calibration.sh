#!/bin/bash
#SBATCH --job-name=huangfu_esc_utest_calibration
#SBATCH --partition=carter-compute
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_esc_torchcnmf_KskillA/Evaluation/Calibration/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_esc_torchcnmf_KskillA/Evaluation/Calibration/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

# Stage 2b: U-test perturbation calibration
# Generates a null distribution of perturbation tests by randomly designating
# subsets of non-targeting guides as "fake" targeting groups, then comparing
# against the real distribution.
#
# Uses --categorical_key sample (single value 'all') so the null is pooled
# across batches, matching the unstratified perturbation association in 2a.
# Run prepare_h5mu_for_eval.sh first.

set -euo pipefail
START_TIME=$(date +%s)

OUT_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result"
RUN_NAME="042926_huangfu_esc_torchcnmf_KskillA"
RUN_DIR="$OUT_DIR/$RUN_NAME"

echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"
echo "Run dir: $RUN_DIR"

mkdir -p "$RUN_DIR/Evaluation/Calibration/logs"

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate
export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

# mdata_guide_path: any of the 8 cNMF_*.h5mu has the same guide info; use K=30
MDATA_GUIDE_PATH="$RUN_DIR/Inference/adata/cNMF_30_2_0.h5mu"

# Guide annotation TSV: built from project's harmonized poolabcd file with
# 'guide_id' column renamed to 'guide_names' (required by U-test fake-test code
# at line 134: guide_target_[args.guide_names_key]).
GUIDE_ANNOTATION_TSV="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Data/guide_annotation.tsv"

# NOTE: --compute_real_perturbation_tests dropped — the eval pipeline already
# wrote per-K perturbation_association_results_all.txt files via the same
# compute_perturbation_association call. Running it again here would just
# duplicate work.

python -u /cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src/Stage2_Evaluation/B_Calibration/Slurm_version/U-test_perturbation_calibration/U-test_perturbation_calibration.py \
    --out_dir "$OUT_DIR" \
    --run_name "$RUN_NAME" \
    --mdata_guide_path "$MDATA_GUIDE_PATH" \
    --guide_annotation_path "$GUIDE_ANNOTATION_TSV" \
    --guide_annotation_key "non-targeting" \
    --data_key "rna" \
    --prog_key "cNMF" \
    --categorical_key "sample" \
    --guide_names_key "guide_names" \
    --guide_targets_key "guide_targets" \
    --guide_assignment_key "guide_assignment" \
    --organism "human" \
    --FDR_method "StoreyQ" \
    --number_run 50 \
    --number_guide 6 \
    --components 30 50 60 80 100 200 250 300 \
    --sel_thresh 2.0 \
    --compute_fake_perturbation_tests

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
printf 'Done @ %s | Elapsed: %dh %dm %ds\n' "$(date)" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
