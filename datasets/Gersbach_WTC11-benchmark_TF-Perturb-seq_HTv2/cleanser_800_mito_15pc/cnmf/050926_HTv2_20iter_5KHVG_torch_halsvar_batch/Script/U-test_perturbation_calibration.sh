#!/bin/bash
#SBATCH --job-name=htv2_utest_calibration
#SBATCH --partition=carter-compute
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/Evaluation/Calibration/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/Evaluation/Calibration/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

# Stage 2b: U-test perturbation calibration on the full benchmark K sweep.
# 50 fake-targeting iterations × 30 K. Pulls guide_annotation.tsv (416 rows:
# 332 targeting + 84 non-targeting [30 NT + 54 OR negative controls]).

set -euo pipefail
START_TIME=$(date +%s)

OUT_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result"
RUN_NAME="050926_HTv2_20iter_5KHVG_torch_halsvar_batch"
RUN_DIR="$OUT_DIR/$RUN_NAME"

echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"
echo "Run dir: $RUN_DIR"

mkdir -p "$RUN_DIR/Evaluation/Calibration/logs"

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate
export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

# Any K's h5mu has the same guide info; use lowest K (=5) for fastest read
MDATA_GUIDE_PATH="$RUN_DIR/Inference/adata/cNMF_5_2_0.h5mu"
GUIDE_ANNOTATION_TSV="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Data/guide_annotation.tsv"

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
    --components 5 6 7 8 9 10 11 12 13 14 15 17 19 21 23 25 27 30 35 40 45 50 55 60 70 80 90 100 150 200 \
    --sel_thresh 2.0 \
    --compute_fake_perturbation_tests

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
printf 'Done @ %s | Elapsed: %dh %dm %ds\n' "$(date)" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
