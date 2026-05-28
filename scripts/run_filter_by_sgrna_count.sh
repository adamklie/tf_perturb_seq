#!/bin/bash
#SBATCH --job-name=filter_sgrna
#SBATCH --mem=300G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --partition=common
#SBATCH --output=/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/post_pipeline_processing/logs/filter_sgrna_%j.out
#SBATCH --error=/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/post_pipeline_processing/logs/filter_sgrna_%j.err

python /hpc/group/gersbachlab/seg95/tf_perturb_seq/scripts/filter_mudata_by_sgrna_count.py \
    --input /hpc/group/gersbachlab/seg95/CRISPR_Pipeline/helen_output_poolabcd_prod_v10/inference_mudata.h5mu \
    --output /hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/post_pipeline_processing/inference_mudata.h5mu \
    --max-guides 15 \
    --log /hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/post_pipeline_processing/logs/filter_sgrna.log