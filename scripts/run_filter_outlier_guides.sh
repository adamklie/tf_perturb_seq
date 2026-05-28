#!/bin/bash
#SBATCH --job-name=filter_outlier_guides
#SBATCH --mem=300G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --partition=common
#SBATCH --output=/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/post_pipeline_processing/logs/filter_outlier_guides_%j.out
#SBATCH --error=/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/post_pipeline_processing/logs/filter_outlier_guides_%j.err

python filter_outlier_guides.py \
    --input  /hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/post_pipeline_processing/inference_mudata.h5mu \
    --output /hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/post_pipeline_processing/inference_mudata_outlier_guide_filt.h5mu \
    --non_targeting_outliers /hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/energy_distance/data_hep/non_targeting_outlier_table.csv \
    --targeting_outliers     /hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/energy_distance/data_hep/targeting_outlier_table.csv \
    --fdr_threshold 0.05