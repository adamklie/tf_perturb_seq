# CHANGE
base_dir=/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq  # Base directory
sample_metadata=$base_dir/2025_10_01_updated_paths_sample_metadata.csv
PIPELINE_PATH=/cellar/users/aklie/opt/CRISPR_Pipeline

# Run it
cd $PIPELINE_PATH
nextflow run main.nf \
    -profile google \
    --input $sample_metadata \
    --outdir gs://igvf-pertub-seq-pipeline-data/adam_test/Hon_WTC11-benchmark_TF-Perturb-seq/2025_10_01/outs/ \
    -with-tower
