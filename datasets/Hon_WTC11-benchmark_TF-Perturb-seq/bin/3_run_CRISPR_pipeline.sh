# CHANGE
base_dir=/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq  # Base directory
sample_metadata=$base_dir/sample_metadata.tsv
PIPELINE_PATH=/cellar/users/aklie/opt/CRISPR_Pipeline

# Run it
cd $PIPELINE_PATH
nextflow run main.nf \
    -profile google \
    --input $sample_metadata \
    --outdir gs://igvf-pertub-seq-pipeline-data/adam_test/Hon_WTC11-benchmark_TF-Perturb-seq/2025_08_14/outs/ \
    -with-tower
