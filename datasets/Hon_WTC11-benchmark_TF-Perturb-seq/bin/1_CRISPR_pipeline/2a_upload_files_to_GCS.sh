# CHANGE
PROJECT=igvf-pertub-seq-pipeline
base_dir=/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq  # Base directory
script=/cellar/users/aklie/opt/CRISPR_Pipeline/download_development/download_igvf.py  # IGVF CRISPR Pipeline script for retrieving fastqs and uploading to GCS

# Authenticate and set project
gcloud auth login
gcloud config set project $PROJECT
gcloud config list

# Run it
python3 $script \
  --sample $base_dir/sample_metadata.tsv \
  --gcs-upload \
  --gcs-bucket igvf-pertub-seq-pipeline-data \
  --gcs-prefix adam_test/Hon_WTC11-benchmark_TF-Perturb-seq/2025_08_14/ \
  --output-dir $base_dir/downloads
  