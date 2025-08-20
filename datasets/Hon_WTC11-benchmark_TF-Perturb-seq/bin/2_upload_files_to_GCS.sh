
# 
gcloud auth login
gcloud config set project igvf-pertub-seq-pipeline
gcloud config list


base_dir=/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq
script=/cellar/users/aklie/opt/CRISPR_Pipeline/download_development/download_igvf.py


python3 $script \
  --sample $base_dir/sample_metadata.tsv \
  --gcs-upload \
  --gcs-bucket igvf-pertub-seq-pipeline-data \
  --gcs-prefix adam_test/Hon_WTC11-benchmark_TF-Perturb-seq/2025_08_14/ \
  --output-dir $base_dir/downloads
  