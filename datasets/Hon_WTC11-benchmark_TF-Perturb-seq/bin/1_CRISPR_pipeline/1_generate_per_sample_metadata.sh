# CHANGE
base_dir=/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq  # Base directory
script=/cellar/users/aklie/opt/CRISPR_Pipeline/download_development/generate_per_sample.py  # IGVF CRISPR pipeline make initial sample sheet script
ACCESSION=IGVFDS4761PYUO  # Analysis set ID in IGVF portal

# Run it
python3 $script \
  --accession $ACCESSION \
  --output $base_dir/sample_metadata.tsv \
  --hash_seqspec $base_dir/bin/1_CRISPR_pipeline/seqspec/yaml_files/hash_seqspec.yml \
  --rna_seqspec $base_dir/bin/1_CRISPR_pipeline/seqspec/yaml_files/rna_seqspec.yml \
  --sgrna_seqspec $base_dir/bin/1_CRISPR_pipeline/seqspec/yaml_files/guide_seqspec.yml
