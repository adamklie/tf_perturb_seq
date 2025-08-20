base_dir=/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq
script=/cellar/users/aklie/opt/CRISPR_Pipeline/download_development/generate_per_sample.py

python3 $script \
  --accession IGVFDS4761PYUO \
  --output $base_dir/sample_metadata.tsv \
  --hash_seqspec $base_dir/yaml_files/hash_seq_spec.yaml \
  --rna_seqspec $base_dir/yaml_files/rna_seq_spec.yaml \
  --sgrna_seqspec $base_dir/yaml_files/sgrna_seq_spec.yaml
  