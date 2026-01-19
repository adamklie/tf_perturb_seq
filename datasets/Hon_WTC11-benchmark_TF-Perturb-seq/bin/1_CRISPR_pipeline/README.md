uv run bash datasets/Hon_WTC11-benchmark_TF-Perturb-seq/bin/1_CRISPR_pipeline/1_generate_per_sample_metadata.sh
DRY_RUN=true uv run bash datasets/Hon_WTC11-benchmark_TF-Perturb-seq/bin/1_CRISPR_pipeline/2_upload_to_gcp.sh
RUN_IN_BACKGROUND=true bash /Users/adamklie/Desktop/projects/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq/bin/1_CRISPR_pipeline/3_run_CRISPR_pipeline.sh