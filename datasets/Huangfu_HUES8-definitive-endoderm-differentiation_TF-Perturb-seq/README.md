# Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq

# Run Pipeline on GCP
```bash
cd $project_root

# Generate samplesheet
uv run bash datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/bin/1_CRISPR_pipeline/1_generate_per_sample_metadata.sh 

# Transfer to GCP -- dry run first
DRY_RUN=true uv run bash datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/bin/2_upload_to_gcp.sh
uv run bash datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/bin/2_upload_to_gcp.sh

# Run pipeline on GCP
RUN_IN_BACKGROUND=true bash datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/bin/1_CRISPR_pipeline/3_run_CRISPR_pipeline.sh
```