# Engreitz_WTC11-benchmark_TF-Perturb-seq

## Dataset

WTC11 benchmark TF Perturb-seq dataset from the Engreitz lab.

## Run Pipeline on GCP

```bash
cd $PROJECT_ROOT

# Generate samplesheet
bash datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/1_generate_per_sample_metadata.sh

# Transfer to GCP -- dry run first
DRY_RUN=true bash datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/2_upload_to_gcp.sh
bash datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/2_upload_to_gcp.sh

# Patch GCP files (decompress .gz seqspecs, onlists, guide designs) -- dry run first
DRY_RUN=true bash datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/patch_files.sh
bash datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/patch_files.sh

# Validate all GCS paths in the patched samplesheet exist
python3 scripts/validate_gcp_paths.py --input datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/sample_metadata_gcp_2026_02_26_patched.csv

# Run pipeline on GCP
RUN_IN_BACKGROUND=true bash datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/3_run_CRISPR_pipeline.sh

# Run QC pipeline
bash datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/4_run_qc_pipeline.sh
```

## Pipeline run reference

cobalt_heron -- downloaded from Synapse (syn73615563) on 2026-02-17
