# Huangfu_WTC11-benchmark_TF-Perturb-seq

## Dataset

WTC11 benchmark TF Perturb-seq dataset from the Huangfu lab.

## Run Pipeline on GCP

```bash
cd $PROJECT_ROOT

# Generate samplesheet
bash datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/1_generate_per_sample_metadata.sh

# Transfer to GCP -- dry run first
DRY_RUN=true bash datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/2_upload_to_gcp.sh
bash datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/2_upload_to_gcp.sh

# Patch .tsv.gz files to unzipped versions (optional)
bash datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/patch_files.sh

# Run pipeline on GCP
RUN_IN_BACKGROUND=true bash datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/3_run_CRISPR_pipeline.sh

# Run QC pipeline
bash datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/4_run_qc_pipeline.sh
```



## Sample Metadata

| Version | File | Notes |
|---------|------|-------|
| Original | `sample_metadata_gcp_2026_01_30.csv` | Per-sample seqspecs from portal (gzipped `.yaml.gz`) |
| Patched | `sample_metadata_gcp_2026_01_30_patched.csv` | Decompressed seqspec/barcode/guide files to `patch/` directory |
| Patched v2 | `sample_metadata_gcp_2026_01_30_patched_v2.csv` | Updated barcode_onlist from `IGVFFI9487JPEN` to `IGVFFI4695IKAL` |
| New seqspecs | `sample_metadata_gcp_2026_02_15.csv` | Replaced seqspec column with new shared seqspecs from `scratch/bioinfolucas/new_03_seqspec/` |

### 2026-02-15: seqspec update
- Based on `sample_metadata_gcp_2026_01_30_patched_v2.csv`
- Replaced all seqspec paths:
  - scRNA: `huangfu_IGVFDS1889TBEY_rna.yaml`
  - gRNA: `huangfu_IGVFDS1889TBEY_guide.yaml`
- All other columns (barcode_onlist, guide_design, etc.) preserved from patched_v2 version
- Uploaded to `gs://igvf-pertub-seq-pipeline-data/Huangfu_WTC11-benchmark_TF-Perturb-seq/`

### 2026-02-13: barcode_onlist update
- Updated `barcode_onlist` from `IGVFFI9487JPEN` to `IGVFFI4695IKAL` across all rows in `sample_metadata.csv`
- Uploaded `IGVFFI4695IKAL` to GCS and decompressed to `patch/IGVFFI4695IKAL.tsv`
- Created `sample_metadata_gcp_2026_01_30_patched_v2.csv` with updated barcode_onlist paths

## Pipeline run reference

friendly_minksy -- my own run with an accompanying tower dashboard at friendly_minsky
pipeline_dashboard -- Gary's latest run (https://www.synapse.org/Synapse:syn72299969)
agitated_hoover -- my own run with an accompanying tower dashboard at agitated_hoover
silver_otter -- downloaded from Synapse (syn72386406) on 2026-02-17