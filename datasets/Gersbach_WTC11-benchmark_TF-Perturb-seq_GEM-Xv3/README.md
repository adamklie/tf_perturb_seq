# Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3

## Dataset

WTC11 benchmark TF Perturb-seq dataset from the Gersbach lab (GEM-X v3 chemistry).

## Run Pipeline on GCP

```bash
cd $PROJECT_ROOT

# Generate samplesheet
bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/1_generate_per_sample_metadata.sh

# Transfer to GCP -- dry run first
DRY_RUN=true bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2_upload_to_gcp.sh
bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2_upload_to_gcp.sh

# Run pipeline on GCP
RUN_IN_BACKGROUND=true bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/3_run_CRISPR_pipeline.sh

#
bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/4_run_qc_pipeline.sh
```


## Sample Metadata

| Version | File | Notes |
|---------|------|-------|
| Original | `sample_metadata_gcp_2026_02_04.csv` | Per-sample seqspecs from portal (gzipped `.yaml.gz`) |
| Patched | `sample_metadata_gcp_2026_02_04_patched.csv` | Decompressed seqspec/barcode/guide files to `patch/` directory |
| New seqspecs | `sample_metadata_gcp_2026_02_15.csv` | Replaced seqspec column with new shared seqspecs from `scratch/bioinfolucas/new_03_seqspec/` |

### 2026-02-15: seqspec update
- Based on `sample_metadata_gcp_2026_02_04_patched.csv`
- Replaced all seqspec paths:
  - scRNA: `alex_IGVFDS6237URFJ_rna_10x_5_gex.yaml`
  - gRNA: `alex_IGVFDS6237URFJ_guide_10x_5_gex.yaml`
- All other columns (barcode_onlist, guide_design, etc.) preserved from patched version
- Uploaded to `gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/`

## Pipeline run reference

frisky_rabbit -- from Synapse (https://www.synapse.org/Synapse:syn73579845) run by Alex Barrera with Sceptre guide inference
jade_dolphin -- downloaded from Synapse (syn73712373) on 2026-02-17
