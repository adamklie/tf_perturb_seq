# Gersbach WTC11 Benchmark TF Perturb-seq (HT v2)

WTC11 benchmark TF Perturb-seq dataset from the Gersbach lab (HT v2 chemistry). Data sourced from IGVF portal (analysis set: IGVFDS6237URFJ).

## Files

| File | Description |
|------|-------------|
| `sample_metadata_gcp_2026_02_15.csv` | Latest samplesheet with updated seqspecs. This is the input to the CRISPR pipeline. |
| `Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2_2026_03_11.config` | Nextflow config for the CRISPR pipeline (supports local, SLURM, and GCP profiles). |
| `_samplesheets/` | Intermediate samplesheets from earlier pipeline steps (`sample_metadata.csv`, `sample_metadata_gcp_2026_02_04.csv`, `sample_metadata_gcp_2026_02_04_patched.csv`, `IGVFDS6237URFJ.updated.v3.tsv`). |
| `runs/` | Current pipeline dashboard outputs. |
| ~~`results/`~~ | Moved to `scratch/2025_11_26/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2_results/`. Stale pipeline outputs from earlier runs. |

## Scripts (run in order)

| Script | Description |
|--------|-------------|
| `1_generate_per_sample_metadata.sh` | Pulls per-sample metadata from the IGVF portal. |
| `2_upload_to_gcp.sh` | Transfers IGVF files to GCS. Supports `DRY_RUN=true`. |
| `3_patch_gcp_files.sh` | Decompresses `.gz` seqspecs, onlists, and guide designs on GCS. |
| `4_run_CRISPR_pipeline.sh` | Runs the CRISPR Nextflow pipeline on GCP. Supports `RUN_IN_BACKGROUND=true`. |
| `5_run_qc_pipeline.sh` | Runs QC on pipeline output (`inference_mudata.h5mu`). Supports `--dry-run`. |
| `6_run_cnmf.sh` | Runs cNMF gene program discovery. Supports `--dry-run`, `--slurm`, `--phase`. |

## Running the pipeline

```bash
cd $PROJECT_ROOT

# 1. Generate samplesheet
bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/1_generate_per_sample_metadata.sh

# 2. Transfer to GCP (dry run first)
DRY_RUN=true bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2_upload_to_gcp.sh
bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2_upload_to_gcp.sh

# 3. Patch GCP files
bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/3_patch_gcp_files.sh

# 4. Run CRISPR pipeline
RUN_IN_BACKGROUND=true bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/4_run_CRISPR_pipeline.sh

# 5. Run QC
bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/5_run_qc_pipeline.sh

# 6. Run cNMF
bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/6_run_cnmf.sh
```

## Samplesheet history

| Version | File | Notes |
|---------|------|-------|
| Original | `sample_metadata_gcp_2026_02_04.csv` | Per-sample seqspecs from portal (gzipped `.yaml.gz`) |
| Patched | `sample_metadata_gcp_2026_02_04_patched.csv` | Decompressed seqspec/barcode/guide files to `patch/` directory |
| Latest | `sample_metadata_gcp_2026_02_15.csv` | Replaced seqspec column with new shared seqspecs (`alex_IGVFDS6673ZFFG_rna_chomium5.yaml` for scRNA, `alex_IGVFDS6673ZFFG_guide_HTV2_chomium5.yaml` for gRNA) |

## Pipeline run reference

- **red_lettuce** -- from Synapse (syn73579418), run by Alex Barrera with Sceptre guide inference
- **amber_sparrow** -- downloaded from Synapse (syn73712262) on 2026-02-17
