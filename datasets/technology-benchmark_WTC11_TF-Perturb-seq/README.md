# TF-Perturb-seq Technology Benchmark

Master tracking document for running the CRISPR Pipeline on WTC11 benchmark datasets.

## Issues

- Gersbach seqspecs not yet validated on portal and may have issues with pipeline?
- Engreitz seqspecs not yet available on portal
- Manual modifications to sample sheet(see below)
The following summarizes what I had to do to the Hon_WTC11-benchmark_TF-Perturb-seq sample sheet to get the pipeline to run successfully on GCP on 2026-01-18.

| Change | Original Value | Patched Value | Reason |
|--------|----------------|---------------|--------|
| **Hash file path** (`barcode_hashtag_map`) | `gs://.../2024/12/20/73fc75aa-.../IGVFFI3028BJFI.tsv.gz` | `gs://.../patch/IGVFFI3028BJFI_with_header.tsv.gz` | Original file missing required `hash_id\tsequence` header |
| **Barcode onlist path** (`barcode_onlist`) | `gs://.../2024/12/05/996cf65e-.../IGVFFI9487JPEN.tsv.gz` | `gs://.../patch/IGVFFI9487JPEN.tsv` | Pipeline renames file to `.txt` without decompressing; `kb count -w` fails on gzipped content |
| **Guide design path** (`guide_design`) | `gs://.../2025/08/04/cc3d00eb-.../IGVFFI5765HMZH.tsv.gz` | `gs://.../patch/IGVFFI5765HMZH.tsv` | Precautionary unzip (see note below) |
| **Lane column** | `3`, `5`, `7`, `8` (various per row) | `` (empty for all rows) | Match working sample sheet structure |
| **Sequencing run column** | `1` | `1` (unchanged) | Already matched working sample sheet |

### Lane Column Details

The lane values were cleared to match the working `sample_sheet.csv` format:

| Row Type | Original Lane | Patched Lane |
|----------|---------------|--------------|
| scRNA rows | `3` | `` (empty) |
| hash rows | `5` | `` (empty) |
| gRNA rows | `7` or `8` | `` (empty) |

### Barcode Onlist Issue Details

The pipeline's `seqSpecParser` module copies the barcode file as-is:
```bash
cp "${barcode_file}" barcode_file.txt
```

This removes the `.gz` extension but doesn't decompress, so `kb count -w barcode_file.txt` fails when trying to read gzipped binary as plain text barcodes.

### Guide Design Note
The guide design file was unzipped out of caution, following the same pattern as the barcode onlist issue.

## Dataset Overview

| Dataset | Lab | Technology | Status |
|---------|-----|------------|---------|
| [Engreitz](../Engreitz_WTC11-benchmark_TF-Perturb-seq/README.md) | Engreitz | CC Perturb-seq | Waiting on seqspecs on portal to start |
| [Gersbach GEM-Xv3](../Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/README.md) | Gersbach | 10x GEM-X 3' | Waiting for seqspecs to validate on portal to start |
| [Gersbach HTv2](../Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/README.md) | Gersbach | 10x HT v2 | Waiting for seqspecs to validate on portal to start |
| [Hon](../Hon_WTC11-benchmark_TF-Perturb-seq/README.md) | Hon | 10x 5' HT v2 w/ HTO | Pipeline run completed with results |
| [Huangfu](../Huangfu_WTC11-benchmark_TF-Perturb-seq/README.md) | Huangfu | 10x 3' v3 | GCP upload completed, need to run pipeline |
---

## Engreitz_WTC11-benchmark_TF-Perturb-seq

- **Analysis Set Accession**: `IGVFDS6673ZFFG`
- **Sample Metadata CSV**: Not yet generated
- **Guide metadata**: `IGVFFI5765HMZH`
- **Seqspecs**: Validating on portal
- **Barcode Onlist**: TODO
- **Status**: Portal missing seqspecs, need to follow up with Engreitz lab
- **Key pipeline parameters:**
   - `ENABLE_DATA_HASHING = false`
   - `Multiplicity_of_infection = 'high'`

## Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3

- **Analysis Set Accession**: `IGVFDS6673ZFFG`
- **Sample Metadata CSV**: TODO
- **Guide Metadata**: TODO
- **Barcode Onlist**: TODO
- **Status**: Waiting for seqspecs to validate on portal to start
- **Key Parameters:**
   - `ENABLE_DATA_HASHING = false`

## Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2

- **Analysis Set Accession**: `IGVFDS6237URFJ`
- **Sample Metadata CSV**: TODO
- **Guide Metadata** TODO
- **Barcode Onlist**: TODO
- **Status**: Waiting for seqspecs to validate on portal to start
- **Key Parameters:**
   - `ENABLE_DATA_HASHING = false`

## Hon_WTC11-benchmark_TF-Perturb-seq

- **Analysis Set Accession**: `IGVFDS4761PYUO`
- **Sample Metadata CSV**: `datasets/Hon_WTC11-benchmark_TF-Perturb-seq/sample_metadata_gcp_2026_01_18_patched.csv`
- **Guide Metadata** Gunzipped `IGVFFI5765HMZH`
- **Seqspecs**: Not from portal, using generic ones that I uploaded to GCS (see .yml files)
- **Barcode Onlist**: Gunzipped `IGVFFI9487JPEN`
- **Hashing Map**: Patched `IGVFFI3028BJFI` with header
- **Status**: Pipeline successfully run to completion with results
- **Key Parameters:**
   - `ENABLE_DATA_HASHING = true`
   - `Multiplicity_of_infection = 'high'`
- **Successful Runs (see configs)**
   - harmless_teacher: initial run with hashing
   - beautiful_tornado: harmless teacher but with cleanser instead of sceptre
   - introverted_frankenstein: disabled hashing, back to sceptre

## Huangfu_WTC11-benchmark_TF-Perturb-seq

- **Analysis Set Accession**: `IGVFDS6673ZFFG`
- **Sample Metadata CSV**: `datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/sample_metadata_gcp_2026_01_30_patched.csv`
- **Guide Metadata** Gunzipped `IGVFFI5765HMZH`
- **Barcode Onlist**: Gunzipped 
- **Modalities**: scRNA, gRNA (no cell hashing)
- **Status**: GCP set-up, need to test running pipeline
- **Key Parameters:**
   - `ENABLE_DATA_HASHING = false`
   - `DUAL_GUIDE = false`
   - `is_10x3v3 = true` (**unique to this dataset!**)
   - `Multiplicity_of_infection = 'high'`
   - `GUIDE_ASSIGNMENT_method = 'sceptre'`
   - `GUIDE_ASSIGNMENT_capture_method = 'CROP-seq'`
   - `use_igvf_reference = true`
