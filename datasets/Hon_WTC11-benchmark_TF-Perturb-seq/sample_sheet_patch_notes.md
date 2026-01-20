# Sample Sheet Patch Summary

## File: `sample_metadata_gcp_2026_01_18_patched.csv`

**Base file:** `sample_metadata_gcp_2026_01_18.csv` (output from download pipeline)

**Date:** 2026-01-20

---

## Changes Made

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

**PRECAUTIONARY CHANGE:** The guide design file was unzipped out of caution, following the same pattern as the barcode onlist issue. This should be investigated to determine:
1. Whether the pipeline properly handles gzipped guide design files
2. If the download pipeline should output unzipped tabular files by default

---

## Files Uploaded to GCP

| File | GCP Location | Description |
|------|--------------|-------------|
| `IGVFFI3028BJFI_with_header.tsv.gz` | `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-benchmark_TF-Perturb-seq/2026_01_18/patch/` | Hash barcode file with header added |
| `IGVFFI9487JPEN.tsv` | `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-benchmark_TF-Perturb-seq/2026_01_18/patch/` | Unzipped barcode onlist file |
| `IGVFFI5765HMZH.tsv` | `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-benchmark_TF-Perturb-seq/2026_01_18/patch/` | Unzipped guide design file (precautionary) |

---

## Pipeline Script Update

| File | Change |
|------|--------|
| `bin/1_CRISPR_pipeline/3_run_CRISPR_pipeline.sh` | `SAMPLE_METADATA` now points to `sample_metadata_gcp_2026_01_18_patched.csv` |

---

## Verification Performed

- `IGVFFI5765HMZH.tsv.gz` (guide file) matches GCP version (`IGVFFI5765HMZH_benchmark_guide_final.tsv`) - identical content, 416 lines, same MD5
- Patched hash file uploaded and verified to have correct header
- Unzipped barcode onlist file uploaded (12 MB)
- Unzipped guide design file uploaded (67.5 KB)

---

## TODO: Investigate

- [ ] Determine if CRISPR pipeline handles gzipped tabular files (guide_design, barcode_onlist, barcode_hashtag_map)
- [ ] Consider updating download pipeline to output unzipped tabular files
- [ ] Fix hash file header issue at source (IGVF portal)
