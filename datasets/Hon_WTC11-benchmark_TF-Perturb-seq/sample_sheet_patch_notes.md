# Sample Sheet Patch Summary

## File: `sample_metadata_gcp_2026_01_18_patched.csv`

**Base file:** `sample_metadata_gcp_2026_01_18.csv` (output from download pipeline)

**Date:** 2026-01-20

---

## Changes Made

| Change | Original Value | Patched Value | Reason |
|--------|----------------|---------------|--------|
| **Hash file path** (`barcode_hashtag_map`) | `gs://.../2024/12/20/73fc75aa-.../IGVFFI3028BJFI.tsv.gz` | `gs://.../patch/IGVFFI3028BJFI_with_header.tsv.gz` | Original file missing required `hash_id\tsequence` header |
| **Lane column** | `3`, `5`, `7`, `8` (various per row) | `` (empty for all rows) | Match working sample sheet structure |
| **Sequencing run column** | `1` | `1` (unchanged) | Already matched working sample sheet |

### Lane Column Details

The lane values were cleared to match the working `sample_sheet.csv` format:

| Row Type | Original Lane | Patched Lane |
|----------|---------------|--------------|
| scRNA rows | `3` | `` (empty) |
| hash rows | `5` | `` (empty) |
| gRNA rows | `7` or `8` | `` (empty) |

---

## Files Uploaded to GCP

| File | GCP Location | Description |
|------|--------------|-------------|
| `IGVFFI3028BJFI_with_header.tsv.gz` | `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-benchmark_TF-Perturb-seq/2026_01_18/patch/` | Hash barcode file with header added |

---

## Pipeline Script Update

| File | Change |
|------|--------|
| `bin/1_CRISPR_pipeline/3_run_CRISPR_pipeline.sh` | `SAMPLE_METADATA` now points to `sample_metadata_gcp_2026_01_18_patched.csv` |

---

## Verification Performed

- `IGVFFI5765HMZH.tsv.gz` (guide file) matches GCP version (`IGVFFI5765HMZH_benchmark_guide_final.tsv`) - identical content, 416 lines, same MD5
- Patched hash file uploaded and verified to have correct header
