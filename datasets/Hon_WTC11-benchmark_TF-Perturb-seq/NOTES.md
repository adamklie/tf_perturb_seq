# Hon_WTC11-benchmark_TF-Perturb-seq - Dataset Notes

## Overview
- **Analysis Set Accession**: (check portal)
- **Sample Metadata CSV**: `sample_metadata_gcp_2026_01_18.csv` (GCS paths)
- **Modalities**: scRNA, gRNA, hash (**includes cell hashing**)
- **Status**: ✅ FULLY CONFIGURED - Pipeline has been run with results

---

## Data Source Files

| File | Description |
|------|-------------|
| `sample_metadata.csv` | Original sample metadata |
| `sample_metadata_gcp_2026_01_18.csv` | GCS-ready CSV with uploaded paths |
| `sample_metadata_gcp_2026_01_18_patched.csv` | Patched version (see `sample_sheet_patch_notes.md`) |

---

## Key Files (Already on GCS)

### 1. Seqspecs
- **Count**: 3
- **Location**: `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-benchmark_TF-Perturb-seq/2026_01_18/`
- **Files**:
  - `rna_seqspec.yml` - for scRNA
  - `guide_seqspec.yml` - for gRNA
  - `hash_seqspec.yml` - for cell hashing
- **Local copies**: `downloads/rna_seqspec.yml`, `downloads/guide_seqspec.yml`, `downloads/hash_seqspec.yml`
- **Status**: ✅ Uploaded to GCS

### 2. Barcode Onlist
- **IGVF Accession**: `IGVFFI9487JPEN`
- **GCS path**: `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-benchmark_TF-Perturb-seq/2026_01_18/.../IGVFFI9487JPEN.tsv.gz`
- **Local copy**: `downloads/IGVFFI9487JPEN.tsv`
- **Status**: ✅ Uploaded to GCS

### 3. Guide Metadata
- **IGVF Accession**: `IGVFFI5765HMZH` (same as Gersbach datasets)
- **GCS path**: `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-benchmark_TF-Perturb-seq/2026_01_18/.../IGVFFI5765HMZH.tsv.gz`
- **Local copies**:
  - `downloads/IGVFFI5765HMZH.tsv.gz`
  - `downloads/IGVFFI5765HMZH.tsv`
- **Status**: ✅ Uploaded to GCS

### 4. Barcode Hashtag Map (for cell hashing)
- **IGVF Accession**: `IGVFFI3028BJFI`
- **GCS path**: `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-benchmark_TF-Perturb-seq/2026_01_18/.../IGVFFI3028BJFI.tsv.gz`
- **Local copy**: `downloads/IGVFFI3028BJFI.tsv.gz`
- **Status**: ✅ Uploaded to GCS

---

## Config Files

### `harmless_teacher.config`
Production config used for GCP pipeline runs.

**Key parameters:**
- `ENABLE_DATA_HASHING = true` (**has cell hashing**)
- `DUAL_GUIDE = false`
- `Multiplicity_of_infection = 'high'`
- `GUIDE_ASSIGNMENT_method = 'sceptre'` (different from Gersbach's cleanser)
- `GUIDE_ASSIGNMENT_capture_method = 'CROP-seq'`
- `reverse_complement_guides = true`
- `spacer_tag = "TTAGCTCTTAAAC"`
- `use_igvf_reference = true`

### Other config files
- `beautiful_tornado.config` - alternative config
- `params.yaml` - YAML params file

---

## Pipeline Results

Results available in `results/1_CRISPR_pipeline/`:
- `2026_01_21/` - Run with dashboard, inference results
- `2026_01_22/` - Additional run

Output files include:
- `pipeline_dashboard/dashboard.html`
- `inference_mudata.h5mu`
- `perturbo_cis_per_guide_output.tsv.gz`
- `perturbo_trans_per_guide_output.tsv.gz`
- Various QC plots and evaluation outputs

---

## Comparison with Gersbach Datasets

| Attribute | Hon | Gersbach GEM-Xv3 | Gersbach HTv2 |
|-----------|-----|------------------|---------------|
| Cell Hashing | ✅ Yes | ❌ No | ❌ No |
| Seqspec approach | 3 files (per modality) | 50 files (per sample) | 68 files (per sample) |
| Barcode Onlist | `IGVFFI9487JPEN` | `IGVFFI1697XAXI` | `IGVFFI9487JPEN` (same) |
| Guide Metadata | `IGVFFI5765HMZH` | `IGVFFI5765HMZH` (same) | `IGVFFI5765HMZH` (same) |
| Guide Assignment | SCEPTRE | CLEANSER | CLEANSER |
| MOI | high | low | low |
| GCS Status | ✅ Uploaded | ❌ Local paths | ❌ Local paths |

---

## Key Differences from Gersbach

1. **Seqspec organization**: Hon uses 3 seqspecs (one per modality) vs Gersbach's per-sample approach
2. **Cell hashing**: Hon has hashing, Gersbach doesn't
3. **Guide assignment**: Hon uses SCEPTRE, Gersbach uses CLEANSER
4. **Already uploaded**: Hon data is already on GCS and pipeline has been run

---

## Session Log

### 2026-01-23
- Dataset is fully configured with GCS paths
- Pipeline has been run with results available
- Uses same guide metadata (`IGVFFI5765HMZH`) as Gersbach datasets
- Uses same barcode onlist as Gersbach HTv2 (`IGVFFI9487JPEN`)
- Created this NOTES.md file
