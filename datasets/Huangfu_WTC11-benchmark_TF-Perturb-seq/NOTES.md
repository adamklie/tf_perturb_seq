# Huangfu_WTC11-benchmark_TF-Perturb-seq - Dataset Notes

## Overview
- **Sample Metadata CSV**: `downloads/2025_11_26_huangfu_cloud_sample_metadata.csv`
- **Modalities**: scRNA, gRNA (no cell hashing)
- **Status**: Partially configured - data on GCS in `beer_test/` location, seqspecs uploaded to new location

---

## Data Source Files

| File | Description |
|------|-------------|
| `downloads/2025_11_26_huangfu_cloud_sample_metadata.csv` | CSV with GCS paths (in `beer_test/` location) |
| `sample_metadata.csv` | Local sample metadata |

---

## Key Files

### 1. Seqspecs
- **Count**: 8 (4 reps × 2 modalities: rna + guide)
- **Pattern**: `michael-beer_Perturb-benchmark-WTC11-[guide-sequencing-]rep{1-4}.yaml`
- **Current GCS location (in CSV)**: `gs://igvf-pertub-seq-pipeline-data/beer_test/huangfu_WTC11_benchmark/yaml_files/`
- **New GCS location**: `gs://igvf-pertub-seq-pipeline-data/Huangfu_WTC11-benchmark_TF-Perturb-seq/seqspecs/` ✅ Uploaded
- **Local copies**: `downloads/*.yaml` (8 files)
- **Status**: ✅ Uploaded to new location, CSV needs path update

**Seqspec files:**
- `michael-beer_Perturb-benchmark-WTC11-rep1.yaml` (scRNA rep1)
- `michael-beer_Perturb-benchmark-WTC11-rep2.yaml` (scRNA rep2)
- `michael-beer_Perturb-benchmark-WTC11-rep3.yaml` (scRNA rep3)
- `michael-beer_Perturb-benchmark-WTC11-rep4.yaml` (scRNA rep4)
- `michael-beer_Perturb-benchmark-WTC11-guide-sequencing-rep1.yaml` (gRNA rep1)
- `michael-beer_Perturb-benchmark-WTC11-guide-sequencing-rep2.yaml` (gRNA rep2)
- `michael-beer_Perturb-benchmark-WTC11-guide-sequencing-rep3.yaml` (gRNA rep3)
- `michael-beer_Perturb-benchmark-WTC11-guide-sequencing-rep4.yaml` (gRNA rep4)

### 2. Barcode Onlist
- **IGVF Accession**: `IGVFFI4695IKAL` (**different from other datasets!**)
- **GCS path**: `gs://igvf-pertub-seq-pipeline-data/beer_test/huangfu_WTC11_benchmark/yaml_files/IGVFFI4695IKAL.tsv`
- **Local copy**: NOT FOUND in downloads/
- **Status**: ⚠️ Needs to be downloaded or moved to new location

### 3. Guide Metadata
- **File**: `benchmark_guide_metadata_v1.tsv` (**different from other datasets which use v3_rc!**)
- **GCS path**: `gs://igvf-pertub-seq-pipeline-data/beer_test/huangfu_WTC11_benchmark/yaml_files/benchmark_guide_metadata_v1.tsv`
- **Local copy**: NOT FOUND in downloads/
- **Status**: ⚠️ Needs to be downloaded or moved to new location

---

## Config Files

### `bin/nextflow.config`

**Key parameters:**
- `ENABLE_DATA_HASHING = false`
- `DUAL_GUIDE = false`
- `is_10x3v3 = true` (**unique to this dataset!**)
- `Multiplicity_of_infection = 'high'`
- `GUIDE_ASSIGNMENT_method = 'sceptre'`
- `GUIDE_ASSIGNMENT_capture_method = 'CROP-seq'`
- `use_igvf_reference = true`

---

## Comparison with Other Datasets

| Attribute | Huangfu | Hon | Gersbach |
|-----------|---------|-----|----------|
| Cell Hashing | ❌ No | ✅ Yes | ❌ No |
| Seqspec approach | 8 files (per rep) | 3 files (per modality) | 50-68 files (per sample) |
| Barcode Onlist | `IGVFFI4695IKAL` | `IGVFFI9487JPEN` | `IGVFFI1697XAXI` / `IGVFFI9487JPEN` |
| Guide Metadata | `v1.tsv` | `IGVFFI5765HMZH` (v3_rc) | `IGVFFI5765HMZH` (v3_rc) |
| Guide Assignment | SCEPTRE | SCEPTRE | CLEANSER |
| MOI | high | high | low |
| 10x 3' v3 | ✅ Yes | ❌ No | ❌ No |
| GCS Status | Partial (`beer_test/`) | ✅ Production location | ❌ Local paths |

---

## Action Items

### Files to Relocate/Download
| File | Current Location | Target Location | Status |
|------|------------------|-----------------|--------|
| Seqspecs (8) | `beer_test/...` | `Huangfu_.../seqspecs/` | ✅ Done |
| Barcode onlist | `beer_test/...` | `Huangfu_.../onlists/` | ❌ TODO |
| Guide metadata | `beer_test/...` | `Huangfu_.../guide_metadata/` | ❌ TODO |

### CSV Updates Needed
Update paths in CSV from:
- `gs://igvf-pertub-seq-pipeline-data/beer_test/huangfu_WTC11_benchmark/yaml_files/`
- → `gs://igvf-pertub-seq-pipeline-data/Huangfu_WTC11-benchmark_TF-Perturb-seq/`

---

## Questions to Resolve

1. **Why different guide metadata?** Uses `v1.tsv` instead of `v3_rc` - is this intentional?
2. **Why different barcode onlist?** Uses `IGVFFI4695IKAL` instead of `IGVFFI9487JPEN` or `IGVFFI1697XAXI`
3. **Should data stay in `beer_test/` or move to production location?**
4. Why is `is_10x3v3 = true` only for this dataset?

---

## Session Log

### 2026-01-23
- Seqspecs uploaded earlier to `gs://igvf-pertub-seq-pipeline-data/Huangfu_WTC11-benchmark_TF-Perturb-seq/seqspecs/`
- Data currently in `beer_test/` GCS location - needs reorganization
- Uses different barcode onlist (`IGVFFI4695IKAL`) and guide metadata (`v1`) from other datasets
- Created this NOTES.md file
