# Huangfu_WTC11-benchmark_TF-Perturb-seq

## Overview

| Attribute | Value |
|-----------|-------|
| **Analysis Set** | (check portal) |
| **Lab** | Huangfu |
| **System** | iPSC (WTC11) |
| **Technology** | 10x 3' v3 |
| **MOI** | High |
| **Cell Hashing** | No |
| **Modalities** | scRNA, gRNA |
| **Pipeline Status** | Partial - data in `beer_test/` location |

---

## Pipeline Inputs Status

| File | Accession/Name | Local | GCS | Action Needed |
|------|----------------|-------|-----|---------------|
| Seqspecs (8) | Named files | ✅ | ✅ New location | Update CSV paths |
| Barcode onlist | `IGVFFI4695IKAL` | ❌ | ✅ `beer_test/` | Move to new location |
| Guide metadata | `benchmark_guide_metadata_v1.tsv` | ❌ | ✅ `beer_test/` | Move to new location |

**WARNING**: Uses different barcode onlist and guide metadata than other datasets!

---

## Key Files

### Sample Metadata
- **Working CSV**: `downloads/2025_11_26_huangfu_cloud_sample_metadata.csv`
- **Current GCS prefix**: `gs://igvf-pertub-seq-pipeline-data/beer_test/huangfu_WTC11_benchmark/`

### Seqspecs (8 files - per replicate)
- **New location**: `gs://igvf-pertub-seq-pipeline-data/Huangfu_WTC11-benchmark_TF-Perturb-seq/seqspecs/` ✅
- **Local**: `downloads/*.yaml`
- **Files**:
  - `michael-beer_Perturb-benchmark-WTC11-rep{1-4}.yaml` (scRNA)
  - `michael-beer_Perturb-benchmark-WTC11-guide-sequencing-rep{1-4}.yaml` (gRNA)

### Barcode Onlist
- **File**: `IGVFFI4695IKAL.tsv` (**different from other datasets!**)
- **Current**: `gs://igvf-pertub-seq-pipeline-data/beer_test/huangfu_WTC11_benchmark/yaml_files/`

### Guide Metadata
- **File**: `benchmark_guide_metadata_v1.tsv` (**v1, not v3_rc like others!**)
- **Current**: `gs://igvf-pertub-seq-pipeline-data/beer_test/huangfu_WTC11_benchmark/yaml_files/`

### Config
- **Primary**: `bin/nextflow.config`
- **Key settings**:
  - `is_10x3v3 = true` (**unique to this dataset!**)
  - `GUIDE_ASSIGNMENT_method = 'sceptre'`
  - `Multiplicity_of_infection = 'high'`
  - `ENABLE_DATA_HASHING = false`

---

## Differences from Other Datasets

| Attribute | Huangfu | Hon | Gersbach |
|-----------|---------|-----|----------|
| Barcode Onlist | `IGVFFI4695IKAL` | `IGVFFI9487JPEN` | `IGVFFI1697XAXI`/`IGVFFI9487JPEN` |
| Guide Metadata | `v1.tsv` | `IGVFFI5765HMZH` (v3_rc) | `IGVFFI5765HMZH` (v3_rc) |
| 10x 3' v3 flag | ✅ Yes | ❌ No | ❌ No |
| Seqspec approach | 8 (per rep) | 3 (per modality) | 50-68 (per sample) |

---

## TODO

- [ ] Move barcode onlist from `beer_test/` to `Huangfu_.../onlists/`
- [ ] Move guide metadata from `beer_test/` to `Huangfu_.../guide_metadata/`
- [ ] Update CSV paths from `beer_test/` to production location
- [ ] Clarify why different guide metadata version (v1 vs v3_rc)
- [ ] Run pipeline on GCP

---

## Questions

1. **Why different guide metadata?** Uses `v1.tsv` instead of `v3_rc` - intentional?
2. **Why different barcode onlist?** Uses `IGVFFI4695IKAL` - is this correct?
3. **Why `is_10x3v3 = true`?** Only dataset with this flag
4. Should data stay in `beer_test/` or move to production location?
