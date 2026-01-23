# Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3

## Overview

| Attribute | Value |
|-----------|-------|
| **Analysis Set** | [IGVFDS6673ZFFG](https://data.igvf.org/analysis-sets/IGVFDS6673ZFFG/) |
| **Lab** | Gersbach (Duke) |
| **System** | iPSC (WTC11) |
| **Technology** | 10x GEM-X 3' v3 |
| **MOI** | Low |
| **Cell Hashing** | No |
| **Modalities** | scRNA, gRNA |
| **Pipeline Status** | Not yet run on GCP |

---

## Pipeline Inputs Status

| File | IGVF Accession | Local | GCS | Action Needed |
|------|----------------|-------|-----|---------------|
| Barcode onlist | `IGVFFI1697XAXI` | ✅ `.gz` | ❌ | Upload to GCS |
| Seqspecs (50) | Various | ❌ | ❌ | Download from portal, confirm usability |
| Guide metadata | `IGVFFI5765HMZH` | ❌ | ❌ | Download from portal |

---

## Key Files

### Sample Metadata
- **Portal export**: `downloads/IGVFDS6673ZFFG.tsv`
- **Working CSV**: `downloads/IGVFDS6673ZFFG.updated.v4.csv` (Duke cluster paths)

### Seqspecs
- **Count**: 50 unique files (one per sample/lane)
- **Pattern**: `IGVFFI####XXXX.formatted.fixed.yaml`
- **Portal availability**: Yes - accessions in portal TSV
- **Note**: Need to confirm if portal seqspecs are usable as-is or need processing

### Barcode Onlist
- **Accession**: `IGVFFI1697XAXI`
- **Local**: `downloads/IGVFFI1697XAXI.tsv.gz` (needs unzipping)

### Guide Metadata
- **Accession**: `IGVFFI5765HMZH`
- **Local file used**: `benchmark_guide_metadata_v3_rc.tsv`
- **Shared with**: HTv2, Hon

### Config
- **Primary**: `downloads/tf_pilot.nextflow.cleanser.config`
- **Alternative**: `downloads/params.config`
- **Key settings**:
  - `GUIDE_ASSIGNMENT_method = 'cleanser'`
  - `Multiplicity_of_infection = 'low'`
  - `ENABLE_DATA_HASHING = false`

---

## TODO

- [ ] Download seqspecs from IGVF portal and confirm pipeline usability
- [ ] Download guide metadata (`IGVFFI5765HMZH`) from portal
- [ ] Unzip barcode onlist (`IGVFFI1697XAXI.tsv.gz`)
- [ ] Upload all files to `gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/`
- [ ] Update CSV paths from Duke cluster to GCS
- [ ] Run pipeline on GCP

---

## Questions

1. Are portal seqspecs usable as-is, or do they need `.formatted.fixed` processing?
2. Should we use CLEANSER (current) or SCEPTRE for guide assignment?
3. Reconcile config differences between `tf_pilot.nextflow.cleanser.config` and `params.config`
