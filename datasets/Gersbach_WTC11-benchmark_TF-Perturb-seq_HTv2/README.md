# Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2

## Overview

| Attribute | Value |
|-----------|-------|
| **Analysis Set** | [IGVFDS6237URFJ](https://data.igvf.org/analysis-sets/IGVFDS6237URFJ/) |
| **Lab** | Gersbach (Duke) |
| **System** | iPSC (WTC11) |
| **Technology** | 10x HT v2 |
| **MOI** | Low |
| **Cell Hashing** | No |
| **Modalities** | scRNA, gRNA |
| **Pipeline Status** | Not yet run on GCP |

**Note**: Same experiment as GEM-Xv3, different chemistry.

---

## Pipeline Inputs Status

| File | IGVF Accession | Local | GCS | Action Needed |
|------|----------------|-------|-----|---------------|
| Barcode onlist | `IGVFFI9487JPEN` | ✅ `.gz` | ❌ | Upload to GCS |
| Seqspecs (68) | Various | ❌ | ❌ | Download from portal, confirm usability |
| Guide metadata | `IGVFFI5765HMZH` | ❌ | ❌ | Download from portal (shared with GEM-Xv3) |

---

## Key Files

### Sample Metadata
- **Portal export**: `downloads/IGVFDS6237URFJ.tsv`
- **Working CSV**: `downloads/IGVFDS6237URFJ.updated.v3.csv` (Duke cluster paths)

### Seqspecs
- **Count**: 68 unique files (one per sample/lane)
- **Pattern**: `IGVFFI####XXXX.simple_formatted.yaml`
- **Note**: Different pattern from GEM-Xv3's `.formatted.fixed.yaml`

### Barcode Onlist
- **Accession**: `IGVFFI9487JPEN`
- **Local**: `downloads/IGVFFI9487JPEN.tsv.gz` (needs unzipping)

### Guide Metadata
- **Accession**: `IGVFFI5765HMZH` (same as GEM-Xv3)
- **Local file used**: `benchmark_guide_metadata_v3_rc.tsv`

### Config
- **Primary**: `bin/tf_pilot.nextflow.cleanser.config`
- **Key settings**: Same as GEM-Xv3

---

## Comparison with GEM-Xv3

| Attribute | HTv2 | GEM-Xv3 |
|-----------|------|---------|
| Analysis Set | IGVFDS6237URFJ | IGVFDS6673ZFFG |
| Seqspec count | 68 | 50 |
| Seqspec pattern | `.simple_formatted.yaml` | `.formatted.fixed.yaml` |
| Barcode onlist | `IGVFFI9487JPEN` | `IGVFFI1697XAXI` |
| Guide metadata | Same (`IGVFFI5765HMZH`) | Same |
| Config | Identical | Identical |

---

## TODO

- [ ] Download seqspecs from IGVF portal and confirm pipeline usability
- [ ] Download guide metadata (`IGVFFI5765HMZH`) from portal
- [ ] Unzip barcode onlist (`IGVFFI9487JPEN.tsv.gz`)
- [ ] Upload all files to `gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/`
- [ ] Update CSV paths from Duke cluster to GCS
- [ ] Run pipeline on GCP

---

## Questions

1. Why different seqspec processing (`.simple_formatted` vs `.formatted.fixed`)?
2. Can guide metadata be shared with GEM-Xv3 in a common location?
