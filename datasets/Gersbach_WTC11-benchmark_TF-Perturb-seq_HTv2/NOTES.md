# Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2 - Dataset Notes

## Overview
- **Analysis Set Accession**: IGVFDS6237URFJ
- **Sample Metadata CSV**: `downloads/IGVFDS6237URFJ.updated.v3.csv`
- **Modalities**: scRNA, gRNA (no cell hashing)
- **Related dataset**: Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3 (same experiment, different chemistry)

---

## Data Source Files

| File | Description |
|------|-------------|
| `downloads/IGVFDS6237URFJ.tsv` | Original portal export with IGVF accessions |
| `downloads/IGVFDS6237URFJ.updated.v3.csv` | Modified CSV with local Duke cluster paths |

---

## Key Files Identified from CSV

### 1. Barcode Onlist
- **IGVF Accession**: `IGVFFI9487JPEN`
- **Current path in CSV**: `/work/aeb84/igvf_crispr/tf_pilot/data/IGVFFI9487JPEN.tsv`
- **Local copy**: `downloads/IGVFFI9487JPEN.tsv.gz` (needs unzipping)
- **Portal TSV note**: `barcode_onlist` column is EMPTY in portal TSV - was added manually to updated CSV
- **Status**: Have local copy (gzipped)

### 2. Seqspecs
- **Count**: 68 unique seqspec files (one per sample/lane)
- **Pattern**: `IGVFFI####XXXX.simple_formatted.yaml` (different from GEM-Xv3's `.formatted.fixed.yaml`)
- **Current path in CSV**: `/work/aeb84/igvf_crispr/tf_pilot/data/IGVFFI*.simple_formatted.yaml`
- **Local copies**: NOT FOUND in downloads/
- **Portal availability**: YES - accessions listed in `downloads/IGVFDS6237URFJ.tsv`
- **Status**: ⏳ Available on IGVF portal, waiting to confirm usability in pipeline

### 3. Guide Metadata
- **IGVF Accession**: `IGVFFI5765HMZH` (same as GEM-Xv3)
- **Local file used**: `benchmark_guide_metadata_v3_rc.tsv`
- **Current path in CSV**: `/work/aeb84/igvf_crispr/tf_pilot/data/benchmark_guide_metadata_v3_rc.tsv`
- **Local copy**: NOT FOUND in downloads/
- **Portal availability**: YES - `IGVFFI5765HMZH`
- **Note**: `_rc` likely means "reverse complement" - shared with GEM-Xv3
- **Status**: ⏳ Available on IGVF portal, waiting to confirm usability in pipeline

---

## Config Files

### `bin/tf_pilot.nextflow.cleanser.config`
Same config as GEM-Xv3 - used on Duke SLURM cluster.

**Key parameters:**
- `ENABLE_DATA_HASHING = false`
- `DUAL_GUIDE = false`
- `Multiplicity_of_infection = 'low'`
- `GUIDE_ASSIGNMENT_method = 'cleanser'`
- `GUIDE_ASSIGNMENT_capture_method = 'crop-seq'`
- `REFERENCE_gencode_gtf = '/hpc/group/gersbachlab/Reference_Data/IGVF/GRCh38/Gencode/v43/IGVFFI7217ZMJZ.gtf'`

---

## Comparison with GEM-Xv3

| Attribute | HTv2 | GEM-Xv3 |
|-----------|------|---------|
| Analysis Set | IGVFDS6237URFJ | IGVFDS6673ZFFG |
| Barcode Onlist | `IGVFFI9487JPEN` | `IGVFFI1697XAXI` |
| Seqspec count | 68 | 50 |
| Seqspec pattern | `.simple_formatted.yaml` | `.formatted.fixed.yaml` |
| Guide metadata | `IGVFFI5765HMZH` (same) | `IGVFFI5765HMZH` (same) |
| Config | Identical | Identical |

---

## Action Items

### Files Status
| File | Local | Portal | Action Needed |
|------|-------|--------|---------------|
| Barcode onlist (`IGVFFI9487JPEN`) | ✅ `.gz` | - | Unzip |
| 68 seqspec YAMLs | ❌ | ✅ | Download & confirm usability |
| Guide metadata (`IGVFFI5765HMZH`) | ❌ | ✅ | Download & confirm usability (shared with GEM-Xv3) |

### Files to Upload to GCS (after confirmation)
Target: `gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/`

- [ ] Seqspecs → `seqspecs/`
- [ ] Barcode onlist → `onlists/`
- [ ] Guide metadata → `guide_metadata/`

### CSV Updates Needed
After uploading, update paths in CSV from:
- `/work/aeb84/igvf_crispr/tf_pilot/data/` → `gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/`

---

## Questions to Resolve
1. **Why different seqspec processing?** HTv2 uses `.simple_formatted.yaml` vs GEM-Xv3's `.formatted.fixed.yaml`
2. Are the portal seqspecs usable as-is, or do they need processing?
3. Can guide metadata be shared between HTv2 and GEM-Xv3 (same accession)?

---

## Session Log

### 2026-01-23
- Initial exploration of dataset structure
- Identified key files from CSV: barcode onlist, seqspecs, guide metadata
- Compared with GEM-Xv3 - same experiment, different chemistry
- Guide metadata is shared (`IGVFFI5765HMZH`)
- Seqspecs and guide metadata available on portal - waiting to confirm pipeline usability
- Created this NOTES.md file
