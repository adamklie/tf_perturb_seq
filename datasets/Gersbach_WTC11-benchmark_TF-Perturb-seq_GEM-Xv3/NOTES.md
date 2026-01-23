# Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3 - Dataset Notes

## Overview
- **Analysis Set Accession**: IGVFDS6673ZFFG
- **Sample Metadata CSV**: `downloads/IGVFDS6673ZFFG.updated.v4.csv`
- **Modalities**: scRNA, gRNA (no cell hashing)

---

## Data Source Files

| File | Description |
|------|-------------|
| `downloads/IGVFDS6673ZFFG.tsv` | Original portal export with IGVF accessions |
| `downloads/IGVFDS6673ZFFG.updated.v4.csv` | Modified CSV with local Duke cluster paths |

---

## Key Files Identified from CSV

### 1. Barcode Onlist
- **IGVF Accession**: `IGVFFI1697XAXI`
- **Current path in CSV**: `/work/aeb84/igvf_crispr/tf_pilot/data/IGVFFI1697XAXI.tsv`
- **Local copy**: `downloads/IGVFFI1697XAXI.tsv.gz` (needs unzipping)
- **Portal TSV note**: `barcode_onlist` column is EMPTY in portal TSV - was added manually to updated CSV
- **Status**: Have local copy (gzipped)

### 2. Seqspecs
- **Count**: 50 unique seqspec files (one per sample/lane)
- **Pattern**: `IGVFFI####XXXX.formatted.fixed.yaml`
- **Current path in CSV**: `/work/aeb84/igvf_crispr/tf_pilot/data/IGVFFI*.formatted.fixed.yaml`
- **Local copies**: NOT FOUND in downloads/
- **Portal availability**: YES - accessions listed in `downloads/IGVFDS6673ZFFG.tsv`
- **Status**: Available on IGVF portal, waiting to confirm usability in pipeline

**Unique seqspec accessions:**
```
IGVFFI0123VGVY, IGVFFI0268RYST, IGVFFI0274PUWA, IGVFFI0499BDIV, IGVFFI0836FIUK,
IGVFFI1335NPXV, IGVFFI1362OCGY, IGVFFI1701AUAR, IGVFFI1806RVDJ, IGVFFI2498HIKH,
IGVFFI2572SITG, IGVFFI2726WCUU, IGVFFI3001OCFQ, IGVFFI3074EJUE, IGVFFI3074VZJW,
IGVFFI3077AIXI, IGVFFI3797VTQZ, IGVFFI4130EZOR, IGVFFI4573QVNI, IGVFFI4678MXPS,
IGVFFI4719TXUG, IGVFFI4987DUPG, IGVFFI5714WXET, IGVFFI5811SVMQ, IGVFFI5820TKFN,
IGVFFI6179UQQW, IGVFFI6467WMOP, IGVFFI6820QRAF, IGVFFI6983XNSI, IGVFFI7139JADC,
IGVFFI7172DLGO, IGVFFI7356TFEN, IGVFFI7379IGPW, IGVFFI7547UAZY, IGVFFI7683SIHY,
IGVFFI7752HNLN, IGVFFI7883PGFF, IGVFFI7933QCZK, IGVFFI7947YZLX, IGVFFI7957XWGJ,
IGVFFI8067TGTN, IGVFFI8449BLVI, IGVFFI8500HGMG, IGVFFI8609HZFG, IGVFFI8791HIOU,
IGVFFI8893TIGY, IGVFFI9161ALQS, IGVFFI9509PDFO, IGVFFI9685JXEU, IGVFFI9818AEJQ
```

### 3. Guide Metadata
- **IGVF Accession**: `IGVFFI5765HMZH` (from portal TSV `guide_design` column for gRNA rows)
- **Local file used**: `benchmark_guide_metadata_v3_rc.tsv`
- **Current path in CSV**: `/work/aeb84/igvf_crispr/tf_pilot/data/benchmark_guide_metadata_v3_rc.tsv`
- **Local copy**: NOT FOUND in downloads/
- **Portal availability**: YES - `IGVFFI5765HMZH`
- **Note**: `_rc` likely means "reverse complement" - relates to config param `reverse_complement_guides = true`
- **Status**: Available on IGVF portal, waiting to confirm usability in pipeline

---

## Config Files

### 1. `tf_pilot.nextflow.cleanser.config`
Main Nextflow config used on Duke SLURM cluster.

**Key parameters:**
- `ENABLE_DATA_HASHING = false`
- `DUAL_GUIDE = false`
- `Multiplicity_of_infection = 'low'`
- `GUIDE_ASSIGNMENT_method = 'cleanser'`
- `GUIDE_ASSIGNMENT_capture_method = 'crop-seq'`
- `REFERENCE_gencode_gtf = '/hpc/group/gersbachlab/Reference_Data/IGVF/GRCh38/Gencode/v43/IGVFFI7217ZMJZ.gtf'`

### 2. `params.config`
Alternative config with different settings.

**Key differences from tf_pilot.nextflow.cleanser.config:**
- `ENABLE_DATA_HASHING = true` (vs false)
- `use_igvf_reference = true`
- `Multiplicity_of_infection = 'high'` (vs low)
- `reverse_complement_guides = true`
- `spacer_tag = "TAGCTCTTAAAC"` (sequence in front of guide)

---

## Action Items

### Files Status
| File | Local | Portal | Action Needed |
|------|-------|--------|---------------|
| Barcode onlist (`IGVFFI1697XAXI`) | ✅ `.gz` | - | Unzip |
| 50 seqspec YAMLs | ❌ | ✅ | Download & confirm usability |
| Guide metadata (`IGVFFI5765HMZH`) | ❌ | ✅ | Download & confirm usability |

### Files to Upload to GCS (after confirmation)
Target: `gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/`

- [ ] Seqspecs → `seqspecs/`
- [ ] Barcode onlist → `onlists/`
- [ ] Guide metadata → `guide_metadata/`

### CSV Updates Needed
After uploading, update paths in CSV from:
- `/work/aeb84/igvf_crispr/tf_pilot/data/` → `gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/`

---

## Questions to Resolve
1. ~~Where are the 50 seqspec YAML files?~~ → Available on IGVF portal
2. ~~Where is the guide metadata file?~~ → Available as `IGVFFI5765HMZH` on portal
3. **Are the portal seqspecs usable as-is, or do they need the `.formatted.fixed` processing?**
4. Which config file should be used for GCP runs - need to reconcile hashing and MOI settings?

---

## Session Log

### 2026-01-23
- Initial exploration of dataset structure
- Identified key files from CSV: barcode onlist, seqspecs, guide metadata
- Reviewed two config files and documented differences
- Discovered portal TSV (`IGVFDS6673ZFFG.tsv`) has IGVF accessions for seqspecs and guide metadata
- Guide metadata accession: `IGVFFI5765HMZH`
- Seqspecs and guide metadata available on portal - waiting to confirm pipeline usability
- Created this NOTES.md file
