# Hon_WTC11-benchmark_TF-Perturb-seq

## Overview

| Attribute | Value |
|-----------|-------|
| **Analysis Set** | [IGVFDS4761PYUO](https://data.igvf.org/analysis-sets/IGVFDS4761PYUO/) |
| **Lab** | Hon |
| **Contact** | Mpathi |
| **System** | iPSC (WTC11) |
| **Technology** | 10x 5' (w/HTO) |
| **MOI** | High |
| **Cell Hashing** | Yes |
| **Modalities** | scRNA, gRNA, hash |
| **Cells Sequenced** | ~750K |
| **Pipeline Status** | ✅ Run completed with results |

---

## Pipeline Inputs Status

| File | IGVF Accession | Local | GCS | Status |
|------|----------------|-------|-----|--------|
| Seqspecs (3) | - | ✅ | ✅ | Done |
| Barcode onlist | `IGVFFI9487JPEN` | ✅ | ✅ | Done |
| Guide metadata | `IGVFFI5765HMZH` | ✅ | ✅ | Done |
| Hashtag map | `IGVFFI3028BJFI` | ✅ | ✅ | Done |

---

## Key Files

### Sample Metadata
- **Production CSV**: `sample_metadata_gcp_2026_01_18.csv`
- **Patched CSV**: `sample_metadata_gcp_2026_01_18_patched.csv`
- **Synapse**: [syn70753575](https://www.synapse.org/Synapse:syn70753575)

### Seqspecs (on GCS)
- `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-benchmark_TF-Perturb-seq/2026_01_18/rna_seqspec.yml`
- `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-benchmark_TF-Perturb-seq/2026_01_18/guide_seqspec.yml`
- `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-benchmark_TF-Perturb-seq/2026_01_18/hash_seqspec.yml`

### Config
- **Production**: `harmless_teacher.config`
- **Key settings**:
  - `GUIDE_ASSIGNMENT_method = 'sceptre'`
  - `Multiplicity_of_infection = 'high'`
  - `ENABLE_DATA_HASHING = true`
  - `reverse_complement_guides = true`
  - `spacer_tag = "TTAGCTCTTAAAC"`

---

## Pipeline Results

Results in `results/1_CRISPR_pipeline/`:
- **2026_01_21/**: Dashboard, inference results
- **2026_01_22/**: Additional run

**Key outputs**:
- `pipeline_dashboard/dashboard.html`
- `inference_mudata.h5mu` ([syn70753405](https://www.synapse.org/Synapse:syn70753405))
- `perturbo_cis_per_guide_output.tsv.gz`
- `perturbo_trans_per_guide_output.tsv.gz`

**QC Analysis**: [Slides](https://docs.google.com/presentation/d/1E5aUUtMdyWWutOi9Uj24nlXjYkLPOfgm7N2XbEqjrtE/edit?slide=id.g39e356c57dd_0_11)

---

## Gene Program Discovery

[syn70776514](https://www.synapse.org/Synapse:syn70776514)

---

## TODO

- [ ] Update guide metadata to latest spec
- [ ] Rerun through pipeline
- [ ] Get CellRanger and PySpade output from Mpathi

---

## Notes

- This dataset serves as the **reference model** for organization (3 seqspecs per modality)
- Uses same guide metadata (`IGVFFI5765HMZH`) as Gersbach datasets
- Uses same barcode onlist as Gersbach HTv2 (`IGVFFI9487JPEN`)
