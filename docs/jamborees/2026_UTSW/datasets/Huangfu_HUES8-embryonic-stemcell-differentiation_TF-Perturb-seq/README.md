# Huangfu HUES8 Embryonic Stem Cell

**Lab**: Huangfu  •  **Cell line**: HUES8  •  **Differentiation**: Embryonic stem cell (undifferentiated)

Production dataset for the [2026 UTSW jamboree](../../README.md).

## Identity

| | |
|---|---|
| Dataset ID | `Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq` |
| IGVF analysis set | `IGVFDS1216AEWT` |
| Construct library | `IGVFDS3299AXST` |
| Guide file | [`IGVFFI8270UPKB`](https://data.igvf.org/tabular-files/IGVFFI8270UPKB/) |
| Guide pools | A-D |
| Perturbation | CRISPRi |
| Assay | 10x 3' v3 |
| Multiplexing | none |
| Measurement sets | 8 |

## Status — all 4 analysis tiers on jamboree Synapse ✅

| Output | Status | Synapse |
|---|---|---|
| CRISPR pipeline | ✅ canonical 3-folder bundle | [`syn74835010`](https://www.synapse.org/Synapse:syn74835010) |
| cNMF | ✅ `042926_huangfu_esc_torchcnmf_KskillA` (selected-k=200, dt=2.0) | [`syn74893846`](https://www.synapse.org/Synapse:syn74893846) |
| Energy distance (Adam) | ✅ complete (⚠ p-value calibration concern) | [`syn74883475`](https://www.synapse.org/Synapse:syn74883475) |
| Energy distance (Sara, comp) | ✅ complete | [`syn74910472`](https://www.synapse.org/Synapse:syn74910472) (`energy_distance_gersbach_comp/`) |
| QC | ✅ mirrored 2026-05-12 | [`syn74918600`](https://www.synapse.org/Synapse:syn74918600) (`qc/`) |

Calibration concern: [Issue 1](../../issues/edistance-calibration.md). Deeper status: [`crispr_pipeline/README.md`](crispr_pipeline/README.md), [`cnmf/README.md`](cnmf/README.md), [`energy_distance/README.md`](energy_distance/README.md).

## Source data

- **HPC**: `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/`
- **GCS canonical**: `gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/2026_04_13/outs/sceptre_v1`
- **Inference MuData**: in `pipeline_dashboard/inference_mudata.h5mu` under [`syn74835010`](https://www.synapse.org/Synapse:syn74835010)

## Pipeline configuration

(from [`reference/experimental_metadata.tsv`](../../reference/experimental_metadata.tsv))

| Param | Value |
|---|---|
| Canonical run label | `sceptre_v1` |
| Canonical config | `Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq_2026_04_13.config` |
| Hashing | false |
| 10x 3' v3 | true |
| Reverse-complement guides | false |
| Spacer tag | `GAGTACATGGGGG` |
| Guide assignment | `sceptre` |
| Capture method | CROP-seq |
| MOI | high |
| Dual guide | false |

## Notes

- Companion to the definitive-endoderm sibling — shares construct library + guide file. The two differ only in differentiation state.
- ESC energy distances are ~5.7× larger than DE (median 527 vs 93). Undifferentiated stem cells appear to have stronger TF-perturbation effects (or higher baseline cell-state variance — see calibration concern in [Issue 1](../../issues/edistance-calibration.md)).
- Pipeline status (per metadata): "Setting up pipeline scripts; portal had blockers (analysis set + seqspecs) audited 2026-03-25."
