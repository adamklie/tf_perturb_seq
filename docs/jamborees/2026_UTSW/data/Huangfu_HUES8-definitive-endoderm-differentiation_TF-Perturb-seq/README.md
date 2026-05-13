# Huangfu HUES8 Definitive Endoderm

**Lab**: Huangfu  •  **Cell line**: HUES8  •  **Differentiation**: Definitive endoderm

Production dataset for the [2026 UTSW jamboree](../../README.md).

## Identity

| | |
|---|---|
| Dataset ID | `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq` |
| IGVF analysis set | `IGVFDS9951KTRR` |
| Construct library | `IGVFDS3299AXST` |
| Guide file | [`IGVFFI8270UPKB`](https://data.igvf.org/tabular-files/IGVFFI8270UPKB/) |
| Guide pools | A-D |
| Perturbation | CRISPRi |
| Assay | 10x 3' v3 |
| Multiplexing | none |
| Measurement sets | 8 |

## Status

| Output | Status | Synapse |
|---|---|---|
| CRISPR pipeline | ✅ canonical 3-folder bundle | [`syn74834952`](https://www.synapse.org/Synapse:syn74834952) |
| cNMF | ✅ complete (k=200, dt=2.0; Config + Inference + Evaluation across 8 K values + selected-k Plot/Interpretation + adata) | [`syn74893844`](https://www.synapse.org/Synapse:syn74893844) |
| Energy distance (Adam) | ✅ complete (⚠ p-value calibration concern) | [`syn74883327`](https://www.synapse.org/Synapse:syn74883327) |
| Energy distance (Sara, comp) | ✅ complete | [`syn74910358`](https://www.synapse.org/Synapse:syn74910358) |
| QC | ✅ **mirrored 2026-05-12** | [`syn74918479`](https://www.synapse.org/Synapse:syn74918479) |
| Calibration | ✅ **2026-05-12** — 4 TSVs (all + cis + direct_target + trans); FDR<0.05 = 14,521 trans tests | [`syn74920615`](https://www.synapse.org/Synapse:syn74920615) |

E-distance p-value calibration concern: [Issue 1](https://github.com/adamklie/tf_perturb_seq/issues/edistance-calibration.md). Deeper status: [`crispr_pipeline/README.md`](crispr_pipeline/README.md), [`cnmf/README.md`](cnmf/README.md), [`energy_distance/README.md`](energy_distance/README.md).

## Source data

- **HPC**: `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/`
- **GCS canonical**: `gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/2026_04_09/outs/muddy_penguin`
- **Inference MuData**: in `pipeline_dashboard/inference_mudata.h5mu` under [`syn74834952`](https://www.synapse.org/Synapse:syn74834952)

## Pipeline configuration

(from [`reference/experimental_metadata.tsv`](../../reference/experimental_metadata.tsv))

| Param | Value |
|---|---|
| Canonical run label | `muddy_penguin` |
| Canonical config | `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq_muddy_penguin.config` |
| Hashing | false |
| 10x 3' v3 | true |
| Reverse-complement guides | false |
| Spacer tag | `GAGTACATGGGGG` (13 bp) |
| Guide assignment | `sceptre` |
| Capture method | CROP-seq |
| MOI | high |
| Dual guide | false |

## Notes

- Two runs were compared on `spacer_tag` length: `entertaining_hamster` (12 bp `GAGTACATGGGG`) vs `muddy_penguin` (13 bp `GAGTACATGGGGG`). `muddy_penguin` is canonical. Production data may carry an extra leading G vs the benchmark.
- Energy distance run completed cleanly but p-values are anti-conservative (all 100 NCs `pval_mean = 0`). Use `distance_mean` as effect-size proxy until the calibration fix lands.
- Shares construct library set + guide file with the embryonic-stem-cell dataset.
