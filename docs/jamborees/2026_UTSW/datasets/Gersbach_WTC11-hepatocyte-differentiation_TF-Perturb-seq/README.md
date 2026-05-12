# Gersbach WTC11 Hepatocyte

**Lab**: Gersbach (Duke)  •  **Cell line**: WTC11  •  **Differentiation**: Hepatocyte (22-day)

Production dataset for the [2026 UTSW jamboree](../../README.md).

## Identity

| | |
|---|---|
| Dataset ID | `Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq` |
| IGVF analysis set | _in progress_ |
| Construct library | `IGVFDS3299AXST` |
| Guide file | [`IGVFFI8270UPKB`](https://data.igvf.org/tabular-files/IGVFFI8270UPKB/) |
| Guide pools | A-D + F |
| Perturbation | CRISPRi |
| Assay | 10x Perturb-seq on NovaSeq X Plus, 25B kit; R1 ~28bp + R2 ~90bp (consistent with 10x 3' v3) |
| Multiplexing | none |
| Measurement sets | 47 (all `in progress` on portal as of 2026-05-07) |

## Status

| Output | Status | Synapse |
|---|---|---|
| CRISPR pipeline | ⚠ Sara's non-canonical layout at `syn70518849` (top-level pipeline outputs flat instead of in `pipeline_outputs/`). Sara to reshape per DE/ESC template. | [`syn70518849`](https://www.synapse.org/Synapse:syn70518849), Sara's full bundle at [`syn74842722`](https://www.synapse.org/Synapse:syn74842722) |
| cNMF | ⏳ pending — Sara to deliver per `schemas/cnmf.json` | — |
| Energy distance | ✓ Sara uploaded into the jamboree (full bundle) | [`syn74902687`](https://www.synapse.org/Synapse:syn74902687) |
| QC (our re-run on Sara's MuData) | ✅ **mirrored 2026-05-12** | [`syn74918946`](https://www.synapse.org/Synapse:syn74918946) (`qc/`) |
| Calibration (our run on Sara's MuData) | ✅ **2026-05-12** — 4 TSVs (all + cis + direct_target + trans) | [`syn74920831`](https://www.synapse.org/Synapse:syn74920831) |

Deeper status: [`crispr_pipeline/README.md`](crispr_pipeline/README.md), [`energy_distance/README.md`](energy_distance/README.md). Open issue: [Gersbach Hep deliverables](../../issues/gersbach-hep-deliverables.md).

## Source data

- **HPC**: `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/` (configs + plots only — no `runs/` dir locally)
- **GCS canonical**: not yet established at the standard `gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/...` path
- **Inference MuData (latest cleanser run)**: Synapse [`syn74728027`](https://www.synapse.org/Synapse:syn74728027) — `gersbach_iPSC_iHep_04.06.2026_cleanser_inference_mudata.h5mu` (30.31 GB)

## Pipeline configuration

(from [`reference/experimental_metadata.tsv`](../../reference/experimental_metadata.tsv))

| Param | Value |
|---|---|
| Canonical run label | _TBD_ |
| Canonical config | `nextflow.config` (lab-specific layout) |
| Hashing | false |
| 10x 3' v3 | false |
| Reverse-complement guides | _?_ |
| Spacer tag | _?_ |
| Guide assignment | `cleanser` |
| Capture method | direct-capture |
| MOI | high |
| Dual guide | false |

## Notes

- **Lab-specific pipeline layout**: Gersbach uses cleanser / direct-capture rather than sceptre / CROP-seq — different inference path than Hon CM and Huangfu. Comparability across datasets has caveats.
- The Synapse parent at `syn70518849` has multiple MuData versions (v1, v2, latest cleanser) and a non-canonical 3-folder layout (has `Perturbo_outputs/`, `cNMF_inputs/`, etc.). We'd rather not derive from this — cleaner to ask Sara for a single canonical run.
- 47 measurement sets, sub-pool naming `10XLane1-8_S{1..47}`.
- Pipeline troubleshooting (Sara). Analysis set is in progress (Gersbach team).
