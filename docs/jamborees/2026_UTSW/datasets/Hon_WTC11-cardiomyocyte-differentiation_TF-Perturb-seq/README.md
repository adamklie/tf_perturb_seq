# Hon WTC11 Cardiomyocyte

**Lab**: Hon  •  **Cell line**: WTC11  •  **Differentiation**: Cardiomyocyte (12-day)

Production dataset for the [2026 UTSW jamboree](../../README.md).

## Identity

| | |
|---|---|
| Dataset ID | `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq` |
| IGVF analysis set | `IGVFDS6332VCTO` |
| Construct library | `IGVFDS3299AXST` |
| Guide file | [`IGVFFI8270UPKB`](https://data.igvf.org/tabular-files/IGVFFI8270UPKB/) |
| Guide pools | A-D + F |
| Perturbation | CRISPRi |
| Assay | 10x 5'-Perturb (Sigma backbone, HT-like chemistry) |
| Multiplexing | HTO |
| Measurement sets | 28 |

## Status

| Output | Status | Synapse |
|---|---|---|
| CRISPR pipeline | ⚠ partial — has `dashboard/` + `pipeline_outputs/`, no `pipeline_info/`. Bug **Weizhou**. | [`syn74520421`](https://www.synapse.org/Synapse:syn74520421) |
| cNMF | ⏳ run pending (gated on full CRISPR bundle) | — |
| Energy distance | ⏳ run pending (gated on full CRISPR bundle) | — |

Deeper status: [`crispr_pipeline/README.md`](crispr_pipeline/README.md). Open issue: [Hon CM CRISPR pipeline gap](../../issues/hon-cm-crispr-bundle.md).

## Source data

- **HPC**: `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/`
- **GCS canonical**: `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/2026_04_15/outs/initial_run` (incomplete — pipeline didn't reach dashboard / outputs stages)
- **Inference MuData (cleaned)**: Synapse [`syn74522725`](https://www.synapse.org/Synapse:syn74522725) — Hon-team-supplied input for downstream cNMF + energy distance

## Pipeline configuration

(from [`reference/experimental_metadata.tsv`](../../reference/experimental_metadata.tsv))

| Param | Value |
|---|---|
| Canonical run label | `initial_run` |
| Canonical config | `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq_initial_run.config` |
| Hashing | true |
| 10x 3' v3 | false |
| Reverse-complement guides | true |
| Spacer tag | `TAGCTCTTAAAC` |
| Guide assignment | `sceptre` |
| Capture method | CROP-seq |
| MOI | high |
| Dual guide | false |

## Notes

- Pipeline troubleshooting (Weizhou) — seqspec hash modality issues. Production yaml's HTO library_spec lacks a clean tag-region position. `seqspec_v2` run with stripped i7/i5 was abandoned.
- Synapse mirror at `syn74520421` is named `2026_04_19_no_spacer` — different from the metadata's `initial_run` canonical label. Reconcile when full bundle arrives.
