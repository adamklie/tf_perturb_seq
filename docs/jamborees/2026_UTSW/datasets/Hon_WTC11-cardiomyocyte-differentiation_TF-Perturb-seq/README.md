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

## Two parallel CRISPR pipeline runs

| Run | Source | Canonical for jamboree? |
|---|---|---|
| **`2026_04_19_no_spacer`** | Weizhou's local CRISPR pipeline run, mirrored to Synapse as [`syn74520421`](https://www.synapse.org/Synapse:syn74520421). MuData at [`syn74522725`](https://www.synapse.org/Synapse:syn74522725). | ✅ **YES** — cNMF + ED are run on this. |
| `seqspec_v3` | Our own Apr–May 2026 GCS run at `gs://.../2026_04_15/outs/seqspec_v3/`. | No — kept as comparison / future canonical. |

## Status

| Output | Status | Synapse |
|---|---|---|
| CRISPR pipeline (Weizhou) | ⏳ partial — `pipeline_dashboard/` + `pipeline_outputs/` uploaded by Weizhou, **missing** `pipeline_info/`. | [`syn74520421`](https://www.synapse.org/Synapse:syn74520421) |
| CRISPR pipeline (seqspec_v3) | ✓ complete on GCS; local pull in flight | n/a — not mirrored to jamboree |
| Energy distance (Adam's run on Weizhou MuData) | ✓ uploaded — partial bundle (results + plots + logs) | [`syn74897350`](https://www.synapse.org/Synapse:syn74897350) |
| Energy distance (Sara's run on Weizhou MuData) | ✓ uploaded — full bundle (h5mu + pickles + results + plots) | [`syn74910330`](https://www.synapse.org/Synapse:syn74910330) (`energy_distance_gersbach_comp/`) |
| QC (Weizhou data, our re-run) | ✓ **mirrored 2026-05-12** | [`syn74917453`](https://www.synapse.org/Synapse:syn74917453) (`qc/`) |
| cNMF | ⏳ Stage 1 against seqspec_v3 **failed** (numba JIT working-dir bug); **Alexandra** is running cNMF against Weizhou's MuData in parallel | — |

Open: [Hon CM CRISPR pipeline gap](../../issues/hon-cm-crispr-bundle.md), [Hon CM cNMF Stage 1](https://github.com/adamklie/tf_perturb_seq/issues/20).

## Local layout (after merge of `2026_04_19_no_spacer` + `weizhou_syn74520421/` into one dir)

```
datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/
├── setup/                              # Adam's input-gen code (scripts/configs/samplesheets/seqspec)
├── synapse_inference_mudata/           # Weizhou's MuData (syn74522725) on HPC (gitignored)
├── HonLabInternal/                     # Hon lab's internal cNMF (gitignored)
├── 2026_04_19_no_spacer/               # Weizhou's CRISPR pipeline run (= syn74520421)
│   ├── qc/                             # our local QC re-run on Weizhou MuData
│   └── energy_distance/                # Adam's ED run on Weizhou MuData
└── seqspec_v3/                         # our own May 2026 GCS run
    ├── crispr_pipeline/                # pulled from GCS 2026-05-12 (excl. h5mu)
    ├── qc/                             # QC on seqspec_v3 output
    └── cnmf/                           # Stage 1 setup; Convert failed; awaiting fix or Alexandra's run
```

## Pipeline configuration (from `reference/experimental_metadata.tsv`)

| Param | Value |
|---|---|
| Canonical run label | `2026_04_19_no_spacer` (Weizhou) |
| Hashing | true |
| 10x 3' v3 | false |
| Reverse-complement guides | true |
| Spacer tag | `TAGCTCTTAAAC` |
| Guide assignment | `sceptre` |
| Capture method | CROP-seq |
| MOI | high |
| Dual guide | false |

## Notes

- Pipeline troubleshooting (Weizhou) — seqspec hash modality issues. Production yaml's HTO library_spec lacks a clean tag-region position. `seqspec_v2` run with stripped i7/i5 was abandoned. `seqspec_v3` is our retry that succeeded on GCS.
- The two ED bundles on Synapse (`energy_distance/` = Adam, `energy_distance_gersbach_comp/` = Sara) were both run against Weizhou's MuData (`2026_04_19_no_spacer`) — kept separate for cross-lab comparison.
