# Hon WTC11 Cardiomyocyte — CRISPR pipeline

Two runs exist; Weizhou's is canonical for the jamboree (cNMF + ED were run on its MuData).

## Weizhou's run (canonical) — `2026_04_19_no_spacer`

| | |
|---|---|
| Synapse | [`syn74520421`](https://www.synapse.org/Synapse:syn74520421) |
| MuData | [`syn74522725`](https://www.synapse.org/Synapse:syn74522725) |
| Status | ⚠ Partial bundle — `dashboard/` + `pipeline_outputs/` uploaded, **missing** `pipeline_info/`. |

```
syn74520421/
├── dashboard/         (syn74526301)   ← named "dashboard", not "pipeline_dashboard"
└── pipeline_outputs/  (syn74520424)   ← perturbo trans-results live here (calibration input)
```

This is the canonical source for downstream analyses on Hon CM: the perturbo trans-results in `pipeline_outputs/` are the input to calibration; the MuData [`syn74522725`](https://www.synapse.org/Synapse:syn74522725) is the input to cNMF (Alexandra) and energy distance (done — Adam + Sara comparison runs).

Outstanding: see [Issue: Hon CM CRISPR pipeline bundle gap](https://github.com/adamklie/tf_perturb_seq/issues/hon-cm-crispr-bundle.md).

## Our `seqspec_v3` run (not canonical)

GCS: `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/2026_04_15/outs/seqspec_v3/`.
Pulled locally (excl. h5mu) to `seqspec_v3/crispr_pipeline/` on 2026-05-12. Not mirrored to the jamboree Synapse folder — kept as a comparison run.

## Schema + walkthrough

- Machine-readable schema: [`schemas/crispr_pipeline.json`](../../../schemas/crispr_pipeline.json)
- Analysis-level walkthrough: [`docs/analysis/CRISPR_PIPELINE_OUTPUTS.md`](../../../../../analysis/CRISPR_PIPELINE_OUTPUTS.md)
