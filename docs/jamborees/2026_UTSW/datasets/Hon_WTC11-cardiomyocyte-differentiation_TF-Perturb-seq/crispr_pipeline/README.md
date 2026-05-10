# Hon WTC11 Cardiomyocyte — CRISPR pipeline

**Status**: ⚠ Partial canonical bundle on Synapse — has `dashboard/` + `pipeline_outputs/`, missing `pipeline_info/`. Awaiting the rest from **Weizhou** (Hon team).

## Synapse

[`syn74520421`](https://www.synapse.org/Synapse:syn74520421) — folder named `2026_04_19_no_spacer`.

```
syn74520421/
├── dashboard/         (syn74526301)   ← named "dashboard", not "pipeline_dashboard"
└── pipeline_outputs/  (syn74520424)
```

`pipeline_info/` (Nextflow params + software versions) is not yet present.

## What's outstanding

See [Issue: Hon CM CRISPR pipeline bundle gap](../../../issues/hon-cm-crispr-bundle.md) for the full ask + reference example.

## Schema + walkthrough

- Machine-readable schema: [`schemas/crispr_pipeline.json`](../../../schemas/crispr_pipeline.json)
- Analysis-level walkthrough: [`docs/analysis/CRISPR_PIPELINE_OUTPUTS.md`](../../../../../analysis/CRISPR_PIPELINE_OUTPUTS.md)
