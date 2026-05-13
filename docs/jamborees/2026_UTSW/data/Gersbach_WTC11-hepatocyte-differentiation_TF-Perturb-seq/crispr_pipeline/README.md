# Gersbach WTC11 Hepatocyte — CRISPR pipeline

**Status**: ⚠ Non-canonical layout exists on Synapse. Awaiting **Sara** (Gersbach team) to deliver a canonical 3-folder bundle.

## Synapse (existing, non-canonical)

[`syn70518849`](https://www.synapse.org/Synapse:syn70518849) — has `Perturbo_outputs/`, `cNMF_inputs/`, multiple MuData files in non-canonical positions. Not in the 3-folder layout we agreed on.

## What's outstanding

A canonical `pipeline_dashboard/` + `pipeline_info/` + `pipeline_outputs/` bundle from a single canonical run, uploaded to `2026_UTSW/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/crispr_pipeline/`.

See [Issue: Gersbach Hep deliverables](https://github.com/adamklie/tf_perturb_seq/issues/gersbach-hep-deliverables.md) for the full ask + reference example (Huangfu DE [`syn74834952`](https://www.synapse.org/Synapse:syn74834952) is the layout target).

## Schema + walkthrough

- Machine-readable schema: [`schemas/crispr_pipeline.json`](../../../schemas/crispr_pipeline.json)
- Analysis-level walkthrough: [`docs/analysis/CRISPR_PIPELINE_OUTPUTS.md`](../../../../../analysis/CRISPR_PIPELINE_OUTPUTS.md)
