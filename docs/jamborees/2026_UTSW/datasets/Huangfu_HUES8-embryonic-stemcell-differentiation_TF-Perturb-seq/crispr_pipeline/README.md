# Huangfu HUES8 Embryonic Stem Cell — CRISPR pipeline

**Status**: ✅ Canonical 3-folder bundle on Synapse.

## Synapse

[`syn74835010`](https://www.synapse.org/Synapse:syn74835010) — canonical layout:

```
syn74835010/
├── pipeline_dashboard/    (~40 GB; dashboard.html, inference_mudata.h5mu, additional_qc/, evaluation_output/, figures/, ...)
├── pipeline_info/         (params JSON + software versions YAML)
└── pipeline_outputs/      (~23 GB; perturbo cis/trans per-element/per-guide TSVs)
```

Mirrored 2026-05-07 from GCS via [`scripts/mirror_pipeline_outputs.py`](../../../scripts/mirror_pipeline_outputs.py) (run from HPC with `--workdir /cellar/users/aklie/scratch/...`).

## Source

`gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/2026_04_13/outs/sceptre_v1/`

## Schema + walkthrough

- Machine-readable schema: [`schemas/crispr_pipeline.json`](../../../schemas/crispr_pipeline.json)
- Analysis-level walkthrough: [`docs/analysis/CRISPR_PIPELINE_OUTPUTS.md`](../../../../../analysis/CRISPR_PIPELINE_OUTPUTS.md)
