# Huangfu HUES8 Definitive Endoderm — CRISPR pipeline

**Status**: ✅ Canonical 3-folder bundle on Synapse.

## Synapse

[`syn74834952`](https://www.synapse.org/Synapse:syn74834952) — canonical layout:

```
syn74834952/
├── pipeline_dashboard/    (~40 GB; dashboard.html, inference_mudata.h5mu, additional_qc/, evaluation_output/, figures/, ...)
├── pipeline_info/         (params JSON + software versions YAML)
└── pipeline_outputs/      (~23 GB; perturbo cis/trans per-element/per-guide TSVs)
```

Mirrored 2026-05-07 from GCS via [`data/mirror_pipeline_outputs.py`](../../../data/mirror_pipeline_outputs.py) (run from HPC with `--workdir /cellar/users/aklie/scratch/...`).

## Source

`gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/2026_04_09/outs/muddy_penguin/`

## Schema + walkthrough

- Machine-readable schema: [`schemas/crispr_pipeline.json`](../../../schemas/crispr_pipeline.json)
- Analysis-level walkthrough: [`docs/analysis/CRISPR_PIPELINE_OUTPUTS.md`](../../../../../analysis/CRISPR_PIPELINE_OUTPUTS.md)
