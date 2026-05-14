# WG6 — Predictive modeling

Topic 3 / Figure 5. The goal is to design and prototype a model trained on the uniformly processed TFP3 outputs that can predict context-specific TF-perturbation impact. WG6 is mostly in discussion mode for this jamboree: end-of-jamboree deliverables are a written task specification (inputs, outputs, validation metrics) and, if time allows, a baseline prototype.

## Questions

- What model architectures best capitalize on TFP3 data as a training source (e.g., predicting context-specific TF perturbation impact)?
- What are the inputs, outputs, task definition, and validation metrics? Can implementation begin during the jamboree?
- Alternatively, which existing tools from member labs can be applied to this data source?

## Data

| Dataset | CRISPR pipeline | cNMF | Energy distance |
|---|:---:|:---:|:---:|
| Hon WTC11 Cardiomyocyte | ready | - | ready |
| Huangfu HUES8 Definitive Endoderm | ready | ready | ready |
| Huangfu HUES8 Embryonic Stem Cell | ready | ready | ready |
| Gersbach WTC11 Hepatocyte | ready | - | - |
| Engreitz WTC11 Endothelial | - | - | - |

Per-dataset cards live under [`../../data/`](../../data/). The planning artifact is [`task_spec_template.md`](task_spec_template.md) — fill it in with the group during the session.

## Issues

- *[FILL IN issue link]*: Task spec

## Working flow

Three steps, in order:

1. **Brainstorm** — agree on the task specification (inputs, outputs, validation metrics); fill in [`task_spec_template.md`](task_spec_template.md). Sketch any baseline figures or model architectures that come up.
2. **Execute** — if time allows, prototype a baseline; commit code (notebooks, scripts, the filled task spec) to this folder on GitHub.
3. **Share** — upload reusable outputs (filled task spec, slides, prototype code, results) to WG6's Synapse folder [`syn74954086`](https://www.synapse.org/Synapse:syn74954086) (mirrored `working_groups/wg6_predictive_modeling/` under [`syn64423137/2026_UTSW/`](https://www.synapse.org/Synapse:syn64423137)). Record the syn ID for each upload in the [Outputs](#outputs) table below so the next person can find it.

## Outputs

| Object | Syn ID | Description | Owner |
|---|---|---|---|
| *[FILL IN as outputs land]* | | | |
