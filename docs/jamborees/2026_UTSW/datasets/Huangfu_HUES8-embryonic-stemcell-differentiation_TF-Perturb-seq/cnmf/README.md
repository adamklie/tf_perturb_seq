# Huangfu HUES8 Embryonic Stem Cell — cNMF

**Status**: ✅ Complete on Synapse [`syn74893846`](https://www.synapse.org/Synapse:syn74893846) — `042926_huangfu_esc_torchcnmf_KskillA`. Selected k = 200, density threshold = 2.0.

## What's being mirrored

| | |
|---|---|
| Run name | `042926_huangfu_esc_torchcnmf_KskillA` |
| Source | HPC: `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_esc_torchcnmf_KskillA/` (27 GB) |
| Selected k | **200** |
| Density threshold | **2.0** (single — this run did not sweep dt=0.2 + 2.0) |
| Mirror target | [`syn74893846`](https://www.synapse.org/Synapse:syn74893846) → `2026_UTSW/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/cnmf/` (no run_name nesting — only one run per dataset for the jamboree) |
| Bundle size | ~9 GB after curation |
| Log on HPC | `docs/jamborees/2026_UTSW/mirror_huangfu_esc_cnmf.log` |

## ⚠ Layout deviations from the schema

Same as the DE sibling — older PerturbNMF tooling produces a layout that differs from the canonical schema. Mirror script auto-detects and uploads as-is:

- File prefix `Inference.` (instead of `<run_name>.`)
- Flat files live in `Inference/` subdir
- `Evaluation/` (instead of `Eval/`)
- `Annotation/` at `Inference/Annotation/`
- Single density threshold `dt_2_0` only

K sweep matches DE: 30, 50, 60, 80, 100, 200, 250, 300.

## Run parameters

Identical to DE — see [`../Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/cnmf/README.md`](../../Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/cnmf/README.md) for the full table comparing the 042926 K-skill-A params to Hon's "guiding light" reference.

## Schema + walkthrough

- Machine-readable schema: [`schemas/cnmf.json`](../../../schemas/cnmf.json)
- Analysis-level walkthrough: [`docs/analysis/cNMF_OUTPUTS.md`](../../../../../analysis/cNMF_OUTPUTS.md), [`docs/analysis/cNMF.md`](../../../../../analysis/cNMF.md)

## Pre-staged Hon-mirroring run (held)

Same as DE — `050926_HuangfuESC_20iter_5KHVG_torch_halsvar_batch.sh` is pre-staged but held since the 042926 run is sufficient.
