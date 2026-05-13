# Gersbach WTC11 benchmark HTv2 (testbed)

**Lab**: Gersbach (Duke)  •  **Cell line**: WTC11  •  **Use**: small testbed for verifying pipeline structure

> ⚠ **Not a production dataset.** HTv2 is mirrored here purely as a small/fast scaffold for shaking down the pipeline structure (CRISPR pipeline mirror, cNMF launch, energy distance schema). It's not part of the 5 production datasets and is not enshrined in the schemas as a "verified reference."

## Status

| Output | Status | Synapse |
|---|---|---|
| CRISPR pipeline | ✅ canonical 3-folder bundle | [`syn74885574`](https://www.synapse.org/Synapse:syn74885574) |
| cNMF | 🔄 testbed run in progress (SLURM job 10577039) | — |
| Energy distance | (used as schema-verification reference; not re-run here) | [`syn74381167`](https://www.synapse.org/Synapse:syn74381167) (Hon-lab benchmark; verified upstream) |

Deeper status: [`crispr_pipeline/README.md`](crispr_pipeline/README.md), [`cnmf/README.md`](cnmf/README.md). Open issue: [HTv2 cNMF testbed verification](https://github.com/adamklie/tf_perturb_seq/issues/htv2-cnmf-testbed.md).

## Source data

- **HPC**: `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/`
- **CRISPR pipeline source (HPC + GCS)**: `runs/cleanser_800_mito_15pc/` for dashboard + outputs, `gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/Benchmark_cleanser_800_mito_15pc/Gersbach_HTV2/pipeline_info/` for params.

## Pipeline configuration

Used the `cleanser_800_mito_15pc` CRISPR pipeline run (cleanser guide assignment, min 800 genes / cell, 15% mito threshold) — this is one of Lucas's benchmarking variants.

For the cNMF testbed run: mirrors Hon WTC11 benchmark's `030726_20iter_5KHVG_torch_halsvar_batch_e7` parameters exactly (halsvar / batch / 20 iter / 5K HVG / density thresholds 0.2 + 2.0 / categorical_key=batch / benchmark k list).

## Notes

- Smaller than production datasets (~37k cells × 13k genes) — runs fast enough to verify pipeline structure end-to-end in hours rather than days.
- HTv2 doesn't get the same canonical-reference treatment as the production datasets: schemas in [`schemas/`](../../schemas/) deliberately stay scoped to production. See [`issues/README.md`](https://github.com/adamklie/tf_perturb_seq/issues/README.md) for the rationale.
