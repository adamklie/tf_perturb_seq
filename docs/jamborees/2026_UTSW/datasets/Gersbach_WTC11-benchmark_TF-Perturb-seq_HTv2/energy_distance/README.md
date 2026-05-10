# Energy distance — Gersbach WTC11 benchmark HTv2 (testbed)

| | |
|---|---|
| Dataset ID | `Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2` |
| Schema | [`../../../schemas/energy_distance.json`](../../../schemas/energy_distance.json) |
| Run label | `cleanser_800_mito_15pc` |
| Source MuData | HPC: `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/runs/cleanser_800_mito_15pc/pipeline_dashboard/inference_mudata.h5mu` (1.1 GB) |
| Status | 🔄 RUNNING — SLURM job 10641084 launched 2026-05-10 (carter-gpu). Smaller dataset (~37k cells) — should finish in a few hours, faster than the production runs. |
| HPC output dir | `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/results/energy_distance/cleanser_800_mito_15pc/` |
| Synapse | _will populate after run completes via `mirror_edistance_outputs.py`_ |

## How to run

```bash
ssh aklie@nrnb-login.ucsd.edu
sbatch /cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/5_run_energy_distance.sh
```

Wraps the shared runner `/cellar/users/aklie/projects/tf_perturb_seq/scripts/run_energy_distance_pipeline.sh` with `--mudata-path` (HPC-local input).

## Mirror outputs to Synapse (after run completes)

```bash
.venv/bin/python docs/jamborees/2026_UTSW/scripts/mirror_edistance_outputs.py \
  --dataset Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2 \
  --source-dir /cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/results/energy_distance/cleanser_800_mito_15pc
```

## Notes

- **Testbed**, not a production dataset. Useful as a self-contained energy-distance reference at small scale to compare against the production Huangfu DE/ESC runs (which had calibration concerns; see [Issue 1](../../../issues/edistance-calibration.md)).
- Source pipeline run is the same `cleanser_800_mito_15pc` that's mirrored to Synapse [`syn74885574`](https://www.synapse.org/Synapse:syn74885574) for the CRISPR pipeline.
- HTv2 has a verified energy distance reference run on Synapse at [`syn74381167`](https://www.synapse.org/Synapse:syn74381167) (Sara's Duke run, used to verify the schema). Our run here is structurally parallel; the `pval_edist_full.csv` should be schema-identical.
