# Energy distance — Gersbach WTC11 benchmark HTv2 (testbed)

| | |
|---|---|
| Dataset ID | `Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2` |
| Schema | [`../../../schemas/energy_distance.json`](../../../schemas/energy_distance.json) |
| Run label | `cleanser_800_mito_15pc` |
| Source MuData | HPC: `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/runs/cleanser_800_mito_15pc/pipeline_dashboard/inference_mudata.h5mu` (1.1 GB) |
| Status | ✅ COMPLETE — SLURM 10641084, 18 min wall time (2026-05-10). All 3 validator layers PASS. **65 targets** (no NCs in run; only positive control + targeting), distance_mean median 19.46, pval_mean median 0.0316 (non-degenerate). Step 3 deferred. |
| HPC output dir | `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/results/energy_distance/cleanser_800_mito_15pc/` |
| Synapse | [`syn74895081`](https://www.synapse.org/Synapse:syn74895081) (mirrored 2026-05-10) |

## Calibration sanity (relevant to [Issue 1](../../../issues/edistance-calibration.md))

The HTv2 testbed run shows **healthier calibration than the production Huangfu DE/ESC runs**:

| | HTv2 testbed (this run) | Huangfu DE | Huangfu ESC |
|---|---:|---:|---:|
| Cells | ~37k | 270k | 190k |
| Targets | 65 | 2267 | 2267 |
| `distance_mean` median (targeting) | **19.46** | 93.21 | 527.33 |
| `pval_mean` median | **0.0316** (non-degenerate) | 0 | 0 |

This supports the "scale-driven anti-conservatism" hypothesis: at smaller cell-count + target-count scale, the permutation null is wider and `pval_mean` doesn't collapse to 0. As the cell count and target count grow (Huangfu DE → ESC), distances inflate and the null tightens — exactly the opposite of what we'd want for calibrated significance testing.

Useful as a calibration reference point for diagnosing the production runs.

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
