# Energy distance — Hon WTC11 Cardiomyocyte

| | |
|---|---|
| Dataset ID | `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq` |
| Schema | [`../../../schemas/energy_distance.json`](../../../schemas/energy_distance.json) |
| Run label | TBD (canonical CRISPR pipeline run is `2026_04_19_no_spacer` per Synapse `syn74520421`) |
| Source MuData | Synapse [`syn74522725`](https://www.synapse.org/Synapse:syn74522725) (Hon Lab upload, 16.65 GB) |
| Status | Not run yet — awaiting `pipeline_info/` from the Hon team to confirm run lineage before kicking off |
| HPC output dir | `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/results/energy_distance/<run>/` (TBD) |
| Synapse | _not yet created_ |

## How to run (when ready)

The dataset doesn't yet have a `5_run_energy_distance.sh` — needs scaffolding once the canonical run label is settled. Template:

```bash
# adapt scripts/run_energy_distance_pipeline.sh to point at the Hon CM MuData (Synapse syn74522725)
sbatch /cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/5_run_energy_distance.sh
```

The runner accepts `--synapse-id syn74522725` if pulling directly from Synapse rather than GCS.

## Notes

- Hon CM is the only dataset where the source MuData lives on Synapse (not GCS), because Hon Lab uploaded it directly. The runner script supports `--synapse-id` for this case.
- The Synapse parent folder for this dataset's CRISPR pipeline is named `2026_04_19_no_spacer`, suggesting the run was done with **no `spacer_tag`** — worth confirming with the Hon team before running e-distance, since spacer config affects guide assignment quality and downstream e-distance interpretability.
