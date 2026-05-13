# Energy distance — Hon WTC11 Cardiomyocyte

Two ED runs on Synapse, **both run on the same source MuData** (Weizhou's, [`syn74522725`](https://www.synapse.org/Synapse:syn74522725)) — kept separate for cross-lab comparison.

| Bundle | Synapse | Owner | Completeness |
|---|---|---|---|
| `energy_distance/` | [`syn74897350`](https://www.synapse.org/Synapse:syn74897350) | Adam | Partial — result CSVs + plots + slurm logs. **Missing**: inference_mudata, pickles. |
| `energy_distance_gersbach_comp/` | [`syn74910330`](https://www.synapse.org/Synapse:syn74910330) | Sara | Full — inference_mudata + pickles + result CSVs + plots. |

Both ran against the `2026_04_19_no_spacer` MuData (Weizhou's CRISPR pipeline output). Run-label on local HPC: `2026_04_19_no_spacer/energy_distance/`.

## How to reproduce

```bash
sbatch datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/2026_04_19_no_spacer/energy_distance/scripts/5_run_energy_distance.sh
```

The runner pulls Weizhou's MuData from Synapse (`--synapse-id syn74522725`).

## Schema + walkthrough

- Machine-readable schema: [`../../../schemas/energy_distance.json`](../../../schemas/energy_distance.json)
- Analysis-level walkthrough: [`docs/analysis/ENERGY_DISTANCE_OUTPUTS.md`](../../../../../analysis/ENERGY_DISTANCE_OUTPUTS.md)

## Notes

- Hon CM is the only production dataset where the source MuData lives on Synapse (not GCS), because Hon Lab uploaded it directly.
- `2026_04_19_no_spacer` is named for the seqspec config used — `no_spacer` means the run was done without `spacer_tag`. Different from our local `seqspec_v3` run on GCS, which kept the spacer.
