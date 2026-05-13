# Energy distance — Huangfu HUES8 Definitive Endoderm

| | |
|---|---|
| Dataset ID | `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq` |
| Schema | [`../../../schemas/energy_distance.json`](../../../schemas/energy_distance.json) |
| Run label | `muddy_penguin` |
| Source MuData | `gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/2026_04_09/outs/muddy_penguin/inference_mudata.h5mu` (18.02 GB) |
| Status | ✅ Steps 1+2+2.1 complete (2026-05-09; SLURM job 10547114, 9h 34m on carter-gpu-02). Step 3 also complete (job 10573244, 6m 29s, filtered to top-200 by distance). |
| HPC output dir | `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/results/energy_distance/muddy_penguin/` |
| Synapse | [`syn74883327`](https://www.synapse.org/Synapse:syn74883327) |
| Validation | All 4 layers PASS (file presence, schema, value ranges, schema-identity vs HTv2 reference). 2267 targets × 46 cols. ⚠ Calibration concern, see below. |
| Step 3 cutoff | `distance_cutoff=105.11` (top-200 by distance; just above NC max 109.9), `pval_cutoff=1.0` (effectively disabled — pvals mis-calibrated). Yields 193 targets in `target_by_target_matrix.csv` (193×193) and TSNE+AffinityPropagation clustering in `edist_embedding_info.csv`. |

## How to run / re-run

From the HPC:

```bash
ssh aklie@nrnb-login.ucsd.edu
sbatch /cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/5_run_energy_distance.sh
```

Wraps the shared runner `/cellar/users/aklie/projects/tf_perturb_seq/scripts/run_energy_distance_pipeline.sh`. Container `edist_pipeline.sif` is shared across datasets; pulled on first run.

## Mirror outputs to Synapse (after run completes)

```bash
# from the HPC, in the project root
.venv/bin/python docs/jamborees/2026_UTSW/data/mirror_edistance_outputs.py \
  --dataset Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq \
  --source-dir /cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/results/energy_distance/muddy_penguin
```

## Notes

- Same source MuData as the Huangfu ESC e-distance run uses for ESC; both share the construct library (`IGVFDS3299AXST`) and guide file (`IGVFFI8270UPKB`).
- `muddy_penguin` is the canonical run because of its 13 bp `spacer_tag` (`GAGTACATGGGGG`); the earlier `entertaining_hamster` (12 bp) had ~6× lower guide UMI capture and is not used.

## ⚠ Calibration concern (2026-05-09)

**The p-values in this run should not be used for significance thresholding.** Diagnosis:
- All 100 negative-control targets have `pval_mean = 0` (should be uniform on [0,1])
- NC distance_mean median = 91.7 vs targeting distance_mean median = 93.2 (ratio 1.02× — practically identical)
- HTv2 verified reference (syn74381183) has distance_mean median **1.78** with frac p<0.05 = 48% (healthy distribution)
- Our distances are 1-2 orders of magnitude larger than HTv2's, and the permutation null is too tight relative to that scale

Likely causes (TBD which combination):
- Used **all genes** with `min_counts≥1` (~9k features) for PCA, vs HTv2 likely used HVG subset
- Library is 14k gRNAs / 600 non-targeting (vs HTv2 60-target benchmark) — much larger non-targeting pool tightens the null
- `use_matched_bg=false` in our config — try `true` if re-running

**Energy distance values are still informative as effect sizes** (targeting range [50, 134] is broader than NC range [76, 110]; some real signal embedded). Preferred path: threshold on raw distance (e.g., distance > NC max ≈ 110 in DE) rather than p-value, until calibration is fixed.

Same calibration issue affects the ESC sibling run; both bundle the same flaw.
