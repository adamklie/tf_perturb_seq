# Energy distance p-value calibration concern — Huangfu DE + ESC

## TL;DR

Both Huangfu production runs (DE + ESC) on 2026-05-09 produce `pval_mean = 0` for all 100 negative-control targets, which is anti-conservative — negative controls should have a uniform p-value distribution. After investigating, the cause appears to be **scale-driven** (the size of the non-targeting pool, ~600 gRNAs in production vs ~30 in the HTv2 benchmark, gives a much tighter permutation null), **not** a pipeline bug, config deviation, or guide-metadata labeling difference. Our pipeline matches Chikara's HTv2 reference run bit-perfect.

We are using raw `distance_mean` as an effect-size proxy in the meantime. **Open question for Chikara: what is the recommended way to recover calibrated p-values at production scale?**

## The phenomenon

| Run | Targets | distance_mean median | NC distance_mean median | NC vs targeting ratio | NC pval_mean=0 |
|---|---|---|---|---|---|
| HTv2 reference (Chikara, [syn74381183](https://www.synapse.org/Synapse:syn74381183)) | 64 | 1.78 | (no explicit NC class) | — | (no NCs) |
| HTv2 our rerun (`cleanser_800_mito_15pc`) | 65 | 19.46 | (no explicit NC class) | — | (no NCs) |
| Our Huangfu DE | 2267 | 93.21 | 91.70 | 1.02× | **100/100** |
| Our Huangfu ESC | 2267 | 527.33 | 528.25 | 1.00× | **100/100** |

Negative controls have essentially the same energy distance as targeting gRNAs in the Huangfu runs, and **all 100 NC targets in each run get pval_mean=0**.

## Plots

- `01_pval_mean_by_type.png` — histogram of pval_mean by type per dataset. NCs piled at 0 in our Huangfu runs.
- `02_distance_mean_by_type.png` — distance_mean histograms. NC and targeting distributions overlap near-perfectly in Huangfu.
- `03_volcano_by_type.png` — distance vs -log10(p) scatter. Huangfu NCs (blue) sit at the top of the volcano next to targeting (red).
- `04_distance_scale_comparison.png` — log-scale comparison of distance_mean across HTv2 reference vs Huangfu. 1-2 orders of magnitude scale difference.

## What we ruled out

### 1. Pipeline / preprocess deviation (ruled out 2026-05-10)

Adam reran the e-distance pipeline on Gersbach HTv2 benchmark with our wrapper (`cleanser_800_mito_15pc`, 65 targets) and compared bit-for-bit against Chikara's HTv2 run on Synapse [syn74895081](https://www.synapse.org/Synapse:syn74895081):

```
65/65 targets overlap, 0 unique to either side
Pearson distance correlation: 1.0000
Pearson pval correlation:     1.0000
Per-target abs diff:          median 0, max 0
```

Identical results: same shape (65 × 46), same distance range [2.71, 79], same pval median (0.0316), same frac p<0.05 (58%). **Our wrapper reproduces Chikara's runs bit-perfect at benchmark scale**, ruling out preprocess / config / submodule / annotation differences as the cause.

### 2. Wrong background source (ruled out 2026-05-10)

Step 2 (`2_e_distance_nontargeting.py`) explicitly uses `non-targeting`-typed gRNAs as the permutation background (line 55: `clear_nt_sgRNA_list = nontargeting_outlier_df[nontargeting_outlier_df["pval_outlier"]>0.05].index.tolist()`). Negative controls are tested as targets, not used as background. Verified.

### 3. Guide-metadata labeling difference (ruled out 2026-05-10)

Both libraries use **OR (olfactory receptor) gene-targeting gRNAs as the de-facto negative control class**, just labeled differently:

| Library | "Negative control" labeling | OR-targeting gRNAs | OR pval_mean median | OR frac p<0.05 |
|---|---|---|---|---|
| HTv2 (416 gRNAs) | Class doesn't exist | 54, labeled `type=targeting` | 0.0654 | 44% |
| Huangfu DE (14k gRNAs) | Explicit `type=negative control` | 592 of 598 NCs match `^OR[digit]` | **0.0** | **100%** |

The labeling difference doesn't affect the test mechanically — the cells go through the same comparison either way. The fact that HTv2 OR-targeting (44% sig) sits below HTv2 real targeting (63% sig) shows calibration works at HTv2's scale even without an explicit NC class. The Huangfu NC class breakdown is a different phenomenon.

## What's left as the cause: non-targeting pool size

| | HTv2 | Huangfu DE/ESC |
|---|---|---|
| Total non-targeting gRNAs | 30 | 600 (20× more) |
| Non-targeting cells fed to permutation null | small pool | 20× larger |
| Permutation null variance | broad | tight (variance ~ 1/√N) |

With a much tighter permutation null in Huangfu, **any observed distance lands in the null tail** — including distances from genuine negative controls that should look indistinguishable from background. This is structural to the test under the current configuration; pre-bottleneck the test by limiting the NT pool, or move to a matched-background scheme (`use_matched_bg=true`?), and calibration may recover.

## Setup details

- Container: `docker.io/takechikara/energy_distance_env:latest` (apptainer .sif)
- Pipeline: `energy_dist_pipeline` @ `5821450a1eacfddb3c83be8c69fd593eaf76a61c` (post-downsampling-fix main)
- Preprocess: matches upstream `preprocess_mudata.py` exactly (modality `gene` → `filter_genes(min_counts=1)` → `normalize_total` → `log1p` → `scale` → `tl.pca(n_comps=50)`); only structural deviation is reading MuData from local disk instead of Synapse
- Step-1/2/2.1 config — identical to upstream `energy_dist_pipeline/config.json` defaults: `threshold_gRNA_num=6`, `combi_count=4`, `total_permute_disco=1000`, `combi_cell_num_max=1000`, `batch_num_basic=120` (filtering); `permute_per_bg=1000`, `num_of_bg=20`, `non_target_pick=2000`, `target_cell_num_max=2000`, `batch_num_basic=200`, `use_matched_bg=false` (permutation)
- Inputs (Huangfu DE/ESC, IGVF `IGVFFI8270UPKB` library pools A-D): 12,934 targeting + 600 non-targeting + 19 positive controls + 598 negative controls (99% OR-targeting). After grouping by `intended_target_promoter`, 2267 unique target regions

## Open questions for Chikara

1. **What's the recommended way to recover calibrated p-values at production scale?** Concretely: should we cap the non-targeting pool size (subsample to ~30 like HTv2)? Switch to `use_matched_bg=true`? Use a variance-correction in the null?
2. Do you have any internal benchmarks at production scale (~2000 targets, ~600 non-targeting gRNAs) where calibration is known good?
3. Is there a known scaling regime where the current default permutation parameters are validated?

## Synapse links

- Huangfu DE bundle: [syn74883327](https://www.synapse.org/Synapse:syn74883327)
- Huangfu ESC bundle: [syn74883475](https://www.synapse.org/Synapse:syn74883475)
- HTv2 reference (Chikara, newer, used for bit-perfect reproducibility check): [syn74895081](https://www.synapse.org/Synapse:syn74895081)
- HTv2 reference (Chikara, older, used for schema cross-check): [syn74381167](https://www.synapse.org/Synapse:syn74381167)
