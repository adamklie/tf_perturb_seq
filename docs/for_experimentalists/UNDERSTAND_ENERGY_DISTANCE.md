# Understanding energy-distance outputs

Energy distance is the headline statistic for "did knocking down this TF move the transcriptome." It collapses thousands of single-cell profiles into one number per perturbation. Companion: [`UNDERSTAND_CRISPR_OUTPUTS.md`](UNDERSTAND_CRISPR_OUTPUTS.md), [`UNDERSTAND_CNMF_OUTPUTS.md`](UNDERSTAND_CNMF_OUTPUTS.md).

For the pipeline reference (run commands, intermediate files), see [`docs/analysis/energy_dist/ENERGY_DISTANCE_OUTPUTS.md`](../analysis/energy_dist/ENERGY_DISTANCE_OUTPUTS.md). This doc is the **interpretation** layer.

---

## What energy distance actually measures

Pick a TF. Look at the cloud of cells with that TF's guides (in PCA / latent space). Look at the cloud of NTC cells. **Energy distance** = a number that says "how different are those two clouds, end-to-end." Big → the perturbation reshaped the transcriptome. Zero → the two clouds are indistinguishable.

Two things make this useful:

1. **It's distribution-aware** — not just "did the mean shift" but "did the *shape* of the cell-state distribution change." Sensitive to subtle effects a t-test would miss.
2. **It has a built-in null** — we run a permutation test by shuffling cell labels thousands of times. The p-value tells us "is the real distance bigger than what shuffling would give."

---

## The files you'll get

A complete energy-distance run produces this layout under the dataset's `energy_distance/` folder:

```
energy_distance/<run_label>/
├── pval_edist_full.csv             # ← THE headline results table
├── targeting_outlier_table.csv     # gRNAs with disagreeing behavior (flagged + dropped)
├── non_targeting_outlier_table.csv # NTC gRNAs with anomalous behavior
├── target_by_target_matrix.csv     # square distance matrix between every TF pair (Step 3, optional)
├── edist_embedding_info.csv        # 2D t-SNE of TFs colored by cluster (Step 3, optional)
├── figures/
│   ├── e-dist_distribution.pdf     # distance distribution + p-value scatter
│   ├── e-dist_cutoff_value.pdf     # heatmap: #sig targets at various thresholds
│   └── e-dist_cutoff_value_NEG_CONTROL.pdf
├── config_step2.json               # the params used for this run
└── logs/
```

The headline file is `pval_edist_full.csv`.

---

## Reading `pval_edist_full.csv`

One row per target (TF). Columns of interest:

| Column | Meaning |
|---|---|
| `target_id` | Format `<ENSG>\|<chr>:<start>-<end>` — joins to TF metadata via the ENSG. |
| `type` | `target` (a TF being knocked down), `positive control` (e.g. AARS — should be huge), `negative control` (safe-harbor region — should be ~baseline). NTCs are *not* in this table; they're the implicit baseline. |
| `cell_count` | Cells used for this target's distance computation. |
| `distance_mean` | Average energy distance across 20 permutation rounds. **The effect-size number.** |
| `pval_mean` | Mean p-value from those 20 rounds. ⚠ See calibration caveat below. |
| `distance_0..distance_19` | Per-permutation values (useful for variance / quality control). |
| `pval_0..pval_19` | Same for p-values. |

`distance_mean_log` and `pval_mean_log` are log-transformed companions for plotting.

---

## How to actually use this

### "Which TFs significantly altered the transcriptome?"

Two complementary thresholds:

| Threshold | Meaning | When to use |
|---|---|---|
| `pval_mean < 0.05` | Permutation-test cutoff | Default — most papers report this. ⚠ See calibration caveat. |
| `distance_mean > (max NC distance)` | Effect-size cutoff: target beats the worst negative control | **Calibration-robust** — works even if p-values are off. |

Recommended workflow:

```python
import pandas as pd
ed = pd.read_csv("pval_edist_full.csv", index_col=0)

targeting = ed[ed["type"] == "target"]
nc = ed[ed["type"] == "negative control"]
nc_max = nc["distance_mean"].max()

sig_by_pval = targeting[targeting["pval_mean"] < 0.05]
sig_by_dist = targeting[targeting["distance_mean"] > nc_max]

print(f"By p-value:        {len(sig_by_pval)}")
print(f"By distance > NCs: {len(sig_by_dist)}")
```

If the two counts disagree wildly, the run has a calibration problem (next section).

### "Which TFs cluster together?"

If `target_by_target_matrix.csv` exists, it's a square TF × TF distance matrix you can hierarchical-cluster:

```python
import seaborn as sns
m = pd.read_csv("target_by_target_matrix.csv", index_col=0)
sns.clustermap(m, cmap="viridis", figsize=(8, 8))
```

TFs in the same cluster moved the transcriptome in similar ways → candidate co-functional / co-regulated pairs.

### "Where do my favorite TFs sit?"

If `edist_embedding_info.csv` exists, it has a 2D t-SNE of TFs:

```python
emb = pd.read_csv("edist_embedding_info.csv")
import matplotlib.pyplot as plt
plt.scatter(emb["x"], emb["y"], c=emb["cluster"], cmap="tab20")
for tf in ["SOX17", "GATA4"]:
    row = emb[emb["index"].str.contains(tf, na=False)]
    if len(row):
        plt.annotate(tf, (row.iloc[0]["x"], row.iloc[0]["y"]))
```

---

## ⚠ The calibration caveat (read this)

Some energy-distance runs produce **anti-conservative** p-values: all of the negative controls (NTC and same-region) come out with `pval_mean == 0`, which is obviously wrong (NTCs *should* be the null). The root cause is usually a mismatch between the PCA basis the test was run on and the cells' actual variance structure (e.g. PCA computed on all genes vs. just HVGs).

**When it's happening:** check `n_NCs_pval_eq_0` in the QC summary. If >0, p-values can't be trusted alone.

**Workaround:** report `distance_mean > NC max` as the calibration-robust significance count and use `distance_mean` as the effect-size estimate. Don't gate on raw `pval_mean < 0.05` until the calibration is fixed.

---

## What to look for (wet-lab takeaways)

1. **The top-distance TFs** for the cell type — these are the master regulators. In endoderm: SOX17, FOXH1. In cardiomyocytes: ISL1, TBX20, MEF2C. In endothelium: ETV2 (when run).
2. **The "always-on" TFs** — TFs significant in *every* lineage. These tend to be basic-cellular-machinery factors (transcription apparatus, telomere maintenance) rather than lineage-specific.
3. **The "lineage-specific" TFs** — TFs significant in only one cell type. These are the lineage-defining regulators and the most interpretable hits.
4. **Where positive controls land** — e.g. AARS should have a *huge* `distance_mean`. If it doesn't, the experiment didn't capture strong perturbations and other numbers should be discounted.
5. **The TF × TF clustering structure** — pairs that cluster tightly are candidate co-regulators or co-functional TFs; pull them for follow-up.

---

## Who to ask

- **Pipeline / method questions** — Chikara Takeuchi (UTSW): <https://github.com/Chikara-Takeuchi/energy_dist_pipeline>.
- **Calibration issues** — Adam Klie (`aklie@ucsd.edu`).
- **Project-specific context** — your dataset's lead in [`docs/TEAM.md`](../TEAM.md).
