# Energy-distance — WG1 analyses

Cross-dataset analyses of the energy-distance pipeline outputs for the five TFP3 production datasets, in preparation for WG1 Figure 1 panels at the [2026 UTSW jamboree](../../../README.md). Source per-dataset tables live under [`docs/jamborees/2026_UTSW/data/<dataset>/energy_distance/`](../../../data/).

## Datasets

| Dataset | ED run | `wg1_significant_tfs.tsv` | Calibration |
|---|---|:---:|---|
| Hon WTC11 Cardiomyocyte | not run on Hon h5mu — *waiting on rerun strategy* | ✅ (from HTv2 testbed pipeline output reshaped) | healthy |
| Huangfu HUES8 Definitive Endoderm | `muddy_penguin` (2026-05-09) | ✅ | anti-conservative |
| Huangfu HUES8 Embryonic Stem Cell | `sceptre_v1` (2026-05-09) | ✅ | anti-conservative |
| Gersbach WTC11 Hepatocyte | not run — awaiting Sara's deliverables | ❌ | — |
| Engreitz WTC11 Endothelial | - | - | — |

See [`data/schemas/energy_distance.json`](../../../data/schemas/energy_distance.json) for the full output schema and the known-issues block.

## Significance criteria

The per-dataset `wg1_significant_tfs.tsv` carries two boolean significance flags:

- **`sig_distance_gt_NC_max`** — `distance_mean` exceeds the max distance over the 100 NTC controls. Non-parametric, calibration-robust. **Use this as the headline criterion** until p-values are recalibrated.
- **`sig_pval_lt_0p05`** — `pval_mean < 0.05` from the permutation test. Currently anti-conservative on the Huangfu runs: all 100 NTCs have `pval_mean = 0`. Tracked in the schema's `known_issues.p_value_calibration_huangfu_runs`; fix is to re-preprocess with HVG selection in `scripts/preprocess_mudata_local.py`. Until then, p-value-based panels carry a hatched-bar caveat.

The brainstorm tasks below also reference an "OR-gene" class. That class only exists in the **HTv2 benchmark** library (54 OR negative controls); the production datasets carry only 100 NTCs + 4 PCs + ~2 k targeting. Tasks framed against OR genes are re-interpreted against NTCs for production data, or flagged as out-of-scope for the jamboree.

## Status

| Task | Description | Status | Output |
|:---:|---|:---:|---|
| 1 | Bar / UpSet of #TFs significantly altering transcriptome per dataset | ✅ | [`results/significant_tf_counts.{tsv,pdf,png}`](results/) |
| 2 | Consistency of guide outlier detection (QC of guides) | ⏳ | — |
| 3 | Jaccard distance of significant-TF sets across datasets | ⏳ | — |
| 4 | Binarize p-values → Venn / UpSet overlaps | ⏳ | — |
| 5 | Cutoff sweeps splitting targeting vs. NTC | ⏳ | (re-framed from "OR genes" — production has no OR class) |
| 6 | Largest cutoff where no NTC is called significant | ⏳ | (re-framed from "OR genes") |
| 7 | Mean-distance scatterplots (pairwise across datasets) | ⏳ | — |
| 8 | Distance heatmap (rank- or quantile-normalize) | ⏳ | — |
| 9 | Repeat per-dataset, then cross-dataset | ⏳ | (umbrella convention for 3-8) |
| 10 | Cluster TF perturbations by energy distance per lineage; flag TFs with sharp cross-lineage differences | ⏳ | — |
| 11 | Spot-based clustering annotation | ⏳ | — |
| 12 | Quantify similarity of clusterings across lineages | ⏳ | — |

## Scripts

Each script does one thing — single input, single output. Run from this folder (`working_groups/wg1_data_qc/edist/`):

| Script | What it does | Input → Output |
|---|---|---|
| [`scripts/count_significant_tfs.py`](scripts/count_significant_tfs.py) | Count significant TFs in one dataset under both criteria | `<wg1_significant_tfs.tsv>` → 1-row TSV |
| [`scripts/combine_count_tables.py`](scripts/combine_count_tables.py) | Concatenate per-dataset count TSVs | directory of 1-row TSVs → combined TSV |
| [`scripts/plot_significant_tf_counts.py`](scripts/plot_significant_tf_counts.py) | Bar plot, two bars per dataset | combined TSV → PDF (or PNG) |

End-to-end for Task 1:

```bash
cd docs/jamborees/2026_UTSW/working_groups/wg1_data_qc/edist

for ds in \
    Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq \
    Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq \
    Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq
do
  uv run python scripts/count_significant_tfs.py \
    --input  ../../../data/$ds/energy_distance/wg1_significant_tfs.tsv \
    --output results/per_dataset/$ds.tsv
done

uv run python scripts/combine_count_tables.py \
  --input  results/per_dataset \
  --output results/significant_tf_counts.tsv

uv run python scripts/plot_significant_tf_counts.py \
  --input  results/significant_tf_counts.tsv \
  --output results/significant_tf_counts.pdf
```

## Output directory

```
results/
  per_dataset/
    <dataset>.tsv                       # 1-row TSV per dataset (count_significant_tfs.py)
  significant_tf_counts.tsv             # combined (combine_count_tables.py)
  significant_tf_counts.{pdf,png}       # bar plot (plot_significant_tf_counts.py)
```

## Headline numbers (Task 1)

![Energy-distance significance per dataset](results/significant_tf_counts.png)

| Dataset | n targeting | sig (distance > NC max) | sig (pval_mean < 0.05) | Calibration |
|---|---:|---:|---:|---|
| Hon CM | 1,932 | **162** | 269 | healthy |
| Huangfu DE | 2,169 | **72** | 2,169 *(all)* | anti-conservative |
| Huangfu ESC | 2,169 | **83** | 2,169 *(all)* | anti-conservative |

The distance-based count is the trustworthy number; the p-value count for the Huangfu runs is degenerate until the HVG-preprocess fix lands.
