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

| Task | Description | Status | Notebook | Output |
|:---:|---|:---:|---|---|
| 1 | Bar of #TFs significantly altering transcriptome per dataset | ✅ | [`significant_tf_counts.ipynb`](significant_tf_counts.ipynb) | [`results/significant_tf_counts/`](results/significant_tf_counts/) |
| 2 | Consistency of guide outlier detection (QC of guides) | ⏳ | [`guide_outlier_jaccard.ipynb`](guide_outlier_jaccard.ipynb) | [`results/guide_outlier_jaccard/`](results/guide_outlier_jaccard/) |
| 5/6 | Cutoff sweeps splitting targeting vs. NTC | ⏳ | [`pval_cutoff_sweep.ipynb`](pval_cutoff_sweep.ipynb) | [`results/pval_cutoff_sweep/`](results/pval_cutoff_sweep/) |
| 7 | Mean-distance scatterplots (pairwise across datasets) | ⏳ | [`pairwise_distance_scatter.ipynb`](pairwise_distance_scatter.ipynb) | [`results/pairwise_distance_scatter/`](results/pairwise_distance_scatter/) |
| 8 | Distance heatmap (rank- or quantile-normalize) | ⏳ | [`distance_heatmap.ipynb`](distance_heatmap.ipynb) | [`results/distance_heatmap/`](results/distance_heatmap/) |
| 3 | Jaccard distance of significant-TF sets across datasets | ⏳ | — | — |
| 4 | Binarize p-values → Venn / UpSet overlaps | ⏳ | — | — |
| 9 | Repeat per-dataset, then cross-dataset | ⏳ | — | (umbrella convention for 3-8) |
| 10 | Cluster TF perturbations by energy distance per lineage; flag TFs with sharp cross-lineage differences | ⏳ | — | — |
| 11 | Spot-based clustering annotation | ⏳ | — | — |
| 12 | Quantify similarity of clusterings across lineages | ⏳ | — | — |

## Notebooks

Each notebook reads from its matching `results/<name>/` subdir (intermediate TSVs are pre-built by the `.py` scripts under [`scripts/`](scripts/)) and writes its plots back into the same subdir.

| Notebook | What it does | Reads | Writes |
|---|---|---|---|
| [`significant_tf_counts.ipynb`](significant_tf_counts.ipynb) | Bar plot, two bars per dataset (`distance > NC max` vs `pval_mean < 0.05`); hatched bars where calibration is anti-conservative. | `results/significant_tf_counts/significant_tf_counts.tsv` | `results/significant_tf_counts/significant_tf_counts.{pdf,png}` |
| [`pval_cutoff_sweep.ipynb`](pval_cutoff_sweep.ipynb) | Per-dataset, sweep p-value cutoffs and split hits into targeting vs NC; right panel = % NCs among hits. | `data/<dataset>/energy_distance/wg1_significant_tfs.tsv` | `results/pval_cutoff_sweep/cutoff_sweep_<short>.pdf` |
| [`pairwise_distance_scatter.ipynb`](pairwise_distance_scatter.ipynb) | Headline pair (Hon CM vs Gersbach Hep) + all-pairs lower-triangle grid of `distance_mean` scatters. | `results/pairwise_distance_scatter/per_target_long.tsv` | `results/pairwise_distance_scatter/distance_scatter_*.pdf` |
| [`distance_heatmap.ipynb`](distance_heatmap.ipynb) | Targets × datasets heatmap (raw / rank / quantile-normalized panels), restricted to targets shared across all 4 datasets. | `results/distance_heatmap/per_target_wide.tsv` | `results/distance_heatmap/distance_heatmap.pdf` |
| [`guide_outlier_jaccard.ipynb`](guide_outlier_jaccard.ipynb) | Pairwise Jaccard heatmap between per-dataset outlier-gRNA sets (`pval_outlier < 0.05`). | `results/guide_outlier_jaccard/guide_outlier_jaccard.tsv` | `results/guide_outlier_jaccard/guide_outlier_jaccard_heatmap.pdf` |

## Scripts

Each script does one thing — single input, single output. Run from this folder (`working_groups/wg1_data_qc/edist/`):

| Script | What it does | Input → Output |
|---|---|---|
| [`scripts/count_significant_tfs.py`](scripts/count_significant_tfs.py) | Count significant TFs in one dataset under both criteria | `<wg1_significant_tfs.tsv>` → 1-row TSV |
| [`scripts/combine_count_tables.py`](scripts/combine_count_tables.py) | Concatenate per-dataset count TSVs | directory of 1-row TSVs → combined TSV |
| [`scripts/build_long_per_target_table.py`](scripts/build_long_per_target_table.py) | Stack per-dataset `wg1_significant_tfs.tsv` into one long table | `data/` dir → long TSV |
| [`scripts/build_wide_per_target_table.py`](scripts/build_wide_per_target_table.py) | Pivot per-dataset `wg1_significant_tfs.tsv` into a wide table (one row per target, dataset stats side-by-side) | `data/` dir → wide TSV |
| [`scripts/compute_guide_outlier_jaccard.py`](scripts/compute_guide_outlier_jaccard.py) | Pairwise Jaccard between per-dataset outlier-gRNA sets | per-dataset `targeting_outlier_table.csv` → Jaccard TSV |

End-to-end for Task 1 (significant TF counts):

```bash
cd docs/jamborees/2026_UTSW/working_groups/wg1_data_qc/edist

for ds in \
    Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq \
    Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq \
    Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq
do
  uv run python scripts/count_significant_tfs.py \
    --input  ../../../data/$ds/energy_distance/wg1_significant_tfs.tsv \
    --output results/significant_tf_counts/per_dataset/$ds.tsv
done

uv run python scripts/combine_count_tables.py \
  --input  results/significant_tf_counts/per_dataset \
  --output results/significant_tf_counts/significant_tf_counts.tsv

# Open significant_tf_counts.ipynb in Jupyter and Run All — writes the {pdf,png}.
```

## Output directory

```
results/
  significant_tf_counts/                # significant_tf_counts.ipynb (+ count + combine .py)
    per_dataset/<dataset>.tsv           #   from count_significant_tfs.py
    significant_tf_counts.tsv           #   from combine_count_tables.py
    significant_tf_counts.{pdf,png}     #   from significant_tf_counts.ipynb
  pval_cutoff_sweep/                    # pval_cutoff_sweep.ipynb (reads from data/, no .py prereq)
    cutoff_sweep_<short>.pdf
  pairwise_distance_scatter/            # pairwise_distance_scatter.ipynb (+ build_long_per_target_table.py)
    per_target_long.tsv
    distance_scatter_HonCM_vs_GersbachHep.pdf
    distance_scatter_all_pairs.pdf
  distance_heatmap/                     # distance_heatmap.ipynb (+ build_wide_per_target_table.py)
    per_target_wide.tsv
    distance_heatmap.pdf
  guide_outlier_jaccard/                # guide_outlier_jaccard.ipynb (+ compute_guide_outlier_jaccard.py)
    guide_outlier_jaccard.tsv
    guide_outlier_jaccard_heatmap.pdf
```

## Headline numbers (Task 1)

![Energy-distance significance per dataset](results/significant_tf_counts/significant_tf_counts.png)

| Dataset | n targeting | sig (distance > NC max) | sig (pval_mean < 0.05) | Calibration |
|---|---:|---:|---:|---|
| Hon CM | 1,932 | **162** | 269 | healthy |
| Huangfu DE | 2,169 | **72** | 2,169 *(all)* | anti-conservative |
| Huangfu ESC | 2,169 | **83** | 2,169 *(all)* | anti-conservative |

The distance-based count is the trustworthy number; the p-value count for the Huangfu runs is degenerate until the HVG-preprocess fix lands.
