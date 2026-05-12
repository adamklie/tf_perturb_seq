# QC outputs

Each row of `samples.tsv` produces three subdirs under `<outdir>/`. All files are prefixed with `<run_name>_<module>_*`.

```
<outdir>/
├── mapping_gene/
│   ├── <run_name>_gene_metrics.tsv
│   ├── <run_name>_gene_knee_plot.png
│   ├── <run_name>_gene_histograms.png
│   ├── <run_name>_gene_histograms_by_batch.png
│   └── <run_name>_gene_cells_per_batch.png
├── mapping_guide/
│   ├── <run_name>_guide_metrics.tsv
│   ├── <run_name>_guide_per_guide_capture.tsv
│   ├── <run_name>_guide_knee_plot.png
│   ├── <run_name>_guide_histograms.png
│   └── <run_name>_guide_histograms_by_batch.png
└── intended_target/
    ├── <run_name>_intended_target_results.tsv
    ├── <run_name>_intended_target_metrics.tsv
    ├── <run_name>_intended_target_volcano.png
    ├── <run_name>_intended_target_log2fc_distribution.png
    └── <run_name>_intended_target_roc_pr_curves.png
```

All `**/qc/**/*.{tsv,png,pdf,h5ad}` are gitignored. They're regenerable from `inference_mudata.h5mu` — don't commit them. Cross-sample comparisons that need the data pull from these TSVs directly or stage a copy under `scratch/`.

## `mapping_gene_metrics.tsv`

Per-batch + overall summary. One row per batch plus one `overall` row.

| Column group | Columns | Source |
|---|---|---|
| Identity | `batch` (`overall` or batch label), `n_cells` | obs |
| UMI | `umi_{median,mean,std,min,max,q25,q75}` | `obs["total_gene_umis"]` |
| Genes | `genes_{median,mean,std,min,max,q25,q75}` | `obs["num_expressed_genes"]` |
| Mito | `mito_{median,mean,std,min,max,q25,q75}` | `obs["percent_mito"]` |

Headline numbers per dataset: `umi_median`, `genes_median`, `mito_median` on the `overall` row.

## `mapping_guide_metrics.tsv`

Per-batch + overall summary.

| Column group | Columns | Source |
|---|---|---|
| Identity | `batch`, `n_cells` | obs |
| Guide UMI | `guide_umi_{median,mean,std,min,max,q25,q75}` | `obs["total_guide_umis"]` |
| Guides/cell | `guides_per_cell_{mean,std,min,max,median}` | computed (binarized `layers["guide_assignment"]`) |
| Assignment | `n_cells_with_guide`, `n_cells_exactly_1_guide`, `frac_cells_with_guide` | computed |
| Cells/guide | `n_guides_total`, `cells_per_guide_{median,mean,std,min,max}` | overall row only |

Headline numbers: `frac_cells_with_guide` (target >0.7 for production datasets), `guides_per_cell_median` (MOI proxy).

## `mapping_guide_per_guide_capture.tsv`

One row per guide. Columns:

- `guide_id` (var index)
- `label` (from `guide.var["label"]` — `targeting` / `non_targeting` / `safe_harbor` etc.)
- `gene_name` (from `guide.var["gene_name"]`)
- `n_cells_detected`, `frac_cells_detected`
- `total_umi`, `mean_umi`, `median_umi`, `std_umi`, `max_umi`

Used to flag dropout guides (`frac_cells_detected` below threshold) and to compute per-target guide redundancy.

## `intended_target_results.tsv`

One row per (targeting guide × its intended target gene). Includes the trans-test result from the pipeline.

Columns:

- `guide_id`
- `gene_name` (intended target)
- `label` (targeting / non_targeting)
- `log2_fc`, `p_value` (from `mdata.uns["trans_per_guide_results"]`)

**Non-targeting guides are excluded from this table** but their NT-vs-target-gene tests are used in the balanced evaluation step inside `intended_target.py` (see `04-module-internals.md`).

## `intended_target_metrics.tsv`

Single-row summary (no per-batch breakdown). Key columns:

| Column | Meaning |
|---|---|
| `n_guides_total` | All `targeting` guides in the dataset |
| `n_guides_tested` | Guides with a trans test for their intended target |
| `n_guides_kd_log2fc` | Guides with `log2_fc <= log2(0.4)` (i.e., ≥60% knockdown) |
| `n_guides_kd_pval` | Guides with `p_value < 0.05` |
| `n_guides_kd_both` | Guides meeting both criteria |
| `frac_*` | Corresponding fractions |
| `auroc` | AUROC for `1 - p_value` discriminating targeting vs non-targeting (balanced) |
| `auprc` | AUPRC for the same scoring |

**AUROC is the headline QC metric** — answers "do the guides knock down their target?" Expectations:

- Production datasets with a working library: AUROC > 0.85
- Benchmark datasets (small TF library, more NTC headroom): AUROC > 0.9
- AUROC < 0.7 = something's wrong (guide assignment, label mismatch, library swap)

## Plots

| Plot | What it shows | When to look |
|---|---|---|
| `gene_knee_plot.png` | log-rank UMI distribution | Cell-barcode filtering sanity check |
| `gene_histograms.png` | 3-panel: UMI, genes, mito with median lines | Per-dataset QC slide |
| `gene_histograms_by_batch.png` | Same 3 panels, density per batch | Identifying outlier lanes |
| `gene_cells_per_batch.png` | Bar chart of cells per batch | Lane-balance check |
| `guide_knee_plot.png` | log-rank guide UMI | Guide capture efficiency |
| `guide_histograms.png` | 3-panel: guide UMI, guides/cell, cells/guide (colored by label) | Whether NTC vs targeting capture diverges |
| `guide_histograms_by_batch.png` | 2-panel: guide UMI, guides/cell per batch | Lane consistency |
| `intended_target_volcano.png` | log2FC vs -log10(p), colored by label, top N labeled | Visual confirmation that targeting log2FCs cluster below 0 |
| `intended_target_log2fc_distribution.png` | Histogram of log2FC, threshold lines at 0 and log2(0.4) | Knockdown strength distribution |
| `intended_target_roc_pr_curves.png` | 2-panel: ROC + PR with AUROC/AUPRC annotated | Discrimination check |

The volcano and ROC/PR plots are the most useful for slides.

## Aggregation across runs

For sweep comparisons, the typical pattern is to pull `intended_target_metrics.tsv` AUROC across runs:

```bash
echo -e "dataset\trun\tauroc\tauprc\tfrac_kd_both" > sweep_summary.tsv
for f in datasets/*/*/qc/intended_target/*_intended_target_metrics.tsv; do
  ds=$(echo "$f" | cut -d/ -f2)
  run=$(echo "$f" | cut -d/ -f3)
  awk -F'\t' -v ds="$ds" -v run="$run" 'NR==2 {print ds"\t"run"\t"$X"\t"$Y"\t"$Z}' "$f"  # replace X/Y/Z with column indices
done >> sweep_summary.tsv
```

For the technology benchmark, the `technology-benchmark_WTC11_TF-Perturb-seq/` dataset directory contains cross-tech comparison notebooks that consume these TSVs.
