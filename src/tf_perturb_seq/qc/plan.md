# QC Module Plan

## Purpose

Take the output of the IGVF CRISPR pipeline (an `inference_mudata.h5mu` file) and generate per-dataset QC summary tables and plots. These outputs are consumed both as standalone per-dataset reports and as inputs to cross-dataset comparison notebooks (e.g. the [technology benchmark](../../../datasets/technology-benchmark_WTC11_TF-Perturb-seq/)).

## Architecture

```
run_qc_pipeline.sh           ← orchestrator (args: --project-root, --input, --outdir, --run-name)
  ├── mapping_gene.py        ← Step 1: Gene expression mapping QC
  ├── mapping_guide.py       ← Step 2: Guide mapping QC
  ├── intended_target.py     ← Step 3: Intended target knockdown QC for guides
```

Each module is a standalone CLI script (`python <module>.py --input ... --outdir ... --prefix ...`). The orchestrator calls them sequentially and prefixes all outputs with `{run-name}_{module}`.

## Input

- `inference_mudata.h5mu` — MuData with at least:
  - `gene` modality: `.obs` has `batch`, `total_gene_umis`, `percent_mito`, and one of `num_expressed_genes` or `log1p_n_genes_by_counts`
  - `guide` modality: `.obs` has `batch`, `total_guide_umis`; `.var` has `label`, `gene_name`; `.layers["guide_assignment"]` is the binary assignment matrix
  - `mdata.uns["trans_per_guide_results"]` — DataFrame with `guide_id`, `gene_id`, `log2_fc`, `p_value`

## Output Directory Layout

```
{outdir}/
  mapping_gene/
    {prefix}_gene_metrics.tsv           ← overall + per-batch summary (n_cells, UMI stats, genes stats, mito stats)
    {prefix}_gene_knee_plot.png
    {prefix}_gene_histograms.png        ← 3-panel: UMI, genes, mito
    {prefix}_gene_histograms_by_batch.png
    {prefix}_gene_cells_per_batch.png
  mapping_guide/
    {prefix}_guide_metrics.tsv          ← overall + per-batch summary (guide UMI stats, guides/cell, cells/guide)
    {prefix}_guide_per_guide_capture.tsv ← per-guide: n_cells_detected, total_umi, mean/median/std/max UMI
    {prefix}_guide_knee_plot.png
    {prefix}_guide_histograms.png       ← 3-panel: guide UMI, guides/cell, cells/guide (colored by label)
    {prefix}_guide_histograms_by_batch.png ← 2-panel: guide UMI, guides/cell by batch
  intended_target/
    {prefix}_intended_target_results.tsv   ← per-guide intended target tests (guide_id, gene_name, label, log2_fc, p_value)
    {prefix}_intended_target_metrics.tsv   ← summary: n_guides, knockdown rates, AUROC, AUPRC
    {prefix}_intended_target_volcano.png
    {prefix}_intended_target_log2fc_distribution.png
    {prefix}_intended_target_roc_pr_curves.png
```

## Module Details

### 1. mapping_gene.py

**Goal:** Summarize gene expression quality per dataset and per batch (lane).

**Metrics (per batch + overall):**
| Column group | Columns | Source |
|---|---|---|
| Cell count | `n_cells` | len(obs) |
| UMI | `umi_{median,mean,std,min,max,q25,q75}` | `obs["total_gene_umis"]` |
| Genes | `genes_{median,mean,std,min,max,q25,q75}` | `obs["num_expressed_genes"]` |
| Mito | `mito_{median,mean,std,min,max,q25,q75}` | `obs["percent_mito"]` |

**Plots:**
- Knee plot — ranked UMI counts (log scale)
- Histograms — 3-panel (UMI, genes, mito) with median lines
- Histograms by batch — same 3 panels, density per batch (step element, colored)
- Cells per batch — bar chart

**Fallback:** If `num_expressed_genes` is missing, derives it from `np.expm1(obs["log1p_n_genes_by_counts"])`. This was needed because the GCP pipeline doesn't always produce `num_expressed_genes`.

### 2. mapping_guide.py

**Goal:** Summarize guide capture quality per dataset and per batch.

**Pre-computation:** `compute_guide_assignment_counts()` — binarizes `layers["guide_assignment"]` to compute `n_guides_per_cell` (obs) and `n_cells_per_guide` (var).

**Metrics (per batch + overall):**
| Column group | Columns | Source |
|---|---|---|
| Cell count | `n_cells` | len(obs) |
| Guide UMI | `guide_umi_{median,mean,std,min,max,q25,q75}` | `obs["total_guide_umis"]` |
| Guides/cell | `guides_per_cell_{mean,std,min,max,median}` | computed `n_guides_per_cell` |
| Assignment | `n_cells_with_guide`, `n_cells_exactly_1_guide`, `frac_cells_with_guide` | computed |
| Cells/guide | `n_guides_total`, `cells_per_guide_{median,mean,std,min,max}` | computed (overall only) |

**Per-guide capture table** (`_per_guide_capture.tsv`):
For each guide in `.var`: `n_cells_detected`, `frac_cells_detected`, `total_umi`, `mean_umi`, `median_umi`, `std_umi`, `max_umi`, plus `label` and `gene_name` from `.var`.

**Plots:**
- Knee plot — ranked guide UMI counts
- Histograms — 3-panel: guide UMI, guides/cell, cells/guide (colored by label)
- Histograms by batch — 2-panel: guide UMI, guides/cell (density per batch)

**Bandaid (Engreitz):** If gene and guide modalities have zero batch ID overlap, copies gene batch labels onto guide obs. This handles the Engreitz dataset where RNA and guide libraries use different analysis set accessions.

### 3. intended_target.py

**Goal:** Evaluate how well guides knock down their intended target genes. Uses `trans_per_guide_results` (not cis) so that non-targeting guides are included as negative controls.

**Workflow:**
1. Build intended target map: `guide.var["intended_target_name"]` → `guide.var["gene_name"]`
2. Load `mdata.uns["trans_per_guide_results"]`
3. Filter to intended targets: keep rows where guide's target gene matches the tested gene
4. Compute knockdown metrics (thresholds: log2FC ≤ log2(0.4), p < 0.05)
5. Build balanced evaluation table: targeting guides (positives) vs non-targeting guides tested against the same genes (negatives), downsampled to match
6. Compute AUROC and AUPRC (score = 1 - p_value)

**Outputs:**
- `_results.tsv` — all intended target test rows (guide_id, gene_name, label, log2_fc, p_value)
- `_metrics.tsv` — summary: n_guides_total, n_guides_tested, knockdown fractions, AUROC, AUPRC

**Plots:**
- Volcano — log2FC vs -log10(p_value), colored by label, top N labeled
- Log2FC distribution — histogram colored by label, threshold lines at 0 and log2(0.4)
- ROC/PR curves — 2-panel with AUROC/AUPRC annotations

## SLURM Execution

For batch processing across datasets, see `bin/2_qc/scripts/qc_array.sh` which reads a manifest TSV (`<run-label>_qc_input.tsv`) and submits one SLURM array task per dataset (4 cores, 16G, 2h, `carter-compute` partition).
