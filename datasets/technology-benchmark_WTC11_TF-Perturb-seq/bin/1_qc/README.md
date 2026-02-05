# Cross-Technology QC Comparison

This directory contains notebooks for comparing QC metrics across the four TF-Perturb-seq technology benchmarks.

## Input Files

All notebooks read from two manifest files in the parent directory:

- **`latest_qc_paths.tsv`**: Paths to QC output files for each dataset
- **`latest_mudata_paths.tsv`**: Paths to MuData files (if needed)

To update the analysis with new runs, simply update these manifest files and re-run the notebooks.

## Output Directory

All outputs are saved to: `results/cross_tech_comparison/`

---

## Notebooks

### 1. `1_metrics_comparison.ipynb` - QC Metrics Comparison

Compares basic gene expression and guide capture metrics across datasets.

#### Plots Generated:

**`gene_metrics_comparison.pdf`**
- 4-panel barplot comparing aggregate gene expression metrics across datasets
- **Number of Cells**: Total cells passing QC in each dataset
- **Median UMI per Cell**: Median total UMI counts per cell (sequencing depth)
- **Median Genes per Cell**: Median number of genes detected per cell (library complexity)
- **Median Mito %**: Median percentage of reads mapping to mitochondrial genes (cell quality indicator; higher = more stressed/dying cells)

**`guide_metrics_comparison.pdf`**
- 4-panel barplot comparing aggregate guide capture metrics
- **Median Guide UMI**: Median total guide UMI counts per cell
- **Median Guides/Cell**: Median number of distinct guides detected per cell
- **Fraction Cells with Guide**: Proportion of cells with at least one guide assigned
- **Median Cells/Guide**: Median number of cells per guide (coverage)

**`metrics_by_lane_comparison.pdf`**
- 2x3 grid showing metrics broken down by individual sequencing lanes/batches
- Each bar represents one lane, colored by dataset
- Value labels on each bar
- Panels:
  - **Singlet scRNA > 0 cells**: Number of cells per lane
  - **Median Tx UMIs/cell**: Median transcriptome UMI per cell
  - **Median Tx Expressed Genes/cell**: Median genes detected per cell
  - **Median Tx % MT UMIs**: Median mitochondrial percentage
  - **Median sgRNA UMIs/cell**: Median guide UMI per cell
  - **Mean assigned sgRNA/cell**: Mean number of guides assigned per cell

#### Intermediate Files:
- `combined_gene_metrics.tsv`: All gene metrics for all datasets/batches
- `combined_guide_metrics.tsv`: All guide metrics for all datasets/batches
- `metrics_summary.tsv`: Summary table of key metrics per dataset
- `metrics_by_lane.tsv`: Per-lane metrics for all datasets

---

### 2. `2_repression_efficiency_comparison.ipynb` - Repression Efficiency

Compares how well guides knock down their intended targets across technologies.

#### Plots Generated:

**`auroc_auprc_comparison.pdf`**
- 2-panel barplot of classification performance metrics
- **auROC** (Area Under ROC Curve): Measures ability to distinguish positive controls from negative controls based on log2FC. 0.5 = random, 1.0 = perfect separation.
- **auPRC** (Area Under Precision-Recall Curve): Similar to auROC but more sensitive to class imbalance. Higher = better.
- *Calculation*: Guides are ranked by log2FC. Positive controls should have more negative log2FC than negative controls. auROC/auPRC quantify this separation.

**`repression_efficiency_boxplot.pdf`**
- Simple boxplot of fold change for positive control guides
- Y-axis: Fold Change (2^log2FC) - values < 1 indicate knockdown
- Median values annotated above each box
- Red line shows median within each box
- *Interpretation*: Lower FC = stronger knockdown. Median around 0.3-0.4 means ~60-70% reduction.

**`log2fc_violin_by_label.pdf`**
- Split violin plot comparing log2FC distributions
- Compares positive_control vs tf_targeting guides
- *Interpretation*: Positive controls should show strong negative log2FC (left shift). TF-targeting guides show variable effects.

**`fc_violin_positive_controls.pdf`**
- Violin plot of fold change for positive controls only
- White dots show medians
- Red dashed line at FC=1 (no effect)

**`log2fc_pairwise_scatter.pdf`**
- Lower-triangular scatter plot matrix
- Each panel: log2FC in dataset X (x-axis) vs dataset Y (y-axis)
- Points colored by Pearson correlation (Reds colormap)
- Diagonal dashed line = perfect agreement
- r value annotated in each panel
- *Interpretation*: Points along diagonal = consistent effects. Scatter = technology-specific variation.

**`log2fc_correlation_matrix.pdf`**
- Lower-triangular scatter matrix with colored dataset labels
- Column labels (top) and row labels (left) show dataset names with colored backgrounds
- Each scatter panel shows guide log2FC correlation between two datasets
- Colorbar shows Pearson r mapping
- *Interpretation*: High correlation (darker red) = technologies agree on guide effects.

#### Intermediate Files:
- `combined_intended_target_metrics.tsv`: auROC, auPRC, and summary stats per dataset
- `combined_intended_target_results.tsv`: Per-guide log2FC and p-values for all datasets

---

### 3. `3_guide_capture_comparison.ipynb` - Guide Capture

Detailed comparison of guide capture efficiency metrics.

#### Plots Generated:

**`cells_per_guide_comparison.pdf`**
- 2-panel barplot
- **Median Cells per Guide**: Median number of cells each guide is detected in
- **Mean Cells per Guide (+/- SE)**: Mean with standard error bars
- *Interpretation*: Higher = better guide representation. Low values may indicate capture bias.

**`guides_per_cell_comparison.pdf`**
- 2-panel barplot of guides per cell statistics
- *Interpretation*: ~1-2 guides/cell is typical. Higher values may indicate multiplet issues or MOI differences.

**`guide_umi_comparison.pdf`**
- 2-panel barplot of guide UMI counts
- *Interpretation*: Higher UMI = stronger guide signal, easier assignment.

**`guide_assignment_comparison.pdf`**
- 2-panel barplot
- **Guide Assignment Rate**: % of cells with at least one guide
- **Single-Guide Cell Rate**: % of cells with exactly one guide
- *Interpretation*: High assignment + high single-guide = clean data.

#### Intermediate Files:
- `guide_capture_summary.tsv`: Summary of all guide capture metrics

---

### 4. `4_deg_comparison.ipynb` - DEG Analysis

Compares genome-wide differential expression (trans effects) across technologies.

#### Plots Generated:

**`trans_metrics_comparison.pdf`**
- 3-panel barplot of trans effect summary statistics
- **Median DEGs/Guide (Targeting)**: Median number of significant DEGs per TF-targeting guide
- **Mean DEGs/Guide (Targeting)**: Mean DEGs per guide
- **Total Significant Tests**: Total number of significant guide-gene associations
- *Calculation*: For each guide, count genes with FDR-adjusted p-value < 0.05.

**`degs_per_guide_violin.pdf`**
- Violin plot of DEG counts per TF-targeting guide
- *Interpretation*: Wide distribution = variable guide activity. High median = strong trans effects.

**`posctrl_degs_barplot.pdf`**
- Grouped barplot of DEG counts for each positive control guide
- Groups by guide, colors by dataset
- *Interpretation*: Positive controls should have consistent DEG patterns across technologies.

**`posctrl_deg_effect_sizes.pdf`**
- Violin plot of log2FC for significant positive control DEGs
- *Interpretation*: Shows the magnitude of trans effects for positive controls.

**`degs_normalized_comparison.pdf`**
- 2-panel barplot
- **Total DEGs**: Raw count of significant DEGs (TF-targeting guides)
- **DEGs per 1000 Cells**: Normalized by cell count to account for dataset size differences
- *Interpretation*: Normalization helps compare datasets with different cell numbers.

#### Intermediate Files:
- `combined_trans_metrics.tsv`: Summary trans statistics per dataset
- `combined_trans_results.tsv.gz`: All per-guide, per-gene trans test results (gzipped due to size)

---

### 5. `5_guide_level_comparison.ipynb` - Guide-Level Comparison

Systematic comparison of individual guide performance across technologies, organized into two phases:

#### Phase 1: Capture Analysis (requires per-guide capture files)

Compares per-guide UMI counts and cell detection rates across technologies. Requires running `extract_per_guide_capture.py` on MuData files first.

**`capture_cells_correlation_matrix.pdf`**
- Lower-triangular scatter matrix with histograms on diagonal
- Each panel: Fraction of cells detecting guide in dataset X vs Y
- Pearson correlation (r) shown in each scatter panel
- *Interpretation*: High r = technologies detect guides in similar proportions of cells.

**`capture_umi_correlation_matrix.pdf`**
- Same format as cells correlation matrix but for mean UMI per guide
- *Interpretation*: High r = technologies capture guides with similar UMI intensities.

**`capture_variability_distribution.pdf`**
- 2-panel figure
- **Left**: Histogram of CV (coefficient of variation) of cell detection across technologies
- **Right**: Scatter of mean capture rate vs CV, colored by concordance category
- *Interpretation*: Low CV = consistent capture across platforms. High CV = technology-specific capture bias.

**`capture_discordant_guides_heatmap.pdf`**
- Heatmap showing top 30 most discordant guides
- Values: Fraction of cells detecting each guide
- *Interpretation*: Shows guides with largest technology-specific capture differences.

**`capture_tech_bias_zscore.pdf`**
- Box plot of z-scores for each technology
- Z-score = deviation from mean capture rate across technologies
- Red dashed lines at z = ±1.5 mark high/low performer thresholds
- *Interpretation*: Technologies with many outliers have technology-specific capture biases.

#### Phase 1 Intermediate Files:
- `guide_cell_detection_matrix.tsv`: Matrix of cell detection fractions with concordance labels
- `guide_umi_matrix.tsv`: Matrix of mean UMI per guide across datasets
- `capture_concordance_summary.tsv`: Summary statistics by concordance category

---

#### Phase 2: Efficacy Analysis

Compares intended target knockdown (log2FC) across technologies.

**`guide_variability_distribution.pdf`**
- 2-panel figure
- **Left**: Histogram of standard deviation of log2FC across datasets per guide
- **Right**: Scatter of mean log2FC vs std dev, colored by guide type
- *Interpretation*: High std = technology-specific performance. Low std = consistent across platforms.

**`guide_pairwise_diff_heatmap.pdf`**
- Heatmap showing top 30 most variable guides
- Columns: Pairwise dataset comparisons (e.g., "Hon_vs_Huangfu")
- Values: Difference in log2FC between datasets
- Blue = better in first dataset, Red = better in second
- *Interpretation*: Identifies guides with largest technology-specific effects.

**`trans_effects_correlation.pdf`** (if trans data available)
- Correlation heatmap of # significant trans DEGs per guide
- *Interpretation*: High correlation = guides have similar trans activity across platforms.

#### Phase 2 Intermediate Files:
- `guide_log2fc_matrix.tsv`: Matrix of log2FC values (guides x datasets) with summary stats
- `guide_detection_matrix.tsv`: Binary matrix showing which guides detected in which datasets

#### Key Analyses:

1. **Capture Concordance**: Classify guides as Concordant, Moderate, or Discordant based on CV of cell detection
2. **Technology-Specific Bias**: Identify guides where specific technologies over/under-perform (z-score analysis)
3. **Guide Coverage**: How many datasets detect each guide (NaN = not detected)
4. **Variability Analysis**: Standard deviation of log2FC across datasets
5. **Exclusive Performers**: Guides that work well in one technology but not others
6. **Missing Guides**: Guides with NaN values indicating poor capture

---

## Supporting Scripts

### `extract_per_guide_capture.py`

Extracts per-guide capture metrics from MuData files for use in Phase 1 of notebook 5.

**Usage:**
```bash
python extract_per_guide_capture.py <mudata_path> <output_path> <dataset_name>
```

**Example:**
```bash
python extract_per_guide_capture.py \
    /path/to/inference_mudata.h5mu \
    /path/to/output/per_guide_capture.tsv \
    Hon_WTC11-benchmark_TF-Perturb-seq
```

**Output columns:**
- `guide_id`: Guide identifier
- `gene_name`: Target gene (if available in MuData)
- `label`: Guide label (positive_control, tf_targeting, etc.)
- `n_cells_detected`: Number of cells with UMI > 0 for this guide
- `frac_cells_detected`: Fraction of total cells with this guide
- `total_umi`: Sum of all UMIs for this guide across cells
- `mean_umi`: Mean UMI among cells with guide detected
- `median_umi`: Median UMI among cells with guide detected
- `std_umi`: Standard deviation of UMI
- `max_umi`: Maximum UMI observed
- `n_cells_assigned`: Number of cells assigned to this guide (if assignment available)
- `total_cells`: Total cells in dataset
- `dataset`: Dataset identifier

### `run_extract_per_guide_capture.sh`

Shell script to run extraction for all datasets listed in `latest_mudata_paths.tsv`.

**Usage:**
```bash
bash run_extract_per_guide_capture.sh
```

**Note:** Run this on the cluster where MuData files are accessible. Output files should be placed in each dataset's QC directory as `per_guide_capture.tsv`.

---

## Color Scheme

All plots use a consistent color scheme defined in `config/colors/technology-benchmark_WTC11_TF-Perturb-seq.yaml`:

| Dataset | Color | Hex |
|---------|-------|-----|
| Hon | Blue | #5B8FF9 |
| Huangfu | Yellow/Gold | #F6BD16 |
| Gersbach GEM-Xv3 | Mint Green | #5AD8A6 |
| Gersbach HTv2 | Coral/Red | #E86452 |

---

## How to Update

1. Run new QC pipelines on your datasets
2. Update `latest_qc_paths.tsv` with new paths:
   ```
   dataset    qc_dir    gene_metrics    guide_metrics    intended_target_results    intended_target_metrics    trans_results    trans_metrics
   ```
3. (Optional) For Phase 1 capture analysis in notebook 5:
   - Update `latest_mudata_paths.tsv` with paths to MuData files
   - Run `extract_per_guide_capture.py` for each dataset
   - Place output as `per_guide_capture.tsv` in each dataset's qc_dir
4. Re-run notebooks in order (1 → 5)
5. New plots and intermediate files will be regenerated in `results/cross_tech_comparison/`

---

## Metric Definitions

### Gene Expression Metrics
- **UMI**: Unique Molecular Identifier - corrected count of RNA molecules
- **Genes detected**: Number of genes with at least 1 UMI
- **Mito %**: Percentage of UMIs from mitochondrial genes (high = cell stress)

### Guide Capture Metrics
- **Guide UMI**: UMI counts for guide RNA sequences
- **Guides/cell**: Number of distinct guides detected per cell
- **Cells/guide**: Number of cells each guide is detected in

### Repression Metrics
- **log2FC**: log2(perturbed / control) - negative = knockdown
- **Fold Change**: 2^log2FC - values < 1 = knockdown
- **auROC**: Area under ROC curve for classifying positive vs negative controls
- **auPRC**: Area under Precision-Recall curve

### Trans Effect Metrics
- **DEG**: Differentially Expressed Gene (FDR < 0.05)
- **Trans effect**: Gene expression change for non-target genes
- **Cis effect**: Expression change for the targeted gene (intended target)
