# Inference Module Plan

## Purpose

Calibrate PerTurbo perturbation inference results using an empirical null to produce reliable per-element effect estimates with well-calibrated p-values.

## Problem Statement

The IGVF CRISPR pipeline (PerTurbo) tests every perturbed element against every detected gene, producing per-element results with:
- `log2_fc`: log2 fold change of gene expression in perturbed vs control cells
- `p_value`: posterior probability from a negative binomial model

These posterior p-values are **not calibrated for frequentist hit calling**. We need to derive empirical p-values by comparing each test statistic against the distribution of test statistics from non-targeting control (NTC) elements.

## Calibration Procedure

1. **Build empirical null** from NTC element results
   - **eCDF** (default): Rank each targeting element's test statistic against the NTC distribution directly
   - **t-fit** (alternative): Convert to z-values, fit a t-distribution to NTC z-values, compute p-values from the fitted distribution
   - Reference: [PerTurbo utils](https://github.com/pinellolab/PerTurbo/blob/main/src/perturbo/utils/utils.py)
2. **Compute empirical p-values** for each targeting element × gene test
3. **Apply BH FDR correction** across all tests in the discovery set

Downstream consumers apply their own FDR cutoffs.

## Inputs

### Pipeline result TSVs (`pipeline_outputs/`)

Per-element and per-guide results from GCP. We use the **per-element TSVs** (not the MuData `.uns` tables, which are per-guide only):

| File | Description |
|---|---|
| `perturbo_trans_per_element_output.tsv.gz` | All trans tests (element × gene pairs) |
| `perturbo_cis_per_element_output.tsv.gz` | Cis-window tests (element × nearby gene pairs) |
| `perturbo_trans_per_guide_output.tsv.gz` | Per-guide trans tests (for reference) |
| `perturbo_cis_per_guide_output.tsv.gz` | Per-guide cis tests (for reference) |
| `inference_mudata.h5mu` | Guide metadata, assignment matrix, gene annotations |

### TF enrichment tables (`tf/benchmark_output/benchmark_tables/`)

| File | Description |
|---|---|
| `enrichment_all.tsv` | ChIP-seq TF enrichment results |
| `tf_order.tsv` | TF ordering for visualization |
| `tf_peak_mapping.tsv` | TF-to-peak mapping |
| `trans_per_element_results_used.tsv` | Trans results used in TF benchmark |

### Derived from MuData

| Field | Source | Description |
|---|---|---|
| `n_cells` per element | `guide.layers["guide_assignment"]` or `guide.var` | Number of cells assigned to each element |
| `intended_target_name` | `guide.var` | Maps element → intended target gene |
| `gene_name` | `guide.var` | Gene symbol of intended target |
| `label` | `guide.var` | Element type: targeting, non_targeting, positive_control, etc. |
| Gene coordinates | `gene.var` | Chromosome, start, end — for cis/trans classification |
| Element coordinates | `guide.var` | Chromosome, start, end — for cis/trans classification |

## Output Tables

All outputs are **per-element level** and **unfiltered** (every test is included). No significance thresholds are applied — consumers choose their own cutoffs.

### a. Calibrated trans results (`{prefix}_calibrated_trans_results.tsv`)

The master table. Every element × gene pair from the trans results, fully annotated.

| Column | Type | Description |
|---|---|---|
| `element_id` | str | Perturbed element identifier (from pipeline) |
| `element_symbol` | str | Gene symbol of the element's intended target |
| `element_label` | str | Element category: `targeting`, `non_targeting`, `positive_control`, etc. |
| `tested_gene_id` | str | Ensembl gene ID of the tested gene |
| `tested_gene_symbol` | str | Gene symbol of the tested gene |
| `n_cells` | int | Number of cells assigned to this element |
| `log2fc` | float | Log2 fold change (negative = downregulation) |
| `log2fc_se` | float | Standard error on log2FC (NaN if unavailable from pipeline) |
| `is_cis` | bool | Element–gene pair within cis window (default: 500 kb) |
| `is_direct_target` | bool | Tested gene is the element's intended target |
| `posterior_pval` | float | Raw posterior probability from PerTurbo |
| `empirical_pval` | float | Calibrated p-value from NTC null distribution |
| `empirical_pval_adj` | float | BH-corrected empirical p-value |

### b. Direct target results (`{prefix}_calibrated_direct_target_results.tsv`)

Subset of (a) where `is_direct_target == True`. Same schema.

### c. Cis results (`{prefix}_calibrated_cis_results.tsv`)

Subset of (a) where `is_cis == True`. Includes direct targets and other cis-window genes. Same schema.

### d. Trans-only results (`{prefix}_calibrated_trans_only_results.tsv`)

Subset of (a) where `is_cis == False`. Same schema.

### e. Pathway enrichment results (`{prefix}_calibrated_pathway_results.tsv`)

Per-element pathway overrepresentation analysis. The consumer provides an FDR cutoff to define the input DEG set from table (a).

| Column | Type | Description |
|---|---|---|
| `element_id` | str | Perturbed element |
| `element_symbol` | str | Gene symbol of the element's intended target |
| `pathway_id` | str | Pathway identifier (GO:XXXXXXX, KEGG:XXXXX, etc.) |
| `pathway_name` | str | Human-readable pathway name |
| `pathway_source` | str | Gene set collection (GO_BP, KEGG, HALLMARK, etc.) |
| `n_degs_in_pathway` | int | DEGs from this element that overlap with this pathway |
| `n_degs_total` | int | Total DEGs for this element (at the consumer's FDR cutoff) |
| `pathway_size` | int | Total genes annotated to this pathway |
| `background_size` | int | Total genes in the background set (all tested genes) |
| `fold_enrichment` | float | Observed/expected ratio |
| `pvalue` | float | Fisher's exact test p-value |
| `pvalue_adj` | float | BH-corrected p-value (across all pathways tested for this element) |

## Implementation Plan

### Phase 1: Data sync update

**File:** `datasets/technology-benchmark_WTC11_TF-Perturb-seq/bin/1_get_data/scripts/sync_benchmark.sh`

Currently syncs `pipeline_dashboard/` or `inference_mudata.h5mu` only. Update to also download:
- `pipeline_outputs/perturbo_*_per_element_output.tsv.gz` — the per-element result TSVs
- `pipeline_outputs/perturbo_*_per_guide_output.tsv.gz` — per-guide results (for reference)
- `tf/benchmark_output/benchmark_tables/` — TF enrichment tables

### Phase 2: Calibration module

**File:** `src/tf_perturb_seq/inference/calibrate.py`

Standalone CLI script. Core calibration logic.

```bash
python calibrate.py \
    --trans-results <perturbo_trans_per_element_output.tsv.gz> \
    --cis-results <perturbo_cis_per_element_output.tsv.gz> \
    --mudata <inference_mudata.h5mu> \
    --outdir <output_dir> \
    --prefix <dataset_name> \
    [--null-method ecdf|t-fit]            # default: t-fit
    [--cis-window 100000]                 # default: 100 kb
    [--non-targeting-label non_targeting]
```

**Internal steps:**
1. Load trans per-element results TSV
2. Load guide metadata + gene coordinates from MuData
3. Identify NTC element results → build empirical null (eCDF or t-fit)
4. Compute empirical p-values for all targeting element tests
5. Apply BH correction across the discovery set
6. Annotate: `is_cis` (from cis results or coordinate-based), `is_direct_target`, `n_cells`, gene symbols
7. Write output tables (a)–(d)

### Phase 3: Pathway analysis module

**File:** `src/tf_perturb_seq/inference/pathways.py`

Lightweight overrepresentation analysis. Dependencies: `scipy.stats` (Fisher's exact / hypergeometric). No heavy frameworks.

```bash
python pathways.py \
    --calibrated-results <calibrated_trans_results.tsv> \
    --gene-sets <pathway.gmt> \
    --outdir <output_dir> \
    --prefix <dataset_name> \
    [--fdr-threshold 0.10]       # FDR cutoff to define DEGs from calibrated results
    [--min-genes 5]              # min DEGs in pathway to test
    [--max-pathway-size 500]     # skip very broad pathways
```

**Internal steps:**
1. Load calibrated results, filter to DEGs at given FDR threshold per element
2. Load gene set collection (GMT format — pre-downloaded MSigDB or GO)
3. For each element with ≥1 DEG: run Fisher's exact test per pathway
4. BH correct across all pathways tested per element
5. Write output table (e)

### Phase 4: Cross-dataset comparison notebooks

**Location:** `datasets/technology-benchmark_WTC11_TF-Perturb-seq/bin/`

Separate step — will be planned after Phases 1–3 are implemented.

## Principles

- Dump source data tables alongside every plot
- Prefer information-rich visualizations (violin over bar, scatter over summary)
- All tables unfiltered; thresholds applied downstream
- Document issues requiring upstream fixes

## Existing Code

| File | Status | Notes |
|---|---|---|
| `trans.py` | Move to /cellar/users/aklie/data/datasets/tf_perturb_seq/scratch/2026_03_18 | Per-guide trans QC using raw posterior p-values + BH. Useful for dataset-level QC summaries. Superseded by `calibrate.py` for hit calling. |
| `mudata_dump.py` | Move to /cellar/users/aklie/data/datasets/tf_perturb_seq/scratch/2026_03_18 | Exports MuData to TSV/MTX. Not used directly by calibration (we read pipeline TSVs), but useful for ad-hoc data extraction. |
