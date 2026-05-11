# WG1 — Overarching data summarization, pipeline & QC cleanup

**Topic 1 / Figure 1.** Goal: establish the primary building blocks for all downstream analysis.

## Questions

From [`../../WORKING_GROUPS.md`](../../WORKING_GROUPS.md):

- Identify which guides are detected and show target repression in each system (cis inference outputs).
- Summarize and compare general statistics (cell counts, gRNA and scRNA MOI, etc.) and data quality metrics (%mito counts) across datasets.
- Extract energy-distance outputs and create a bar plot (or UpSet plot) of the number of TFs per production dataset that significantly alter the transcriptome.
- UpSet plot of shared TFs with significant effects across lineages.
- Cluster TF perturbations by their energy distance in each lineage; label cases where TFs have very different energy distances across lineages.
- Identify inferred downstream trans targets of each guide/TF in each system and look at the distribution of overlap.
- Develop and apply strategies to reduce the impact of technical differences between datasets.
- Discussion of ways to integrate clonal expansion, doublet removal, and DEG calibration into upstream analysis pipeline.

## Artifacts in this folder

| ID | File | Status | What it answers | Source data |
|---|---|---|---|---|
| WG1-A | `qc_summary.tsv` | ✅ landed | Per-dataset cell counts, library QC, knockdown stats — at-a-glance | `cross_dataset_pipeline_summary.tsv` + `experimental_metadata_simplified.tsv` |
| WG1-B | `edistance_significance.tsv` | 🟡 partial (2/5 datasets) | n TFs called significant per dataset, under multiple thresholds incl. the calibration-robust `distance_mean > NC max` proxy | per-dataset `pval_edist_full.csv` on Synapse |
| WG1-C | `shared_tfs_upset.tsv` | 🟡 partial (need ≥3 datasets) | UpSet-ready table: per-TF boolean membership across datasets | WG1-B per-dataset significance lists |
| WG1-D | `tf_distance_similarity_long.tsv` | 🟡 partial | Per-TF cross-lineage `distance_mean` comparison + divergence classification | per-dataset `pval_edist_full.csv` |

> **⚠ Calibration caveat for everything ED-based**: the Huangfu DE/ESC `pval_mean` values are anti-conservative (see [`../../issues/edistance-calibration.md`](../../issues/edistance-calibration.md)). All WG1 derivatives above use `distance_mean > NC max` as the calibration-robust significance proxy. Raw `pval_mean<0.05` counts are included but flagged.

## Per-dataset companions (under `datasets/<dataset>/<analysis>/`)

| Artifact | Path pattern | Status |
|---|---|---|
| WG1-E Trans-target counts (per perturbation) | `datasets/<dataset>/crispr_pipeline/wg1_trans_target_counts.tsv` | ✅ ready to build for the 3 datasets with canonical CRISPR bundles (Huangfu DE/ESC, Hon CM via syn74520421) |
