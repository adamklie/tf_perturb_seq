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
| WG1-B | `edistance_summary.tsv` | 🟡 partial (2/5 datasets) | n TFs called significant per dataset, under multiple thresholds incl. the calibration-robust `distance_mean > NC max` proxy. Wide format, slide-deck-friendly. | `reference/cross_dataset_edistance_summary.tsv` + per-dataset placeholders |
| WG1-C | `shared_tfs_upset.tsv` | 🟡 partial (need ≥3 datasets) | UpSet-ready table: per-TF boolean membership across datasets | WG1-B per-dataset significance lists |
| WG1-D | `tf_cross_lineage.tsv` | ✅ landed (2 datasets so far; widens as more land) | Per-TF, wide format. One row per `target_id`; per-dataset columns (`distance_mean_<ds>`, `distance_rank_<ds>`, `sig_dist_gt_NC_max_<ds>`) auto-widen as new datasets get a `wg1_significant_tfs.tsv`. Cross-dataset `classification` column (convergent_significant / `<ds>`-specific / discordant_partial / convergent_nonsignificant). | per-dataset `wg1_significant_tfs.tsv` (joined by `target_id`) |

> **⚠ Calibration caveat for everything ED-based**: the Huangfu DE/ESC `pval_mean` values are anti-conservative (see [`../../issues/edistance-calibration/`](../../issues/edistance-calibration/)). All WG1 derivatives above use `distance_mean > NC max` as the calibration-robust significance proxy. Raw `pval_mean<0.05` counts are included but flagged.

### WG1-D first snapshot (2 datasets: Huangfu DE × Huangfu ESC, 2026-05-11)

| classification | n_TFs |
|---|---:|
| convergent_significant (sig in both lineages) | **3** |
| HuangfuDE-specific | 82 |
| HuangfuESC-specific | 80 |
| convergent_nonsignificant | 2132 |

The 3 convergent TFs are **TERF2** (Myb/SANT), **GTF2B**, **ZNF574** (C2H2 ZF). Out of 73 + 83 = 156 per-dataset significant calls, only 3 overlap — high cross-lineage discordance worth flagging in the WG1 calibration-debug discussion. Refresh this table once Hon CM ED + Gersbach Hep ED are wired in; convergent set may grow.

## Per-dataset companions (under `datasets/<dataset>/<analysis>/`)

| Artifact | Path pattern | Status |
|---|---|---|
| WG1-B detail Per-TF significance + ranking | `datasets/<dataset>/energy_distance/wg1_significant_tfs.tsv` | ✅ ready for Huangfu DE + ESC (2273 rows × 16 cols each; joined with TF metadata) |
| WG1-E Trans-target counts (per perturbation) | `datasets/<dataset>/crispr_pipeline/wg1_trans_target_counts.tsv` | ✅ landed for Huangfu DE (1,741 perturbations, median 5 sig trans targets, max 4,361 = SOX17) + ESC (1,452 perturbations, median 2 sig trans, max 938 = STRAP). Hon CM blocked on CRISPR pipeline mirror; Gersbach/Engreitz blocked upstream. |

### WG1-E first snapshot (2 datasets: Huangfu DE × Huangfu ESC, 2026-05-11)

Per-perturbation count of significant trans-target genes at FDR<0.05 (per-TF BH; see WG4-A for the underlying edge list). Top-of-the-distribution TFs match the expected lineage masters:

- **Huangfu DE**: SOX17 (4,361 sig trans targets — DE master), FOXH1 (3,862), SETDB1 (2,405), SMARCC1 (1,677), ARID1A (1,305), SOX11 (1,084), DBX1 (951), SOX4 (783), SMAD3 (765)
- **Huangfu ESC**: STRAP (938), SETDB1 (712), RCOR2 (508), SALL4 (461), GRHL2 (365), POU5F1 / OCT4 (270), DNMT1 (186), KAT2A (198)

41 TFs in DE and 15 in ESC drive >100 significant trans targets — these are the candidate "master regulator" perturbations for lineage-specific deep-dives. Tail is long: median is 5 (DE) / 2 (ESC) trans targets per perturbation, so most TFs have small trans footprints.
