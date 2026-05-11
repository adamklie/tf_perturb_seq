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

### WG1-D snapshot (3 datasets: Hon CM × Huangfu DE × Huangfu ESC, 2026-05-11)

| classification | n_TFs |
|---|---:|
| **convergent_significant** (sig in all 3 lineages) | **3** |
| discordant_partial (sig in 2 of 3 lineages) | 13 |
| HonCM-specific | 154 |
| HuangfuDE-specific | 123 |
| HuangfuESC-specific | 68 |
| convergent_nonsignificant | 2032 |

**The 3 convergent_significant TFs survive the addition of Hon CM**: TERF2, GTF2B, ZNF574 are sig in *all three* lineages. Absolute distances vary wildly between lineages — e.g. TERF2 is 6.6 in Hon CM vs 638 in Huangfu ESC — but each clears its own dataset's NC_max threshold. These are the candidate "always wired" TFs across the perturb-seq library.

**Discordant_partial pattern is striking**: 10 of the 13 discordant TFs are sig in Hon CM + Huangfu ESC but *not* Huangfu DE (NOC3L, RPF1, ZMAT2, ATF5, SRF, TAF1A, TAF11, SALL4, …). Huangfu DE is the selective lineage; Hon CM and Huangfu ESC respond to a broader set of TF perturbations. Worth digging into during the calibration-debug discussion (Issue 1).

Hon CM dominates the lineage-specific list (154) — also has the highest per-lineage sig count overall (164 vs 73 DE / 83 ESC). Refresh again once Gersbach Hep ED + Engreitz Endo ED land; convergent set may shrink.

## Per-dataset companions (under `datasets/<dataset>/<analysis>/`)

| Artifact | Path pattern | Status |
|---|---|---|
| WG1-B detail Per-TF significance + ranking | `datasets/<dataset>/energy_distance/wg1_significant_tfs.tsv` | ✅ ready for Hon CM (2,036 rows; 164 sig) + Huangfu DE (2,273 rows; 73 sig) + Huangfu ESC (2,273 rows; 83 sig). Joined with TF metadata. |
| WG1-E Trans-target counts (per perturbation) | `datasets/<dataset>/crispr_pipeline/wg1_trans_target_counts.tsv` | ✅ landed for Huangfu DE (1,741 perturbations, median 5 sig trans targets, max 4,361 = SOX17) + ESC (1,452 perturbations, median 2 sig trans, max 938 = STRAP). Hon CM blocked on CRISPR pipeline mirror; Gersbach/Engreitz blocked upstream. |

### WG1-E first snapshot (2 datasets: Huangfu DE × Huangfu ESC, 2026-05-11)

Per-perturbation count of significant trans-target genes at FDR<0.05 (per-TF BH; see WG4-A for the underlying edge list). Top-of-the-distribution TFs match the expected lineage masters:

- **Huangfu DE**: SOX17 (4,361 sig trans targets — DE master), FOXH1 (3,862), SETDB1 (2,405), SMARCC1 (1,677), ARID1A (1,305), SOX11 (1,084), DBX1 (951), SOX4 (783), SMAD3 (765)
- **Huangfu ESC**: STRAP (938), SETDB1 (712), RCOR2 (508), SALL4 (461), GRHL2 (365), POU5F1 / OCT4 (270), DNMT1 (186), KAT2A (198)

41 TFs in DE and 15 in ESC drive >100 significant trans targets — these are the candidate "master regulator" perturbations for lineage-specific deep-dives. Tail is long: median is 5 (DE) / 2 (ESC) trans targets per perturbation, so most TFs have small trans footprints.
