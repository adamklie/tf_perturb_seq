# WG1 — Overarching data summarization, pipeline and QC cleanup

Topic 1 / Figure 1. Goal: establish the primary building blocks for all downstream analysis.

Files under [`examples/`](examples/) are illustrative starting points, not finished deliverables.

## Questions

- Identify which guides are detected and show target repression in each system (cis inference outputs).
- Summarize and compare general statistics (cell counts, gRNA and scRNA MOI, etc.) and data quality metrics (%mito counts) across datasets.
- Extract energy-distance outputs and create a bar plot (or UpSet plot) of the number of TFs per production dataset that significantly alter the transcriptome.
- UpSet plot of shared TFs with significant effects across lineages.
- Cluster TF perturbations by their energy distance in each lineage; label cases where TFs have very different energy distances across lineages.
- Identify inferred downstream trans targets of each guide/TF in each system and look at the distribution of overlap.
- Develop and apply strategies to reduce the impact of technical differences between datasets.
- Discussion of ways to integrate clonal expansion, doublet removal, and DEG calibration into upstream analysis pipeline.

## Artifacts in `examples/`

| ID | File | Status | What it answers | Source data |
|---|---|---|---|---|
| WG1-A | [`examples/qc_summary.tsv`](examples/qc_summary.tsv) | ready | Per-dataset cell counts, library QC, knockdown stats — at-a-glance | `cross_dataset_pipeline_summary.tsv` + `experimental_metadata_simplified.tsv` |
| WG1-B | [`examples/edistance_summary.tsv`](examples/edistance_summary.tsv) | caveat (2/5 datasets) | n TFs called significant per dataset, under multiple thresholds incl. the calibration-robust `distance_mean > NC max` proxy. Wide format, slide-deck-friendly. | `reference/cross_dataset_edistance_summary.tsv` + per-dataset placeholders |
| WG1-C | `examples/shared_tfs_upset.tsv` | caveat (need ≥3 datasets) | UpSet-ready table: per-TF boolean membership across datasets | WG1-B per-dataset significance lists |
| WG1-D | [`examples/tf_cross_lineage.tsv`](examples/tf_cross_lineage.tsv) | ready (auto-widens as more datasets land) | Per-TF, wide format. One row per `target_id`; per-dataset columns auto-widen as new datasets get a `wg1_significant_tfs.tsv`. Cross-dataset `classification` column (convergent_significant / `<ds>`-specific / discordant_partial / convergent_nonsignificant). | per-dataset `wg1_significant_tfs.tsv` (joined by `target_id`) |

Status legend: `ready` (built and refreshable), `caveat` (partial coverage; gets stronger as datasets land), `blocked` (upstream data missing).

## Issues

- *[FILL IN issue link]*: Huangfu DE/ESC `pval_mean` values are anti-conservative. All WG1 derivatives use `distance_mean > NC max` as the calibration-robust significance proxy. Raw `pval_mean<0.05` counts are included but flagged.

### WG1-D snapshot (3 datasets: Hon CM × Huangfu DE × Huangfu ESC, 2026-05-11)

| classification | n_TFs |
|---|---:|
| convergent_significant (sig in all 3 lineages) | 3 |
| discordant_partial (sig in 2 of 3 lineages) | 13 |
| HonCM-specific | 154 |
| HuangfuDE-specific | 123 |
| HuangfuESC-specific | 68 |
| convergent_nonsignificant | 2032 |

The 3 convergent_significant TFs (TERF2, GTF2B, ZNF574) clear each dataset's NC_max threshold despite varying absolute distances across lineages.

Of 13 discordant_partial TFs, 10 are sig in Hon CM + Huangfu ESC but not Huangfu DE (NOC3L, RPF1, ZMAT2, ATF5, SRF, TAF1A, TAF11, SALL4, ...). Hon CM dominates the lineage-specific list (154) and the per-lineage sig count (164 vs 73 DE / 83 ESC). Refresh again once Gersbach Hep ED + Engreitz Endo ED land.

## Per-dataset companions (under `../../data/<dataset>/<analysis>/`)

| Artifact | Path pattern | Status |
|---|---|---|
| Per-TF significance + ranking (WG1-B detail) | `data/<dataset>/energy_distance/wg1_significant_tfs.tsv` | ready for Hon CM (2,036 rows; 164 sig) + Huangfu DE (2,273 rows; 73 sig) + Huangfu ESC (2,273 rows; 83 sig). Joined with TF metadata. |
| Trans-target counts per perturbation (WG1-E) | `data/<dataset>/crispr_pipeline/wg1_trans_target_counts.tsv` | ready for Hon CM (2,065 perts, median 49 sig trans), Huangfu DE (1,741 perts, median 5), Huangfu ESC (1,452 perts, median 2). Gersbach/Engreitz blocked upstream. |

### WG1-E snapshot (3 datasets, 2026-05-11)

Per-perturbation count of significant trans-target genes at FDR<0.05 (per-TF BH; see WG4-A for the underlying edge list). Top perturbations by trans-target count:

- Hon CM: ISL1 (2,769), SOX11, TADA2B, SOX4, TBX20 (1,859), ZNF787, MEF2C (1,743), RCOR2, HAND1 (1,638), CHAMP1.
- Huangfu DE: SOX17 (4,361), FOXH1 (3,862), SETDB1 (2,405), SMARCC1 (1,677), ARID1A (1,305), SOX11 (1,084), DBX1 (951), SOX4 (783), SMAD3 (765).
- Huangfu ESC: STRAP (938), SETDB1 (712), RCOR2 (508), SALL4 (461), GRHL2 (365), POU5F1 (270), DNMT1 (186), KAT2A (198).

TFs with >100 sig trans targets: 537 in Hon CM (vs 41 DE / 15 ESC) — driven partly by Hon CM's denser trans-effect signal (median 49 vs 5/2). Whether the density reflects biology or technical factors (newer `seqspec_v3` pipeline, WTC11 vs HUES8) is open — see the harmonization caveat in WG4-A.

## Run the examples

[`examples/examples.py`](examples/examples.py) loads the artifacts above and prints quick views:

- §1 — QC summary at-a-glance (cell counts, UMI medians, AUROC).
- §2 — Energy-distance significance per dataset (calibration-robust + raw p-value).
- §3 — Cross-lineage classification counts + the convergent_significant TF list.
- §4 — UpSet-ready boolean matrix from `tf_cross_lineage.tsv`.
- §5 — Per-dataset top-10 targets by `distance_mean`, joined with TF metadata.
- §6 — Optional plot snippet (commented).

Run end-to-end: `uv run python working_groups/wg1_data_qc/examples/examples.py` from the jamboree folder root.
