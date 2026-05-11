# Working-group outputs — 2026 UTSW jamboree

One folder per working group. Each folder holds the cross-dataset rolled-up artifacts (TSVs / markdown / figures) that answer that group's questions — sized for slide decks and exploratory notebooks, not schema-completeness.

For per-dataset companion artifacts (e.g., per-dataset trans-target counts, per-dataset top-N-genes-per-program), see the matching `datasets/<dataset>/<analysis>/` folder. Each WG README points at its per-dataset companions.

Status legend:

- ✅ **landed** — built, in place, refreshable
- 🟡 **partial** — built from data we have; gets stronger as more datasets / runs land
- 🔴 **blocked** — upstream data missing (cNMF runs / collaborator deliveries / external resource)
- 🤔 **discussion** — analysis choice needs human decision

| WG | Folder | Topic / Figure | Lead questions |
|---|---|---|---|
| 1 | [`wg1_data_qc/`](wg1_data_qc/) | Topic 1 / Fig 1 | Guide detection & repression; transcriptome-wide significance; cross-lineage shared TFs; trans-target overlap; technical harmonization |
| 2 | [`wg2_gene_programs/`](wg2_gene_programs/) | Topic 2.1 / Fig 2 | Cross-lineage program similarity; lineage-shared vs lineage-specific programs; regulators per program; bifurcation TFs |
| 3 | [`wg3_disease_gwas/`](wg3_disease_gwas/) | Topic 2.2 / Fig 2 | Disease-gene TFs; convergent vs divergent activity across lineages; GWAS variants near TFs |
| 4 | [`wg4_grn_inference/`](wg4_grn_inference/) | Topic 2.3 / Fig 3 | TF→gene network edge lists; cross-lineage network structure; multiome integration |
| 5 | [`wg5_tf_family_case_studies/`](wg5_tf_family_case_studies/) | Topic 2.4 / Fig 4 | TF family activity scorecard; per-family pathway enrichment; family-specific case studies |
| 6 | [`wg6_predictive_modeling/`](wg6_predictive_modeling/) | Topic 3 / Fig 5 | Model task spec; baseline model results |

Scientific scope: [`../TOPICS.md`](../TOPICS.md), [`../WORKING_GROUPS.md`](../WORKING_GROUPS.md). Onboarding: [`../GETTING_STARTED.md`](../GETTING_STARTED.md).

## Quick-reference: all pre-computed WG artifacts (2026-05-11)

For loading the underlying Synapse-mirrored bundles (CRISPR pipeline / cNMF / energy distance), see [`../GETTING_STARTED.md`](../GETTING_STARTED.md). The artifacts below are derivative roll-ups — open them directly in `pandas` / a spreadsheet to skip the upstream pull-and-process work.

### Cross-dataset roll-ups (sized for slide decks)

| Artifact | Status | Top finding | Load |
|---|---|---|---|
| [`wg1_data_qc/qc_summary.tsv`](wg1_data_qc/qc_summary.tsv) | ✅ | 5-row × ~18-col QC matrix for slide decks | `pd.read_csv(p, sep="\t")` |
| [`wg1_data_qc/edistance_summary.tsv`](wg1_data_qc/edistance_summary.tsv) | 🟡 (2/5) | Per-dataset ED significance counts with multi-threshold + calibration note | `pd.read_csv(p, sep="\t")` |
| [`wg1_data_qc/tf_cross_lineage.tsv`](wg1_data_qc/tf_cross_lineage.tsv) | ✅ (3 lineages, auto-widens) | One row per TF, per-dataset distance/rank/sig + classification. **TERF2/GTF2B/ZNF574** convergent across all 3 lineages (Hon CM + Huangfu DE + Huangfu ESC) | `pd.read_csv(p, sep="\t")` |
| [`wg3_disease_gwas/disease_tf_activity.tsv`](wg3_disease_gwas/disease_tf_activity.tsv) | ✅ (3 lineages) | 512 disease-flagged TFs × per-dataset ED; Hon CM-specific count (39) leads | `pd.read_csv(p, sep="\t")` |
| [`wg3_disease_gwas/tf_convergence_scorecard.tsv`](wg3_disease_gwas/tf_convergence_scorecard.tsv) | ✅ (3 lineages) | Refined classification (convergent_high/moderate/low, divergent_\<lineage\>, divergent_partial); top ESC-on/DE-off: MEF2A, DNAJC21, RB1 | `pd.read_csv(p, sep="\t")` |
| [`wg4_grn_inference/network_structure_by_lineage.tsv`](wg4_grn_inference/network_structure_by_lineage.tsv) | ✅ (3 lineages) | **Hon CM is 5× DE's edge count, 18× ESC's** (220k / 41k / 12k); all cross-lineage jaccards 0.4–3% (massive rewiring everywhere) | `pd.read_csv(p, sep="\t")` |
| [`wg5_tf_family_case_studies/family_activity_scorecard.tsv`](wg5_tf_family_case_studies/family_activity_scorecard.tsv) | ✅ (3 lineages) | 84 families ≥3 members; **23** candidates_for_deepdive (Hon CM adds 10). Top: C2H2 ZF, Myb/SANT (4→11 with Hon CM), HMG/Sox (new) | `pd.read_csv(p, sep="\t")` |

### Per-dataset companions (under `../datasets/<dataset>/<analysis>/`)

| Artifact pattern | Status | What it gives you |
|---|---|---|
| `energy_distance/wg1_significant_tfs.tsv` | ✅ Hon CM (2,036 / 164 sig) + Huangfu DE (2,273 / 73 sig) + Huangfu ESC (2,273 / 83 sig) | Per-target ED + ranking + significance flags joined with TF metadata |
| `crispr_pipeline/wg4_tf_gene_edges_FDR05.tsv` | ✅ Hon CM (219,587) + Huangfu DE (41,418) + ESC (12,366) | TF→gene edges at per-TF BH FDR<0.05; joined with TF identity (symbol/family/DBD) |
| `crispr_pipeline/wg1_trans_target_counts.tsv` | ✅ Hon CM (2,065 perts) + Huangfu DE (1,741) + ESC (1,452) | Per-perturbation trans-target counts; top targets list. Hon CM top: ISL1 (2,769). DE top: SOX17 (4,361). ESC top: STRAP (938) |
| `crispr_pipeline/<prefix>_calibrated_*_results.tsv` (calibrated DE tables) | 🔴 deferred | Per-dataset calibrated DE tables (NTC-based empirical p-values + BH; cis / direct_target / trans-only splits). Plan + draft code at [`src/tf_perturb_seq/inference/calibrate.py`](../../../../src/tf_perturb_seq/inference/calibrate.py). Will run once the Synapse mirror is filled out; output schema documented in `src/tf_perturb_seq/inference/plan.md`. Drives [Issue #11](https://github.com/adamklie/tf_perturb_seq/issues/11). |

### What you actually do with these

- **WG1**: open `tf_cross_lineage.tsv`, filter `classification == "convergent_significant"` for cross-lineage hits, or `*-specific` for lineage-distinguishing TFs.
- **WG2**: blocked on cNMF until Huangfu DE/ESC cNMF mirrors mature; check [`wg2_gene_programs/README.md`](wg2_gene_programs/README.md) for the current state.
- **WG3**: open `tf_convergence_scorecard.tsv`, sort by `convergence_class` then `max_distance_across_datasets` for deep-dive selection. Mind the calibration-driven `divergent_HuangfuDE` caveat (most "DE-specific" entries have higher absolute distance in ESC).
- **WG4**: open `network_structure_by_lineage.tsv` for the per-lineage stats; load `datasets/<id>/crispr_pipeline/wg4_tf_gene_edges_FDR05.tsv` for the edge lists themselves (NetworkX-ready: `nx.from_pandas_edgelist(df, "intended_target_name", "gene_id")`).
- **WG5**: open `family_activity_scorecard.tsv`, filter `candidate_for_deepdive == True`, pick a family, then pull its members' edges from `wg4_tf_gene_edges_FDR05.tsv` for the family-specific target gene list.

> **⚠ Universal caveat for ED-based artifacts**: `pval_mean` is anti-conservative for the Huangfu runs (see [`../issues/edistance-calibration/`](../issues/edistance-calibration/)). All cross-lineage classifications above use `distance_mean > NC max` as the calibration-robust significance proxy. Raw `pval_mean<0.05` counts are reported but flagged with ⚠ in the per-WG READMEs.

## Refresh / extend

Each TSV has a corresponding `scripts/build_<artifact>.py` that re-runs the build off the upstream sources. Add a new dataset by:

1. Mirror its CRISPR pipeline / ED output to Synapse (via the corresponding `scripts/mirror_*.py`).
2. Drop the per-dataset derivative TSVs under `datasets/<new_id>/<analysis>/` (the per-dataset build scripts handle this).
3. Re-run the cross-dataset scripts (`build_wg1_tf_cross_lineage.py`, `build_wg3_disease_tf_activity.py`, `build_wg3_tf_convergence_scorecard.py`, `build_wg4_network_structure_by_lineage.py`, `build_wg5_tf_family_scorecard.py`). They all auto-discover datasets by filesystem scan — no code edits needed when new datasets land.
