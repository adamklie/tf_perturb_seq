# Working-group outputs — 2026 UTSW jamboree

One folder per working group. Each folder holds a per-WG `README.md` (scope, questions, status) and an `examples/` subdirectory with illustrative build scripts, runnable examples, and small roll-up TSVs sized for slide decks.

Everything under `<wg>/examples/` is illustrative — starting points, not finished analyses or deliverables. Participants are expected to extend, replace, or supersede them with the analyses they actually run during the jamboree.

For per-dataset companion artifacts (e.g., per-dataset trans-target counts, per-dataset top-N-genes-per-program), see the matching `data/<dataset>/<analysis>/` folder. Each WG README points at its per-dataset companions.

## Working groups

| WG | Folder | Topic / Figure | Lead questions | Runnable examples |
|---|---|---|---|---|
| 1 | [`wg1_data_qc/`](wg1_data_qc/) | Topic 1 / Fig 1 | Guide detection & repression; transcriptome-wide significance; cross-lineage shared TFs; trans-target overlap; technical harmonization | [`wg1_data_qc/examples/examples.py`](wg1_data_qc/examples/examples.py) |
| 2 | [`wg2_gene_programs/`](wg2_gene_programs/) | Topic 2.1 / Fig 2 | Cross-lineage program similarity; lineage-shared vs lineage-specific programs; regulators per program; bifurcation TFs | [`wg2_gene_programs/examples/examples.py`](wg2_gene_programs/examples/examples.py) (gated on cNMF) |
| 3 | [`wg3_disease_gwas/`](wg3_disease_gwas/) | Topic 2.2 / Fig 2 | Disease-gene TFs; convergent vs divergent activity across lineages; GWAS variants near TFs | [`wg3_disease_gwas/examples/examples.py`](wg3_disease_gwas/examples/examples.py) |
| 4 | [`wg4_grn_inference/`](wg4_grn_inference/) | Topic 2.3 / Fig 3 | TF→gene network edge lists; cross-lineage network structure; multiome integration | [`wg4_grn_inference/examples/examples.py`](wg4_grn_inference/examples/examples.py) |
| 5 | [`wg5_tf_family_case_studies/`](wg5_tf_family_case_studies/) | Topic 2.4 / Fig 4 | TF family activity scorecard; per-family pathway enrichment; family-specific case studies | [`wg5_tf_family_case_studies/examples/examples.py`](wg5_tf_family_case_studies/examples/examples.py) |
| 6 | [`wg6_predictive_modeling/`](wg6_predictive_modeling/) | Topic 3 / Fig 5 | Model task spec; baseline model results | [`wg6_predictive_modeling/task_spec_template.md`](wg6_predictive_modeling/task_spec_template.md) (discussion) |

Onboarding: [`../README.md`](../README.md). Each WG's `examples/examples.py` is self-contained — dataset paths and small loader helpers are inlined at the top of each script (no shared library import). Run any example end-to-end from the jamboree folder root with `uv run python working_groups/wg<N>_*/examples/examples.py`.

## Pre-computed WG artifacts

For loading the underlying Synapse-mirrored sources (CRISPR pipeline / cNMF / energy distance), see [`../README.md`](../README.md). The artifacts below are derivative roll-ups — open them directly in `pandas` / a spreadsheet to skip the upstream pull-and-process work.

### Cross-dataset roll-ups (sized for slide decks)

| Artifact | Status | Description |
|---|---|---|
| [`wg1_data_qc/examples/qc_summary.tsv`](wg1_data_qc/examples/qc_summary.tsv) | ready | 5-row QC matrix sized for slide decks |
| [`wg1_data_qc/examples/edistance_summary.tsv`](wg1_data_qc/examples/edistance_summary.tsv) | caveat (2/5 datasets) | Per-dataset ED significance counts with multi-threshold + calibration note |
| [`wg1_data_qc/examples/tf_cross_lineage.tsv`](wg1_data_qc/examples/tf_cross_lineage.tsv) | ready (3 lineages, auto-widens) | One row per TF, per-dataset distance/rank/sig + classification |
| [`wg3_disease_gwas/examples/disease_tf_activity.tsv`](wg3_disease_gwas/examples/disease_tf_activity.tsv) | ready (3 lineages) | Disease-flagged TFs × per-dataset ED |
| [`wg3_disease_gwas/examples/tf_convergence_scorecard.tsv`](wg3_disease_gwas/examples/tf_convergence_scorecard.tsv) | ready (3 lineages) | Refined classification (convergent_high/moderate/low, divergent_\<lineage\>, divergent_partial) |
| [`wg4_grn_inference/examples/network_structure_by_lineage.tsv`](wg4_grn_inference/examples/network_structure_by_lineage.tsv) | ready (3 lineages) | Per-lineage edge count, degree distribution, cross-lineage edge/TF jaccards |
| [`wg5_tf_family_case_studies/examples/family_activity_scorecard.tsv`](wg5_tf_family_case_studies/examples/family_activity_scorecard.tsv) | ready (3 lineages) | 84 families with ≥3 members, candidate-for-deepdive flag |

Status legend: `ready` (built and refreshable), `caveat` (built from partial data; gets stronger as more datasets land), `blocked` (upstream data missing), `discussion` (analysis choice needs human decision).

Load with `pd.read_csv(p, sep="\t")`.

### Per-dataset companions (under `../data/<dataset>/<analysis>/`)

| Artifact pattern | Status | Description |
|---|---|---|
| `energy_distance/wg1_significant_tfs.tsv` | ready (Hon CM, Huangfu DE, Huangfu ESC) | Per-target ED + ranking + significance flags joined with TF metadata |
| `crispr_pipeline/wg4_tf_gene_edges_FDR05.tsv` | ready (Hon CM, Huangfu DE, Huangfu ESC) | TF→gene edges at per-TF BH FDR<0.05; joined with TF identity (symbol/family/DBD) |
| `crispr_pipeline/wg1_trans_target_counts.tsv` | ready (Hon CM, Huangfu DE, Huangfu ESC) | Per-perturbation trans-target counts; top targets list |
| `crispr_pipeline/<prefix>_calibrated_*_results.tsv` | discussion (deferred) | Per-dataset calibrated DE tables. Plan + draft code at [`src/tf_perturb_seq/inference/calibrate.py`](../../../../src/tf_perturb_seq/inference/calibrate.py). Drives [Issue #11](https://github.com/adamklie/tf_perturb_seq/issues/11). |

### How to use these

- WG1: open `tf_cross_lineage.tsv`, filter `classification == "convergent_significant"` for cross-lineage hits, or `*-specific` for lineage-distinguishing TFs.
- WG2: blocked on cNMF until production mirrors mature; check [`wg2_gene_programs/README.md`](wg2_gene_programs/README.md) for the current state.
- WG3: open `tf_convergence_scorecard.tsv`, sort by `convergence_class` then `max_distance_across_datasets` for deep-dive selection. Mind the calibration-driven `divergent_HuangfuDE` caveat.
- WG4: open `network_structure_by_lineage.tsv` for the per-lineage stats; load `data/<id>/crispr_pipeline/wg4_tf_gene_edges_FDR05.tsv` for the edge lists themselves (NetworkX-ready: `nx.from_pandas_edgelist(df, "intended_target_name", "gene_id")`).
- WG5: open `family_activity_scorecard.tsv`, filter `candidate_for_deepdive == True`, pick a family, then pull its members' edges from `wg4_tf_gene_edges_FDR05.tsv` for the family-specific target gene list.

## Issues

- *[FILL IN issue link]*: `pval_mean` is anti-conservative for the Huangfu runs. All cross-lineage classifications above use `distance_mean > NC max` as the calibration-robust significance proxy. Raw `pval_mean<0.05` counts are reported but flagged in the per-WG READMEs.

## Refresh / extend

Each TSV has a corresponding per-WG `build_*.py` (under each `working_groups/wg<N>/examples/`) that re-runs the build off the upstream sources. Add a new dataset by:

1. Mirror its CRISPR pipeline / ED output to Synapse (via the corresponding `data/scripts/mirror_*.py`).
2. Drop the per-dataset derivative TSVs under `data/<new_id>/<analysis>/` (the per-dataset build scripts handle this).
3. Re-run the cross-dataset scripts under `working_groups/wg<N>_*/examples/`. They all auto-discover datasets by filesystem scan — no code edits needed when new datasets land.
