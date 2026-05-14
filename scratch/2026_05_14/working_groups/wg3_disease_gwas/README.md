# WG3 — Disease and GWAS

Topic 2.2 / Figure 2. Goal: connect TF regulatory activity to human disease.

Files under [`examples/`](examples/) are illustrative starting points, not finished deliverables.

## Questions

- Identify which TFs regulate disease/GWAS genes in each lineage.
- For TFs that are disease genes in multiple lineages, assess whether their activity is convergent or divergent.
- Identify and annotate GWAS variants near important TFs or the regulatory elements upstream of their downstream targets.

## Artifacts in `examples/`

| ID | File | Status | What it answers | Source data |
|---|---|---|---|---|
| WG3-A | [`examples/disease_tf_activity.tsv`](examples/disease_tf_activity.tsv) | ready (3 lineages so far) | Per-TF: is it a disease gene (Mondo/OMIM via HPO), and how strongly does its perturbation alter the transcriptome in each lineage? Auto-widens as more datasets land. | `tf_metadata.tsv` + per-dataset `wg1_significant_tfs.tsv` + HPO genes_to_disease (MONDO + OMIM) cached at `reference/gene_disease_associations.tsv` |
| WG3-B | [`examples/tf_convergence_scorecard.tsv`](examples/tf_convergence_scorecard.tsv) | ready (3 lineages so far) | Disease TFs with ED data in ≥2 lineages, refined classification (convergent_high / convergent_moderate / convergent_low / divergent_\<lineage\> / divergent_partial), plus distance range + max/min ratio | WG3-A |
| WG3-C | `examples/gwas_variants_near_tfs.tsv` (planned) | blocked | Per-TF GWAS variants in its locus or upstream regulatory elements | External GWAS catalog (not yet in repo) |

Status legend: `ready` (built and refreshable), `caveat` (partial coverage), `blocked` (upstream missing).

## Source choices

Gene-disease associations come from HPO's `genes_to_disease.txt` (https://hpo.jax.org/, released alongside the HPO ontology). It provides per-gene disease IDs across MONDO + OMIM (the agreed Mondo-based sources). Pulled fresh via [`reference/scripts/fetch_hpo_gene_disease.py`](../../reference/scripts/fetch_hpo_gene_disease.py); cached at [`../../reference/gene_disease_associations.tsv`](../../reference/gene_disease_associations.tsv). 5,090 unique gene symbols, 512 overlap with the TF library.

## Issues

- *[FILL IN issue link]*: `distance_mean` is the trustworthy signal for ED; `pval_mean` is anti-conservative for the Huangfu runs.

## WG3-A snapshot (3 lineages: Hon CM × Huangfu DE × Huangfu ESC, 2026-05-11)

512 of 1,983 TFs in the library have HPO disease associations (~26%). Cross-lineage classification within disease-TFs:

| classification | n_TFs |
|---|---:|
| convergent_significant (sig in all 3 lineages) | 0 |
| discordant_partial (sig in 2 of 3 lineages) | 3 |
| HonCM-specific | 39 |
| HuangfuESC-specific | 19 |
| HuangfuDE-specific | 14 |
| no_data (gene absent from ED runs) | 5 |
| convergent_nonsignificant | 432 |

Hon CM contributes the most lineage-specific disease TFs (39 vs 19 ESC, 14 DE) — consistent with Hon CM having the highest overall TF hit rate in WG1-D (164 sig total). Top ESC-specific disease TFs by rank shift: MEF2A (rank 1926 → 4), DNAJC21, RB1, SMARCB1, MYCN, CREB3L3, TNXB, GZF1, KMT2B. None of the WG1-D convergent_significant trio (TERF2 / GTF2B / ZNF574) are HPO-disease-flagged — worth a manual scan against ClinVar / GWAS Catalog when WG3 meets.

## WG3-B snapshot (507 of 512 WG3-A TFs have data in ≥2 lineages, 2026-05-11)

Refines the WG3-A classification by combining sig-or-not with magnitude (75th-percentile cutoff on max_distance_across_datasets = 558.3 — within the disease-TF cohort, not the full TF library).

| refined class | n_TFs | meaning |
|---|---:|---|
| convergent_high | 0 | Sig in all lineages with data. Still empty after Hon CM addition — consistent with WG3-A and the calibration concern. |
| convergent_moderate | 98 | Not sig in any lineage with data, but ≥75th-percentile max_distance. Borderline cases where calibration may be the limiting factor. |
| convergent_low | 334 | Not sig in any lineage, low magnitude across the board. Quiet TFs across lineages. |
| divergent_HonCM | 38 | Sig in Hon CM only. Highest count. |
| divergent_HuangfuESC | 19 | Sig in ESC only. |
| divergent_HuangfuDE | 15 | Sig in DE only. Caveat: most of these have higher absolute distance in ESC and pass the DE threshold mainly because DE's NC_max is lower. |
| divergent_partial | 3 | Sig in 2 of 3 lineages. |

Calibration-driven artifact to watch: "divergent_HuangfuDE" TFs typically have higher absolute distance in ESC than DE; they count as DE-specific by the calibration-robust threshold, but treat as starting points for hand-review.

Top divergent_HuangfuESC by max_distance: MEF2A, DNAJC21, RB1, PPP1R13L, IRF8, SMARCB1, MYCN, CREB3L3, TNXB, GZF1. Top convergent_moderate by max_distance (high-magnitude borderline — flag for re-check after calibration resolves): SPEN, ZNF750, PAX9, ATOH7, DLX4, FOXE1, ZIC2, HR, ALX3, MSX2.

## Run the examples

[`examples/examples.py`](examples/examples.py) walks the WG3 workflow:

- §1 — Top-20 disease-flagged TFs by max distance across lineages.
- §2 — Convergence-class counts.
- §3 — Deep-dive on the top convergent TF: per-lineage downstream targets from `wg4_tf_gene_edges_FDR05.tsv`.
- §4 — Disease-association lookup for that TF (from `reference/gene_disease_associations.tsv`).
- §5 — Divergent TFs sorted by distance range.

Run end-to-end: `uv run python working_groups/wg3_disease_gwas/examples/examples.py` from the jamboree folder root.
