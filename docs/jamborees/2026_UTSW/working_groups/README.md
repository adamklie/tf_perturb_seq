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
