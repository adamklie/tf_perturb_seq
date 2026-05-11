# WG3 — Disease & GWAS

**Topic 2.2 / Figure 2.** Goal: connect TF regulatory activity to human disease.

## Questions

From [`../../WORKING_GROUPS.md`](../../WORKING_GROUPS.md):

- Identify which TFs regulate disease/GWAS genes in each lineage.
- For TFs that are disease genes in multiple lineages, assess whether their activity is convergent or divergent.
- Identify and annotate GWAS variants near important TFs or the regulatory elements upstream of their downstream targets.

## Artifacts in this folder

| ID | File | Status | What it answers | Source data |
|---|---|---|---|---|
| WG3-A | `disease_tf_activity.tsv` | ✅ ready (Mondo) | Per-TF: is it a disease gene (Mondo), and how strongly does its perturbation alter the transcriptome in each lineage? | `tf_metadata.tsv` + per-dataset `pval_edist_full.csv` + Mondo Disease Ontology |
| WG3-B | `tf_convergence_scorecard.tsv` | 🟡 partial | For disease-relevant TFs with data in ≥2 lineages: convergent_high / convergent_low / divergent classification | WG3-A + cross-lineage ED join |
| WG3-C | `gwas_variants_near_tfs.tsv` | 🔴 blocked | Per-TF GWAS variants in its locus or upstream regulatory elements | External GWAS catalog (not yet in repo) |

## Source choices

- **Disease-gene list**: [Mondo Disease Ontology](https://mondo.monarchinitiative.org/) (open, comprehensive). Pulled fresh from the public OBO/JSON release whenever the script runs. WG3-A's `mondo_disease_ids` column carries the matched Mondo IDs per TF for downstream filtering.

> **⚠ Calibration caveat for ED-based fields**: same as WG1 — `distance_mean` is the trustworthy signal; `pval_mean` is anti-conservative for the Huangfu runs (see [`../../issues/edistance-calibration.md`](../../issues/edistance-calibration.md)).
