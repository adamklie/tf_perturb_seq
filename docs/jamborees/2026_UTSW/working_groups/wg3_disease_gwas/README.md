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
| WG3-A | `disease_tf_activity.tsv` | ✅ landed (2 lineages so far) | Per-TF: is it a disease gene (Mondo/OMIM via HPO), and how strongly does its perturbation alter the transcriptome in each lineage? Auto-widens as more datasets land. | `tf_metadata.tsv` + per-dataset `wg1_significant_tfs.tsv` + HPO genes_to_disease (MONDO + OMIM) cached at `reference/gene_disease_associations.tsv` |
| WG3-B | `tf_convergence_scorecard.tsv` | ✅ landed (2 lineages so far) | Disease TFs with ED data in ≥2 lineages, refined classification (convergent_high / convergent_moderate / convergent_low / divergent_\<lineage\>), plus distance range + max/min ratio | WG3-A |
| WG3-C | `gwas_variants_near_tfs.tsv` | 🔴 blocked | Per-TF GWAS variants in its locus or upstream regulatory elements | External GWAS catalog (not yet in repo) |

## Source choices

- **Gene-disease associations**: HPO's `genes_to_disease.txt` (https://hpo.jax.org/, released alongside the HPO ontology). Provides per-gene disease IDs across MONDO + OMIM (the agreed Mondo-based sources). Pulled fresh via [`scripts/fetch_hpo_gene_disease.py`](../../../../scripts/fetch_hpo_gene_disease.py); cached at [`../../reference/gene_disease_associations.tsv`](../../reference/gene_disease_associations.tsv). 5,090 unique gene symbols → 512 overlap with our TF library.

> **⚠ Calibration caveat for ED-based fields**: same as WG1 — `distance_mean` is the trustworthy signal; `pval_mean` is anti-conservative for the Huangfu runs (see [`../../issues/edistance-calibration/`](../../issues/edistance-calibration/)).

## WG3-A first snapshot (2 lineages: Huangfu DE × Huangfu ESC, 2026-05-11)

- **512 of 1,983 TFs** in the library have HPO disease associations (~26%)
- Cross-lineage classification within disease-TFs:

| classification | n_TFs |
|---|---:|
| convergent_significant (sig in both lineages) | **0** |
| HuangfuESC-specific | 21 |
| HuangfuDE-specific | 16 |
| no_data (gene absent from ED runs) | 5 |
| convergent_nonsignificant | 470 |

Top ESC-specific disease TFs (sorted by distance): **MEF2A** (rank 1926 → rank 4 across lineages — striking shift), **DNAJC21**, **RB1**, **SMARCB1**, **MYCN**, **CREB3L3**, **TNXB**, **GZF1**, **KMT2B**. Top DE-specific: covered in the TSV. Notably none of the WG1-D convergent_significant trio (TERF2 / GTF2B / ZNF574) are HPO-disease-flagged — they're real TFs but lack curated disease associations in the HPO release; worth a manual scan against ClinVar / GWAS Catalog when WG3 meets.

## WG3-B snapshot (507 of 512 WG3-A TFs have data in ≥2 lineages, 2026-05-11)

Refines the WG3-A classification by combining sig-or-not with magnitude (75th-percentile cutoff on max_distance_across_datasets = 558.3). 75th percentile chosen within the disease-TF cohort, not the full TF library.

| refined class | n_TFs | meaning |
|---|---:|---|
| **convergent_high** | **0** | Sig in all lineages with data. **Empty**: consistent with WG3-A and the calibration concern (Issue 1). |
| convergent_moderate | 105 | Not sig in any lineage with data, but ≥75th-percentile max_distance. Borderline cases where calibration may be the limiting factor — re-check after Issue 1 resolves. |
| convergent_low | 365 | Not sig in any lineage, low magnitude across the board. Genuinely quiet TFs. |
| divergent_HuangfuESC | 21 | Sig in ESC only. Real ESC-specific responsive disease TFs. |
| divergent_HuangfuDE | 16 | Sig in DE only. **Caveat**: most of these still have higher absolute distance in ESC — they pass the DE threshold mainly because DE's NC_max is lower. |

**Calibration-driven artifact to watch**: the "divergent_HuangfuDE" TFs typically have `distance_mean_HuangfuESC` ≈ 500–560 vs `distance_mean_HuangfuDE` ≈ 110. They count as DE-specific by the calibration-robust threshold, but the absolute distance pattern says ESC is at least as responsive. Treat WG3-B's `divergent_*` labels as starting points for hand-review, not as final calls.

**Top divergent_HuangfuESC by max_distance** (real "ESC-on, DE-off" disease TFs to inspect): MEF2A, DNAJC21, RB1, PPP1R13L, IRF8, SMARCB1, MYCN, CREB3L3, TNXB, GZF1.

**Top convergent_moderate by max_distance** (high-magnitude borderline disease TFs — flag for re-check after calibration resolves): SPEN, ZNF750, PAX9, ATOH7, DLX4, FOXE1, ZIC2, HR, ALX3, MSX2.
