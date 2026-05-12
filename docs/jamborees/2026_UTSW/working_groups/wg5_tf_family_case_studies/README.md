# WG5 — TF family case studies

**Topic 2.4 / Figure 4.** Goal: deep-dive analysis of specific TFs / TF families with newly implicated roles in lineage differentiation (e.g., ZNF factors).

## Questions

From [`../../WORKING_GROUPS.md`](../../WORKING_GROUPS.md):

- Identify interesting TF families using outputs from Working Groups 1–3.
- Pathway analysis of top gene programs for these TFs, with DEG fold-change overlays.
- Use ATAC-seq data to assess differential accessibility of binding motifs upstream of top regulated genes.
- Highlight GWAS SNPs in or near these TFs or their regulatory elements linked to disease in the relevant cell type.

## Artifacts in this folder

| ID | File | Status | What it answers | Source data |
|---|---|---|---|---|
| WG5-A | `family_activity_scorecard.tsv` | ✅ landed (3 lineages so far; widens as more land) | Per-family: n members, n significant in any lineage, mean / max distance per dataset, disease fraction, cross-lineage rollups, candidate_for_deepdive flag | `tf_metadata.tsv` (`jaspar_tf_family` w/ fallback to `lambert_2018_dbd`) + per-dataset `wg1_significant_tfs.tsv` + HPO MONDO+OMIM disease flags from `reference/gene_disease_associations.tsv` |
| WG5-B | `family_<family>_enrichment.tsv` (one file per family selected) | 🔴 blocked | KEGG / GO pathway enrichment for the family's regulated genes | cNMF programs + target gene lists (cNMF not yet run for all datasets) |

## Source choices

- **Family taxonomy**: prefer `jaspar_tf_family` (clean, curated TF family taxonomy). For TFs not in JASPAR core, fall back to `lambert_2018_dbd` prefixed with `DBD:`. TFs without either get grouped under `unannotated` (a coverage gap, never marked deepdive-candidate).
- **Min family size**: 3 — sub-3 families are case-studies, not statistical units.
- **Candidate-for-deepdive heuristic**: family has ≥2 significant members in any lineage (`distance_mean > NC_max`, calibration-robust) AND ≥1 disease-gene member. Adjustable; meant as the WG5 starting shortlist, not a final filter.

### WG5-A snapshot (3 lineages: Hon CM × Huangfu DE × Huangfu ESC, 2026-05-11)

84 families ≥ 3 members; **23 candidates_for_deepdive** (up from 13 with the addition of Hon CM). Top by n_sig_in_any_lineage (excluding `unannotated`):

| family | n_members | n_disease | n_sig_any_lineage | notes |
|---|---:|---:|---:|---|
| C2H2 ZF (DBD) | 483 | 47 | 46 | Nearly doubled with Hon CM (26→46). Dominates lineage-specific WG1-D ZNF callouts. |
| More than 3 adjacent zinc fingers (JASPAR) | 177 | 18 | 14 | Subset of C2H2 ZF as defined by JASPAR class. |
| Homeodomain (DBD) | 52 | 18 | 11 | 35% disease-gene fraction — strongest disease-density family in our library. |
| **Myb/SANT (DBD)** | 28 | 4 | **11** | **Jumped from 4 → 11 with Hon CM**. Hosts TERF2 (WG1-D convergent_significant). |
| bHLH (DBD) | 49 | 9 | 7 | Includes MyoD/HAND/NeuroD-family lineage drivers. |
| Multiple dispersed zinc fingers (JASPAR) | 47 | 11 | 7 | |
| HMG/Sox (DBD) | 30 | 10 | 6 | **New entry with Hon CM addition** — SOX-family is a cardiomyocyte story. |
| Three-zinc finger Kruppel-related (JASPAR) | 28 | 5 | 6 | 21% hit rate — small focused family, high payoff per perturbation. |
| HOX (JASPAR) | 49 | 12 | 4 | Body-plan TFs — relevant cross-lineage. |
| Ets-related | 27 | 6 | 4 | |
| Paired-related HD factors | 28 | 16 | 3 | Highest disease-fraction (57%) among small families. |
| FOX | 27 | 10 | 3 | FOXH1 is a DE master (4k trans targets in WG1-E). |
| MBD (DBD) | 8 | 1 | 3 | **New entry with Hon CM** — methyl-CpG binding domain proteins. |
| POU domain factors | 18 | 8 | 2 | **New entry** — includes POU5F1 (OCT4). |
| Tal-related | 21 | 8 | 2 | |
| HD-LIM | 8 | 2 | 2 | |

**`unannotated` flagged as a coverage gap (unchanged)**: 359 TFs lack both JASPAR and Lambert annotation. Worth a manual curation pass before WG5 meets.

**Notable Hon CM-driven shifts**: Myb/SANT (4→11 sig), HMG/Sox (new top family), C2H2 ZF (26→46). These suggest cardiomyocyte differentiation engages a broader set of TF families than the ESC↔DE comparison alone shows. Three new candidate-for-deepdive families enter on Hon CM's contribution: HMG/Sox, MBD, POU domain factors.

## Run the examples

[`examples.py`](examples.py) walks the WG5 workflow:

- §1 — `candidate_for_deepdive` families ranked by `n_sig_in_any_lineage`.
- §2 — Pull family members from `reference/tf_metadata.tsv` (default pick: C2H2 ZF).
- §3 — Per-landed-dataset edge counts for those members from `wg4_tf_gene_edges_FDR05.tsv`.
- §4 — JASPAR motif lookup against `reference/jaspar_core_tf_metadata.tsv`.
- §5 — GO-enrichment stub (gseapy).

Run end-to-end: `uv run python working_groups/wg5_tf_family_case_studies/examples.py` from the jamboree folder root.

## Out of scope here

- **ATAC-seq accessibility + motif analysis**: multiome stream (same as WG4-C).
- **GWAS SNP overlay**: external resource (same as WG3-C).
