# WG4 — GRN inference

**Topic 2.3 / Figure 3.** Goal: build causal / mechanistic gene regulatory networks from TFP3 perturbation data, integrated with multiome.

## Questions

From [`../../WORKING_GROUPS.md`](../../WORKING_GROUPS.md):

- Integrate multiome data (E2G linking, ChromBPNet models) with TF-gene regulatory networks from Perturb-seq.
- Compare TF importance inferred from multiome vs Perturb-seq; identify which TFs act through direct binding vs indirect mechanisms.
- Assess common themes in how disease genes are regulated.
- Characterize how the structure of TF networks changes across lineages.

## Artifacts in this folder

| ID | File | Status | What it answers | Source data |
|---|---|---|---|---|
| WG4-B | `network_structure_by_lineage.tsv` | ✅ landed (2 lineages so far; widens as more land) | Per-lineage network summary stats (degree distributions, TF→TF fraction, cross-lineage edge/TF overlap) | WG4-A per-dataset edge lists |

## Per-dataset companions (under `datasets/<dataset>/crispr_pipeline/`)

| Artifact | Path pattern | Status |
|---|---|---|
| WG4-A TF→gene edge list (FDR 0.05) | `datasets/<dataset>/crispr_pipeline/wg4_tf_gene_edges_FDR05.tsv` | ✅ landed for Huangfu DE + ESC (Synapse syn74834952 / syn74835010). Hon CM via syn74520421. Gersbach Hep blocked on Sara. |

### WG4-A first snapshot (2 datasets: Huangfu DE × Huangfu ESC, 2026-05-11)

Per-dataset edge list, filtered at per-TF Benjamini-Hochberg FDR<0.05. The per-TF BH approach treats each TF's genome-wide tests as a family (lenient/standard), vs. genome-wide BH which is harsh and biases against polygenic TFs.

| dataset | n_edges | n_TFs_with_sig_edges | median_edges_per_sig_TF | TFs with >100 sig edges |
|---|---:|---:|---:|---:|
| Huangfu DE  | 41,418 | 1,741 | 5 | 41 |
| Huangfu ESC | 12,366 | 1,452 | 2 | 15 |

DE has ~3× more significant edges than ESC, consistent with a more transcriptionally responsive cellular state during commitment to definitive endoderm vs. the relatively buffered pluripotent ESC baseline. Top-out-degree TFs match lineage biology (SOX17 / FOXH1 in DE; POU5F1 / SALL4 in ESC — see WG1-E for the per-TF counts).

> **⚠ Cis vs trans**: this is `perturbo_trans_per_element_output.tsv.gz` — trans (off-target locus) effects only. Cis (on-target knockdown of the perturbed gene itself) is reported separately in `perturbo_cis_*`. Use these edges for downstream-gene network analysis; don't double-count the TF's own gene.

### WG4-B first snapshot (2 lineages: Huangfu DE × Huangfu ESC, 2026-05-11)

| metric | Huangfu DE | Huangfu ESC |
|---|---:|---:|
| n_TFs with ≥1 sig edge | 1,741 | 1,452 |
| n_target_genes | 8,410 | 5,064 |
| n_edges | 41,418 | 12,366 |
| mean TF out-degree | 23.8 | 8.5 |
| median TF out-degree | 5 | 2.5 |
| max TF out-degree (top "hub" perturbation) | 4,361 (SOX17) | 938 (STRAP) |
| 95th-percentile TF out-degree | 44 | 21 |
| median gene in-degree | 3 | 1 |
| max gene in-degree (most-targeted gene) | 94 | 33 |
| n TF→TF edges (target is another perturbed TF) | 4,018 | 1,281 |
| fraction TF→TF | 9.7% | 10.4% |
| median \|log2_fc\| of sig edges | 0.28 | 0.34 |

**Cross-lineage rewiring (DE ↔ ESC)**:
- 1,228 TFs are perturbed and significant in *both* lineages (~70% of ESC's set; ~71% of DE's).
- But only **1,533 edges** are shared (TF → same target gene in both lineages).
- **Jaccard edge overlap = 0.029** — meaning ~3% of the union of significant edges is convergent. The downstream targets are almost entirely lineage-specific even when the same TF is wired in.

Implication for WG4 Q4 ("how does network structure change across lineages?"): substantial rewiring is the headline. The DE network is ~3.4× denser than ESC's, with much heavier tails (4,361 vs 938 out-degree max). The TF→TF backbone fraction is similar (~10%) in both, suggesting the architectural backbone is lineage-invariant even though the leaf targets are not.

## Out of scope here

- **Multiome integration (WG4-C)**: E2G linking + ChromBPNet outputs are a separate analysis stream from TFP3. Pull those in once they exist; this folder won't track them.
