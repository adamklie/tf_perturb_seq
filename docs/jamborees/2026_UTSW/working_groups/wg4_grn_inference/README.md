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
| WG4-B | `network_structure_by_lineage.tsv` | ✅ landed (3 lineages so far; widens as more land) | Per-lineage network summary stats (degree distributions, TF→TF fraction, cross-lineage edge/TF overlap) | WG4-A per-dataset edge lists |

## Per-dataset companions (under `datasets/<dataset>/crispr_pipeline/`)

| Artifact | Path pattern | Status |
|---|---|---|
| WG4-A TF→gene edge list (FDR 0.05) | `datasets/<dataset>/crispr_pipeline/wg4_tf_gene_edges_FDR05.tsv` | ✅ landed for **Hon CM** (GCS: `gs://igvf-pertub-seq-pipeline-data/.../2026_04_15/outs/seqspec_v3/pipeline_outputs/`) + Huangfu DE (Synapse syn74834952) + Huangfu ESC (syn74835010). Gersbach Hep + Engreitz Endo blocked. |
| **Trans-DE results table** (folks-ready, IGVF-shaped) | `datasets/<dataset>/crispr_pipeline/trans_de_results.tsv.gz` | ✅ landed for all 3 datasets with full gene-symbol + chr/start/end/strand annotation. Modeled on IGVF's `global differential expression` files (e.g. IGVFFI5989UAVX). Hon CM: 219,587 rows (10.6 MB). Huangfu DE: 41,418 rows (2.1 MB). Huangfu ESC: 12,366 rows (644 KB). |

### WG4-A snapshot (3 datasets: Hon CM × Huangfu DE × Huangfu ESC, 2026-05-11)

Per-dataset edge list, filtered at per-TF Benjamini-Hochberg FDR<0.05. The per-TF BH approach treats each TF's genome-wide tests as a family (lenient/standard), vs. genome-wide BH which is harsh and biases against polygenic TFs.

| dataset | n_edges | n_TFs_with_sig_edges | median_edges_per_sig_TF | TFs with >100 sig edges |
|---|---:|---:|---:|---:|
| **Hon CM**   | **219,587** | **2,065** | **49** | **537** |
| Huangfu DE   | 41,418 | 1,741 | 5 | 41 |
| Huangfu ESC  | 12,366 | 1,452 | 2 | 15 |

**Hon CM is the most responsive system by a wide margin**: 5× DE's edge count, 18× ESC's. Median 49 sig edges per TF (vs 5 DE / 2 ESC); 537 TFs have >100 sig edges (vs 41 / 15). Note that Hon CM's ED-significant count (164) is also higher than the Huangfu runs (73 / 83), so the two readouts agree on the direction even if trans-effect magnitudes are larger than ED would predict.

> **⚠ Trans-DEG count discrepancy across datasets** ([issue #11](https://github.com/adamklie/tf_perturb_seq/issues/11)): The same perturbation in matched WTC11 iPSC benchmarks shows trans DEG counts varying 3× across technologies (Engreitz ~1.4K vs Huangfu ~4.2K at FDR<0.1), while direct-target and cis hits are consistent. Candidate causes: reads/cell (11K–147K), cells/element, pipeline QC, calibration sensitivity, guide-assignment differences. Our 3-lineage spread is the same phenomenon but across production lineages — Hon CM's 5–18× density vs the Huangfu runs may partly reflect WTC11 vs HUES8 + newer `seqspec_v3` pipeline. Treat absolute counts cautiously; effect-size correlation across datasets (per-gene log2_fc) is the more comparable signal.

DE has ~3× ESC's edge count, consistent with definitive-endoderm being a more transcriptionally responsive committed state vs. the buffered pluripotent ESC baseline. Top-out-degree TFs match lineage biology (SOX17 / FOXH1 in DE; POU5F1 / SALL4 in ESC — see WG1-E for per-TF counts).

> **⚠ Cis vs trans**: this is `perturbo_trans_per_element_output.tsv.gz` — trans (off-target locus) effects only. Cis (on-target knockdown of the perturbed gene itself) is reported separately in `perturbo_cis_*`. Use these edges for downstream-gene network analysis; don't double-count the TF's own gene.

### WG4-B snapshot (3 lineages: Hon CM × Huangfu DE × Huangfu ESC, 2026-05-11)

| metric | **Hon CM** | Huangfu DE | Huangfu ESC |
|---|---:|---:|---:|
| n_TFs with ≥1 sig edge | **2,065** | 1,741 | 1,452 |
| n_target_genes | 7,257 | 8,410 | 5,064 |
| n_edges | **219,587** | 41,418 | 12,366 |
| mean TF out-degree | **106.3** | 23.8 | 8.5 |
| median TF out-degree | **49** | 5 | 2.5 |
| max TF out-degree (top "hub" perturbation) | 2,769 | 4,361 (SOX17) | 938 (STRAP) |
| 95th-percentile TF out-degree | 387.8 | 44 | 21 |
| median gene in-degree | **16** | 3 | 1 |
| max gene in-degree | 448 | 94 | 33 |
| n TF→TF edges | 17,278 | 4,018 | 1,281 |
| fraction TF→TF | 7.9% | 9.7% | 10.4% |
| median \|log2_fc\| of sig edges | 0.26 | 0.28 | 0.34 |

**Cross-lineage edge overlap (jaccard, all pairs)**:

| pair | shared edges | shared TFs | jaccard edges |
|---|---:|---:|---:|
| Huangfu DE ↔ ESC  | 1,533 | 1,228 | **0.029** |
| Hon CM ↔ Huangfu DE  | 2,645 | 1,676 | **0.010** |
| Hon CM ↔ Huangfu ESC | 919 | 1,381 | **0.004** |

**Headline for WG4 Q4**: massive cross-lineage rewiring is universal. DE↔ESC has the highest overlap (still only ~3% of the union) — they share lab + cell line (Huangfu/HUES8). Hon CM ↔ either Huangfu is much lower (0.4–1.0% of edges shared), partly because Hon CM is WTC11 + a different pipeline version. The "same TF perturbed in two lineages → mostly different downstream targets" story holds across all 3 pairs.

The TF→TF backbone fraction is similar (~8–10%) across all 3, suggesting the architectural backbone (TFs regulating other TFs) is more conserved than the leaf-gene targets. This is the calibration-robust takeaway for the figure.

**Per-lineage degree-distribution shape**: Hon CM has the densest network (mean 106 out-degree, median 49) with the heaviest median gene in-degree (16). DE has the longest tail (max out-degree 4,361 = SOX17). ESC is sparse end-to-end. Hon CM's density may partly reflect the newer pipeline (`seqspec_v3`) — flag for harmonization.

## Out of scope here

- **Multiome integration (WG4-C)**: E2G linking + ChromBPNet outputs are a separate analysis stream from TFP3. Pull those in once they exist; this folder won't track them.
