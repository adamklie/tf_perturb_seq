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
| WG4-B | `network_structure_by_lineage.tsv` | 🟡 partial | Per-lineage network summary stats (n_TFs, n_edges, mean TF outdegree, median gene indegree, cross-lineage edge overlap) | WG4-A per-dataset edge lists |

## Per-dataset companions (under `datasets/<dataset>/crispr_pipeline/`)

| Artifact | Path pattern | Status |
|---|---|---|
| WG4-A TF→gene edge list (FDR 0.05) | `datasets/<dataset>/crispr_pipeline/wg4_tf_gene_edges_FDR0p05.tsv` | ✅ ready for Huangfu DE/ESC (Synapse syn74834952 / syn74835010). Hon CM via syn74520421. Gersbach Hep blocked on Sara. |

## Out of scope here

- **Multiome integration (WG4-C)**: E2G linking + ChromBPNet outputs are a separate analysis stream from TFP3. Pull those in once they exist; this folder won't track them.
