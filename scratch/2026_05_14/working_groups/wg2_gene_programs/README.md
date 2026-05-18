# WG2 — Gene program annotation and interpretation

Topic 2.1 / Figure 2. Goal: biological questions about TF activity across lineages — which programs are conserved, which are lineage-specific, who regulates them.

Files under [`examples/`](examples/) are illustrative starting points, not finished deliverables.

## Result location

All results from the analyses below live under:

```
/oak/stanford/groups/engreitz/Users/ymo/Testing_Scripts/Meta_program/Result/
```

Referred to below as `$META_RESULT`. Results are NOT mirrored into this repo — only paths are recorded here. Analysis scripts that produced them now live in [`scripts/`](scripts/).

## Artifacts

Coverage caveat for all C/D/E rows: built from 3 of 5 datasets (Hon_CM, Huangfu_definitive, Huangfu_embryonic). Gersbach_Hep and Engreitz_Endo are still upstream. Meta-clustering ran at `meta_k ∈ {20,30,50,60,70,80,100}`; **k=60 is the canonical k** (only k where Shared_GO is computed). Other k values live under `$META_RESULT/meta_consensus/k_<k>/`.

| ID | File | Status | What it answers | Source data |
|---|---|---|---|---|
| WG2-C | `$META_RESULT/meta_consensus/k_60/cluster_summary.k_60.tsv` (+ `cluster_assignments.k_60.txt`, `meta_spectra.k_60.median.txt`, `clustering.k_60.png`) | caveat (3/5 datasets, k=60) | Meta-consensus clustering of programs across datasets: KMeans on l2-normalized stacked `gene_spectra_score` (rows = programs across datasets, k=60 meta-clusters). One row per meta-cluster with n_programs, n_datasets, member programs, category. | per-dataset `gene_spectra_score.k_<sel>.dt_2_0.txt` (Hon_CM k=80, Huangfu_definitive k=50, Huangfu_embryonic k=50) |
| WG2-D | `$META_RESULT/meta_consensus/k_60/shared_vs_specific.k_60.tsv` | caveat (3/5 datasets, k=60) | Each program classified as `shared_all` / `shared_partial` / `specific_<dataset>` from its meta-cluster membership | WG2-C cluster assignments |
| WG2-E | `$META_RESULT/Shared_GO/Shared_GO_summary.k_60.tsv` (also `.xlsx`, `_long.k_60.tsv`) | discussion (3/5 datasets, k=60) | Per-meta-cluster worksheet: top-10 genes + top-10 GO per member program, cross-dataset gene/GO Jaccard, intersections | WG2-C clusters + per-program top genes + Enrichr GO |

Status legend: `ready` (built), `caveat` (partial coverage), `blocked` (upstream missing), `discussion` (needs human decision).

Note: WG2-C was re-framed from "all-vs-all cosine similarity matrix" (original spec) to the meta-consensus KMeans clustering output — the clustering is the downstream artifact a similarity matrix would have been used to produce.

## Analyses

Four analyses produced the results in `$META_RESULT`. Each section lists what it does, the input it consumes, and what files it writes.

### Meta-consensus clustering (`$META_RESULT/meta_consensus/k_<k>/`)

Stacks `gene_spectra_score` matrices from each per-dataset cNMF run into a single (programs × genes) matrix, l2-normalizes rows, then runs KMeans (n_init=10, random_state=1; mirrors `torch_cnmf/cnmf.py:1588`) at every `k ∈ {20,30,50,60,70,80,100}`. Each program is assigned to a meta-cluster; clusters with members from multiple datasets are `shared_all`/`shared_partial`, single-dataset clusters are `specific_<dataset>`. Inputs: `Inference.gene_spectra_score.k_<sel>.dt_2_0.txt` for Hon_CM (k=80), Huangfu_embryonic (k=50), Huangfu_definitive (k=50). Outputs per k:

- `cluster_assignments.k_<k>.txt` — `<dataset>|<program_id>` → `meta_cluster_id`
- `cluster_summary.k_<k>.tsv` — one row per meta-cluster: n_programs, n_datasets, per-dataset counts, category, member list
- `shared_vs_specific.k_<k>.tsv` — per-program (dataset, program_id, meta_cluster_id, category) — **WG2-D**
- `meta_spectra.k_<k>.median.txt` — median consensus spectra (meta-clusters × genes)
- `category_counts.k_<k>.txt`, `silhouette.k_<k>.txt` — quality summaries
- `clustering.k_<k>.png`, `composition.k_<k>.png`, `top_genes.k_<k>.png`, `category_breakdown.k_<k>.png` — visualizations
- `cluster_<id>_zoom.png` — focused view of clusters listed under `zoom_clusters` in config

### Shared GO worksheet (`$META_RESULT/Shared_GO/`)

For each meta-cluster (k=60 only), extracts top-N genes per member program from the consensus spectra, queries Enrichr GO_Biological_Process_2023, and computes pairwise gene-list + GO-term Jaccard within and across datasets. Outputs:

- `Shared_GO_summary.k_60.tsv` (+ `.xlsx`) — one row per meta-cluster: meta top-10 genes, cross-dataset GO intersection, cross-dataset gene intersection, per-member top-10 genes + top-10 GO terms — **WG2-E**
- `Shared_GO_long.k_60.tsv` — long-form per-(meta_cluster, member_program) for downstream filtering

### PCA / UMAP of programs (`$META_RESULT/PCA/k_<k>/`, `$META_RESULT/UMAP/k_<k>/`)

2D embeddings of the same stacked-and-l2-normalized program matrix used for meta-consensus, with points colored by meta-cluster assignment from KMeans. Sanity-checks that the KMeans clusters separate in low-dim and identifies outlier programs. Outputs per k (k ∈ {20,30,50,80,100} — k=60 not yet rendered):

- `pca_programs_by_cluster.png`
- `umap_programs_by_cluster.png`

### CM × DE × SC overlap (`$META_RESULT/CM_ED_SC_analysis/`)

3-way set overlap across Hon_CM, Huangfu_DE, and Huangfu_ESC of (a) expressed genes, (b) all TFs targeted by guides, and (c) expressed TFs targeted by guides. Frames whether shared meta-clusters reflect shared regulators or just shared underlying biology. Outputs:

- `summary_counts.tsv`, `set_sizes.tsv` — per-category region counts
- `regions/{expressed_genes,regulators,expressed_regulators}__<region>.tsv` — 7 region files per category (common_all_three, CM_only, DE_only, SC_only, CM_DE_not_SC, CM_SC_not_DE, DE_SC_not_CM)
- `venn_{expressed_genes,regulators,expressed_regulators}.{pdf,png}` — Venn diagrams
- `CM_ED_SC_overlap.xlsx` — workbook summary

## Scripts

The analysis pipeline lives in [`scripts/`](scripts/) (moved here from `Testing_Scripts/Meta_program/scripts/`). Results still write to `$META_RESULT` per each script's config.

| Subdir | What it does | Entry point |
|---|---|---|
| [`scripts/meta_consensus/`](scripts/meta_consensus/) | KMeans meta-consensus across per-dataset cNMF runs (WG2-C, WG2-D). Sweeps `meta_k ∈ {20,30,50,60,70,80,100}`. | `meta_consensus.py` (driver: `run_meta_consensus.py`, config: `config.yaml`) |
| [`scripts/Shared_GO/`](scripts/Shared_GO/) | Per-meta-cluster top-gene + GO Jaccard worksheet (WG2-E). | `shared_GO_analysis.py` |
| [`scripts/PCA/`](scripts/PCA/) | PCA visualization of programs colored by meta-cluster. | `run_pca.py`, `pca_plots.py` |
| [`scripts/UMAP/`](scripts/UMAP/) | UMAP visualization of programs colored by meta-cluster. | `run_umap.py`, `umap_plots.py` |
| [`scripts/CM_ED_SC_analysis/`](scripts/CM_ED_SC_analysis/) | 3-way regulator/expressed-gene overlap (CM × DE × ESC), Venn diagrams. | `Script/CM_ED_SC_overlap.py` (SLURM: `Script/run_CM_ED_SC_overlap.sh`) |

Conda env: `NMF_Benchmarking`. The meta_consensus config (`scripts/meta_consensus/config.yaml`) hardcodes per-dataset Inference dirs and the output_dir — update those paths when onboarding new datasets.

## Questions

- Extract gene programs for each production dataset; create a heatmap of loading similarities across programs, labeled by lineage of origin.
- Identify programs highly similar across lineages (expected basic cellular processes) and annotate them.
- Identify lineage-specific programs and assess whether they correspond to lineage-specific biological processes.
- (Stretch) Assess whether gene programs can be mapped to in-vivo counterparts.
- For programs shared across lineages: which perturbations alter program usage in the same direction? Which TFs contribute to lineage-agnostic vs lineage-specific programs?
- Within lineages: identify regulators of each program (barplot or heatmap of TF perturbations associated with each program).
- For shared programs: do they share regulators, or have distinct ones? Pull out interesting case studies.
- Assess whether any TFs contribute to lineage bifurcations.
- Apply Percoder to assess perturbation sensitivity of programs across lineages.

## Per-dataset companions (under `../../data/<dataset>/cnmf/<run_name>/`)

| Artifact | Path pattern | Status |
|---|---|---|
| WG2-A Top-N genes per program | `data/<dataset>/cnmf/<run_name>/wg2_top_genes_per_program.tsv` | ready for Huangfu DE/ESC (Synapse [`syn74895462`](https://www.synapse.org/Synapse:syn74895462) / [`syn74895475`](https://www.synapse.org/Synapse:syn74895475)) |
| WG2-B Regulators per program | `data/<dataset>/cnmf/<run_name>/wg2_regulators_per_program.tsv` | ready for Huangfu DE/ESC (perturbation_association_results files exist) |

## Issues

- *[FILL IN issue link]*: Hon CM cNMF — not run yet — production launch pre-staged but held; gated on validating the HTv2 testbed.
- *[FILL IN issue link]*: Gersbach Hep cNMF — awaiting Sara's deliverables.
- *[FILL IN issue link]*: Engreitz cNMF — blocked on portal data.

## Run the examples

[`examples/examples.py`](examples/examples.py) is currently gated — it prints the recipes for each step (load `gene_spectra_score`, build similarity matrix, classify programs, extract top-N genes, regulators per program) and points at the code to run once cNMF sources are mirrored. Replace the `print(...)` blocks with the indented code as data lands.

Run as-is: `uv run python working_groups/wg2_gene_programs/examples/examples.py` — it will print the recipes without computing anything.
