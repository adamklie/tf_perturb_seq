# WG2 — Gene program annotation and interpretation

Topic 2.1 / Figure 2. Goal: biological questions about TF activity across lineages — which programs are conserved, which are lineage-specific, who regulates them.

Files under [`examples/`](examples/) are illustrative starting points, not finished deliverables.

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

## Artifacts

Coverage caveat for all C/D/E rows: built from 3 of 5 datasets (Hon_CM, Huangfu_definitive, Huangfu_embryonic). Gersbach_Hep and Engreitz_Endo are still upstream. Meta-clustering ran at `meta_k ∈ {20,30,50,60,70,80,100}`; **k=60 is the canonical k** (only k where Shared_GO is computed). Other k values live under the same root.

Results are NOT mirrored into this directory — paths below point at the live outputs under `/oak/stanford/groups/engreitz/Users/ymo/Testing_Scripts/Meta_program/Result/` (referred to below as `$META_RESULT`).

| ID | File | Status | What it answers | Source data |
|---|---|---|---|---|
| WG2-C | `$META_RESULT/meta_consensus/k_60/cluster_summary.k_60.tsv` (+ `cluster_assignments.k_60.txt`, `meta_spectra.k_60.median.txt`, `clustering.k_60.png`) | caveat (3/5 datasets, k=60) | Meta-consensus clustering of programs across datasets: KMeans on l2-normalized stacked `gene_spectra_score` (rows = programs across datasets, k=60 meta-clusters). One row per meta-cluster with n_programs, n_datasets, member programs, category. | per-dataset `gene_spectra_score.k_<sel>.dt_2_0.txt` (Hon_CM k=80, Huangfu_definitive k=50, Huangfu_embryonic k=50) |
| WG2-D | `$META_RESULT/meta_consensus/k_60/shared_vs_specific.k_60.tsv` | caveat (3/5 datasets, k=60) | Each program classified as `shared_all` / `shared_partial` / `specific_<dataset>` from its meta-cluster membership | WG2-C cluster assignments |
| WG2-E | `$META_RESULT/Shared_GO/Shared_GO_summary.k_60.tsv` (also `.xlsx`, `_long.k_60.tsv`) | discussion (3/5 datasets, k=60) | Per-meta-cluster worksheet: top-10 genes + top-10 GO per member program, cross-dataset gene/GO Jaccard, intersections | WG2-C clusters + per-program top genes + Enrichr GO |

Status legend: `ready` (built), `caveat` (partial coverage), `blocked` (upstream missing), `discussion` (needs human decision).

Note: WG2-C was re-framed from "all-vs-all cosine similarity matrix" (original spec) to the meta-consensus KMeans clustering output — the clustering is the downstream artifact a similarity matrix would have been used to produce.

## Scripts

The analysis pipeline now lives in [`scripts/`](scripts/) (moved here from `Testing_Scripts/Meta_program/scripts/`). Results still write to `$META_RESULT` per each script's config.

| Subdir | What it does | Entry point |
|---|---|---|
| [`scripts/meta_consensus/`](scripts/meta_consensus/) | KMeans meta-consensus across per-dataset cNMF runs (WG2-C, WG2-D). Sweeps `meta_k ∈ {20,30,50,60,70,80,100}`. | `meta_consensus.py` (driver: `run_meta_consensus.py`, config: `config.yaml`) |
| [`scripts/Shared_GO/`](scripts/Shared_GO/) | Per-meta-cluster top-gene + GO Jaccard worksheet (WG2-E). | `shared_GO_analysis.py` |
| [`scripts/PCA/`](scripts/PCA/) | PCA visualization of programs colored by meta-cluster. | `run_pca.py`, `pca_plots.py` |
| [`scripts/UMAP/`](scripts/UMAP/) | UMAP visualization of programs colored by meta-cluster. | `run_umap.py`, `umap_plots.py` |
| [`scripts/CM_ED_SC_analysis/`](scripts/CM_ED_SC_analysis/) | 3-way regulator/expressed-gene overlap (CM × DE × ESC), Venn diagrams. | `Script/CM_ED_SC_overlap.py` (SLURM: `Script/run_CM_ED_SC_overlap.sh`) |

Conda env: `NMF_Benchmarking`. The meta_consensus config (`scripts/meta_consensus/config.yaml`) hardcodes per-dataset Inference dirs and the output_dir — update those paths when onboarding new datasets.

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
