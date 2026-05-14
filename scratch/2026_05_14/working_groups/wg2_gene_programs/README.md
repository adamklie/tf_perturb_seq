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

## Artifacts in `examples/`

| ID | File | Status | What it answers | Source data |
|---|---|---|---|---|
| WG2-C | `examples/program_similarity_matrix.tsv` (planned) | caveat (DE × ESC only) | All-vs-all cosine similarity of program loadings across datasets | per-dataset `gene_spectra_score.k_<sel>.dt_2_0.txt` |
| WG2-D | `examples/program_classification.tsv` (planned) | caveat | Each program classified as lineage_shared / lineage_specific / cell_state_program | WG2-C similarity matrix |
| WG2-E | `examples/program_annotation_worksheet.tsv` (planned) | discussion | Pre-filled worksheet for group hand-curation of program names and biology | WG2-A top-genes + GO enrichment |

Status legend: `ready` (built), `caveat` (partial coverage), `blocked` (upstream missing), `discussion` (needs human decision).

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
