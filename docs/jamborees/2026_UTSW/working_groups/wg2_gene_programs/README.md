# WG2 — Gene program annotation and interpretation

Topic 2.1 / Figure 2. The goal is to use cNMF gene programs to ask biological questions about TF activity across lineages: which programs are conserved, which are lineage-specific, and who regulates them. End-of-jamboree deliverables include a cross-lineage program similarity panel, a shortlist of annotated lineage-shared and lineage-specific programs, and per-program regulator summaries for the lineages with cNMF available.

## Questions

- Extract gene programs per dataset and cluster them by loading similarity across lineages.
- Identify programs that are highly similar across lineages (basic cellular processes) and annotate them.
- Identify lineage-specific programs and ask whether they correspond to lineage-specific biology.
- For shared programs, which perturbations alter program usage in the same direction? Which TFs contribute to lineage-agnostic vs lineage-specific programs?
- Within each lineage, which TFs are the strongest regulators of each program?
- Do shared programs share regulators or have distinct ones?
- Do any TFs contribute to lineage bifurcations?
- (Stretch) Map programs to in-vivo counterparts; apply Percoder to assess perturbation sensitivity of programs across lineages.

## Data

WG2 is gated on cNMF availability.

| Dataset | cNMF |
|---|:---:|
| Hon WTC11 Cardiomyocyte | - |
| Huangfu HUES8 Definitive Endoderm | ready |
| Huangfu HUES8 Embryonic Stem Cell | ready |
| Gersbach WTC11 Hepatocyte | - |
| Engreitz WTC11 Endothelial | - |

Per-dataset cards live under [`../../data/`](../../data/).

## Issues

- *[FILL IN issue link]*: Hon CM cNMF not yet run
- *[FILL IN issue link]*: Gersbach Hep cNMF not yet run
- *[FILL IN issue link]*: Engreitz Endo cNMF — blocked on portal data.

## Working flow

Three steps, in order:

1. **Brainstorm** — sketch example figures, summary tables, and pseudocode that answer the questions above. Capture this in a notebook, doc, or notes in this folder.
2. **Execute** — run the analyses; commit code (notebooks, scripts, supporting docs) to this folder on GitHub.
3. **Share** — upload reusable outputs (figures, tables, intermediate data) to WG2's Synapse folder [`syn74954081`](https://www.synapse.org/Synapse:syn74954081) (mirrored `working_groups/wg2_gene_programs/` under [`syn64423137/2026_UTSW/`](https://www.synapse.org/Synapse:syn64423137)). Record the syn ID for each upload in the [Outputs](#outputs) table below so the next person can find it.

## Outputs

| Object | Syn ID | Description | Owner |
|---|---|---|---|
| *[FILL IN as outputs land]* | | | |
