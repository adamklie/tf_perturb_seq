# WG4 — GRN inference

Topic 2.3 / Figure 3. The goal is to build causal and mechanistic gene regulatory networks from TFP3 perturbation data, then integrate them with multiome (E2G linking, ChromBPNet) to separate direct from indirect regulation. End-of-jamboree deliverables include per-lineage edge lists, a cross-lineage network-structure panel for Figure 3, and a written plan for the multiome-integration follow-up.

## Questions

- Integrate multiome data (E2G linking, ChromBPNet) with TF-gene networks from Perturb-seq.
- Compare TF importance inferred from multiome vs Perturb-seq; identify which TFs act through direct binding vs indirect mechanisms.
- What are the common themes in how disease genes are regulated?
- How does the structure of TF networks change across lineages?

## Data

| Dataset | CRISPR pipeline |
|---|:---:|
| Hon WTC11 Cardiomyocyte | ready |
| Huangfu HUES8 Definitive Endoderm | ready |
| Huangfu HUES8 Embryonic Stem Cell | ready |
| Gersbach WTC11 Hepatocyte | ready |
| Engreitz WTC11 Endothelial | - |

Per-dataset cards live under [`../../data/`](../../data/).

## Issues

- [Issue #11](https://github.com/adamklie/tf_perturb_seq/issues/11): cross-dataset DEG count discrepancy. Matched WTC11 iPSC benchmarks show trans DEG counts varying several-fold across technologies while direct-target and cis hits stay consistent. Candidate causes include reads per cell, cells per element, pipeline QC, calibration sensitivity, and guide-assignment differences. Treat absolute trans-edge counts cautiously across datasets; per-gene effect-size correlation is the more comparable signal.
- *[FILL IN issue link]*: Calibrated DE tables (per dataset) are deferred — plan and draft implementation live at [`../../../../../src/tf_perturb_seq/inference/calibrate.py`](../../../../../src/tf_perturb_seq/inference/calibrate.py).

## Working flow

Three steps, in order:

1. **Brainstorm** — sketch example figures, summary tables, and pseudocode that answer the questions above. Capture this in a notebook, doc, or notes in this folder.
2. **Execute** — run the analyses; commit code (notebooks, scripts, supporting docs) to this folder on GitHub.
3. **Share** — upload reusable outputs (figures, tables, intermediate data) to WG4's Synapse folder [`syn74954084`](https://www.synapse.org/Synapse:syn74954084) (mirrored `working_groups/wg4_grn_inference/` under [`syn64423137/2026_UTSW/`](https://www.synapse.org/Synapse:syn64423137)). Record the syn ID for each upload in the [Outputs](#outputs) table below so the next person can find it.

## Outputs

| Object | Syn ID | Description | Owner |
|---|---|---|---|
| *[FILL IN as outputs land]* | | | |
