# For experimentalists — reading our CRISPRi Perturb-seq data

Hi! You've been handed (or asked to look at) results from a CRISPRi Perturb-seq experiment from this project. This folder is a 1-stop reference for understanding what's in the bundle and how to read each output, **without writing any code**.

> This is the **generic** guide — it does not assume any specific dataset or jamboree. If you're a participant in a specific event (e.g. the 2026 UTSW jamboree), there will also be an event-specific dashboard (`docs/jamborees/<event>/communications/FOR_EXPERIMENTALISTS.md`) that lists the status of *that event's* datasets. Read this guide first; then go to the dashboard.

---

## What is CRISPRi Perturb-seq?

Briefly: we knock down ~2,000 transcription factors (TFs) one at a time across human cell types, then do single-cell RNA-seq to read out the transcriptome of each cell. We can then ask, for any TF: *what genes did its knockdown move?* and *what gene programs did it disrupt?* — and compare those answers across cell types.

The biology we look for follows the framing in two recent papers ([citations](CITATIONS.md)):

- **Lineage-master TFs with broad trans effects** — a small number of TFs that, when knocked down, perturb thousands of downstream genes (e.g. SOX17 in endoderm, ISL1 in cardiomyocytes).
- **Convergent regulatory modules across cell types** — TFs that change the same gene programs in different lineages are likely in shared complexes/pathways.
- **Lineage-specific TFs** — TFs that matter only in one cell type often correspond to that lineage's identity.
- **Dosage-sensitive hubs** — TFs whose graded knockdown produces graded phenotypes; therapeutic candidates.
- **Disease-gene TFs with unexpected regulators** — patient-mutation genes that cluster among targets of the same upstream TF; the data deorphanizes TFs by guilt-by-association with known disease genes.

Keep these in the back of your head as you read the outputs.

---

## What's in this folder

| File | What it covers |
|---|---|
| [`GLOSSARY.md`](GLOSSARY.md) | Plain-English definitions: cis, trans, perturbo, sceptre, NTC, MOI, knockdown, FDR, calibration, energy distance, cNMF program. |
| [`DATA_FORMATS.md`](DATA_FORMATS.md) | What a `.h5mu` / `.h5ad` / `.tsv.gz` is and how to open one. What "AnnData" / "MuData" mean. Synapse / GCS / IGVF portal basics. |
| [`QUICK_START.md`](QUICK_START.md) | Step-by-step: where to download, how to peek at the first table without code, and how to load one in Python or R if you want to. |
| [`UNDERSTAND_CRISPR_OUTPUTS.md`](UNDERSTAND_CRISPR_OUTPUTS.md) | The IGVF CRISPR pipeline outputs: the dashboard HTML, inference MuData, perturbo cis/trans TSVs. What columns mean and when to trust a result. |
| [`UNDERSTAND_ENERGY_DISTANCE.md`](UNDERSTAND_ENERGY_DISTANCE.md) | Per-target energy distance + permutation p-values. How to read the volcano / cutoff plots. The calibration caveat in one paragraph. |
| [`UNDERSTAND_CNMF_OUTPUTS.md`](UNDERSTAND_CNMF_OUTPUTS.md) | cNMF gene programs. What a "program" is, how to find which programs your TF altered, the k-selection plot, top-loaded genes per program. |
| [`COMPLETE_DATASET_CONTENTS.md`](COMPLETE_DATASET_CONTENTS.md) | The canonical file inventory: what a "complete" Perturb-seq bundle should contain, cross-referenced against the IGVF DACC submission spec and a reference IGVF analysis set ([IGVFDS6332VCTO](https://data.igvf.org/analysis-sets/IGVFDS6332VCTO/)). |
| [`FAQ.md`](FAQ.md) | Common questions: "why are some p-values exactly zero?", "how many TFs should be significant?", "what's the difference between cis and trans?" |
| [`CITATIONS.md`](CITATIONS.md) | Methods papers + key Perturb-seq references; quick 1-line summary of what each one contributes. |

---

## How to use this guide

1. **First time?** Read [`QUICK_START.md`](QUICK_START.md) — 10 minutes, then you can poke at a result file without code.
2. **Looking at a specific output?** Open the matching `UNDERSTAND_*` file — they each have a "What you'll see" section followed by a "How to read it" section.
3. **Stuck on a term?** [`GLOSSARY.md`](GLOSSARY.md) is alphabetical.
4. **Wondering if a file is missing?** Check [`COMPLETE_DATASET_CONTENTS.md`](COMPLETE_DATASET_CONTENTS.md) — it lists what a fully-baked dataset should have.
5. **Need to talk to someone?** Each `UNDERSTAND_*` doc ends with a "Who to ask" section. If in doubt, check [`docs/REFERENCES.md`](../REFERENCES.md) for Slack + GitHub issue pointers.

If you find anything confusing, that's a doc bug — please flag it in a GitHub issue so we can fix the explanation for the next reader.
