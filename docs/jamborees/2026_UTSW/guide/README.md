# Analysis guide — reading our CRISPRi Perturb-seq data

A 1-stop reference for what's in a CRISPRi Perturb-seq dataset directory and how to read each output. Written so anyone walking into the jamboree — whether you live in JupyterLab or at the bench — can navigate the outputs and pull the numbers they need.

> Pair this with the per-dataset directory at [`../data/`](../data/), which lists what's mirrored to Synapse for each production dataset in this event.

---

## What is CRISPRi Perturb-seq?

Briefly: we knock down ~2,000 transcription factors (TFs) one at a time across human cell types, then do single-cell RNA-seq to read out the transcriptome of each cell. We can then ask, for any TF: *what genes did its knockdown move?* and *what gene programs did it disrupt?* — and compare those answers across cell types.

Patterns to look for in the outputs:

- **Lineage-defining TFs with broad trans effects** — a small number of TFs that, when knocked down, perturb thousands of downstream genes.
- **Convergent regulatory modules across cell types** — TFs that change the same gene programs in different lineages are candidates for shared complexes or pathways.
- **Lineage-specific TFs** — TFs that matter only in one cell type often correspond to that lineage's identity.
- **Dosage-sensitive hubs** — TFs whose graded knockdown produces graded phenotypes; therapeutic candidates.
- **Disease-gene TFs with unexpected regulators** — patient-mutation genes that cluster among targets of the same upstream TF; the data can deorphanize TFs by association with known disease genes.

Keep these in the back of your head as you read the outputs.

---

## What's in this folder

| File | What it covers |
|---|---|
| [`GLOSSARY.md`](GLOSSARY.md) | Plain-English definitions: cis, trans, perturbo, sceptre, NTC, MOI, knockdown, FDR, calibration, energy distance, cNMF program. |
| [`DATA_FORMATS.md`](DATA_FORMATS.md) | What a `.h5mu` / `.h5ad` / `.tsv.gz` is and how to open one. What "AnnData" / "MuData" mean. Synapse / GCS / IGVF portal basics. |
| [`CRISPR.md`](CRISPR.md) | The IGVF CRISPR pipeline outputs: the dashboard HTML, inference MuData, perturbo cis/trans TSVs. What columns mean and when to trust a result. |
| [`ENERGY_DISTANCE.md`](ENERGY_DISTANCE.md) | Per-target energy distance + permutation p-values. How to read the volcano / cutoff plots. The calibration caveat in one paragraph. |
| [`CNMF.md`](CNMF.md) | cNMF gene programs. What a "program" is, how to find which programs your TF altered, the k-selection plot, top-loaded genes per program. |
| [`COMPLETE_DATASET_CONTENTS.md`](COMPLETE_DATASET_CONTENTS.md) | The shared file inventory: what a fully-baked Perturb-seq dataset should contain, cross-referenced against the IGVF DACC submission spec and a reference IGVF analysis set ([IGVFDS6332VCTO](https://data.igvf.org/analysis-sets/IGVFDS6332VCTO/)). |
| [`FAQ.md`](FAQ.md) | Common questions: "why are some p-values exactly zero?", "how many TFs should be significant?", "what's the difference between cis and trans?" |

---

## How to use this guide

1. **First time?** The [jamboree homebase](../README.md) has the env-setup and "load my first dataset" snippet — start there.
2. **Looking at a specific output?** Open the matching interpretation file (`CRISPR.md`, `ENERGY_DISTANCE.md`, `CNMF.md`) — they each have a "What you'll see" section followed by a "How to read it" section.
3. **Stuck on a term?** [`GLOSSARY.md`](GLOSSARY.md) is alphabetical.
4. **Wondering if a file is missing?** Check [`COMPLETE_DATASET_CONTENTS.md`](COMPLETE_DATASET_CONTENTS.md) — it lists what a fully-baked dataset should have.
5. **Need to talk to someone?** Each interpretation doc ends with a "Who to ask" section. If in doubt, file a GitHub issue.

If you find anything confusing, that's a doc bug — please flag it in a GitHub issue so we can fix the explanation for the next reader.
