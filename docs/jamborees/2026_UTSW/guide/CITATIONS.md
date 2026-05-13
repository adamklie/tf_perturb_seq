# Citations

One-line summary + canonical URL for each tool / paper that produced or framed the outputs in this project.

---

## Methods used in our pipeline

### CRISPR data processing (IGVF CRISPR pipeline)
**Repository:** <https://github.com/IGVF/CRISPR_Pipeline>
A Nextflow pipeline that takes raw fastqs + a guide library, maps reads, assigns guides to cells, calls cis + trans differential expression, and emits a QC dashboard. Output: `pipeline_dashboard/`, `pipeline_outputs/`, `pipeline_info/`.

### perturbo — primary differential expression
A Bayesian mixed-effects model for CRISPR Perturb-seq DE calling, used in our pipeline as the canonical trans-DE caller. Run in both per-element (one row per perturbation) and per-guide (one row per guide) modes.
**Method:** https://github.com/insitro/perturbo

### sceptre — alternative differential expression
A robust DE caller for single-cell CRISPR screens that uses conditional randomization to handle confounders.
**Method:** Barry et al., *Genome Biology* 2024. https://katsevich-lab.github.io/sceptre/

### Energy distance method (transcriptome-wide perturbation effect)
A distribution-based test for whether two single-cell populations differ. Pipeline implemented by Chikara Takeuchi (UTSW).
**Pipeline:** <https://github.com/Chikara-Takeuchi/energy_dist_pipeline>
**Method:** Peidli et al., *scPerturb*, *Nature Methods* 2024. https://www.nature.com/articles/s41592-023-02144-y

### cNMF — gene-program discovery
Consensus non-negative matrix factorization for discovering gene programs from single-cell data.
**Method:** Kotliar et al., *eLife* 2019. https://elifesciences.org/articles/43803
**Implementation used here:** torch-cNMF wrapped by PerturbNMF. https://github.com/EngreitzLab/PerturbNMF

---

## Reference annotations

### Lambert 2018 TF list
The canonical human TF inventory with DNA-binding domain annotations.
**Lambert et al., *Cell* 2018.** https://doi.org/10.1016/j.cell.2018.01.029

### JASPAR
TF binding motif database (we use the CORE collection for human).
**Castro-Mondragon et al., JASPAR 2024.** https://jaspar.elixir.no/

### GENCODE V43
Gene model and annotation set used as the reference universe for all our analyses.
https://www.gencodegenes.org/

---

## Reference Perturb-seq papers (the "look for these things" framing)

The kinds of biology + figure types our docs lean on:

### Paper 1 — Comprehensive perturbation of TFs in human cardiomyocytes
TF Perturb-seq in iPSC-CMs revealing the regulatory architecture of congenital heart disease.
**Key takeaways:** lineage-master TFs with broad trans effects (TBX5, MEF2 family); convergent regulatory modules; deorphanization of disease-gene TFs via guilt-by-association.
**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12724611/

### Paper 2 — Benchmarking and optimizing Perturb-seq in differentiating hPSCs
Methods + QC best practices for multi-lineage Perturb-seq, including how clonal-selection artifacts surface (e.g. TP53 perturbation expanding clones).
**Key takeaways:** sustained-perturbation phenotypes; how to detect clonal artifacts during the 8–12-day differentiation; progressive regulatory-network assembly.
**URL:** https://www.cell.com/stem-cell-reports/fulltext/S2213-6711(25)00317-0

### Other foundational Perturb-seq references
- **Dixit et al., *Cell* 2016** — the original Perturb-seq paper. https://doi.org/10.1016/j.cell.2016.11.038
- **Replogle et al., *Cell* 2022** — single-cell CRISPRi at scale (Weissman lab). https://doi.org/10.1016/j.cell.2022.05.013
- **Joung et al., *Cell* 2023** — ORF + CRISPRi atlas of ~3,500 TF isoforms. https://doi.org/10.1016/j.cell.2023.05.005

---

## ENCODE / IGVF context

### IGVF consortium
The umbrella project this data lives under.
**URL:** https://www.igvf.org/

### IGVF data portal
Where raw + processed data are released.
**URL:** https://data.igvf.org/

### ENCODE CRISPR benchmarking
Established the QC framework adopted by IGVF.
**Citation:** https://www.nature.com/articles/s41592-024-02404-5

---

## Project-specific resources

See [`docs/REFERENCES.md`](https://github.com/adamklie/tf_perturb_seq/blob/main/docs/REFERENCES.md) for the full project resource map (Slack, Synapse, internal Google Docs, GitHub).
