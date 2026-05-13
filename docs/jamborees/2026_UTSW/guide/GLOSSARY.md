# Glossary

Plain-language definitions. Italicized cross-references like *[NTC]* point at other entries below.

---

**AnnData** — A Python data structure for single-cell data. A 2-D matrix (cells × genes) plus a `.obs` table (per-cell metadata, like cell-type or guide-identity) and a `.var` table (per-gene metadata). Saved as `.h5ad`. See [`DATA_FORMATS.md`](DATA_FORMATS.md).

**BH FDR / Benjamini–Hochberg false-discovery rate** — A method for adjusting p-values when you're testing thousands of hypotheses at once (e.g. one test per gene). "FDR < 0.05" means: among the hits we call significant, ≤5% are expected to be false alarms.

**Calibration** — Whether a statistical test gives the right number of false positives when nothing is happening. A *[NTC]* is the natural null: if you compare two random groups of NTC cells, you should see ~5% of "significant" results at p < 0.05. If you see way more or way fewer, the test is mis-calibrated.

**cis effect** — When a guide knocks down its **intended target** gene's transcript level. We measure cis effect to confirm the guide is doing what it should ("did the knockdown work?"). Synonym: *on-target* effect.

**cNMF (consensus non-negative matrix factorization)** — A method that decomposes a cell × gene matrix into "gene programs" (groups of co-varying genes) and "cell usages" (how much each cell uses each program). Think of programs as biological modules — cell cycle, mitochondrial respiration, lineage identity, etc. The output of cNMF tells you which programs are active in which cells, and which TFs disrupted which programs.

**CRISPRi (CRISPR interference)** — A version of CRISPR that *suppresses* a gene's expression without cutting the DNA. Uses a catalytically-dead Cas9 ("dCas9") tethered to a transcriptional repressor (e.g. KRAB) and a guide RNA targeting the gene's promoter. Knockdown is typically 70–95%.

**dataset bundle** — All the files we package for one experiment (one cell type × one perturbation library × one lab). See [`COMPLETE_DATASET_CONTENTS.md`](COMPLETE_DATASET_CONTENTS.md).

**energy distance** — A statistic that measures how different two distributions of single-cell profiles are. "How different are the cells with TF X knocked down vs. the non-targeting baseline?" — answered in one number per TF. Bigger = the perturbation moved the transcriptome more. The associated p-value tests "is this distance significantly larger than what you'd see by chance?" See [`UNDERSTAND_ENERGY_DISTANCE.md`](UNDERSTAND_ENERGY_DISTANCE.md).

**ENSG / Ensembl gene ID** — A stable identifier for a human gene, like `ENSG00000164853`. Different gene-name databases use different symbols (HGNC, GENCODE, UCSC), but ENSGs are stable across them.

**FDR** — see *[BH FDR]*.

**gene program** — A group of genes that vary together across cells. Outputs of *[cNMF]*. Examples: a "cell cycle program," a "stress-response program," a "cardiomyocyte-identity program." Each program is defined by a "loading" score per gene (how strongly that gene belongs to the program).

**GRN (gene regulatory network)** — A directed graph: TFs → downstream target genes. Built from trans-DE results: "TF X knockdown moved gene Y" → "TF X regulates gene Y." Comparing GRNs across cell types shows how regulatory wiring rewires across lineages.

**guide / sgRNA / gRNA** — A small RNA that directs Cas9 (or dCas9) to a specific genomic location. Each TF in our library typically has 4–6 guides; multiple guides per TF lets us check that effects are reproducible across guides (not artifacts of one bad guide).

**h5ad** — File extension for a saved *[AnnData]*.

**h5mu / MuData** — File extension for a saved **MuData** — a multi-modal extension of AnnData that holds *both* the RNA modality (cells × genes) *and* the guide modality (cells × guides) in one object. This is the canonical CRISPR-pipeline output. See [`DATA_FORMATS.md`](DATA_FORMATS.md).

**HVG (highly variable genes)** — The subset of genes whose expression varies the most across cells in a dataset. cNMF runs on the top ~2,000 HVGs (other genes are too noisy or uninformative).

**inference MuData** — The final MuData a CRISPR pipeline emits, containing the QC-filtered RNA + guide modalities with one cell per row, ready for downstream analysis. The canonical file is `inference_mudata.h5mu`.

**intended target** — The gene a guide is supposed to knock down. Stored in the guide annotation column `intended_target_name` (a symbol like `SOX17`) or `intended_target_ensembl` (the ENSG).

**knockdown efficiency** — How much a guide actually suppresses its target's transcript (typically expressed as % suppression vs. NTC cells). Computed as part of *[QC]*.

**library** — The full set of guides used in an experiment. Our TF guide library has ~14,000 guides covering ~2,000 TFs + non-targeting + positive controls + negative controls. Pools A-D (and sometimes F) refer to physically distinct library subsets.

**log2FC / log2 fold change** — How much a gene changed expression, on log scale. `log2FC = -1` means halved; `log2FC = 1` means doubled. We report log2FC per perturbation × gene.

**MOI (multiplicity of infection)** — How many guides each cell got on average. Higher MOI = more cells get useful guide-cell pairings, but harder to attribute effects to a single guide. Typically targeted at 1.5–3.

**MuData** — see *[h5mu]*.

**negative control** — A guide that targets a region not expected to affect the cell (e.g. a safe-harbor gene like *AAVS1*). Distinct from *[NTC]* (non-targeting) — negative controls *target* a region but a benign one; NTCs don't target anything in the genome.

**NTC / non-targeting control** — A guide whose sequence doesn't match anywhere in the genome (or matches no protein-coding gene). Cells with an NTC are the "baseline" — they got the CRISPR machinery but no perturbation. We use NTCs to define the null distribution for all our significance tests.

**perturbo** — The differential-expression caller used at the end of the IGVF CRISPR pipeline. Reports per-guide and per-element (per-intended-target) log2FC + p-values for every gene, in both cis and trans flavors.

**positive control** — A guide that targets a gene expected to have a strong effect (e.g. an essential gene like *AARS*). Used to confirm the experiment "works" before trusting other results.

**Perturb-seq** — Single-cell RNA-seq where each cell also carries a guide RNA. After sequencing, each cell has both a gene-expression profile and a "which guide did it have" label.

**sceptre** — An alternative differential-expression caller used in the IGVF CRISPR pipeline (alongside perturbo). The pipeline runs both; perturbo trans-DE is the canonical output we report.

**spacer / protospacer** — The 19–20-nt sequence within a guide RNA that matches the target genomic site.

**Synapse** — A data-sharing platform (synapse.org) where we mirror final outputs for collaborators. Each file gets a stable ID like `syn74834952`. See [`DATA_FORMATS.md`](DATA_FORMATS.md).

**trans effect** — When a guide knockdown changes a gene other than its intended target. Trans effects capture *downstream* regulatory consequences — the biology Perturb-seq is built to find.

**TSS / transcription start site** — The genomic coordinate where transcription of a gene begins. Many guide libraries target the promoter near the TSS to maximize CRISPRi knockdown.

**TF Universe** — A DACC-spec file listing every TF whose promoter is targeted in this experiment. See [`COMPLETE_DATASET_CONTENTS.md`](COMPLETE_DATASET_CONTENTS.md).

**Element Universe** — A DACC-spec BED file listing every genomic element (promoter, enhancer) targeted by a guide in this experiment.

**UMAP** — A 2-D dimensionality-reduction plot of single cells. Used for visual QC and to color cells by guide identity / cell state.

---

If a term you need isn't here, ping us — adding entries is cheap.
