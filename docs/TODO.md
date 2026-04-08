# TODO

## Key deadlines

| Date | Milestone |
|------|-----------|
| **Mar 12** | Benchmark datasets through CRISPR, cNMF, E-distance pipelines |
| **Mar 26** | >= 2 production datasets through CRISPR pipeline |
| **Apr 9** | Iteration on / additional production datasets through CRISPR pipeline |
| **Apr 23** | >= 2 production datasets through cNMF and E-distance pipelines |
| **May 7** | Preparation for jamboree |
| **May 13-16** | Jamboree at UTSW, Dallas |

---

## Pre-jamboree: Must have

### Benchmark datasets through all pipelines (by Mar 12)
- [ ] Complete CRISPR pipeline runs across all 5 benchmark datasets
- [ ] Run energy distance on all benchmark datasets
- [ ] Resolve promoter-level target name issue in E-distance (use coordinates to differentiate)
- [ ] Run cNMF k-selection and program inference on benchmark datasets
- [ ] Update cNMF input path for Huangfu benchmark after pipeline run completes

### Get >= 2 production datasets through CRISPR pipeline (by Mar 26)
- [ ] Gersbach hepatocyte production dataset — in progress
- [ ] Hon cardiomyocyte production dataset — starting soon
- [ ] Run QC on pipeline outputs for both

### Production datasets through cNMF + E-distance (by Apr 23)
- [ ] Run energy distance on >= 2 production datasets
- [ ] Run cNMF on >= 2 production datasets
- [ ] Prepare clean, lightweight input files for jamboree (avoid environment/loading bottlenecks)

## Pre-jamboree: Nice to have

- [ ] E2G links from multiome data in non-perturbed state
- [ ] ChromBPNet models from multiome data in non-perturbed state (Revant)
- [ ] Additional datasets beyond the initial two through the pipeline
- [ ] Touch base with non-coding FG on existing resources

---

## Jamboree working groups and analysis goals

### Working Group 1: Data summarization + QC (Fig 1)
**Point person:** Chikara

- [ ] Extract energy distance outputs; create bar/UpSet plot of # TFs with significant transcriptome effects per lineage
- [ ] Cluster TF perturbations by energy distance in each lineage; label cases with very different effects across lineages
- [ ] Create plots summarizing general statistics (cell counts, gRNA/scRNA MOI, %mito) across datasets
- [ ] Kick off remaining CRISPR pipeline runs for production datasets and analyze QC

### Working Group 2: Gene program annotation + interpretation (Fig 2)
- [ ] Extract gene programs for each production dataset; create heatmap of loading similarities across programs, labeled by lineage
- [ ] Identify programs highly similar across lineages (expect basic cellular processes) vs. lineage-specific programs
- [ ] For shared programs: do the same perturbations change usage in the same way across lineages?
- [ ] Within lineages: identify regulators of each program (which TF perturbations change usage the most)
- [ ] For shared programs: identify shared vs. unique regulators; pull out interesting case studies
- [ ] For case studies: look at top genes and check for enriched TF binding motifs upstream
- [ ] Apply Percoder (Chikara)

### Working Group 3: GRN inference (Fig 3)
- [ ] Apply GRN methods on TFP3 datasets (details TBD)

---

## Key analysis questions to answer at jamboree

1. Which guides are detected and have significant target repression in each system?
2. Which guides/TFs have significant transcriptome-wide effects in each system?
3. What are the inferred downstream targets of each guide/TF in each system?
4. What are the affected gene programs in each system?
5. How do programs compare across lineages — which are lineage-specific vs. lineage-agnostic?
6. Can we build GRN models from these data, and what do they tell us?

---

## Repo organization and documentation
- [ ] Add READMEs for datasets missing documentation (7 datasets)
- [ ] Fill in Stage 1 and Stage 5 documentation in main README.md
- [ ] Clean up stray slurm-*.out files from repo root
- [ ] Add slurm-*.out pattern to .gitignore

## Action items from team meetings
- [ ] **Lucas**: Update pipeline to output TF enrichment automatically if specified
- [ ] **Logan**: Output confidence/uncertainty estimates and cell numbers with Perturbo outputs
- [ ] **Sushama**: Curate list of validated TF-target pairs (internal + published, with direction and source)
- [ ] **Adam**: Implement QC metrics documented in TFP3 Planning doc

---

## Completed
- [x] First CRISPR pipeline runs across benchmark datasets started (Mar 2026)
- [x] Guide metadata uploaded to portal (Mpathi)
- [x] QC standards defined (Feb 2025 team meeting)
- [x] Sub-aims defined and assigned (Feb 2025)
- [x] Repository created and shared (Mar 2025)
- [x] Project documentation reorganized into docs/ (Mar 2026)
