Working Groups
Working Group 1: Overarching data summarization, pipeline & QC cleanup (Topic 1 / Fig. 1) - full group, prior to splitting into subgroups
Goal: Focus on establishing the primary building blocks for all downstream analysis.
Guide detection & target repression
Identify which guides are detected and show target repression in each system (i.e. cis inference outputs)
Summarize and compare general statistics (cell counts, gRNA and scRNA MOI, etc.) and data quality metrics (%mito counts) across datasets
Transcriptome-wide significance
Extract energy distance outputs and create a bar plot (or UpSet plot) of the number of TFs per production dataset that significantly alter the transcriptome
UpSet plot of shared TFs with significant effects across lineages
Cluster TF perturbations by their energy distance in each lineage; label cases where TFs have very different energy distances across lineages
Downstream target inference
Identify inferred downstream trans targets of each guide/TF in each system and look at the distribution of overlap
Technical harmonization
Develop and apply strategies to reduce the impact of technical differences between datasets
Discussion of ways to integrate clonal expansion, doublet removal, and DEG calibration into upstream analysis pipeline

Working Group 2: Gene program annotation & interpretation (Topic 2 / Fig. 2) 
Potential participants: Sara Geraghty, Alexandra Mo, Hazel (Jiahui) Zhao, Ashlesha Gogate
Goal: Focus on biological questions about TF activity across lineages and disease.
Gene program inference and cross-lineage mapping (Topic 2.1)
Extract gene programs for each production-level dataset; create a heatmap of loading similarities across programs, labeled by lineage of origin
Identify programs that are highly similar across lineages (expected to be basic cellular processes) and annotate them
Identify lineage-specific programs and assess whether they correspond to lineage-specific biological processes
(Stretch) Assess whether gene programs can be mapped to in vivo counterparts
TF activity: lineage-specific vs. lineage-agnostic (Topic 2.1a–c)
For programs shared across lineages: which perturbations alter program usage in the same direction? Which TFs contribute to lineage-agnostic vs. lineage-specific programs?
Within lineages: identify regulators of each program (barplot or heatmap of TF perturbations associated with each program)
For shared programs: do they share regulators, or have distinct ones? Pull out interesting case studies
Assess whether any TFs contribute to lineage bifurcations
TF dosage sensitivity (Topic 2.1d)
Applying Percoder to assess perturbation sensitivity of programs across lineages

Working Group 3: Disease & GWAS (Topic 2.2 / Fig. 2)
Potential participants: Sara Geraghty, Sushama Sivakumar, Gary Yang, Hazel Zhao, Ashlesha Gogate 
Goal: Connect TF regulatory activity to human disease.
Identify which TFs regulate disease/GWAS genes in each lineage
For TFs that are disease genes in multiple lineages, assess whether their activity is convergent or divergent
Identify and annotate GWAS variants near important TFs or the regulatory elements upstream of their downstream targets

Working Group 4: GRN inference (Topic 2.3 / Fig. 3)
Potential participants: Adam Klie, Taosha Gao, Weizhou Qian, Sushama Sivakumar, Gary Yang, Denis Torre 
Goal: Apply GRN methods to TFP3 datasets to build causal and mechanistic gene regulatory networks.
Integrate multiome data (E2G linking, ChromBPNet models) with TF-gene regulatory networks from Perturb-seq
Compare TF importance inferred from multiome vs. Perturb-seq; identify which TFs act through direct binding vs. indirect mechanisms
Assess common themes in how disease genes are regulated
Characterize how the structure of TF networks changes across lineages

Working Group 5: TF family case studies (Topic 2.4 / Fig. 4) 
Potential participants: [TBD]
Goal: Deep-dive analysis of specific TFs or TF families with newly implicated roles in lineage differentiation (e.g., ZNF factors).
Identify interesting TF families using outputs from Working Groups 1–3
Pathway analysis of top gene programs for these TFs, with DEG fold-change overlays
Use ATAC-seq data to assess differential accessibility of binding motifs upstream of top regulated genes
Highlight GWAS SNPs in or near these TFs or their regulatory elements linked to disease in the relevant cell type

[Optional] Working Group 6: Predictive modeling (Topic 3 / Fig. 5)
Potential participants: Weizhou Qian, Sid Raghavan
Goal: Design and begin implementation of a simple predictive model trained on uniformly processed outputs from five lineages
Brainstorm model architectures that capitalize on the TFP3 data as a training source, that (for example) are able to more effectively predict context-specific TF perturbation impact
Specifically define inputs, outputs, task, and validation metrics, and potentially begin implementation
Alternatively, apply existing tools from member labs to this data source

