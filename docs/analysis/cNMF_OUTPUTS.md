Torch cNMF minimal data requirement 
Inference for torch-cNMF 
Minimal required input 
adata object
├── adata/ 
├── X (gene expression count matrix) 
├── obs/
├── cell_ID (obs_names)
├── sample
├── var/
├── gene_names (var_names)
├── uns/
├── guide_names 
├── guide_targets
├── obsm/
├── guide_assignment 
├── X_PCA
├── X_umap

Explanation 
X 
must be unnormalized data with raw mRNA count. It can be filtered (removed cells or genes with low quality). This is because cNMF has its normalization step (TPM normalization). 
obs
cell_ID
Name of cells 
sample
Conditions on cells; if only has 1 condition, still need this column
var
gene_names
Name of genes (make sure to include the names you want to use for plot, ig. Ensembl ID or gene symbol)  
uns
guide_names
Name of the guide, length must be the same as guide_targets
guide_targets
Name of the gene, the guide is targeting to; name here must be the same names as in gene_names in var; length must be the same as guide_names
obsm 
guide_assignment
Cells x guide matrix; guide length must be the same as guide_targets
When set shuffle_cells=True, must have guide_assignment to make sure cells are shuffled with the correct guide assignment per cell
X_PCA, X_umap
Embedding for running UMAP plot
When set shuffle_cells=True, must have X_PCA, X_umap to make sure cells are shuffled with the correct index per cell


Output (one mdata object per k and per density threshold)
├── mdata/ 
├── rna/ 
├── X (gene expression count matrix) 
├── obs/
├── cell_ID (obs_names)
├── sample
├── var/
├── gene_names (var_names)
├── uns/
├── guide_names 
├── guide_targets 
├── obsm/
├── guide_assignment 
├── X_PCA
├── X_umap
├── cNMF/ 
├── X/ (cell x program loading matrix) 
├── obs/ (direct copy from rna.obs)
├── cell_ID (obs_names)
├── sample
├── obsm/ (direct copy from rna.obsm)
├── guide_assignment 
├── X_PCA  
├── X_umap  
├── uns / (direct copy from rna.uns)
├── guide_names 
├── guide_targets 
├── var_names (direct copy from rna.var_names)	
├── varm/
├── loadings (program x gene matrix)
Evaluation
Minimal required input 
mdata object( if inference step is working properly, this is a direct input from inference output) 
├── mdata/ 
├── rna/ 
├── X (gene expression count matrix) 
├── obs/
├── sample
├── uns/
├── guide_names 
├── guide_targets (must match the same naming in var)
├── obsm/
├── guide_assignment 
├── X_PCA
├── X_umap
├── var/
├── gene_names	
├── cNMF/ 
├── X (cell x program loading matrix) 
├── obs/
├── sample
├── obsm/
├── guide_assignment 
├── X_PCA (direct copy from rna mod)
├── X_umap (direct copy from rna mod)
├── uns /
├── guide_names 
├── guide_targets (must match the same naming in var)
├── gene_names	
├── varm/
├── loadings (program x gene matrix)

Required resource files
X_normalized_file: TPM normalized count matrix (cNMF’s output)
gwas_file: GWAS file for trait enrichment test
Loci_file: contains loci for scanning motif in .txt (e.g.: output from scE2G)
Motif_file: motif file in .meme (e.g.: hocomoco_meme.meme)
Seq_file: sequence file in .fa (e.g.:hg38.fa)


Optional resource files
Guide annotation: 
“guide_names” for each guide id 
“guide_targets” for the intended targeting gene
“targeting” indicating True for targeting/False for non-targeting 

Output (those files will be stored in the following example structure when K=30 with cell conditions from D0 to D3)
├── Eval/ 
├── 30_0.4/ 
├── 30_categorical_association_posthoc.csv 
├── 30_categorical_association_results.csv 
├── 30_explained_variance_summary.csv 
├── 30_explained_variance.csv 
├── 30_geneset_enrichment.csv 
├── 30_go_term_enrichment.csv 
├── 30_perturbation_association_d0.csv 
├── 30_perturbation_association_d1.csv 
├── 30_perturbation_association_d2.csv 
├── 30_perturbation_association_d3.csv 
└── 30_trait_enrichment.csv




Plotting + Excel Summarization 
Minimal required input
mdata object( if inference step is working properly, this is a direct input from inference output) 
├── mdata/ 
├── rna/ 
├── X/ (gene expression count matrix) 
├── obs/
├── sample
├── uns/
├── guide_names 
├── guide_targets (must match the same naming in var)
├── obsm/
├── guide_assignment 
├── X_PCA
├── X_umap
├── var/
├── var_names/	
├── cNMF/ 
├── X/ (cell x program loading matrix) 
├── obs/
├── sample
├── obsm/
├── guide_assignment 
├── X_PCA
├── X_umap
├── uns /
├── guide_names 
├── guide_targets (must match the same naming in var)
├── var_names/	
├── varm/
├── loadings (program x gene matrix)

Resources( if evaluation step is working properly, this is a direct input from evaluation output) 
├── Eval/ 
├── 30_0.4/ 
├── 30_categorical_association_posthoc.csv 
├── 30_categorical_association_results.csv 
├── 30_explained_variance_summary.csv 
├── 30_explained_variance.csv 
├── 30_geneset_enrichment.csv 
├── 30_go_term_enrichment.csv 
├── 30_perturbation_association_d0.csv 
├── 30_perturbation_association_d1.csv 
├── 30_perturbation_association_d2.csv 
├── 30_perturbation_association_d3.csv 
└── 30_trait_enrichment.csv

Output( 4 plotting options + 1 excel option)  

K-selection plots: Stability &Error, GO/Genesets/Trait enrichment, perturbation sensitivity, explained variances, program dot plot by conditions 

Compare model plots (with same K): clustermap and boxplots for shared gene, GO/Genesets/Trait enrichment, perturbation sensitivity; coefficient of variances  

Program QC plots: program UMAP, program violin plot, program loading correlations, top GO term plot, top loading genes, volcano plot + dot plot + waterfall plot + bar plot for regulated programs per condition of cells example for one gene

Perturbed-gene plots: gene UMAP, guide UMAP, gene dotplot, gene loading correlations, top loading programs, volcano plot + dot plot + waterfall plot + bar plot for regulated programs per condition of cells  example for one program

Excel summarization: Integrate mdata + evaluation results information together  example for K=30, density threshold = 0.4
