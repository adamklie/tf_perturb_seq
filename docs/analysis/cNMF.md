cNMF pipeline
Gene Network Evaluation Dashboard: https://github.com/EngreitzLab/gene_network_evaluation
Hon Lab cNMF repo: https://github.com/adamklie/tf_perturb_seq/tree/main/scripts/cNMF-TF-perturbseq-CMs-Honlab
Engreitz Lab GPU-accelerated cNMF repo: https://github.com/EngreitzLab/cNMF_benchmarking/tree/main/cNMF_benchmarking_pipeline/Inference/torch-cNMF

Torch-cNMF Handbook: cNMF pipeline handbook

Docker container for torch-cNMF: https://hub.docker.com/repository/docker/igvf/torch-cnmf/general 

Guidelines for uniformly running gene expression program (GEP) inference

Overview:
Discover GEPs in the gene expression matrix from each TF Perturb-seq screen
Preferred discovery method is consensus non-negative matrix factorization (cNMF) using the above Docker container

Running torch-cNMF:
Pull Docker container (Singularity instructions)
singularity pull docker://igvf/torch-cnmf:v01
Convert inference_mudata.h5mu to .h5ad AnnData object with this format and give them informative names. Note that .X must be unnormalized.
>>>import anndata as ad
>>>import mudata as md
>>>from scipy.sparse import csr_matrix, issparse
>>>import scanpy as sc
>>>
>>>mdata = md.read_h5mu("inference_mudata.h5mu")
>>>gene = mdata.mod['gene']
>>>guide = mdata.mod['guide']
>>>adata = ad.AnnData(
    X=gene.X,
    obs=gene.obs,
    var=gene.var
)
>>># Add guide assignment
>>>GA = guide.layers["guide_assignment"]
>>># Has to be dense
>>>if issparse(GA):
    GA = GA.toarray()
>>>adata.obsm["guide_assignment"] = GA
>>>guide_names = guide.var["guide_id"].astype(str).to_numpy()
>>>guide_targets = guide.var["intended_target_name"].astype(str).to_numpy()
>>>
>>>for k in guide.uns.keys():
    adata.uns[k] = guide.uns[k]
>>>for k in mdata.uns.keys():
    adata.uns[k] = mdata.uns[k]


>>>adata.uns["guide_names"]   = guide_names
>>>adata.uns["guide_targets"] = guide_targets
>>>
>>>sc.tl.pca(adata, n_comps=50)
>>>sc.pp.neighbors(adata)
>>>sc.tl.umap(adata)
>>>
>>>adata.write_h5ad("gersbach_ipsc_bridge_10xHTv2_GEX_forcNMF.h5ad")
(Option 1) Within Singularity shell:
singularity shell -B ~/my/path:/mnt torch-cnmf_v01.sif
Singularity> cd /mnt/igvf_cnmf/torch_based_cNMF_docker
Singularity> ls
Singularity> python
(Option 2) Launch an interactive Jupyter notebook session (example on Duke DCC, be sure to change paths):

^Creates a new kernel using the .sif image and allows you to select it in an OnDemand Apptainer on your HPC
Download https://github.com/EngreitzLab/cNMF_benchmarking/blob/main/cNMF_benchmarking_pipeline/Inference/torch-cNMF/JupterNote_Version/torch-cNMF_inference_pipeline.ipynb 
(Option 3) Run cNMF using SLURM script on HPC
Download https://github.com/EngreitzLab/cNMF_benchmarking/tree/main/cNMF_benchmarking_pipeline/Inference/torch-cNMF/Slurm_Version and adapt parameters for your HPC
Across all possibilities, run cNMF with the following parameters: 
Number of iterations: 100
Number of HVGs: 1000
Values of k: 
Benchmark: [5,6,7,8,9,10,11,12,13,14,15,17,19,21,23,25,27,30,35,40,45,50,55,60,70,80,90,100,150,200]
Production: [5,6,7,8,9,10,11,12,13,14,15,17,19,21,23,25,27,30,35,40,45,50,55,60,70,80,90,100,150,200,250,500]
Seed: 123
Species: human
use_GPU: True
*Other initial parameters are left as defaults
Density_thresholds: [0.01, 0.05, 2.0] (in the run_cnmf_consensus step)
NOTE: Run with multiple density thresholds, but set sel_thresh = 2.0 for all datasets.
Update run_name to an informative name (e.g. "ipsc_bridge_HTv2")

Naming and storing results files:
All output files and figures should be stored in Synapse (https://www.synapse.org/Synapse:syn64423137). Locate the folder with the name of the dataset you are handling, and locate/create a directory called “outputs”. Within this directory, create a sub-directory called “cNMF_outputs”. Create subdirectories for “plots” and “files”.
Output figures:
Stability Error Plot (e.g. ipsc_bridge_HTv2.k_selection.png)
Clustergrams and local density histograms for each value of k (e.g. ipsc_bridge_HTv2.clustering.k_<X>.dt_2_0.png, where <X> is replaced with the k value used to generate that plot)
Output files:
A README.txt file with (1) the density threshold you used, (2) the k values you selected (below), and (3) any other notes.
The generated directory with your run_name with its full contents, which includes 3 sub-directories (adata, loading, and prog_data), as well as gene_spectra_score, gene_spectra_tpm, spectra_consensus, usages_consensus, and starcat files. 
IMPORTANT: Link the Synapse page with these outputs into the “Gene Program Discovery” section of this spreadsheet.

Evaluation Step: 
Create a local copy of the following repo: https://github.com/EngreitzLab/cNMF_benchmarking/tree/main/cNMF_benchmarking_pipeline/
> git clone https://github.com/EngreitzLab/cNMF_benchmarking.git
> cd ./cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation
Jupyter version:
Open cNMF_evaluation_pipeline.ipynb as a Jupyter notebook and select torch-cNMF as your kernel (see above).
Update “components” with the values of k run previously, and “sel_thresholds” to match what you selected above.
Update main path, out_dir, and run_name to match your dataset and paths. 
Download guide metadata file and set path to this file as guide_annotation_path. 
Benchmark:benchmark_guide_metadata_v4
Production: production_guide_metadata_v6


Plotting/ K-selection:
Examine the Stability-Error plot from cNMF run. As a group, select the top 10-15 k values that appear to maximize stability and minimize error.
Jupyter version:
Open cNMF_k_selection.ipynb as a Jupyter notebook and select torch-cNMF as your kernel (see above).
Update “K” with the 10-15 selected values of k. 
Update output_directory, run_name, save_folder_name, and eval_folder_name to match your dataset and paths.
Update sel_threhs to match what you selected above.
Run the notebook.
Naming and storing results files:
Save output from plot_stability_error() as error_stability.png.
Save output from plot_enrichment(), which consists of GO by k, genesets by k, and traits by k plots, as enrichment.png.
Save the output of plot_perturbation(), which consists of the number of unique regulators by k and total number of regulators by k plots, as perturbation.png.
Save the output of plot_explained_variance() as explained_variance_by_k.png.
Add all these plots to Synapse as described above.

NOTE: For all datasets, decide on k-selection as a group jointly (clinical review board style)

Additional plot generation (optional)
Program QC plots: program UMAP, program violin plot, program loading correlations, top GO term plot, top loading genes, volcano plot + dot plot + waterfall plot + bar plot for regulated programs per condition of cells, example for one gene
Perturbed-gene plots: gene UMAP, guide UMAP, gene dotplot, gene loading correlations, top loading programs, volcano plot + dot plot + waterfall plot + bar plot for regulated programs per condition of cells,  example for one program
Excel summarization: Integrate mdata + evaluation results information together,  example for K=30, density threshold = 0.4