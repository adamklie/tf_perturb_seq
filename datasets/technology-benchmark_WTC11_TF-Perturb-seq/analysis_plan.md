Analyses

Color code: `/Users/adamklie/Desktop/projects/tf_perturb_seq/config/colors`

1. Hon lab comparison (Sceptre and CLEANSER vs internal)
    1. Metrics
    2. Guide UMIs and assignment
    3. Repression efficiency (Sceptre and CLEANSER)

1. Cross technology comparison

1. Metrics
    1. Barplots something like this:
    
    ![image.png](Update/image.png)
    
2. Guide capture (NEED SOME IDEAS HERE)
    1. Cells per guide comparisons?
3. Repression efficiency —> plot this a couple different ways
    1. Barplot of auROC and auPRC: need to grab results from output directories
    2. Simple FC violin plot with medians
    
    ![image.png](Update/image%201.png)
    
    c. Scatter/heatmap of log2FCs: `datasets/technology-benchmark_WTC11_TF-Perturb-seq/bin/1_qc/intended_target_comparison.ipynb` 
    
4. DEGs
    1. Barplot of number of DEGs for positive controls
    2. Normalized by number of cells?
    4. Look at effect sizes for positive control DEGs pairwise comparison