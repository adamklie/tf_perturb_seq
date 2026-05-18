"""Interactive test for the meta-consensus pipeline.

Run as a script or step through cells (#%%) in VS Code / Jupyter.
All config values are inlined here — edit them directly to experiment.
"""

#%% Imports + path
import os
import sys
import numpy as np
# Make the local modules importable when running from anywhere
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else \
       "/oak/stanford/groups/engreitz/Users/ymo/Testing_Scripts/Meta_program/scripts/meta_consensus"
sys.path.insert(0, HERE)

import pandas as pd

from meta_consensus import (
    DatasetSpec,
    build_shared_vs_specific,
    categorize_clusters,
    intersect_genes,
    load_spectra,
    normalize_rows,
    run_kmeans_consensus,
    stack_programs,
    write_category_counts,
)
from plotting import (
    plot_category_breakdown,
    plot_clustergram,
    plot_composition_bar,
    plot_k_selection,
    plot_top_genes_heatmap,
)


#%% Config — edit these to test
DATASETS = [
    DatasetSpec(
        name="Hon_CM",
        inference_dir="/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Hon_CM/Result/050626_1M_CM_torch_halsvar_dataloader/Inference",
        k=100,
        dt="2_0",
    ),
    DatasetSpec(
        name="Huangfu_embryonic",
        inference_dir="/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_embryonic-stemcell/Result/Adam_run/Inference",
        k=100,
        dt="2_0",
    ),
    DatasetSpec(
        name="Huangfu_definitive",
        inference_dir="/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_definitive-endoderm/Result/Adam_run/Inference",
        k=100,
        dt="2_0",
    ),
]

INPUT_FILE_TYPE = "gene_spectra_score"   # gene_spectra_score | spectra | gene_spectra_tpm
NORMALIZATION  = "l2"                    # l2 | row_zscore | row_sum | none
META_K_LIST    = [20, 30, 50, 80, 100]
OUTPUT_DIR     = "/oak/stanford/groups/engreitz/Users/ymo/Testing_Scripts/Meta_program/Result/meta_consensus_test"
KMEANS_N_INIT  = 10
KMEANS_RANDOM_STATE = 1


#%% 1. Load each dataset's gene_spectra file
dfs = []
for spec in DATASETS:
    df = load_spectra(spec, INPUT_FILE_TYPE)
    print(f"{spec.name}: shape={df.shape}  k={spec.k} dt={spec.dt}")
    dfs.append(df)


#%% 2. Intersect genes across datasets
common = intersect_genes(dfs)
print(f"\nCommon genes: {len(common)}")
print(f"  per-dataset: {[d.shape[1] for d in dfs]}")


#%% 3. Stack into one (P_total x G_common) matrix
S = stack_programs(dfs, common)
print(f"Stacked: {S.shape}  (rows = '<dataset>|<program_id>')")
print("Example rows:", S.index[:5].tolist(), "...", S.index[-3:].tolist())


#%% 4. Row-normalize (default L2)
N = normalize_rows(S, NORMALIZATION)
print(f"Normalization: {NORMALIZATION}")
print(f"Row norms (first 3): {(N.values ** 2).sum(axis=1)[:3]}  (should be ~1.0 for l2)")


#%% 5. Run meta-consensus for one k (try k=50 first) — inspect interactively
K = 50
labels, meta_spectra, sil = run_kmeans_consensus(
    N, K, n_init=KMEANS_N_INIT, random_state=KMEANS_RANDOM_STATE
)
print(f"k={K}  silhouette={sil:.4f}")
print(f"meta_spectra shape: {meta_spectra.shape}")
print("labels head:")
print(labels.head())


#%% 6. Categorize clusters (shared_all / shared_partial / specific_<ds>)
dataset_names = [d.name for d in DATASETS]
cluster_summary = categorize_clusters(labels, dataset_names)
shared_vs_specific = build_shared_vs_specific(labels, cluster_summary)

print("Category counts:")
print(cluster_summary["category"].value_counts())
print("\ncluster_summary head:")
print(cluster_summary.head())


#%% 7. Full sweep over META_K_LIST: write outputs + plots
os.makedirs(OUTPUT_DIR, exist_ok=True)
silhouettes = {}

for k in META_K_LIST:
    print(f"\n# k={k}")
    k_dir = os.path.join(OUTPUT_DIR, f"k_{k}")
    os.makedirs(k_dir, exist_ok=True)

    labels, meta_spectra, sil = run_kmeans_consensus(
        N, k, n_init=KMEANS_N_INIT, random_state=KMEANS_RANDOM_STATE
    )
    silhouettes[k] = sil
    print(f"  silhouette = {sil:.4f}")

    cluster_summary = categorize_clusters(labels, dataset_names)
    shared_vs_specific = build_shared_vs_specific(labels, cluster_summary)

    meta_spectra.to_csv(os.path.join(k_dir, f"meta_spectra.k_{k}.median.txt"), sep="\t")
    labels.to_csv(os.path.join(k_dir, f"cluster_assignments.k_{k}.txt"), sep="\t", header=True)
    with open(os.path.join(k_dir, f"silhouette.k_{k}.txt"), "w") as f:
        f.write(f"{sil:.6f}\n")
    cluster_summary.to_csv(os.path.join(k_dir, f"cluster_summary.k_{k}.tsv"), sep="\t", index=False)
    shared_vs_specific.to_csv(os.path.join(k_dir, f"shared_vs_specific.k_{k}.tsv"), sep="\t", index=False)
    write_category_counts(cluster_summary, dataset_names,
                          os.path.join(k_dir, f"category_counts.k_{k}.txt"))

    plot_clustergram(N, labels,
                     os.path.join(k_dir, f"clustering.k_{k}.png"),
                     title=f"Meta-consensus clustergram (k={k}, sil={sil:.3f})")
    plot_composition_bar(cluster_summary, dataset_names,
                         os.path.join(k_dir, f"composition.k_{k}.png"),
                         title=f"Dataset composition per meta-cluster (k={k})")
    plot_category_breakdown(cluster_summary, dataset_names,
                            os.path.join(k_dir, f"category_breakdown.k_{k}.png"),
                            title=f"Meta-cluster categories (k={k})")
    plot_top_genes_heatmap(meta_spectra,
                           os.path.join(k_dir, f"top_genes.k_{k}.png"),
                           n_top=20,
                           title=f"Top genes per meta-program (k={k})")


#%% 8. k-selection curve across all k
pd.Series(silhouettes).sort_index().to_csv(
    os.path.join(OUTPUT_DIR, "k_selection.txt"), sep="\t", header=False
)
plot_k_selection(silhouettes, os.path.join(OUTPUT_DIR, "k_selection.png"))
print(f"\nDone. Outputs in {OUTPUT_DIR}")
