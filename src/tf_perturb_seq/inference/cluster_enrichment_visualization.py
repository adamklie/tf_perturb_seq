#!/usr/bin/env python3
"""
cluster_enrichment_visualization.py

Generates UMAP plots from the outputs of cluster_enrichment.py:

  1. UMAP colored by Leiden cluster
  2. UMAP colored by guide target (each TF shown as a separate color,
     NT cells shown in grey) -- split into pages of up to 12 TFs each
  3. For each cluster: UMAP highlighting cells in that cluster vs rest

Inputs (from cluster_enrichment.py output directory):
  cluster_labels.csv  -- barcode, cluster
  pca_coords.npz      -- X_pca (n_cells x n_pcs), barcodes

UMAP is recomputed from the saved PCA coordinates. No expression data
needs to be loaded.

Usage:
    python cluster_enrichment_visualization.py \\
        --clustering_dir /path/to/clustering_output \\
        --output_dir /path/to/output \\
        --mudata_path /path/to/inference_mudata.h5mu
"""

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

warnings.filterwarnings("ignore")

# --- Default paths ------------------------------------------------------------
_DEFAULT_MUDATA = (
    "/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/"
    "helen_output_poolabcd_prod_v9/pipeline_outputs/inference_mudata.h5mu"
)
_DEFAULT_GUIDE_META = (
    "/hpc/group/gersbachlab/seg95/crispr-pipeline-personal/"
    "example-data/production_guide_metadata_v6.tsv"
)
ENSG_TO_SYM = "/hpc/home/seg95/lab-storage/ref/ensembl_to_symbol.tsv"


# --- CLI ----------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--clustering_dir", required=True,
        help="Directory containing cluster_labels.csv and pca_coords.npz "
             "(output of cluster_enrichment.py)",
    )
    p.add_argument(
        "--output_dir", default=None,
        help="Output directory for plots (default: same as clustering_dir)",
    )
    p.add_argument(
        "--mudata_path", default=_DEFAULT_MUDATA,
        help="Path to inference_mudata.h5mu (needed for guide assignments)",
    )
    p.add_argument(
        "--guide_metadata", default=_DEFAULT_GUIDE_META,
        help="Path to guide metadata TSV",
    )
    p.add_argument(
        "--n_neighbors", type=int, default=15,
        help="k-NN neighbors for UMAP graph (default: 15)",
    )
    p.add_argument(
        "--min_dist", type=float, default=0.3,
        help="UMAP min_dist parameter (default: 0.3)",
    )
    p.add_argument(
        "--downsample", type=int, default=200_000,
        help="Cells to plot (balanced across clusters, default: 200,000)",
    )
    p.add_argument(
        "--rare_thresh", type=int, default=50,
        help="Clusters with fewer cells than this are kept whole during "
             "downsampling (default: 50)",
    )
    p.add_argument(
        "--tfs_per_page", type=int, default=12,
        help="TFs per page in the per-TF UMAP grid (default: 12)",
    )
    p.add_argument(
        "--min_cells_tf", type=int, default=50,
        help="Skip TFs with fewer than this many cells (default: 50)",
    )
    p.add_argument(
        "--random_seed", type=int, default=42,
    )
    p.add_argument(
        "--skip_umap", action="store_true", default=False,
        help="Load precomputed UMAP from umap_coords.npz in clustering_dir",
    )
    return p.parse_args()


# --- Helpers ------------------------------------------------------------------
def build_ensg_to_symbol():
    ensg2sym = {}
    if not os.path.exists(ENSG_TO_SYM):
        return ensg2sym
    ref = pd.read_csv(ENSG_TO_SYM, sep="\t", header=None,
                      names=["versioned_id", "symbol"])
    ref["ensg"] = ref["versioned_id"].str.split(".").str[0]
    return dict(zip(ref["ensg"], ref["symbol"]))


def balanced_downsample(cell_indices, group_labels, target_total, rare_thresh, rng):
    """
    Downsample cell_indices to target_total, stratified by group_labels.
    Rare groups (<= rare_thresh cells) are kept in full.
    Remaining budget is split equally across abundant groups.
    """
    unique_groups = np.unique(group_labels)
    group_to_idx  = {g: cell_indices[group_labels == g] for g in unique_groups}

    rare     = {g: idx for g, idx in group_to_idx.items() if len(idx) <= rare_thresh}
    abundant = {g: idx for g, idx in group_to_idx.items() if len(idx) > rare_thresh}

    n_rare   = sum(len(v) for v in rare.values())
    budget   = target_total - n_rare
    n_abund  = len(abundant)

    selected = [idx for idx in rare.values()]

    if budget <= 0 or n_abund == 0:
        return np.sort(np.concatenate(selected)) if selected else np.array([], dtype=int)

    per_group = budget // n_abund
    sample_sizes = {}
    for g, idx in abundant.items():
        n = min(len(idx), per_group)
        chosen = rng.choice(idx, size=n, replace=False)
        selected.append(chosen)
        sample_sizes[g] = n

    # Redistribute shortfall round-robin
    shortfall   = target_total - sum(len(s) for s in selected)
    selected_set = set(np.concatenate(selected).tolist())
    abund_list   = list(abundant.keys())
    while shortfall > 0:
        added_any = False
        for g in abund_list:
            if shortfall <= 0:
                break
            available = [i for i in abundant[g] if i not in selected_set]
            if available:
                extra = rng.choice(available, size=1)[0]
                selected.append(np.array([extra]))
                selected_set.add(extra)
                shortfall -= 1
                added_any = True
        if not added_any:
            break

    return np.sort(np.concatenate(selected).astype(int))


# --- Step 1: Load cluster labels + PCA ----------------------------------------
def load_clustering(clustering_dir):
    cluster_csv = os.path.join(clustering_dir, "cluster_labels.csv")
    pca_npz     = os.path.join(clustering_dir, "pca_coords.npz")

    if not os.path.exists(cluster_csv):
        raise FileNotFoundError(f"cluster_labels.csv not found in {clustering_dir}")
    if not os.path.exists(pca_npz):
        raise FileNotFoundError(
            f"pca_coords.npz not found in {clustering_dir}. "
            f"Re-run cluster_enrichment.py (without --skip_clustering) "
            f"to regenerate it."
        )

    print(f"  Loading cluster labels: {cluster_csv}")
    cl_df    = pd.read_csv(cluster_csv)
    barcodes = cl_df["barcode"].tolist()
    clusters = cl_df["cluster"].astype(str).values
    print(f"  {len(barcodes):,} cells, {len(np.unique(clusters))} clusters")

    print(f"  Loading PCA coords: {pca_npz}")
    pca_data   = np.load(pca_npz, allow_pickle=True)
    X_pca      = pca_data["X_pca"]
    pca_barcodes = pca_data["barcodes"].astype(str).tolist()

    # Align order (should already match, but be safe)
    if pca_barcodes != barcodes:
        print("  Aligning PCA barcode order to cluster_labels order ...")
        pca_bc_idx = {bc: i for i, bc in enumerate(pca_barcodes)}
        order  = [pca_bc_idx[bc] for bc in barcodes if bc in pca_bc_idx]
        X_pca  = X_pca[order, :]

    return barcodes, clusters, X_pca


# --- Step 2: Compute UMAP -----------------------------------------------------
def compute_umap(X_pca, args, clustering_dir):
    umap_path = os.path.join(clustering_dir, "umap_coords.npz")

    if args.skip_umap and os.path.exists(umap_path):
        print(f"  Loading precomputed UMAP: {umap_path}")
        data = np.load(umap_path)
        return data["X_umap"]

    try:
        import umap as umap_lib
        print(f"  Computing UMAP with umap-learn "
              f"(n_neighbors={args.n_neighbors}, min_dist={args.min_dist}) ...")
        reducer = umap_lib.UMAP(
            n_neighbors=args.n_neighbors,
            min_dist=args.min_dist,
            random_state=args.random_seed,
            verbose=True,
        )
        X_umap = reducer.fit_transform(X_pca)
    except ImportError:
        print("  umap-learn not found, using scanpy UMAP ...")
        import scanpy as sc
        import anndata as ad
        adata_tmp = ad.AnnData(X=np.zeros((X_pca.shape[0], 1)))
        adata_tmp.obsm["X_pca"] = X_pca
        sc.pp.neighbors(adata_tmp, n_neighbors=args.n_neighbors,
                        use_rep="X_pca", random_state=args.random_seed)
        sc.tl.umap(adata_tmp, min_dist=args.min_dist,
                   random_state=args.random_seed)
        X_umap = adata_tmp.obsm["X_umap"]

    np.savez_compressed(umap_path, X_umap=X_umap)
    print(f"  UMAP coords saved: {umap_path}")
    return X_umap


# --- Step 3: Load guide assignments -------------------------------------------
def load_guide_assignments(barcodes, args):
    """Returns cell_target array (one label per cell, parallel to barcodes)."""
    import mudata as md
    import scipy.sparse as sp

    print(f"  Loading guide assignments from {args.mudata_path} ...")
    mdata     = md.read_h5mu(args.mudata_path, backed="r")
    guide_mod = mdata["guide"]
    guide_var = guide_mod.var
    layer     = guide_mod.layers["guide_assignment"]
    guide_barcodes = list(guide_mod.obs_names)
    guide_names    = list(guide_mod.var_names)

    layer_csc = layer.tocsc() if sp.issparse(layer) else sp.csc_matrix(layer)
    mdata.file.close()

    if "targeting" in guide_var.columns:
        raw    = guide_var["targeting"]
        is_tgt = (raw.str.upper().str.strip() == "TRUE") \
                 if raw.dtype == object else raw.astype(bool)
    else:
        is_tgt = pd.Series(True, index=guide_var.index)

    tgt_cols = [i for i, v in enumerate(is_tgt) if v]

    if "intended_target_name" in guide_var.columns:
        col_to_target = {i: guide_var["intended_target_name"].iloc[i]
                         for i in tgt_cols}
    else:
        col_to_target = {i: guide_names[i].split("#")[0] for i in tgt_cols}

    obs_bc_to_pos = {bc: i for i, bc in enumerate(barcodes)}
    n_cells       = len(barcodes)
    cell_target   = np.full(n_cells, "non-targeting", dtype=object)

    for col_i in tgt_cols:
        raw_indices = layer_csc[:, col_i].nonzero()[0]
        target      = col_to_target[col_i]
        for raw_i in raw_indices:
            pos = obs_bc_to_pos.get(guide_barcodes[raw_i], -1)
            if pos >= 0:
                cell_target[pos] = target

    ensg2sym     = build_ensg_to_symbol()
    target_syms  = np.array([
        ensg2sym.get(t, t) if t != "non-targeting" else "non-targeting"
        for t in cell_target
    ], dtype=object)

    n_tgt = (cell_target != "non-targeting").sum()
    n_nt  = (cell_target == "non-targeting").sum()
    print(f"  Targeting cells: {n_tgt:,}  |  NT/unassigned: {n_nt:,}")

    return target_syms


# --- Plot helpers -------------------------------------------------------------
def scatter(ax, xy, c, cmap=None, vmin=None, vmax=None,
            color=None, s=1.5, alpha=0.5, rasterized=True, **kw):
    if color is not None:
        ax.scatter(xy[:, 0], xy[:, 1], c=color, s=s, alpha=alpha,
                   linewidths=0, rasterized=rasterized, **kw)
    else:
        ax.scatter(xy[:, 0], xy[:, 1], c=c, cmap=cmap, vmin=vmin, vmax=vmax,
                   s=s, alpha=alpha, linewidths=0, rasterized=rasterized, **kw)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top","right","left","bottom"]].set_visible(False)


def cluster_colormap(n):
    """Return n distinct colors for cluster labels."""
    if n <= 10:
        return plt.cm.tab10.colors[:n]
    elif n <= 20:
        return plt.cm.tab20.colors[:n]
    else:
        return [plt.cm.hsv(i / n) for i in range(n)]


# --- Analysis 1: UMAP by cluster ----------------------------------------------
def plot_umap_by_cluster(X_umap, clusters, sampled_idx, out_dir):
    print("\n  Plotting UMAP colored by Leiden cluster ...")
    unique_cls = sorted(np.unique(clusters), key=int)
    colors     = cluster_colormap(len(unique_cls))
    cl_to_col  = {cl: colors[i] for i, cl in enumerate(unique_cls)}

    xy   = X_umap[sampled_idx]
    cl_s = clusters[sampled_idx]

    fig, ax = plt.subplots(figsize=(10, 8))
    for cl in unique_cls:
        mask = cl_s == cl
        ax.scatter(xy[mask, 0], xy[mask, 1],
                   c=[cl_to_col[cl]], s=1.5, alpha=0.5,
                   linewidths=0, rasterized=True, label=f"C{cl} (n={mask.sum():,})")
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top","right","left","bottom"]].set_visible(False)
    ax.set_title(f"Leiden clusters  (n={len(sampled_idx):,} cells, "
                 f"{len(unique_cls)} clusters)",
                 fontsize=13, fontweight="bold")

    # Legend outside -- two columns if many clusters
    ncol = 2 if len(unique_cls) > 10 else 1
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1),
              fontsize=7, frameon=False, ncol=ncol,
              markerscale=4, handlelength=0.8)

    plt.tight_layout()
    out = os.path.join(out_dir, "umap_leiden_clusters.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# --- Analysis 2: Per-cluster highlight UMAPs ----------------------------------
def plot_umap_per_cluster(X_umap, clusters, sampled_idx, out_dir):
    print("\n  Plotting per-cluster highlight UMAPs ...")
    unique_cls = sorted(np.unique(clusters), key=int)
    colors     = cluster_colormap(len(unique_cls))
    cl_to_col  = {cl: colors[i] for i, cl in enumerate(unique_cls)}

    xy   = X_umap[sampled_idx]
    cl_s = clusters[sampled_idx]

    ncols = min(4, len(unique_cls))
    nrows = int(np.ceil(len(unique_cls) / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(5 * ncols, 4 * nrows),
                             squeeze=False)

    for i, cl in enumerate(unique_cls):
        ax   = axes[i // ncols][i % ncols]
        mask = cl_s == cl
        # Background
        ax.scatter(xy[~mask, 0], xy[~mask, 1],
                   c=["#dddddd"], s=0.8, alpha=0.3,
                   linewidths=0, rasterized=True)
        # Highlighted cluster
        ax.scatter(xy[mask, 0], xy[mask, 1],
                   c=[cl_to_col[cl]], s=1.8, alpha=0.7,
                   linewidths=0, rasterized=True)
        ax.set_title(f"Cluster {cl}  (n={mask.sum():,})", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        ax.spines[["top","right","left","bottom"]].set_visible(False)

    for j in range(len(unique_cls), nrows * ncols):
        axes[j // ncols][j % ncols].set_visible(False)

    plt.suptitle("Per-cluster UMAP highlights", fontsize=13, fontweight="bold")
    plt.tight_layout()
    out = os.path.join(out_dir, "umap_per_cluster_highlights.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# --- Analysis 3: Per-TF UMAPs (paginated) -------------------------------------
def plot_umap_per_tf(X_umap, clusters, target_syms, sampled_idx,
                     tfs_per_page, min_cells_tf, out_dir):
    print("\n  Plotting per-TF UMAP pages ...")
    unique_cls = sorted(np.unique(clusters), key=int)
    colors     = cluster_colormap(len(unique_cls))
    cl_to_col  = {cl: colors[i] for i, cl in enumerate(unique_cls)}

    xy      = X_umap[sampled_idx]
    cl_s    = clusters[sampled_idx]
    tgt_s   = target_syms[sampled_idx]

    # Count cells per TF in the full (unsampled) dataset
    all_tfs = sorted([t for t in np.unique(target_syms)
                      if t != "non-targeting"
                      and (target_syms == t).sum() >= min_cells_tf])
    print(f"  {len(all_tfs):,} TFs with >= {min_cells_tf} cells")

    n_pages = int(np.ceil(len(all_tfs) / tfs_per_page))
    ncols   = min(4, tfs_per_page)
    nrows   = int(np.ceil(tfs_per_page / ncols))

    for page in range(n_pages):
        page_tfs = all_tfs[page * tfs_per_page: (page + 1) * tfs_per_page]
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(5 * ncols, 4 * nrows),
                                 squeeze=False)

        for i, tf in enumerate(page_tfs):
            ax   = axes[i // ncols][i % ncols]
            mask = tgt_s == tf

            # Background in grey
            ax.scatter(xy[~mask, 0], xy[~mask, 1],
                       c=["#dddddd"], s=0.8, alpha=0.25,
                       linewidths=0, rasterized=True)
            # TF cells colored by their cluster
            for cl in unique_cls:
                cl_mask = mask & (cl_s == cl)
                if cl_mask.sum() == 0:
                    continue
                ax.scatter(xy[cl_mask, 0], xy[cl_mask, 1],
                           c=[cl_to_col[cl]], s=2.5, alpha=0.8,
                           linewidths=0, rasterized=True,
                           label=f"C{cl}")
            n_shown = mask.sum()
            ax.set_title(f"{tf}  (n={n_shown:,})", fontsize=8, fontweight="bold")
            ax.set_xticks([]); ax.set_yticks([])
            ax.spines[["top","right","left","bottom"]].set_visible(False)

        for j in range(len(page_tfs), nrows * ncols):
            axes[j // ncols][j % ncols].set_visible(False)

        plt.suptitle(
            f"TF perturbation cells on UMAP (colored by cluster)\n"
            f"Page {page + 1} of {n_pages}",
            fontsize=11, fontweight="bold",
        )
        plt.tight_layout()
        out = os.path.join(out_dir, f"umap_per_tf_page{page+1:03d}.png")
        plt.savefig(out, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {out}  ({len(page_tfs)} TFs)")


# --- Main ---------------------------------------------------------------------
def main():
    args       = parse_args()
    clust_dir  = args.clustering_dir
    out_dir    = args.output_dir or clust_dir
    os.makedirs(out_dir, exist_ok=True)
    rng        = np.random.default_rng(args.random_seed)

    print(f"Clustering dir   : {clust_dir}")
    print(f"Output dir       : {out_dir}")
    print(f"Downsample target: {args.downsample:,}")

    # --- Load clustering outputs ----------------------------------------------
    print(f"\n{'='*60}")
    print("STEP 1: Loading cluster labels and PCA coords")
    print(f"{'='*60}")
    barcodes, clusters, X_pca = load_clustering(clust_dir)

    # --- Compute UMAP from PCA ------------------------------------------------
    print(f"\n{'='*60}")
    print("STEP 2: Computing UMAP")
    print(f"{'='*60}")
    X_umap = compute_umap(X_pca, args, clust_dir)
    del X_pca  # free memory

    # --- Load guide assignments -----------------------------------------------
    print(f"\n{'='*60}")
    print("STEP 3: Loading guide assignments")
    print(f"{'='*60}")
    target_syms = load_guide_assignments(barcodes, args)

    # --- Balanced downsample for plotting -------------------------------------
    n_cells     = len(barcodes)
    all_idx     = np.arange(n_cells)
    target_total = min(args.downsample, n_cells)

    print(f"\n  Balanced downsample: {n_cells:,} -> {target_total:,} cells "
          f"(stratified by cluster) ...")
    sampled_idx = balanced_downsample(
        all_idx, clusters, target_total, args.rare_thresh, rng
    )
    print(f"  Sampled {len(sampled_idx):,} cells across "
          f"{len(np.unique(clusters[sampled_idx]))} clusters")

    # --- Plots ----------------------------------------------------------------
    print(f"\n{'='*60}")
    print("STEP 4: Generating plots")
    print(f"{'='*60}")

    plot_umap_by_cluster(X_umap, clusters, sampled_idx, out_dir)
    plot_umap_per_cluster(X_umap, clusters, sampled_idx, out_dir)
    plot_umap_per_tf(X_umap, clusters, target_syms, sampled_idx,
                     args.tfs_per_page, args.min_cells_tf, out_dir)

    print("\n=== Complete ===")


if __name__ == "__main__":
    main()