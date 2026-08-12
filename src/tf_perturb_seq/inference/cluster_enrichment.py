#!/usr/bin/env python3
"""
cluster_enrichment.py

Leiden-clusters all ~1M cells from the TF perturb-seq screen using the gene
expression modality, then tests whether cells assigned to each TF perturbation
are enriched or depleted in particular clusters compared to a background.

Workflow:
  1. Load gene modality from MuData (backed mode)
  2. Select 2,000 HVGs, scale, PCA (50 components)
  3. Build k-NN graph and run Leiden clustering
  4. Load guide assignments from guide modality
  5. For each TF with >= min_cells cells, test enrichment per cluster using
     Fisher's exact test (one-sided, separately for enrichment and depletion)
  6. BH correction across all TF-cluster pairs
  7. Output heatmap of enrichment scores + ranked tables

Background modes (--background):
  "nt"  : compare TF-assigned cells to non-targeting cells only
           (tests cell-intrinsic effects of each perturbation)
  "all" : compare TF-assigned cells to all other cells
           (tests whether TF cells are over/underrepresented vs global frequency)

Memory strategy:
  - Gene matrix is loaded in backed mode; only HVG columns are extracted
    into memory (n_cells x 2000 floats ~ 16GB for 1M cells)
  - Neighbor graph and Leiden clustering work in embedding space (50-dim PCA)
  - Guide assignment matrix is CSC-sparse; loaded once and kept sparse

Usage:
    python cluster_enrichment.py \\
        --output_dir /path/to/output \\
        --background nt            # or "all"

    # With custom parameters:
    python cluster_enrichment.py \\
        --output_dir /path/to/output \\
        --background nt \\
        --resolution 0.5 \\
        --n_hvgs 2000 \\
        --n_pcs 50 \\
        --n_neighbors 15 \\
        --min_cells 50 \\
        --fdr_thresh 0.05
"""

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

warnings.filterwarnings("ignore")

# --- Paths --------------------------------------------------------------------
MUDATA_PATH = (
    "/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/"
    "helen_output_poolabcd_prod_v9/pipeline_outputs/inference_mudata.h5mu"
)
ENSG_TO_SYM = "/hpc/home/seg95/lab-storage/ref/ensembl_to_symbol.tsv"


# --- CLI ----------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--output_dir",  required=True)
    p.add_argument("--background",  choices=["nt", "all"], default="nt",
                   help="Background for enrichment test: 'nt' (non-targeting cells) "
                        "or 'all' (all other cells). Default: nt")
    p.add_argument("--resolution",  type=float, default=0.5,
                   help="Leiden resolution (higher = more clusters). Default: 0.5")
    p.add_argument("--n_hvgs",      type=int,   default=2000,
                   help="Number of highly variable genes for PCA. Default: 2000")
    p.add_argument("--n_pcs",       type=int,   default=50,
                   help="Number of PCA components. Default: 50")
    p.add_argument("--n_neighbors", type=int,   default=15,
                   help="k-NN neighbors for graph construction. Default: 15")
    p.add_argument("--min_cells",   type=int,   default=50,
                   help="Skip TFs with fewer than this many cells. Default: 50")
    p.add_argument("--fdr_thresh",  type=float, default=0.05,
                   help="BH FDR threshold for enrichment. Default: 0.05")
    p.add_argument("--top_n_tfs",   type=int,   default=40,
                   help="TFs to show in heatmap (ranked by max enrichment). "
                        "Default: 40")
    p.add_argument("--random_seed", type=int,   default=42)
    p.add_argument(
        "--skip_clustering", action="store_true", default=False,
        help="Skip PCA/clustering and load cluster labels from "
             "cluster_labels.csv in output_dir (for re-running enrichment "
             "with different parameters without recomputing clusters)."
    )
    return p.parse_args()


# --- Helpers ------------------------------------------------------------------
def build_ensg_to_symbol():
    ensg2sym = {}
    if not os.path.exists(ENSG_TO_SYM):
        print(f"  WARNING: {ENSG_TO_SYM} not found")
        return ensg2sym
    ref = pd.read_csv(ENSG_TO_SYM, sep="\t", header=None,
                      names=["versioned_id", "symbol"])
    ref["ensg"] = ref["versioned_id"].str.split(".").str[0]
    return dict(zip(ref["ensg"], ref["symbol"]))


def bh_correct(pvals):
    pvals  = np.asarray(pvals, dtype=float)
    finite = np.isfinite(pvals)
    padj   = np.full(len(pvals), np.nan)
    pf     = pvals[finite]
    m      = finite.sum()
    if m == 0:
        return padj
    try:
        from statsmodels.stats.multitest import multipletests
        _, adj, _, _ = multipletests(pf, method="fdr_bh")
    except ImportError:
        from scipy.stats import rankdata
        ranks = rankdata(pf, method="ordinal")
        adj   = np.minimum(
            np.minimum.accumulate((pf * m / ranks)[::-1])[::-1], 1.0
        )
    padj[finite] = adj
    return padj


# --- Step 1: Cluster cells ----------------------------------------------------
def cluster_cells(args):
    """
    Load gene modality, select HVGs, run PCA + Leiden clustering.
    Returns (obs_names, cluster_labels) as parallel arrays.
    """
    import mudata as md
    import scanpy as sc
    import anndata as ad

    print(f"\n{'='*60}")
    print("STEP 1: Clustering cells")
    print(f"{'='*60}")
    print(f"  Loading gene modality (backed) from {MUDATA_PATH} ...")
    mdata = md.read_h5mu(MUDATA_PATH, backed="r")
    gene  = mdata["gene"]
    n_obs = gene.n_obs
    n_vars = gene.n_vars
    print(f"  Gene modality: {n_obs:,} cells x {n_vars:,} genes (backed)")

    # --- Select HVGs on a subsample, then apply to full dataset ---------------
    # Loading 1M x 11k into memory is very slow. Instead:
    #   1. Sample 100k cells to estimate HVGs (seurat_v3 flavor)
    #   2. Extract only the HVG columns for all 1M cells
    # This avoids holding the full matrix in memory at any point.
    hvg_subsample = min(100_000, n_obs)
    print(f"  Estimating HVGs on {hvg_subsample:,}-cell subsample ...")
    rng         = np.random.default_rng(args.random_seed)
    sub_idx     = np.sort(rng.choice(n_obs, size=hvg_subsample, replace=False))
    X_sub       = gene.X[sub_idx, :]
    if not sp.issparse(X_sub):
        X_sub = sp.csr_matrix(X_sub)

    # Filter zero-count genes in subsample
    gene_counts    = np.asarray(X_sub.sum(axis=0)).flatten()
    expressed_mask = gene_counts > 0
    n_expressed    = expressed_mask.sum()
    print(f"  {n_expressed:,} / {n_vars:,} genes expressed in subsample")

    adata_sub = ad.AnnData(
        X   = X_sub[:, expressed_mask],
        var = gene.var.iloc[expressed_mask].copy(),
    )
    del X_sub

    n_hvgs = min(args.n_hvgs, n_expressed)
    sc.pp.highly_variable_genes(adata_sub, n_top_genes=n_hvgs, flavor="seurat_v3")
    hvg_mask_in_expressed = adata_sub.var["highly_variable"].values
    # Map back to original gene indices
    expressed_indices = np.where(expressed_mask)[0]
    hvg_indices       = expressed_indices[hvg_mask_in_expressed]
    print(f"  Selected {len(hvg_indices):,} HVGs")
    del adata_sub

    # --- Load only HVG columns for all cells ----------------------------------
    print(f"  Loading {len(hvg_indices):,} HVG columns for all {n_obs:,} cells ...")
    # Extract column-by-column to avoid peak memory of full matrix
    # For a CSR matrix backed on disk, column slicing triggers row iteration --
    # more efficient to load in chunks of columns.
    CHUNK = 200
    hvg_chunks = []
    for start in range(0, len(hvg_indices), CHUNK):
        cols = hvg_indices[start:start + CHUNK]
        chunk = gene.X[:, cols]
        if not sp.issparse(chunk):
            chunk = sp.csr_matrix(chunk)
        hvg_chunks.append(chunk)
        if (start // CHUNK) % 5 == 0:
            print(f"    Loaded {min(start + CHUNK, len(hvg_indices)):,} / "
                  f"{len(hvg_indices):,} HVG columns ...")
    X_hvg = sp.hstack(hvg_chunks, format="csr")
    del hvg_chunks
    print(f"  HVG matrix shape: {X_hvg.shape}")

    mdata.file.close()

    hvg_var = gene.var.iloc[hvg_indices].copy()
    adata_hvg = ad.AnnData(
        X   = X_hvg,
        obs = gene.obs.copy(),
        var = hvg_var,
    )
    del X_hvg

    # --- PCA ------------------------------------------------------------------
    # Use randomized SVD -- much faster than arpack at 1M cells
    print(f"  Scaling and running PCA ({args.n_pcs} components, randomized SVD) ...")
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=args.n_pcs, svd_solver="randomized",
              random_state=args.random_seed)
    print(f"  PCA complete")

    # --- k-NN graph + Leiden --------------------------------------------------
    # Use pynndescent for approximate k-NN -- much faster than exact at 1M cells
    print(f"  Building k-NN graph (k={args.n_neighbors}, transformer=pynndescent) ...")
    try:
        sc.pp.neighbors(adata_hvg, n_neighbors=args.n_neighbors,
                        n_pcs=args.n_pcs, method="pynndescent",
                        random_state=args.random_seed)
        print(f"  Used pynndescent (approximate k-NN)")
    except Exception as e:
        print(f"  pynndescent failed ({e}), falling back to exact k-NN ...")
        sc.pp.neighbors(adata_hvg, n_neighbors=args.n_neighbors,
                        n_pcs=args.n_pcs, random_state=args.random_seed)

    print(f"  Running Leiden clustering (resolution={args.resolution}) ...")
    sc.tl.leiden(adata_hvg, resolution=args.resolution,
                 random_state=args.random_seed)

    cluster_labels = adata_hvg.obs["leiden"].values.astype(str)
    obs_names      = list(adata_hvg.obs_names)
    n_clusters     = len(np.unique(cluster_labels))
    print(f"  Found {n_clusters} clusters")
    for cl in sorted(np.unique(cluster_labels), key=int):
        n = (cluster_labels == cl).sum()
        print(f"    Cluster {cl}: {n:,} cells ({100*n/len(cluster_labels):.1f}%)")

    # Save PCA embedding so visualization script can compute UMAP without
    # re-running the full clustering pipeline.
    pca_coords = adata_hvg.obsm["X_pca"]
    return obs_names, cluster_labels, pca_coords


# --- Step 2: Load guide assignments ------------------------------------------
def load_guide_assignments(obs_names):
    """
    Load guide assignment matrix and build:
      - cell_to_target: array of intended_target_name per cell (or "non-targeting")
      - target_to_cells: dict of target -> list of cell indices (in obs_names order)
      - nt_cell_indices: array of cell indices assigned to NT guides
    """
    import mudata as md

    print(f"\n{'='*60}")
    print("STEP 2: Loading guide assignments")
    print(f"{'='*60}")

    mdata     = md.read_h5mu(MUDATA_PATH, backed="r")
    guide_mod = mdata["guide"]
    guide_var = guide_mod.var

    layer          = guide_mod.layers["guide_assignment"]
    guide_barcodes = list(guide_mod.obs_names)
    guide_names    = list(guide_mod.var_names)

    print(f"  Guide matrix: {layer.shape[0]:,} cells x {layer.shape[1]:,} guides")

    # Convert to CSC for fast per-guide column slicing
    if not sp.issparse(layer):
        layer_csc = sp.csc_matrix(layer)
    else:
        layer_csc = layer.tocsc()
    mdata.file.close()

    # Targeting flag
    if "targeting" in guide_var.columns:
        raw    = guide_var["targeting"]
        is_tgt = (raw.str.upper().str.strip() == "TRUE") \
                 if raw.dtype == object else raw.astype(bool)
    else:
        is_tgt = pd.Series(True, index=guide_var.index)

    tgt_cols = [i for i, v in enumerate(is_tgt) if v]
    nt_cols  = [i for i, v in enumerate(is_tgt) if not v]

    # intended_target_name per targeting guide
    if "intended_target_name" in guide_var.columns:
        col_to_target = {i: guide_var["intended_target_name"].iloc[i]
                         for i in tgt_cols}
    else:
        col_to_target = {i: guide_names[i].split("#")[0] for i in tgt_cols}

    # Align to obs_names order from clustering
    bc_to_idx = {bc: i for i, bc in enumerate(guide_barcodes)}
    obs_order = [bc_to_idx[bc] for bc in obs_names if bc in bc_to_idx]
    n_cells   = len(obs_names)

    # Map each cell to its guide target
    cell_target = np.full(n_cells, "non-targeting", dtype=object)
    for col_i in tgt_cols:
        assigned = layer_csc[:, col_i].nonzero()[0]
        # Remap to obs_order position
        for raw_idx in assigned:
            bc = guide_barcodes[raw_idx]
            if bc in bc_to_idx:
                obs_i = obs_names.index(bc) if bc in obs_names else -1
                # Use a faster lookup
        # Rebuild using obs_order
    # Faster approach: build obs_names index
    obs_bc_to_pos = {bc: i for i, bc in enumerate(obs_names)}

    cell_target = np.full(n_cells, "non-targeting", dtype=object)
    for col_i in tgt_cols:
        col   = layer_csc[:, col_i]
        raw_indices = col.nonzero()[0]
        target = col_to_target[col_i]
        for raw_i in raw_indices:
            bc = guide_barcodes[raw_i]
            pos = obs_bc_to_pos.get(bc, -1)
            if pos >= 0:
                cell_target[pos] = target

    # NT cells
    nt_mask = np.zeros(n_cells, dtype=bool)
    for col_i in nt_cols:
        raw_indices = layer_csc[:, col_i].nonzero()[0]
        for raw_i in raw_indices:
            bc = guide_barcodes[raw_i]
            pos = obs_bc_to_pos.get(bc, -1)
            if pos >= 0:
                nt_mask[pos] = True
    nt_cell_indices = np.where(nt_mask)[0]

    # Build target -> cell index list
    target_to_cells = {}
    for i, tgt in enumerate(cell_target):
        target_to_cells.setdefault(tgt, []).append(i)

    n_targets = len([t for t in target_to_cells if t != "non-targeting"])
    print(f"  {n_targets:,} unique targeting elements")
    print(f"  NT cells: {len(nt_cell_indices):,}")
    print(f"  Unassigned / non-targeting cells: "
          f"{(cell_target == 'non-targeting').sum():,}")

    return cell_target, target_to_cells, nt_cell_indices


# --- Step 3: Enrichment testing -----------------------------------------------
def test_enrichment(cluster_labels, cell_target, target_to_cells,
                    nt_cell_indices, args, ensg2sym):
    """
    For each TF with >= min_cells cells, test enrichment/depletion in each
    cluster using Fisher's exact test.

    Returns a DataFrame with columns:
      target_symbol, cluster, odds_ratio, pval_enrich, pval_deplete,
      padj_enrich, padj_deplete, n_tf_in_cluster, n_tf_total,
      n_bg_in_cluster, n_bg_total
    """
    from scipy.stats import fisher_exact

    print(f"\n{'='*60}")
    print("STEP 3: Enrichment testing")
    print(f"{'='*60}")
    print(f"  Background mode: {args.background}")

    clusters     = sorted(np.unique(cluster_labels), key=int)
    n_clusters   = len(clusters)
    cluster_idx  = {cl: np.where(cluster_labels == cl)[0] for cl in clusters}

    # Background cell indices
    if args.background == "nt":
        bg_indices = nt_cell_indices
        bg_label   = "NT cells"
    else:
        bg_indices = np.arange(len(cluster_labels))
        bg_label   = "all cells"
    print(f"  Background: {len(bg_indices):,} {bg_label}")

    rows = []
    targets = [t for t in target_to_cells
               if t != "non-targeting"
               and len(target_to_cells[t]) >= args.min_cells]
    print(f"  Testing {len(targets):,} TFs (>= {args.min_cells} cells) "
          f"x {n_clusters} clusters ...")

    for tgt in targets:
        tf_cells = np.array(target_to_cells[tgt])
        n_tf     = len(tf_cells)
        tf_set   = set(tf_cells.tolist())

        # Background: exclude TF cells from background
        if args.background == "nt":
            bg_set = set(bg_indices.tolist())
        else:
            bg_set = set(bg_indices.tolist()) - tf_set
        n_bg = len(bg_set)

        sym = ensg2sym.get(tgt, tgt)

        for cl in clusters:
            cl_set = set(cluster_idx[cl].tolist())

            a = len(tf_set  & cl_set)           # TF in cluster
            b = len(tf_set  - cl_set)           # TF not in cluster
            c = len(bg_set  & cl_set)           # BG in cluster
            d = len(bg_set  - cl_set)           # BG not in cluster

            if a + b == 0 or c + d == 0:
                continue

            _, pval_enrich  = fisher_exact([[a, b], [c, d]],
                                           alternative="greater")
            _, pval_deplete = fisher_exact([[a, b], [c, d]],
                                           alternative="less")

            # Log2 odds ratio (with pseudocount)
            p_tf = (a + 0.5) / (n_tf + 1)
            p_bg = (c + 0.5) / (n_bg + 1)
            lor  = np.log2(p_tf / p_bg)

            rows.append({
                "intended_target_name": tgt,
                "target_symbol":        sym,
                "cluster":              cl,
                "log2_odds_ratio":      lor,
                "pval_enrich":          pval_enrich,
                "pval_deplete":         pval_deplete,
                "n_tf_in_cluster":      a,
                "n_tf_total":           n_tf,
                "n_bg_in_cluster":      c,
                "n_bg_total":           n_bg,
                "pct_tf_in_cluster":    100 * a / n_tf if n_tf > 0 else 0,
                "pct_bg_in_cluster":    100 * c / n_bg if n_bg > 0 else 0,
            })

    df = pd.DataFrame(rows)
    print(f"  {len(df):,} TF-cluster pairs tested")

    # BH correction across all pairs (separately for enrichment and depletion)
    df["padj_enrich"]  = bh_correct(df["pval_enrich"].values)
    df["padj_deplete"] = bh_correct(df["pval_deplete"].values)

    n_sig_e = (df["padj_enrich"]  < args.fdr_thresh).sum()
    n_sig_d = (df["padj_deplete"] < args.fdr_thresh).sum()
    print(f"  Significant enrichments:  {n_sig_e:,} (BH FDR < {args.fdr_thresh})")
    print(f"  Significant depletions:   {n_sig_d:,} (BH FDR < {args.fdr_thresh})")

    return df, clusters


# --- Step 4: Output -----------------------------------------------------------
def save_results(df, clusters, args, out_dir):
    print(f"\n{'='*60}")
    print("STEP 4: Saving results")
    print(f"{'='*60}")

    # Full results table
    csv_path = os.path.join(out_dir, "cluster_enrichment_results.csv")
    df.sort_values(["target_symbol", "cluster"]).to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    # Significant hits only
    sig = df[(df["padj_enrich"]  < args.fdr_thresh) |
             (df["padj_deplete"] < args.fdr_thresh)].copy()
    sig_csv = os.path.join(out_dir, "cluster_enrichment_significant.csv")
    sig.sort_values("padj_enrich").to_csv(sig_csv, index=False)
    print(f"  Saved: {sig_csv}")

    # --- Heatmap: log2 odds ratio, top TFs by max |LOR| ----------------------
    # Pivot to TF x cluster matrix
    lor_mat = df.pivot_table(
        index="target_symbol", columns="cluster",
        values="log2_odds_ratio", aggfunc="mean",
    )
    # Sort clusters numerically
    lor_mat = lor_mat[[c for c in sorted(clusters, key=int)
                       if c in lor_mat.columns]]

    # Select top TFs by maximum absolute LOR across any cluster
    max_lor = lor_mat.abs().max(axis=1)
    top_tfs = max_lor.nlargest(args.top_n_tfs).index
    lor_top = lor_mat.loc[lor_mat.index.isin(top_tfs)]

    # Also mask cells that are not significant (grey them out)
    sig_mask_e = df.pivot_table(
        index="target_symbol", columns="cluster",
        values="padj_enrich", aggfunc="min",
    ).reindex(index=lor_top.index, columns=lor_top.columns)
    sig_mask_d = df.pivot_table(
        index="target_symbol", columns="cluster",
        values="padj_deplete", aggfunc="min",
    ).reindex(index=lor_top.index, columns=lor_top.columns)
    is_sig = (sig_mask_e < args.fdr_thresh) | (sig_mask_d < args.fdr_thresh)

    # Apply NaN to non-significant cells so they render as grey
    lor_plot = lor_top.copy().astype(float)
    lor_plot[~is_sig] = np.nan

    vmax   = np.nanpercentile(lor_top.values.astype(float), 95)
    vmax   = max(vmax, 0.5)
    n_tfs  = len(lor_top)
    n_cl   = len(lor_top.columns)

    fig, ax = plt.subplots(figsize=(max(8, n_cl * 0.8), max(6, n_tfs * 0.3)))
    cmap    = plt.cm.RdBu_r.copy()
    cmap.set_bad("#dddddd")   # grey for non-significant

    im = ax.imshow(
        lor_plot.values.astype(float),
        aspect="auto", cmap=cmap,
        vmin=-vmax, vmax=vmax,
        interpolation="nearest",
    )
    plt.colorbar(im, ax=ax, label="log2 odds ratio", shrink=0.6, pad=0.02)

    ax.set_xticks(range(n_cl))
    ax.set_xticklabels([f"C{c}" for c in lor_top.columns], fontsize=9)
    ax.set_yticks(range(n_tfs))
    ax.set_yticklabels(lor_top.index, fontsize=8)
    ax.set_xlabel("Leiden cluster", fontsize=11)
    ax.set_ylabel("TF perturbation", fontsize=11)
    ax.set_title(
        f"TF cluster enrichment  (background: {args.background}, "
        f"resolution: {args.resolution})\n"
        f"Top {args.top_n_tfs} TFs by |log2 OR|  |  "
        f"grey = not significant (BH FDR >= {args.fdr_thresh})",
        fontsize=11, fontweight="bold",
    )
    plt.tight_layout()
    heatmap_path = os.path.join(out_dir, "cluster_enrichment_heatmap.png")
    plt.savefig(heatmap_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {heatmap_path}")

    # --- Per-cluster barplot: which TFs are most enriched in each cluster -----
    n_cols = min(4, len(clusters))
    n_rows = int(np.ceil(len(clusters) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(6 * n_cols, 4 * n_rows),
                             squeeze=False)
    for i, cl in enumerate(sorted(clusters, key=int)):
        ax  = axes[i // n_cols][i % n_cols]
        sub = df[(df["cluster"] == cl) &
                 (df["padj_enrich"] < args.fdr_thresh)].copy()
        sub = sub.nlargest(15, "log2_odds_ratio")
        if len(sub) == 0:
            ax.set_title(f"Cluster {cl}\n(no significant enrichments)", fontsize=9)
            ax.axis("off")
            continue
        sub_s = sub.sort_values("log2_odds_ratio", ascending=True)
        ax.barh(sub_s["target_symbol"], sub_s["log2_odds_ratio"],
                color="#d62728", edgecolor="white", linewidth=0.3)
        ax.set_title(f"Cluster {cl}  (n={len(sub):,} sig.)", fontsize=9)
        ax.set_xlabel("log2 odds ratio", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=7)

    # Hide unused axes
    for j in range(len(clusters), n_rows * n_cols):
        axes[j // n_cols][j % n_cols].set_visible(False)

    plt.suptitle(
        f"Top enriched TFs per cluster  (BH FDR < {args.fdr_thresh}, "
        f"background: {args.background})",
        fontsize=12, fontweight="bold",
    )
    plt.tight_layout()
    bar_path = os.path.join(out_dir, "cluster_enrichment_per_cluster.png")
    plt.savefig(bar_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {bar_path}")


# --- Main ---------------------------------------------------------------------
def main():
    args    = parse_args()
    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)

    print(f"Output directory : {out_dir}")
    print(f"Background mode  : {args.background}")
    print(f"Leiden resolution: {args.resolution}")
    print(f"HVGs             : {args.n_hvgs}")
    print(f"PCs              : {args.n_pcs}")
    print(f"k-NN             : {args.n_neighbors}")
    print(f"Min cells / TF   : {args.min_cells}")
    print(f"FDR threshold    : {args.fdr_thresh}")

    ensg2sym = build_ensg_to_symbol()

    # --- Step 1: Cluster ------------------------------------------------------
    cluster_csv = os.path.join(out_dir, "cluster_labels.csv")
    if args.skip_clustering and os.path.exists(cluster_csv):
        print(f"\n  Loading precomputed cluster labels from {cluster_csv} ...")
        cl_df          = pd.read_csv(cluster_csv)
        obs_names      = cl_df["barcode"].tolist()
        cluster_labels = cl_df["cluster"].astype(str).values
        print(f"  {len(obs_names):,} cells, "
              f"{len(np.unique(cluster_labels))} clusters")
    else:
        obs_names, cluster_labels, pca_coords = cluster_cells(args)
        cl_df = pd.DataFrame({"barcode": obs_names, "cluster": cluster_labels})
        cl_df.to_csv(cluster_csv, index=False)
        print(f"  Cluster labels saved: {cluster_csv}")
        # Save PCA coordinates for downstream UMAP visualization
        pca_path = os.path.join(out_dir, "pca_coords.npz")
        np.savez_compressed(pca_path,
                            X_pca=pca_coords,
                            barcodes=np.array(obs_names))
        print(f"  PCA coords saved:     {pca_path}")

    # --- Step 2: Guide assignments --------------------------------------------
    cell_target, target_to_cells, nt_cell_indices = load_guide_assignments(
        obs_names
    )

    # --- Step 3: Enrichment ---------------------------------------------------
    df, clusters = test_enrichment(
        cluster_labels, cell_target, target_to_cells,
        nt_cell_indices, args, ensg2sym,
    )

    # --- Step 4: Save ---------------------------------------------------------
    save_results(df, clusters, args, out_dir)

    print("\n=== Complete ===")


if __name__ == "__main__":
    main()