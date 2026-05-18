"""Plotting for meta-consensus.

Clustergram is a stripped-down port of torch_cnmf/cnmf.py:1666-1759
(no local-density histogram, with a dataset-of-origin side bar).
"""

from __future__ import annotations

from typing import Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import gridspec
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics.pairwise import euclidean_distances


def _order_within_clusters(topics_dist: np.ndarray, labels: pd.Series) -> List[int]:
    """Within each cluster, sort rows by hierarchical leaves_list."""
    order: List[int] = []
    for cl in sorted(set(labels)):
        mask = (labels == cl).values
        idx = np.where(mask)[0]
        if mask.sum() > 1:
            sub = topics_dist[mask, :][:, mask]
            cl_dist = squareform(sub, checks=False)
            cl_dist[cl_dist < 0] = 0
            link = linkage(cl_dist, "average")
            leaves = leaves_list(link)
            order += list(idx[leaves])
        else:
            order += list(idx)
    return order


def plot_clustergram(
    N: pd.DataFrame,
    labels: pd.Series,
    out_path: str,
    title: str = "",
    draw_cluster_boxes: bool = True,
    box_color: str = "red",
    box_linewidth: float = 1.0,
) -> None:
    """Pairwise distance heatmap with cluster + dataset side bars.

    If ``draw_cluster_boxes`` is True, draws a rectangular outline around
    each KMeans cluster's diagonal block on the distance heatmap.
    """
    topics_dist = euclidean_distances(N.values)
    spectra_order = _order_within_clusters(topics_dist, labels)

    datasets = np.array([p.split("|", 1)[0] for p in N.index])
    unique_ds = sorted(set(datasets))
    ds_to_int = {d: i for i, d in enumerate(unique_ds)}
    ds_codes = np.array([ds_to_int[d] for d in datasets])

    width_ratios = [0.5, 0.5, 9, 1]
    height_ratios = [0.5, 0.5, 9]
    fig = plt.figure(figsize=(sum(width_ratios), sum(height_ratios)))
    gs = gridspec.GridSpec(
        len(height_ratios), len(width_ratios), fig,
        0.02, 0.02, 0.98, 0.94,
        height_ratios=height_ratios, width_ratios=width_ratios,
        wspace=0.02, hspace=0.02,
    )

    # main distance heatmap
    dist_ax = fig.add_subplot(gs[2, 2], xticks=[], yticks=[])
    D = topics_dist[spectra_order, :][:, spectra_order]
    dist_im = dist_ax.imshow(D, interpolation="none", cmap="viridis",
                             aspect="auto", rasterized=True)

    if draw_cluster_boxes:
        from matplotlib.patches import Rectangle
        ordered_labels = labels.values[spectra_order]
        n = len(ordered_labels)
        start = 0
        while start < n:
            cl = ordered_labels[start]
            end = start
            while end < n and ordered_labels[end] == cl:
                end += 1
            # imshow extent: each cell spans [i-0.5, i+0.5]
            size = end - start
            rect = Rectangle(
                (start - 0.5, start - 0.5), size, size,
                fill=False, edgecolor=box_color, linewidth=box_linewidth,
            )
            dist_ax.add_patch(rect)
            start = end

    # left side bar: dataset of origin
    ds_left = fig.add_subplot(gs[2, 0], xticks=[], yticks=[])
    ds_left.imshow(ds_codes[spectra_order].reshape(-1, 1),
                   interpolation="none", cmap="tab10", aspect="auto")
    ds_left.set_ylabel("dataset", fontsize=8)

    # left side bar: cluster id
    cl_left = fig.add_subplot(gs[2, 1], xticks=[], yticks=[])
    cl_left.imshow(labels.values[spectra_order].reshape(-1, 1),
                   interpolation="none", cmap="Spectral", aspect="auto")

    # top side bars (mirror of left)
    ds_top = fig.add_subplot(gs[0, 2], xticks=[], yticks=[])
    ds_top.imshow(ds_codes[spectra_order].reshape(1, -1),
                  interpolation="none", cmap="tab10", aspect="auto")
    cl_top = fig.add_subplot(gs[1, 2], xticks=[], yticks=[])
    cl_top.imshow(labels.values[spectra_order].reshape(1, -1),
                  interpolation="none", cmap="Spectral", aspect="auto")

    # colorbar
    cbar_ax = fig.add_subplot(gs[2, 3])
    cbar_ax.set_axis_off()
    inner = gridspec.GridSpecFromSubplotSpec(8, 1, subplot_spec=gs[2, 3])
    cb = fig.add_subplot(inner[3, 0])
    fig.colorbar(dist_im, cax=cb, orientation="horizontal")
    cb.set_title("Euclidean distance", fontsize=8)

    # dataset legend
    leg_ax = fig.add_subplot(inner[6, 0])
    leg_ax.set_axis_off()
    handles = [plt.Line2D([0], [0], marker="s", linestyle="",
                           color=plt.cm.tab10(i / max(len(unique_ds) - 1, 1)),
                           label=d)
               for i, d in enumerate(unique_ds)]
    leg_ax.legend(handles=handles, loc="center", fontsize=8, frameon=False)

    if title:
        fig.suptitle(title, fontsize=10)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_composition_bar(
    cluster_summary: pd.DataFrame,
    dataset_names: List[str],
    out_path: str,
    title: str = "",
) -> None:
    """Stacked bar: programs per dataset per cluster.

    Sorted so shared clusters appear left, dataset-specific on right.
    """
    df = cluster_summary.copy()
    category_order = {"shared_all": 0, "shared_partial": 1}
    df["sort_key"] = df["category"].map(category_order).fillna(2).astype(int)
    df = df.sort_values(["sort_key", "n_programs"], ascending=[True, False])

    fig, ax = plt.subplots(figsize=(max(8, len(df) * 0.25), 5))
    bottom = np.zeros(len(df))
    cmap = plt.get_cmap("tab10")
    for i, ds in enumerate(dataset_names):
        col = f"n_{ds}"
        vals = df[col].values
        ax.bar(np.arange(len(df)), vals, bottom=bottom,
               label=ds, color=cmap(i / max(len(dataset_names) - 1, 1)))
        bottom = bottom + vals

    ax.set_xticks(np.arange(len(df)))
    ax.set_xticklabels(df["meta_cluster_id"].astype(str), rotation=90, fontsize=7)
    ax.set_xlabel("meta-cluster (sorted: shared → specific)")
    ax.set_ylabel("# input programs")
    ax.set_title(title or "Dataset composition per meta-cluster")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_category_breakdown(
    cluster_summary: pd.DataFrame,
    dataset_names: List[str],
    out_path: str,
    title: str = "",
) -> None:
    """Bar chart of category counts."""
    cats = ["shared_all", "shared_partial"] + [f"specific_{d}" for d in dataset_names]
    counts = [int((cluster_summary["category"] == c).sum()) for c in cats]

    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(cats, counts, color="steelblue")
    ax.set_ylabel("# meta-clusters")
    ax.set_title(title or "Meta-cluster category breakdown")
    ax.tick_params(axis="x", rotation=30)
    for b, c in zip(bars, counts):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(),
                str(c), ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_top_genes_heatmap(
    meta_spectra: pd.DataFrame,
    out_path: str,
    n_top: int = 20,
    title: str = "",
) -> None:
    """Heatmap of top-n genes per meta-program (union across programs)."""
    top_genes_set = set()
    for cid in meta_spectra.index:
        top = meta_spectra.loc[cid].sort_values(ascending=False).head(n_top).index
        top_genes_set.update(top.tolist())
    top_genes = sorted(top_genes_set)
    sub = meta_spectra.loc[:, top_genes]

    fig, ax = plt.subplots(figsize=(max(8, len(top_genes) * 0.12),
                                     max(4, len(meta_spectra) * 0.2)))
    im = ax.imshow(sub.values, aspect="auto", cmap="viridis")
    ax.set_yticks(np.arange(len(meta_spectra)))
    ax.set_yticklabels(meta_spectra.index, fontsize=7)
    ax.set_xticks(np.arange(len(top_genes)))
    ax.set_xticklabels(top_genes, rotation=90, fontsize=6)
    ax.set_title(title or f"Top {n_top} genes per meta-program (union)")
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_k_selection(silhouettes: Dict[int, float], out_path: str) -> None:
    """Silhouette score vs k."""
    ks = sorted(silhouettes.keys())
    vals = [silhouettes[k] for k in ks]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(ks, vals, "o-")
    ax.set_xlabel("meta-consensus k")
    ax.set_ylabel("silhouette score (euclidean)")
    ax.set_title("k selection")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
