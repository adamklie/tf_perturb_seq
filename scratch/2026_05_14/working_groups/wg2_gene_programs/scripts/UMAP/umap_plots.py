"""UMAP visualization of raw input programs.

UMAP is fit once on the stacked + normalized program matrix N
(rows = individual programs across all datasets). Plots reuse that same
embedding, varying only the coloring / annotations.
"""

from __future__ import annotations

from typing import Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap


def fit_umap(
    N: pd.DataFrame,
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    metric: str = "euclidean",
    random_state: int = 0,
) -> Tuple[np.ndarray, "umap.UMAP"]:
    """Fit UMAP on N (rows = programs, cols = genes) and return (embedding, model)."""
    model = umap.UMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=random_state,
    )
    embedding = model.fit_transform(N.values)
    return embedding, model


def plot_umap_programs_by_dataset(
    embedding: np.ndarray,
    index: pd.Index,
    out_path: str,
    title: str = "",
) -> None:
    """Scatter UMAP1 vs UMAP2 of input programs, colored by dataset of origin."""
    datasets = np.array([p.split("|", 1)[0] for p in index])
    unique_ds = sorted(set(datasets))
    cmap = plt.get_cmap("tab10")

    fig, ax = plt.subplots(figsize=(7, 6))
    for i, ds in enumerate(unique_ds):
        mask = datasets == ds
        ax.scatter(embedding[mask, 0], embedding[mask, 1],
                   s=24, alpha=0.7,
                   color=cmap(i / max(len(unique_ds) - 1, 1)),
                   label=f"{ds} (n={int(mask.sum())})", edgecolors="none")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(title or "UMAP of input programs (colored by dataset)")
    ax.legend(loc="best", fontsize=8, frameon=True)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_umap_programs_by_cluster(
    embedding: np.ndarray,
    index: pd.Index,
    labels: pd.Series,
    out_path: str,
    title: str = "",
    annotate_centroids: bool = False,
) -> None:
    """Scatter UMAP1 vs UMAP2 of input programs.

    Points are colored by **dataset of origin** (same palette as the global
    umap_programs_by_dataset plot, for visual consistency). If
    `annotate_centroids=True`, the meta-cluster id is overlaid at each
    cluster's UMAP centroid. Default is False (clean plot).
    """
    cluster_ids = labels.reindex(index).values.astype(int)
    datasets = np.array([p.split("|", 1)[0] for p in index])
    unique_ds = sorted(set(datasets))
    cmap = plt.get_cmap("tab10")

    fig, ax = plt.subplots(figsize=(7, 6))
    for i, ds in enumerate(unique_ds):
        mask = datasets == ds
        ax.scatter(embedding[mask, 0], embedding[mask, 1],
                   s=24, alpha=0.75,
                   color=cmap(i / max(len(unique_ds) - 1, 1)),
                   label=f"{ds} (n={int(mask.sum())})",
                   edgecolors="none")

    if annotate_centroids:
        for cl in sorted(set(cluster_ids)):
            mask = cluster_ids == cl
            if mask.sum() == 0:
                continue
            cx, cy = embedding[mask, 0].mean(), embedding[mask, 1].mean()
            ax.text(cx, cy, str(cl), fontsize=7, ha="center", va="center",
                    color="black",
                    bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none",
                              alpha=0.7))

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(title or "UMAP of input programs (color=dataset, label=meta-cluster)")
    ax.legend(loc="best", fontsize=8, frameon=True)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
