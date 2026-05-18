"""PCA visualization of raw input programs.

PCA is fit once on the stacked + normalized program matrix N
(rows = individual programs across all datasets, not meta-consensus medians).
Plots reuse that same PC1/PC2 projection, varying only the coloring.
"""

from __future__ import annotations

from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def fit_pca(N: pd.DataFrame, n_components: int = 2) -> tuple[np.ndarray, PCA]:
    """Fit PCA on N (rows = programs, cols = genes) and return (scores, model)."""
    pca = PCA(n_components=n_components, random_state=0)
    scores = pca.fit_transform(N.values)
    return scores, pca


def _axis_label(pca: PCA, i: int) -> str:
    return f"PC{i + 1} ({pca.explained_variance_ratio_[i] * 100:.1f}% var)"


def plot_pca_programs_by_dataset(
    scores: np.ndarray,
    pca: PCA,
    index: pd.Index,
    out_path: str,
    title: str = "",
) -> None:
    """Scatter PC1 vs PC2 of input programs, colored by dataset of origin."""
    datasets = np.array([p.split("|", 1)[0] for p in index])
    unique_ds = sorted(set(datasets))
    cmap = plt.get_cmap("tab10")

    fig, ax = plt.subplots(figsize=(7, 6))
    for i, ds in enumerate(unique_ds):
        mask = datasets == ds
        ax.scatter(scores[mask, 0], scores[mask, 1],
                   s=24, alpha=0.7, color=cmap(i / max(len(unique_ds) - 1, 1)),
                   label=f"{ds} (n={int(mask.sum())})", edgecolors="none")
    ax.set_xlabel(_axis_label(pca, 0))
    ax.set_ylabel(_axis_label(pca, 1))
    ax.set_title(title or "PCA of input programs (colored by dataset)")
    ax.legend(loc="best", fontsize=8, frameon=True)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_pca_programs_by_cluster(
    scores: np.ndarray,
    pca: PCA,
    index: pd.Index,
    labels: pd.Series,
    out_path: str,
    title: str = "",
    annotate_centroids: bool = False,
) -> None:
    """Scatter PC1 vs PC2 of input programs.

    Points are colored by **dataset of origin** (same palette as the global
    pca_programs_by_dataset plot, for visual consistency). If
    `annotate_centroids=True`, the meta-cluster id is overlaid at each
    cluster's PC1/PC2 centroid. Default is False (clean plot).
    """
    cluster_ids = labels.reindex(index).values.astype(int)
    datasets = np.array([p.split("|", 1)[0] for p in index])
    unique_ds = sorted(set(datasets))
    cmap = plt.get_cmap("tab10")

    fig, ax = plt.subplots(figsize=(7, 6))
    for i, ds in enumerate(unique_ds):
        mask = datasets == ds
        ax.scatter(scores[mask, 0], scores[mask, 1],
                   s=24, alpha=0.75,
                   color=cmap(i / max(len(unique_ds) - 1, 1)),
                   label=f"{ds} (n={int(mask.sum())})",
                   edgecolors="none")

    if annotate_centroids:
        for cl in sorted(set(cluster_ids)):
            mask = cluster_ids == cl
            if mask.sum() == 0:
                continue
            cx, cy = scores[mask, 0].mean(), scores[mask, 1].mean()
            ax.text(cx, cy, str(cl), fontsize=7, ha="center", va="center",
                    color="black",
                    bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none",
                              alpha=0.7))

    ax.set_xlabel(_axis_label(pca, 0))
    ax.set_ylabel(_axis_label(pca, 1))
    ax.set_title(title or "PCA of input programs (color=dataset, label=meta-cluster)")
    ax.legend(loc="best", fontsize=8, frameon=True)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
