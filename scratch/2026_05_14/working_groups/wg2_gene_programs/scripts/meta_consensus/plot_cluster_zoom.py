"""Zoomed cluster figure for the meta-consensus clustergram.

Reproduces the main pairwise-distance heatmap (same data + ordering as
``plotting.plot_clustergram``) on the left, with one meta-cluster boxed and
labeled, and a large zoom panel on the right showing the cluster's 4xN
sub-distance heatmap plus a GO-term annotation block (shared term across
datasets + dataset-specific top terms).

Standalone CLI; also importable via ``plot_cluster_zoom``.

Usage:
    python plot_cluster_zoom.py --config config.yaml --k 60 --cluster 52
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib import gridspec
from matplotlib.patches import Rectangle
from matplotlib.patheffects import withStroke
from sklearn.metrics.pairwise import euclidean_distances

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from meta_consensus import (  # noqa: E402
    DatasetSpec,
    intersect_genes,
    load_spectra,
    normalize_rows,
    stack_programs,
)
from plotting import _order_within_clusters  # noqa: E402


def _load_normalized(cfg: dict) -> pd.DataFrame:
    datasets = [DatasetSpec(**d) for d in cfg["datasets"]]
    input_file_type = cfg.get("input_file_type", "gene_spectra_score")
    normalization = cfg.get("normalization", "l2")
    dfs = [load_spectra(d, input_file_type) for d in datasets]
    common = intersect_genes(dfs)
    return normalize_rows(stack_programs(dfs, common), normalization)


def _load_labels(k_dir: str, k: int) -> pd.Series:
    path = os.path.join(k_dir, f"cluster_assignments.k_{k}.txt")
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df["meta_cluster_id"]


def _go_terms_for_member(go_long: pd.DataFrame, cluster_id: int, dataset: str) -> List[str]:
    """Union of top10 GO terms across that dataset's members in this cluster."""
    sub = go_long[
        (go_long["meta_cluster_id"] == cluster_id) & (go_long["dataset"] == dataset)
    ]
    terms: List[str] = []
    for raw in sub["top10_GO"].dropna():
        for t in str(raw).split(";"):
            t = t.strip()
            if t and t not in terms:
                terms.append(t)
    return terms


def _dataset_specific_terms(
    go_long: pd.DataFrame,
    cluster_id: int,
    target_ds: str,
    all_datasets: List[str],
    top_n: int = 4,
) -> List[str]:
    """Top terms enriched in target_ds members but absent from other datasets' members."""
    own = _go_terms_for_member(go_long, cluster_id, target_ds)
    others: set = set()
    for ds in all_datasets:
        if ds == target_ds:
            continue
        others.update(_go_terms_for_member(go_long, cluster_id, ds))
    unique = [t for t in own if t not in others]
    return unique[:top_n]


def _build_dataset_codes(programs: List[str]) -> Tuple[np.ndarray, List[str], Dict[str, int]]:
    datasets = np.array([p.split("|", 1)[0] for p in programs])
    unique_ds = sorted(set(datasets))
    ds_to_int = {d: i for i, d in enumerate(unique_ds)}
    codes = np.array([ds_to_int[d] for d in datasets])
    return codes, unique_ds, ds_to_int


def plot(
    N: pd.DataFrame,
    labels: pd.Series,
    cluster_id: int,
    cluster_summary: pd.DataFrame,
    go_summary: pd.DataFrame,
    go_long: pd.DataFrame,
    out_path: str,
    title: str = "",
) -> None:
    """Render the combined main + zoom + GO figure."""
    if cluster_id not in set(labels.values):
        raise ValueError(f"cluster_id={cluster_id} not present in labels")

    # Distance matrix and ordering (matches plotting.plot_clustergram).
    topics_dist = euclidean_distances(N.values)
    spectra_order = _order_within_clusters(topics_dist, labels)
    ordered_labels = labels.values[spectra_order]
    ordered_programs = [N.index[i] for i in spectra_order]

    ds_codes, unique_ds, ds_to_int = _build_dataset_codes(list(N.index))
    n_ds = len(unique_ds)

    # Identify the on-diagonal block for the highlighted cluster.
    mask = ordered_labels == cluster_id
    block_idx = np.where(mask)[0]
    block_start = int(block_idx[0])
    block_end = int(block_idx[-1]) + 1
    block_size = block_end - block_start
    member_programs = [ordered_programs[i] for i in block_idx]
    member_datasets = [p.split("|", 1)[0] for p in member_programs]

    # ---------- figure & layout ----------
    fig = plt.figure(figsize=(20, 11))
    outer = gridspec.GridSpec(
        1, 2, figure=fig,
        left=0.04, right=0.97, top=0.93, bottom=0.05,
        width_ratios=[9, 5], wspace=0.12,
    )

    # ===== LEFT: main clustergram (mirrors plotting.plot_clustergram structure) =====
    left = gridspec.GridSpecFromSubplotSpec(
        3, 4, subplot_spec=outer[0, 0],
        height_ratios=[0.5, 0.5, 9],
        width_ratios=[0.5, 0.5, 9, 1],
        hspace=0.02, wspace=0.02,
    )

    dist_ax = fig.add_subplot(left[2, 2], xticks=[], yticks=[])
    D = topics_dist[spectra_order, :][:, spectra_order]
    dist_im = dist_ax.imshow(
        D, interpolation="none", cmap="viridis", aspect="auto", rasterized=True
    )

    # red boxes for every cluster
    n = len(ordered_labels)
    start = 0
    while start < n:
        cl = ordered_labels[start]
        end = start
        while end < n and ordered_labels[end] == cl:
            end += 1
        size = end - start
        rect = Rectangle(
            (start - 0.5, start - 0.5), size, size,
            fill=False, edgecolor="red", linewidth=0.8,
        )
        dist_ax.add_patch(rect)
        start = end

    # highlight box for the focal cluster (thicker, white-with-black-stroke)
    pad = max(1.5, block_size * 0.2)
    highlight = Rectangle(
        (block_start - 0.5 - pad, block_start - 0.5 - pad),
        block_size + 2 * pad, block_size + 2 * pad,
        fill=False, edgecolor="white", linewidth=2.2,
    )
    highlight.set_path_effects([withStroke(linewidth=4.0, foreground="black")])
    dist_ax.add_patch(highlight)

    # label tag near the focal cluster
    cx = block_start + block_size / 2.0
    tag_x = cx + max(8, block_size + 6)
    tag_y = cx - max(8, block_size + 6)
    tag_y = max(tag_y, 1)
    txt = dist_ax.text(
        tag_x, tag_y, f"#{cluster_id}",
        color="white", fontsize=14, fontweight="bold",
        ha="left", va="center",
    )
    txt.set_path_effects([withStroke(linewidth=3.0, foreground="black")])
    # connector line from label to box
    dist_ax.annotate(
        "", xy=(block_start + block_size - 0.5, block_start - 0.5),
        xytext=(tag_x - 0.5, tag_y),
        arrowprops=dict(arrowstyle="-", color="white", lw=1.2,
                        shrinkA=2, shrinkB=2,
                        path_effects=[withStroke(linewidth=2.5, foreground="black")]),
    )

    # left side bars: dataset of origin + cluster id
    ds_left = fig.add_subplot(left[2, 0], xticks=[], yticks=[])
    ds_left.imshow(
        ds_codes[spectra_order].reshape(-1, 1),
        interpolation="none", cmap="tab10", aspect="auto",
    )
    ds_left.set_ylabel("dataset", fontsize=8)

    cl_left = fig.add_subplot(left[2, 1], xticks=[], yticks=[])
    cl_left.imshow(
        labels.values[spectra_order].reshape(-1, 1),
        interpolation="none", cmap="Spectral", aspect="auto",
    )

    ds_top = fig.add_subplot(left[0, 2], xticks=[], yticks=[])
    ds_top.imshow(
        ds_codes[spectra_order].reshape(1, -1),
        interpolation="none", cmap="tab10", aspect="auto",
    )
    cl_top = fig.add_subplot(left[1, 2], xticks=[], yticks=[])
    cl_top.imshow(
        labels.values[spectra_order].reshape(1, -1),
        interpolation="none", cmap="Spectral", aspect="auto",
    )

    # shared colorbar (right of left panel)
    cbar_holder = fig.add_subplot(left[2, 3])
    cbar_holder.set_axis_off()
    inner = gridspec.GridSpecFromSubplotSpec(8, 1, subplot_spec=left[2, 3])
    cb_ax = fig.add_subplot(inner[3, 0])
    fig.colorbar(dist_im, cax=cb_ax, orientation="horizontal")
    cb_ax.set_title("Euclidean distance", fontsize=8)

    leg_ax = fig.add_subplot(inner[6, 0])
    leg_ax.set_axis_off()
    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="",
                   color=plt.cm.tab10(i / max(n_ds - 1, 1)), label=d)
        for i, d in enumerate(unique_ds)
    ]
    leg_ax.legend(handles=handles, loc="center", fontsize=8, frameon=False)

    # ===== RIGHT: zoom panel + GO annotations =====
    right = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer[0, 1],
        height_ratios=[5.5, 5.0], hspace=0.18,
    )

    # zoom heatmap with its own dataset side bar
    zoom_block = gridspec.GridSpecFromSubplotSpec(
        2, 3, subplot_spec=right[0, 0],
        height_ratios=[0.35, 6.0],
        width_ratios=[0.35, 6.0, 0.6],
        hspace=0.03, wspace=0.03,
    )

    sub_dist = D[block_start:block_end, block_start:block_end]

    # zoom imshow
    zoom_ax = fig.add_subplot(zoom_block[1, 1])
    zoom_im = zoom_ax.imshow(
        sub_dist, interpolation="none", cmap="viridis", aspect="auto",
        rasterized=True,
    )
    zoom_ax.set_xticks(np.arange(block_size))
    zoom_ax.set_yticks(np.arange(block_size))
    zoom_ax.set_xticklabels(member_programs, rotation=45, ha="right", fontsize=9)
    zoom_ax.set_yticklabels(member_programs, fontsize=9)
    zoom_ax.tick_params(axis="both", length=0)
    zoom_ax.set_title(
        f"Cluster #{cluster_id} (n={block_size}, "
        f"{len(set(member_datasets))} datasets)",
        fontsize=12, pad=8,
    )

    # zoom side bars (dataset codes)
    member_codes = np.array([ds_to_int[d] for d in member_datasets])
    z_left = fig.add_subplot(zoom_block[1, 0], xticks=[], yticks=[])
    z_left.imshow(member_codes.reshape(-1, 1), interpolation="none",
                  cmap="tab10", aspect="auto", vmin=0, vmax=max(n_ds - 1, 1))
    z_top = fig.add_subplot(zoom_block[0, 1], xticks=[], yticks=[])
    z_top.imshow(member_codes.reshape(1, -1), interpolation="none",
                 cmap="tab10", aspect="auto", vmin=0, vmax=max(n_ds - 1, 1))

    # zoom colorbar
    z_cb_holder = fig.add_subplot(zoom_block[1, 2])
    z_cb_holder.set_axis_off()
    z_cb_inner = gridspec.GridSpecFromSubplotSpec(5, 1, subplot_spec=zoom_block[1, 2])
    z_cb = fig.add_subplot(z_cb_inner[2, 0])
    fig.colorbar(zoom_im, cax=z_cb, orientation="vertical")
    z_cb.tick_params(labelsize=7)
    z_cb.set_title("dist", fontsize=8, pad=4)

    # ===== GO annotation block =====
    go_ax = fig.add_subplot(right[1, 0])
    go_ax.set_axis_off()

    # shared GO term(s)
    cs_row = go_summary[go_summary["meta_cluster_id"] == cluster_id]
    if len(cs_row) == 0:
        shared = "(no GO summary row)"
    else:
        raw = cs_row.iloc[0].get("GO_intersect_across_datasets", "")
        if pd.isna(raw) or str(raw).strip() == "":
            shared = "(none)"
        else:
            shared = "; ".join(t.strip() for t in str(raw).split(";") if t.strip())

    # category + member count
    cs_summary = cluster_summary[cluster_summary["meta_cluster_id"] == cluster_id]
    if len(cs_summary):
        category = cs_summary.iloc[0]["category"]
        datasets_present = sorted({d for d in member_datasets})
    else:
        category = "?"
        datasets_present = sorted(set(member_datasets))

    # dataset-specific top terms
    ds_specific: Dict[str, List[str]] = {}
    for ds in datasets_present:
        ds_specific[ds] = _dataset_specific_terms(
            go_long, cluster_id, ds, datasets_present, top_n=4
        )

    # render text block
    y = 0.98
    line_h = 0.062
    go_ax.text(
        0.0, y,
        f"Cluster #{cluster_id} — {category} ({len(member_programs)} programs, "
        f"{len(datasets_present)}/{n_ds} datasets)",
        fontsize=12, fontweight="bold", transform=go_ax.transAxes, va="top",
    )
    y -= line_h * 1.6

    go_ax.text(
        0.0, y, "Shared GO term(s) across all datasets:",
        fontsize=10, fontweight="bold", color="#1f6feb",
        transform=go_ax.transAxes, va="top",
    )
    y -= line_h
    go_ax.text(
        0.02, y, shared if shared else "(none)",
        fontsize=10, transform=go_ax.transAxes, va="top", wrap=True,
    )
    y -= line_h * 1.6

    palette = {d: plt.cm.tab10(i / max(n_ds - 1, 1)) for i, d in enumerate(unique_ds)}
    for ds in datasets_present:
        go_ax.text(
            0.0, y, f"{ds}-specific top terms:",
            fontsize=10, fontweight="bold", color=palette[ds],
            transform=go_ax.transAxes, va="top",
        )
        y -= line_h
        terms = ds_specific[ds] or ["(no unique terms)"]
        for t in terms:
            go_ax.text(
                0.04, y, f"• {t}",
                fontsize=9, transform=go_ax.transAxes, va="top",
            )
            y -= line_h * 0.85
        y -= line_h * 0.4

    # top suptitle
    if title:
        fig.suptitle(title, fontsize=14, y=0.985)

    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="Path to config.yaml")
    ap.add_argument("--k", type=int, default=60, help="Meta-consensus k to use")
    ap.add_argument("--cluster", type=int, default=52, help="Meta-cluster id to zoom into")
    ap.add_argument("--out", default=None, help="Output PNG path (default: <k_dir>/cluster_<id>_zoom.png)")
    ap.add_argument(
        "--go_summary",
        default=None,
        help="Path to Shared_GO_summary.k_<k>.tsv "
             "(default: <output_dir>/../Shared_GO/Shared_GO_summary.k_<k>.tsv)",
    )
    ap.add_argument(
        "--go_long",
        default=None,
        help="Path to Shared_GO_long.k_<k>.tsv "
             "(default: <output_dir>/../Shared_GO/Shared_GO_long.k_<k>.tsv)",
    )
    args = ap.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    out_root = cfg["output_dir"]
    k_dir = os.path.join(out_root, f"k_{args.k}")
    if not os.path.isdir(k_dir):
        raise FileNotFoundError(f"k directory not found: {k_dir}")

    out_path = args.out or os.path.join(k_dir, f"cluster_{args.cluster}_zoom.png")

    go_dir_default = os.path.join(os.path.dirname(out_root.rstrip("/")), "Shared_GO")
    go_summary_path = args.go_summary or os.path.join(
        go_dir_default, f"Shared_GO_summary.k_{args.k}.tsv"
    )
    go_long_path = args.go_long or os.path.join(
        go_dir_default, f"Shared_GO_long.k_{args.k}.tsv"
    )

    print(f"# Loading normalized stacked spectra (config={args.config})")
    N = _load_normalized(cfg)
    print(f"  shape: {N.shape}")

    print(f"# Loading cluster assignments for k={args.k}")
    labels = _load_labels(k_dir, args.k)

    cluster_summary = pd.read_csv(
        os.path.join(k_dir, f"cluster_summary.k_{args.k}.tsv"), sep="\t"
    )
    go_summary = pd.read_csv(go_summary_path, sep="\t")
    go_long = pd.read_csv(go_long_path, sep="\t")

    print(f"# Plotting cluster #{args.cluster} zoom -> {out_path}")
    plot(
        N=N,
        labels=labels,
        cluster_id=args.cluster,
        cluster_summary=cluster_summary,
        go_summary=go_summary,
        go_long=go_long,
        out_path=out_path,
        title=f"Meta-consensus k={args.k} — cluster #{args.cluster} zoom",
    )
    print("# Done.")


if __name__ == "__main__":
    main()
