"""Per-lane and aggregated mapping figures for sgRNA and scRNA modalities.

Reads results/cross_production_qc/upstream/per_lane_mapping_summary.tsv and emits,
for modality in {scRNA, Guide}:

  per_lane_mapping_<mod>.{pdf,png}    one bar per measurement set, colored by dataset
  aggregated_mapping_<mod>.{pdf,png}  one bar per dataset (sum reads, weighted-mean
                                       alignment%, sum barcodes)

Each figure is 5 horizontally arranged panels:
  Total Reads, Mapped Reads, Alignment %, Detected Barcodes, % Reads in Onlist
The 5th panel is NaN for runs from non-kallisto mappers (e.g. Gersbach Hep) — its
bars are drawn with hatching to flag missing data.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REPO = ROOT.parents[5]

COLORS_YAML = REPO / "config" / "colors" / "production_TF-Perturb-seq.yaml"
MAPPING_TSV = ROOT / "results" / "cross_production_qc" / "upstream" / "per_lane_mapping_summary.tsv"
OUT_DIR = ROOT / "results" / "cross_production_qc" / "upstream"

PANELS = [
    ("total_reads", "Total Reads", "Million reads", 1e6, "{:.0f}"),
    ("paired_reads_mapped", "Mapped Reads", "Million reads", 1e6, "{:.0f}"),
    ("alignment_pct", "Alignment %", "%", 1, "{:.1f}"),
    ("detected_barcodes", "Detected Barcodes", "Thousand barcodes", 1e3, "{:.0f}"),
    ("pct_reads_in_onlist", "% Reads in Onlist", "%", 1, "{:.1f}"),
]


def style_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", which="both", length=0)


def plot_per_lane(df, modality, colors, order, out_prefix):
    sub = df[df["modality"] == modality].copy()
    if sub.empty:
        print(f"  [skip] no rows for modality={modality}")
        return

    sub["__order"] = sub["dataset"].map({d: i for i, d in enumerate(order)})
    sub = sub.sort_values(["__order", "measurement_set"]).reset_index(drop=True)
    n = len(sub)
    x = np.arange(n)
    short_to_dataset = dict(zip(sub["short_name"], sub["dataset"]))

    fig, axes = plt.subplots(1, len(PANELS), figsize=(4.0 * len(PANELS), 4.2))
    fig.suptitle(
        {"scRNA": "scRNA Metrics", "Guide": "sgRNA Metrics", "Hashing": "HTO Metrics"}.get(modality, modality)
        + " — per measurement set",
        fontsize=13,
        fontweight="bold",
    )

    for ax, (key, title, ylabel, scale, _fmt) in zip(axes, PANELS):
        vals = sub[key].astype(float) / scale
        bar_colors = [colors.get(short_to_dataset[s], "#888") for s in sub["short_name"]]
        is_nan = vals.isna()
        ax.bar(x, vals.fillna(0), color=bar_colors, width=0.85, edgecolor="white", linewidth=0.5)
        if is_nan.any():
            ax.bar(
                x[is_nan],
                np.full(is_nan.sum(), ax.get_ylim()[1] if not vals.dropna().empty else 1.0),
                color="none",
                hatch="///",
                edgecolor="#bbb",
                linewidth=0.5,
                width=0.85,
            )
        ax.set_title(title, fontsize=11)
        ax.set_ylabel(ylabel)
        ax.set_xticks([])
        style_axes(ax)

    short_arr = sub["short_name"].values
    spans = []
    for short in pd.unique(short_arr):
        positions = np.where(short_arr == short)[0]
        spans.append((short, positions[0], positions[-1], colors[short_to_dataset[short]]))

    underline_y = -0.025
    label_y = -0.055
    for ax in axes:
        ax.set_xlim(-0.5, n - 0.5)
        for short, i0, i1, color in spans:
            ax.annotate(
                "",
                xy=(i1 + 0.4, underline_y),
                xytext=(i0 - 0.4, underline_y),
                xycoords=("data", "axes fraction"),
                arrowprops=dict(arrowstyle="-", color=color, lw=2.5),
                annotation_clip=False,
            )
            ax.annotate(
                short,
                xy=((i0 + i1) / 2, label_y),
                xycoords=("data", "axes fraction"),
                ha="right",
                va="top",
                rotation=30,
                rotation_mode="anchor",
                fontsize=9,
                color=color,
                fontweight="bold",
                annotation_clip=False,
            )

    plt.tight_layout(rect=(0, 0.10, 1, 0.94))
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"{out_prefix}.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"  wrote {out.relative_to(ROOT)}")
    plt.close(fig)


def plot_aggregated(df, modality, colors, order, out_prefix):
    sub = df[df["modality"] == modality].copy()
    if sub.empty:
        print(f"  [skip] no rows for modality={modality}")
        return

    def _agg(group: pd.DataFrame) -> pd.Series:
        total = group["total_reads"].sum()
        return pd.Series(
            {
                "total_reads": total,
                "paired_reads_mapped": group["paired_reads_mapped"].sum(),
                "alignment_pct": (
                    100 * group["paired_reads_mapped"].sum() / total if total else np.nan
                ),
                "detected_barcodes": group["detected_barcodes"].sum(),
                "pct_reads_in_onlist": (
                    (group["pct_reads_in_onlist"] * group["total_reads"]).sum() / total
                    if total and group["pct_reads_in_onlist"].notna().any()
                    else np.nan
                ),
            }
        )

    agg = sub.groupby(["dataset", "short_name"], dropna=False).apply(_agg, include_groups=False).reset_index()
    agg["__order"] = agg["dataset"].map({d: i for i, d in enumerate(order)})
    agg = agg.sort_values("__order").reset_index(drop=True)

    fig, axes = plt.subplots(1, len(PANELS), figsize=(3.0 * len(PANELS), 3.8))
    fig.suptitle(
        {"scRNA": "scRNA Metrics", "Guide": "sgRNA Metrics", "Hashing": "HTO Metrics"}.get(modality, modality)
        + " — aggregated across lanes",
        fontsize=13,
        fontweight="bold",
    )

    x = np.arange(len(agg))
    bar_colors = [colors.get(d, "#888") for d in agg["dataset"]]

    for ax, (key, title, ylabel, scale, fmt) in zip(axes, PANELS):
        vals = agg[key].astype(float) / scale
        ax.bar(x, vals.fillna(0), color=bar_colors, width=0.7, edgecolor="white", linewidth=0.5)
        for i, (v, raw) in enumerate(zip(vals, agg[key])):
            if pd.isna(raw):
                ax.text(
                    i,
                    (ax.get_ylim()[1] or 1) * 0.5,
                    "n/a",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="#999",
                    fontweight="bold",
                )
            else:
                ax.text(i, v, fmt.format(v), ha="center", va="bottom", fontsize=8)
        ax.set_title(title, fontsize=11)
        ax.set_ylabel(ylabel)
        ax.set_xticks(x)
        ax.set_xticklabels(agg["short_name"], rotation=45, ha="right", fontsize=9)
        style_axes(ax)

    plt.tight_layout(rect=(0, 0, 1, 0.92))
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"{out_prefix}.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"  wrote {out.relative_to(ROOT)}")
    plt.close(fig)


def main() -> int:
    df = pd.read_csv(MAPPING_TSV, sep="\t")
    with open(COLORS_YAML) as f:
        palette = yaml.safe_load(f)
    colors, order = palette["dataset_colors"], palette["dataset_order"]

    for modality, mod_label in [("scRNA", "scrna"), ("Guide", "guide")]:
        print(f"== {modality} ==")
        plot_per_lane(df, modality, colors, order, f"per_lane_mapping_{mod_label}")
        plot_aggregated(df, modality, colors, order, f"aggregated_mapping_{mod_label}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
