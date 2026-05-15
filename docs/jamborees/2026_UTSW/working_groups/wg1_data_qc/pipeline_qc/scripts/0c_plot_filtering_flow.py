"""Cell-filtering funnel diagram across the 4 production datasets.

Reads results/cross_production_qc/upstream/filtering_funnel.tsv and emits:
  results/cross_production_qc/upstream/filtering_flow.{pdf,png}      — linear-width funnel
  results/cross_production_qc/upstream/filtering_flow_log.{pdf,png}  — log-width funnel

One column per dataset, each step a colored trapezoid whose width is proportional
to the cell count, with the count label inside and the stage name beneath.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.patches import Polygon

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REPO = HERE.parents[6]  # scripts -> pipeline_qc -> wg1_data_qc -> working_groups -> 2026_UTSW -> jamborees -> docs -> repo
COLORS_YAML = REPO / "config" / "colors" / "production_TF-Perturb-seq.yaml"
FUNNEL_TSV = ROOT / "results" / "cross_production_qc" / "upstream" / "filtering_funnel.tsv"
OUT_DIR = ROOT / "results" / "cross_production_qc" / "upstream"


STAGE_LABELS = {
    "Concatenated scRNA anndata (unfiltered)": "Concatenated\n(unfiltered)",
    "Barcode filter (knee)": "Barcode filter\n(knee)",
    "Barcode filter (knee2)": "Barcode filter\n(knee)",
    "Min genes per cell (>= 2500)": "Min genes\n(≥ 2500)",
    "Intersection with hashing barcodes": "Hashing\nintersection",
    "Demultiplex filter (HTO)": "Demultiplex\nfilter (HTO)",
    "Intersection with guide barcodes": "Guide\nintersection",
    "Final cells in MuData": "Final\n(MuData)",
}


def short_stage(raw: str) -> str:
    if raw in STAGE_LABELS:
        return STAGE_LABELS[raw]
    if raw.startswith("Mito filter"):
        pct = raw.split("<")[-1].strip(" )%")
        try:
            return f"Mito filter\n(< {float(pct):.0f}%)"
        except ValueError:
            return "Mito filter"
    return raw


def fmt_count(n: float) -> str:
    n = float(n)
    if n >= 1e9:
        return f"{n/1e9:.1f}B"
    if n >= 1e6:
        v = n / 1e6
        return f"{v:.1f}M" if v < 10 else f"{v:.0f}M"
    if n >= 1e3:
        return f"{n/1e3:.0f}K"
    return f"{int(n)}"


def load_colors() -> tuple[dict, list]:
    with open(COLORS_YAML) as f:
        cfg = yaml.safe_load(f)
    return cfg["dataset_colors"], cfg["dataset_order"]


def draw_funnel(
    ax: plt.Axes,
    short_name: str,
    color: str,
    stages: pd.DataFrame,
    scale: str,
    global_max: float,
    n_stages_max: int,
) -> None:
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-(n_stages_max + 1.0) * 1.2, 0.6)
    ax.set_aspect("auto")
    ax.axis("off")
    ax.text(0, 0.25, short_name, ha="center", va="bottom", fontsize=11, fontweight="bold", color=color)

    if scale == "log":
        widths = np.log10(np.clip(stages["cells"].values, 1, None)) / np.log10(global_max)
    else:
        widths = stages["cells"].values / global_max
    widths = np.clip(widths, 0.05, 1.0)

    box_h = 0.6
    gap = 0.55

    for i, (_, row) in enumerate(stages.iterrows()):
        w_top = widths[i]
        w_bot = widths[i + 1] if i + 1 < len(widths) else w_top
        y_top = -(i * (box_h + gap))
        y_bot = y_top - box_h

        ax.add_patch(
            Polygon(
                [
                    (-w_top / 2, y_top),
                    (w_top / 2, y_top),
                    (w_bot / 2, y_bot),
                    (-w_bot / 2, y_bot),
                ],
                facecolor=color,
                edgecolor="white",
                linewidth=1.0,
                alpha=0.85,
            )
        )
        ax.text(
            0,
            (y_top + y_bot) / 2,
            fmt_count(row["cells"]),
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
            color="white",
        )
        ax.text(
            0,
            y_bot - gap * 0.35,
            short_stage(row["stage"]),
            ha="center",
            va="top",
            fontsize=7,
            color="#444",
        )


def main() -> int:
    funnel = pd.read_csv(FUNNEL_TSV, sep="\t")
    if funnel.empty:
        print("filtering_funnel.tsv is empty", file=sys.stderr)
        return 1
    colors, order = load_colors()

    datasets = [d for d in order if d in funnel["dataset"].unique()]
    n_datasets = len(datasets)
    n_stages_max = funnel.groupby("dataset").size().max()

    fig_w = 2.4 * n_datasets
    fig_h = 0.85 * n_stages_max + 1.0

    for scale in ("linear", "log"):
        fig, axes = plt.subplots(1, n_datasets, figsize=(fig_w, fig_h), sharey=True)
        if n_datasets == 1:
            axes = [axes]
        fig.suptitle("Cell Filtering Flow", fontsize=13, fontweight="bold", y=0.97)

        global_max_raw = funnel["cells"].max()
        for ax, ds in zip(axes, datasets):
            stages = funnel[funnel["dataset"] == ds].sort_values("stage_index").reset_index(drop=True)
            draw_funnel(
                ax,
                short_name=stages["short_name"].iloc[0],
                color=colors[ds],
                stages=stages,
                scale=scale,
                global_max=global_max_raw,
                n_stages_max=n_stages_max,
            )

        plt.tight_layout(rect=(0, 0, 1, 0.94))
        for ext in ("pdf", "png"):
            suffix = "_log" if scale == "log" else ""
            out = OUT_DIR / f"filtering_flow{suffix}.{ext}"
            fig.savefig(out, dpi=200, bbox_inches="tight")
            print(f"wrote {out.relative_to(ROOT)}")
        plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
