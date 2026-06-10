"""Cell yield by filter config vs expected.

For each benchmark dataset, plot how many cells survive each scRNA filter
config (200, 500, 800, 2000 UMI, knee2) compared to the expected number of
cells per the experimental design (`dataset_overview.tsv`).

Outputs:
  cells_vs_expected.pdf   — small multiples, one panel per dataset
  cells_vs_expected_overlaid.pdf  — single-panel, all 5 datasets overlaid
  cells_vs_expected.tsv   — observed n_cells, expected, ratio per (dataset, config)

Colors are dataset-specific (per `config/colors/technology-benchmark_WTC11_TF-Perturb-seq.yaml`).
"""
from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

PROJECT_ROOT = Path("/cellar/users/aklie/projects/tf_perturb_seq")
BASE = PROJECT_ROOT / "datasets" / "technology-benchmark_WTC11_TF-Perturb-seq"
RESULTS = BASE / "results" / "cross_tech_comparison"
OUT = RESULTS / "cells_vs_expected"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.append(str(PROJECT_ROOT / "config"))
from loader import load_colors
cfg = load_colors("technology-benchmark_WTC11_TF-Perturb-seq")
DATASET_COLORS = cfg["dataset_colors"]
DATASET_ORDER = cfg["dataset_order"]

RUNS = [
    ("cleanser_extremes_200_mito_15pc", "200 UMI"),
    ("cleanser_500_mito_15pc", "500 UMI"),
    ("cleanser_800_mito_15pc", "800 UMI"),
    ("cleanser_extremes_2000_mito_15pc", "2000 UMI"),
    ("cleanser_knee2_mito_15pc", "Knee2"),
]
RUN_LABELS = [lbl for _, lbl in RUNS]

SHORT = {
    "Hon_WTC11-benchmark_TF-Perturb-seq": "Hon",
    "Huangfu_WTC11-benchmark_TF-Perturb-seq": "Huangfu",
    "Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3": "Gersbach GEM-Xv3",
    "Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2": "Gersbach HTv2",
    "Engreitz_WTC11-benchmark_TF-Perturb-seq": "Engreitz",
}

overview = pd.read_csv(BASE / "dataset_overview.tsv", sep="\t")
expected = dict(zip(overview["dataset"], overview["expected_cells"]))

# Load metrics_summary.tsv from each run
frames = []
for run, label in RUNS:
    f = RESULTS / run / "metrics_summary.tsv"
    if not f.exists():
        print(f"  skip {label}: {f} missing")
        continue
    df = pd.read_csv(f, sep="\t")[["dataset", "n_cells"]].copy()
    df["run"] = run
    df["run_label"] = label
    frames.append(df)
metrics = pd.concat(frames, ignore_index=True)

metrics["expected_cells"] = metrics["dataset"].map(expected)
metrics["ratio_observed_to_expected"] = metrics["n_cells"] / metrics["expected_cells"]
metrics.to_csv(OUT / "cells_vs_expected.tsv", sep="\t", index=False)
print(f"Saved table: cells_vs_expected.tsv")
print(metrics.pivot_table(index="dataset", columns="run_label", values="n_cells", aggfunc="first")[RUN_LABELS].round(0).to_string())


def fmt_k(x, _):
    return f"{x/1000:.0f}K"


# Small multiples — one panel per dataset
ds_in = [d for d in DATASET_ORDER if d in metrics["dataset"].values]
n_ds = len(ds_in)
ncols = min(3, n_ds)
nrows = int(np.ceil(n_ds / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), sharey=False)
axes = np.atleast_2d(axes).flatten()

for ax, ds in zip(axes, ds_in):
    sub = metrics[metrics["dataset"] == ds].set_index("run_label").reindex(RUN_LABELS)
    color = DATASET_COLORS[ds]
    x = np.arange(len(RUN_LABELS))
    bars = ax.bar(x, sub["n_cells"].values, color=color, alpha=0.85, edgecolor="white")
    exp = expected.get(ds, np.nan)
    if not np.isnan(exp):
        ax.axhline(exp, color="black", linestyle="--", linewidth=1.5, label=f"Expected ({exp/1000:.0f}K)")
    # annotate ratio above each bar
    for xi, val in zip(x, sub["n_cells"].values):
        if pd.notna(val) and not np.isnan(exp):
            ratio = val / exp
            ax.text(xi, val * 1.02, f"{ratio:.1f}×", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(RUN_LABELS, rotation=30, ha="right", fontsize=9)
    ax.set_title(SHORT.get(ds, ds), fontsize=11, color=color, fontweight="bold")
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_k))
    ax.set_ylabel("n_cells")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=8, loc="upper right")

for ax in axes[n_ds:]:
    ax.set_visible(False)

fig.suptitle("Cells per filter config vs experimental expectation", fontsize=13, fontweight="bold", y=1.00)
plt.tight_layout()
plt.savefig(OUT / "cells_vs_expected.pdf", dpi=300, bbox_inches="tight")
plt.savefig(OUT / "cells_vs_expected.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved: cells_vs_expected.pdf / .png")

# Overlaid view — log y, all 5 datasets together, expected as horizontal lines
fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(len(RUN_LABELS))
width = 0.16
n_ds_in = len(ds_in)

for j, ds in enumerate(ds_in):
    sub = metrics[metrics["dataset"] == ds].set_index("run_label").reindex(RUN_LABELS)
    color = DATASET_COLORS[ds]
    offset = (j - (n_ds_in - 1) / 2) * width
    ax.bar(x + offset, sub["n_cells"].values, width, color=color, alpha=0.85, edgecolor="white",
           label=SHORT.get(ds, ds))
    exp = expected.get(ds, np.nan)
    if not np.isnan(exp):
        ax.hlines(exp, x[0] + offset - width / 2, x[-1] + offset + width / 2,
                  colors=color, linestyles="--", linewidth=1.5, alpha=0.9)

ax.set_xticks(x)
ax.set_xticklabels(RUN_LABELS, fontsize=10)
ax.set_ylabel("n_cells (log scale)")
ax.set_yscale("log")
ax.yaxis.set_major_formatter(FuncFormatter(fmt_k))
ax.set_title("Cells per filter config (bars) vs expected (dashed lines)", fontsize=12, fontweight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), title="Dataset", fontsize=9)
plt.tight_layout()
plt.savefig(OUT / "cells_vs_expected_overlaid.pdf", dpi=300, bbox_inches="tight")
plt.savefig(OUT / "cells_vs_expected_overlaid.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved: cells_vs_expected_overlaid.pdf / .png")

# Ratio plot — per dataset, observed / expected as a function of filter config
fig, ax = plt.subplots(figsize=(10, 5))
for ds in ds_in:
    sub = metrics[metrics["dataset"] == ds].set_index("run_label").reindex(RUN_LABELS)
    color = DATASET_COLORS[ds]
    ax.plot(RUN_LABELS, sub["ratio_observed_to_expected"].values, "o-", color=color,
            label=SHORT.get(ds, ds), linewidth=2, markersize=8)

ax.axhline(1.0, color="black", linestyle="--", linewidth=1, alpha=0.5, label="1× (matches expected)")
ax.set_ylabel("observed / expected cells")
ax.set_title("Cell-yield ratio vs experimental expectation", fontsize=12, fontweight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), title="Dataset", fontsize=9)
plt.tight_layout()
plt.savefig(OUT / "cells_vs_expected_ratio.pdf", dpi=300, bbox_inches="tight")
plt.savefig(OUT / "cells_vs_expected_ratio.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved: cells_vs_expected_ratio.pdf / .png")

print(f"\nDone. Outputs in: {OUT}")
