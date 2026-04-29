"""Per-lane mapping-rate comparison for the two DE runs.

Reads /tmp/de_compare/mapping_summary.tsv (compiled from each run's
mappingscrna/<MS>/run_info.json and mappingguide/<MS>/run_info.json) and
produces per-lane plots of:
  - p_pseudoaligned (mapping rate)
  - n_pseudoaligned vs n_processed (read counts)

Stratifies by modality (scRNA vs guide). Only spacer_tag differs between runs,
so scRNA mapping is identical and guide mapping is the dramatic story.
"""
from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

PROJECT_ROOT = Path("/cellar/users/aklie/projects/tf_perturb_seq")
DATASET = "Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq"
OUT = PROJECT_ROOT / "datasets" / DATASET / "results" / "de_run_comparison"
OUT.mkdir(parents=True, exist_ok=True)

COLORS = {"hamster": "#9aa0a6", "penguin": "#1a73e8"}
LABELS = {"hamster": "hamster (12 bp)", "penguin": "penguin (13 bp)"}

df = pd.read_csv("/tmp/de_compare/mapping_summary.tsv", sep="\t")
# Clean MS name
df["ms_short"] = df["measurement_set"].str.replace("_ks_guide_out", "", regex=False)
print(f"Loaded {len(df)} rows")
print(df.head())

ms_order = sorted(df["ms_short"].unique())


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def fmt_M(x, _):
    return f"{x/1e6:.0f}M" if x >= 1e6 else f"{x/1e3:.0f}K"


# ============ Plot 5: per-lane mapping rate ============
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
modalities = [("mappingscrna", "scRNA mapping rate"), ("mappingguide", "Guide mapping rate")]

for ax, (mod, title) in zip(axes, modalities):
    sub = df[df["modality"] == mod]
    x = np.arange(len(ms_order))
    width = 0.4
    for j, run in enumerate(["hamster", "penguin"]):
        vals = []
        for ms in ms_order:
            r = sub[(sub["run"] == run) & (sub["ms_short"] == ms)]
            vals.append(r["p_pseudoaligned"].iloc[0] if len(r) else np.nan)
        ax.bar(x + (j - 0.5) * width, vals, width, color=COLORS[run], alpha=0.9,
               edgecolor="white", label=LABELS[run])
    ax.set_xticks(x)
    ax.set_xticklabels(ms_order, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("% reads pseudoaligned")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_ylim(0, 100)
    ax.axhline(50, color="black", linestyle="--", linewidth=0.5, alpha=0.3)
    ax.legend(fontsize=9, loc="upper right" if mod == "mappingscrna" else "upper left")
    style(ax)

fig.suptitle("DE production: per-lane mapping rate (hamster 12 bp vs penguin 13 bp)",
             fontsize=12, fontweight="bold", y=1.00)
plt.tight_layout()
plt.savefig(OUT / "05_per_lane_mapping_rate.pdf", dpi=300, bbox_inches="tight")
plt.savefig(OUT / "05_per_lane_mapping_rate.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved: 05_per_lane_mapping_rate.pdf / .png")


# ============ Plot 6: read counts (processed + pseudoaligned) ============
fig, axes = plt.subplots(2, 2, figsize=(14, 9))

for col, (mod, title) in enumerate(modalities):
    sub = df[df["modality"] == mod]
    # Top row: n_processed (same per MS but plot for visual completeness)
    ax = axes[0, col]
    x = np.arange(len(ms_order))
    width = 0.4
    for j, run in enumerate(["hamster", "penguin"]):
        vals = []
        for ms in ms_order:
            r = sub[(sub["run"] == run) & (sub["ms_short"] == ms)]
            vals.append(r["n_processed"].iloc[0] if len(r) else np.nan)
        ax.bar(x + (j - 0.5) * width, vals, width, color=COLORS[run], alpha=0.9,
               edgecolor="white", label=LABELS[run])
    ax.set_xticks(x)
    ax.set_xticklabels(ms_order, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Reads processed")
    ax.set_title(f"{mod}: n_processed (input)", fontsize=10)
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_M))
    if col == 0:
        ax.legend(fontsize=9, loc="upper right")
    style(ax)

    # Bottom row: n_pseudoaligned (the key difference for guide)
    ax = axes[1, col]
    for j, run in enumerate(["hamster", "penguin"]):
        vals = []
        for ms in ms_order:
            r = sub[(sub["run"] == run) & (sub["ms_short"] == ms)]
            vals.append(r["n_pseudoaligned"].iloc[0] if len(r) else np.nan)
        ax.bar(x + (j - 0.5) * width, vals, width, color=COLORS[run], alpha=0.9,
               edgecolor="white", label=LABELS[run])
    ax.set_xticks(x)
    ax.set_xticklabels(ms_order, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Reads pseudoaligned")
    ax.set_title(f"{mod}: n_pseudoaligned (mapped)", fontsize=10)
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_M))
    style(ax)

fig.suptitle("DE production: per-lane read counts (input vs mapped)",
             fontsize=12, fontweight="bold", y=1.00)
plt.tight_layout()
plt.savefig(OUT / "06_per_lane_read_counts.pdf", dpi=300, bbox_inches="tight")
plt.savefig(OUT / "06_per_lane_read_counts.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved: 06_per_lane_read_counts.pdf / .png")


# ============ Aggregate summary ============
agg = df.groupby(["run", "modality"]).agg(
    n_processed=("n_processed", "sum"),
    n_pseudoaligned=("n_pseudoaligned", "sum"),
    mean_p_pseudoaligned=("p_pseudoaligned", "mean"),
    median_p_pseudoaligned=("p_pseudoaligned", "median"),
).reset_index()
agg["aggregate_p_pseudoaligned"] = 100 * agg["n_pseudoaligned"] / agg["n_processed"]
agg.to_csv(OUT / "mapping_aggregate.tsv", sep="\t", index=False)
print()
print("=== Aggregate mapping summary ===")
print(agg.round(2).to_string(index=False))
print(f"\nDone. Outputs in: {OUT}")
