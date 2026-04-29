"""Careful comparison of the two Huangfu DE production runs.

Both at gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/2026_04_09/outs/

Verified from pipeline_info/params_*.json (everything else identical):
    entertaining_hamster: spacer_tag = "GAGTACATGGGG"   (12 bp)
    muddy_penguin:        spacer_tag = "GAGTACATGGGGG"  (13 bp)

Pulls the additional_qc/ TSVs from each run's pipeline_dashboard and produces a
side-by-side comparison covering cell yield, guide capture, knockdown sensitivity,
trans signal, and per-class log2FC. Saves figures + a numeric summary TSV.

Outputs:
    results/de_run_comparison/
      summary_metrics.tsv         — pivot of headline numbers
      01_cell_and_guide_yield.pdf — n_cells, UMI, guide UMI, frac with guide, guides/cell
      02_knockdown_sensitivity.pdf — AUROC, AUPRC, frac_strong, frac_significant, log2FC
      03_trans_signal.pdf          — significant tests, hits per guide
      04_per_class_log2fc.pdf      — intended-target + trans log2FC by guide class
      README.md                    — narrative summary
"""
from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

PROJECT_ROOT = Path("/cellar/users/aklie/projects/tf_perturb_seq")
DATASET = "Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq"
BASE = PROJECT_ROOT / "datasets" / DATASET
OUT = BASE / "results" / "de_run_comparison"
OUT.mkdir(parents=True, exist_ok=True)

# Local paths (synced yesterday from gs://...muddy_penguin/pipeline_dashboard/additional_qc and entertaining_hamster/...)
QC_PATHS = {
    "hamster (12 bp)": Path("/tmp/de_compare/hamster"),
    "penguin (13 bp)": Path("/tmp/de_compare/penguin"),
}
COLORS = {
    "hamster (12 bp)": "#9aa0a6",  # gray (older)
    "penguin (13 bp)": "#1a73e8",  # blue (newer / canonical)
}
RUN_ORDER = list(QC_PATHS.keys())

POSITIVE_CONTROL_GENES = {"CD151", "CD55", "CD81", "NGFRAP1", "TFRC"}


def classify_guides(gene_name, targeting):
    is_target = targeting == True
    is_pc = gene_name.isin(POSITIVE_CONTROL_GENES) & is_target
    is_or = gene_name.str.match(r"^OR[0-9]", na=False) & is_target
    return np.where(
        ~is_target, "non_targeting",
        np.where(is_pc, "positive_control",
                 np.where(is_or, "negative_control_OR", "tf_targeting")))


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def fmt_k(x, _):
    return f"{x/1000:.0f}K" if x >= 1000 else f"{x:.0f}"


# === Load ===
def read_batch_all(path):
    df = pd.read_csv(path, sep="\t")
    return df[df["batch"] == "all"].iloc[0]


rows = []
for run, qc in QC_PATHS.items():
    g = read_batch_all(qc / "gene" / "gene_metrics.tsv")
    gd = read_batch_all(qc / "guide" / "guide_metrics.tsv")
    it = pd.read_csv(qc / "intended_target" / "intended_target_metrics.tsv", sep="\t").iloc[0]
    tr = pd.read_csv(qc / "trans" / "trans_metrics.tsv", sep="\t").iloc[0]
    rows.append({
        "run": run,
        # gene
        "n_cells": g["n_cells"],
        "umi_median_per_cell": g["umi_median"],
        "umi_mean_per_cell": g["umi_mean"],
        "mito_median": g["mito_median"],
        # guide
        "guide_umi_median": gd["guide_umi_median"],
        "guide_umi_mean": gd["guide_umi_mean"],
        "guides_per_cell_mean": gd["guides_per_cell_mean"],
        "guides_per_cell_median": gd["guides_per_cell_median"],
        "frac_cells_with_guide": gd["frac_cells_with_guide"],
        "n_cells_with_guide": gd["n_cells_with_guide"],
        "n_cells_exactly_1_guide": gd["n_cells_exactly_1_guide"],
        "cells_per_guide_median": gd["cells_per_guide_median"],
        # intended target
        "n_guides_tested": it["n_guides_tested"],
        "auroc": it["auroc"],
        "auprc": it["auprc"],
        "median_log2fc_intended": it["median_log2fc"],
        "frac_strong_knockdowns": it["frac_strong_knockdowns"],
        "frac_significant": it["frac_significant"],
        "n_strong_and_significant": it["n_strong_and_significant"],
        # trans
        "trans_total_significant": tr["total_significant_tests"],
        "trans_mean_significant_per_targeting": tr["mean_significant_per_guide_targeting"],
        "trans_median_significant_per_targeting": tr["median_significant_per_guide_targeting"],
        "trans_mean_significant_per_nt": tr["mean_significant_per_guide_nt"],
        "trans_median_genome_log2fc_targeting": tr["median_genome_log2fc_targeting"],
        "trans_median_genome_log2fc_nt": tr["median_genome_log2fc_nt"],
    })

summary = pd.DataFrame(rows).set_index("run").reindex(RUN_ORDER)
# Compute deltas (penguin minus hamster)
deltas = summary.loc["penguin (13 bp)"] - summary.loc["hamster (12 bp)"]
ratios = summary.loc["penguin (13 bp)"] / summary.loc["hamster (12 bp)"]
summary.loc["Δ (penguin − hamster)"] = deltas
summary.loc["× (penguin / hamster)"] = ratios
summary.to_csv(OUT / "summary_metrics.tsv", sep="\t")
print(f"Saved: summary_metrics.tsv")
print()
print("=== Headline numbers (batch=all) ===")
key_cols = ["n_cells", "umi_median_per_cell", "guide_umi_median", "guides_per_cell_mean",
            "frac_cells_with_guide", "auroc", "auprc", "frac_significant",
            "trans_total_significant"]
print(summary[key_cols].round(3).T.to_string())


# === Plot 1: cell + guide yield ===
panels1 = [
    ("n_cells", "Total cells", "k", None),
    ("umi_median_per_cell", "Median UMI / cell", None, None),
    ("guide_umi_median", "Median guide UMI / cell", None, None),
    ("frac_cells_with_guide", "Fraction cells with guide", None, None),
    ("guides_per_cell_mean", "Mean guides / cell", None, None),
    ("cells_per_guide_median", "Median cells / guide", None, None),
]
fig, axes = plt.subplots(2, 3, figsize=(13, 7))
axes = axes.flatten()
for ax, (col, title, fmt, _) in zip(axes, panels1):
    vals = [summary.loc[r, col] for r in RUN_ORDER]
    colors = [COLORS[r] for r in RUN_ORDER]
    bars = ax.bar(RUN_ORDER, vals, color=colors, alpha=0.9, edgecolor="white")
    for r, v in zip(RUN_ORDER, vals):
        txt = f"{v:,.0f}" if v > 100 else f"{v:.3f}"
        ax.text(RUN_ORDER.index(r), v, txt, ha="center", va="bottom", fontsize=9)
    if fmt == "k":
        ax.yaxis.set_major_formatter(FuncFormatter(fmt_k))
    ax.set_title(title, fontsize=10)
    ax.tick_params(axis="x", labelsize=8)
    style(ax)
fig.suptitle("DE production runs: cell + guide yield (only spacer_tag differs)", fontsize=12, fontweight="bold", y=1.00)
plt.tight_layout()
plt.savefig(OUT / "01_cell_and_guide_yield.pdf", dpi=300, bbox_inches="tight")
plt.savefig(OUT / "01_cell_and_guide_yield.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved: 01_cell_and_guide_yield.pdf / .png")


# === Plot 2: knockdown sensitivity ===
panels2 = [
    ("auroc", "auROC (intended target)", 0.5),
    ("auprc", "auPRC (intended target)", None),
    ("frac_strong_knockdowns", "Frac strong knockdowns\n(log2FC < log2(0.4))", None),
    ("frac_significant", "Frac significant\n(p<0.05)", None),
    ("n_strong_and_significant", "N strong + significant", None),
    ("median_log2fc_intended", "Median log2FC\n(intended target)", 0),
]
fig, axes = plt.subplots(2, 3, figsize=(13, 7))
axes = axes.flatten()
for ax, (col, title, hline) in zip(axes, panels2):
    vals = [summary.loc[r, col] for r in RUN_ORDER]
    colors = [COLORS[r] for r in RUN_ORDER]
    ax.bar(RUN_ORDER, vals, color=colors, alpha=0.9, edgecolor="white")
    for r, v in zip(RUN_ORDER, vals):
        txt = f"{v:,.0f}" if abs(v) > 100 else f"{v:.3f}"
        ax.text(RUN_ORDER.index(r), v if v >= 0 else v - 0.005, txt,
                ha="center", va="bottom" if v >= 0 else "top", fontsize=9)
    if hline is not None:
        ax.axhline(hline, color="black", linestyle="--", linewidth=0.5, alpha=0.4)
    ax.set_title(title, fontsize=10)
    ax.tick_params(axis="x", labelsize=8)
    style(ax)
fig.suptitle("DE production runs: knockdown sensitivity", fontsize=12, fontweight="bold", y=1.00)
plt.tight_layout()
plt.savefig(OUT / "02_knockdown_sensitivity.pdf", dpi=300, bbox_inches="tight")
plt.savefig(OUT / "02_knockdown_sensitivity.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved: 02_knockdown_sensitivity.pdf / .png")


# === Plot 3: trans signal ===
panels3 = [
    ("trans_total_significant", "Total significant\ntrans tests (FDR<0.05)", None),
    ("trans_mean_significant_per_targeting", "Mean significant trans hits\nper targeting guide", None),
    ("trans_median_significant_per_targeting", "Median significant trans hits\nper targeting guide", None),
    ("trans_mean_significant_per_nt", "Mean significant trans hits\nper non-targeting guide (null)", None),
    ("trans_median_genome_log2fc_targeting", "Median trans log2FC\n(targeting)", 0),
    ("trans_median_genome_log2fc_nt", "Median trans log2FC\n(non-targeting null)", 0),
]
fig, axes = plt.subplots(2, 3, figsize=(13, 7))
axes = axes.flatten()
for ax, (col, title, hline) in zip(axes, panels3):
    vals = [summary.loc[r, col] for r in RUN_ORDER]
    colors = [COLORS[r] for r in RUN_ORDER]
    ax.bar(RUN_ORDER, vals, color=colors, alpha=0.9, edgecolor="white")
    for r, v in zip(RUN_ORDER, vals):
        txt = f"{v:,.0f}" if abs(v) > 100 else f"{v:.4f}"
        ax.text(RUN_ORDER.index(r), v if v >= 0 else v * 0.95, txt,
                ha="center", va="bottom" if v >= 0 else "top", fontsize=9)
    if hline is not None:
        ax.axhline(hline, color="black", linestyle="--", linewidth=0.5, alpha=0.4)
    ax.set_title(title, fontsize=10)
    ax.tick_params(axis="x", labelsize=8)
    style(ax)
fig.suptitle("DE production runs: trans signal", fontsize=12, fontweight="bold", y=1.00)
plt.tight_layout()
plt.savefig(OUT / "03_trans_signal.pdf", dpi=300, bbox_inches="tight")
plt.savefig(OUT / "03_trans_signal.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved: 03_trans_signal.pdf / .png")


# === Plot 4: per-class log2FC ===
print("\nLoading per-class trans log2FC...")
class_rows = []
for run, qc in QC_PATHS.items():
    pg = pd.read_csv(qc / "trans" / "trans_per_guide_summary.tsv", sep="\t")
    pg["guide_class"] = classify_guides(pg["gene_name"], pg["targeting"])
    for cls in ["non_targeting", "negative_control_OR", "positive_control", "tf_targeting"]:
        sub = pg[pg["guide_class"] == cls]
        class_rows.append({
            "run": run, "guide_class": cls,
            "n_guides": len(sub),
            "trans_median_log2fc": sub["median_log2fc"].median() if len(sub) else np.nan,
            "trans_frac_significant": sub["frac_significant"].median() if len(sub) else np.nan,
        })
class_df = pd.DataFrame(class_rows)
print(class_df.round(4).to_string(index=False))

# Intended-target per-class median log2FC (from intended_target_results.tsv)
it_rows = []
for run, qc in QC_PATHS.items():
    df = pd.read_csv(qc / "intended_target" / "intended_target_results.tsv", sep="\t")
    df["guide_class"] = classify_guides(df["gene_name"], df["targeting"])
    for cls in ["positive_control", "tf_targeting"]:  # NT/OR have no intended target
        sub = df[df["guide_class"] == cls]
        it_rows.append({
            "run": run, "guide_class": cls,
            "intended_median_log2fc": sub["log2_fc"].median() if len(sub) else np.nan,
            "n_guides": len(sub),
        })
it_df = pd.DataFrame(it_rows)
print("\nIntended-target per-class median log2FC:")
print(it_df.round(4).to_string(index=False))

# Plot — 2x4 grid: two flavors x four classes
class_labels = ["non_targeting", "negative_control_OR", "positive_control", "tf_targeting"]
fig, axes = plt.subplots(2, 4, figsize=(15, 7))
# Top row: trans median log2FC
for j, cls in enumerate(class_labels):
    ax = axes[0, j]
    sub = class_df[class_df["guide_class"] == cls]
    vals = [sub[sub["run"] == r]["trans_median_log2fc"].iloc[0] if len(sub[sub["run"] == r]) else np.nan for r in RUN_ORDER]
    n = sub.iloc[0]["n_guides"] if len(sub) else 0
    colors = [COLORS[r] for r in RUN_ORDER]
    ax.bar(RUN_ORDER, vals, color=colors, alpha=0.9, edgecolor="white")
    for r, v in zip(RUN_ORDER, vals):
        if not np.isnan(v):
            ax.text(RUN_ORDER.index(r), v, f"{v:.4f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=8)
    ax.axhline(0, color="black", linestyle="--", linewidth=0.5, alpha=0.4)
    ax.set_title(f"trans log2FC\n{cls} (n={n})", fontsize=9)
    ax.tick_params(axis="x", labelsize=7)
    style(ax)
# Bottom row: intended-target log2FC (only for PC + TF)
for j, cls in enumerate(class_labels):
    ax = axes[1, j]
    if cls in ["positive_control", "tf_targeting"]:
        sub = it_df[it_df["guide_class"] == cls]
        vals = [sub[sub["run"] == r]["intended_median_log2fc"].iloc[0] if len(sub[sub["run"] == r]) else np.nan for r in RUN_ORDER]
        n = sub.iloc[0]["n_guides"] if len(sub) else 0
        colors = [COLORS[r] for r in RUN_ORDER]
        ax.bar(RUN_ORDER, vals, color=colors, alpha=0.9, edgecolor="white")
        for r, v in zip(RUN_ORDER, vals):
            if not np.isnan(v):
                ax.text(RUN_ORDER.index(r), v, f"{v:.3f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=8)
        ax.axhline(0, color="black", linestyle="--", linewidth=0.5, alpha=0.4)
        ax.set_title(f"intended log2FC\n{cls} (n={n})", fontsize=9)
    else:
        ax.text(0.5, 0.5, "N/A\n(no intended target)", ha="center", va="center",
                fontsize=10, color="gray", transform=ax.transAxes)
        ax.set_title(f"intended log2FC\n{cls}", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=7)
    style(ax)
fig.suptitle("DE production runs: per-class log2FC (top: trans median per guide; bottom: intended-target)",
             fontsize=11, fontweight="bold", y=1.00)
plt.tight_layout()
plt.savefig(OUT / "04_per_class_log2fc.pdf", dpi=300, bbox_inches="tight")
plt.savefig(OUT / "04_per_class_log2fc.png", dpi=300, bbox_inches="tight")
plt.close()
print("\nSaved: 04_per_class_log2fc.pdf / .png")
print(f"\nDone. Outputs in: {OUT}")
