"""Parameter-sweep small multiples per dataset.

Emits a multi-page PDF (all figures) plus individual PNGs per figure.
One panel per dataset; x = filtering parameter (500 / 800 / Knee);
lines colored by dataset; solid = primary, dashed = secondary.

Also emits per-lane variants for metrics that are available by lane.
"""
from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

PROJECT_ROOT = Path("/cellar/users/aklie/projects/tf_perturb_seq")
BASE = PROJECT_ROOT / "datasets" / "technology-benchmark_WTC11_TF-Perturb-seq"
RESULTS = BASE / "results" / "cross_tech_comparison"
OUT = RESULTS / "cross_parameter"
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

ENGREITZ_CAVEAT = "Engreitz Knee drop reflects lane-splitting issue (not biology)"


def load_combined(filename):
    frames = []
    for run, label in RUNS:
        f = RESULTS / run / filename
        if not f.exists():
            continue
        df = pd.read_csv(f, sep="\t")
        df["run"], df["run_label"] = run, label
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


metrics = load_combined("metrics_summary.tsv")
targets = load_combined("combined_intended_target_metrics.tsv")
guide_metrics = load_combined("combined_guide_metrics.tsv")
guide_all = guide_metrics[guide_metrics["batch"] == "all"].copy() if len(guide_metrics) else guide_metrics
lane_metrics = load_combined("metrics_by_lane.tsv")


def _fmt_axis(ax, col):
    if col == "n_cells":
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x/1000:.0f}K"))


def plot_across_params(df, y_cols, title, ylabel, out_stem,
                       bounded_01=False, add_engreitz_caveat=False,
                       twin_second=False):
    """Small multiples per dataset, across filtering params.

    y_cols: list of 1 or 2 column names (solid, dashed).
    """
    datasets = [d for d in DATASET_ORDER if d in df["dataset"].values]
    n = len(datasets)
    fig, axes = plt.subplots(1, n, figsize=(3.0 * n, 4.0), sharey=(not twin_second))
    if n == 1:
        axes = [axes]

    x = np.arange(len(RUN_LABELS))
    for ax, ds in zip(axes, datasets):
        color = DATASET_COLORS.get(ds, "#444444")
        sub = df[df["dataset"] == ds].set_index("run_label").reindex(RUN_LABELS)

        ax.plot(x, sub[y_cols[0]], marker="o", color=color, linewidth=2,
                label=y_cols[0])

        if len(y_cols) == 2:
            if twin_second:
                ax2 = ax.twinx()
                ax2.plot(x, sub[y_cols[1]], marker="s", color=color, linewidth=2,
                         linestyle="--", alpha=0.8, label=y_cols[1])
                ax2.tick_params(axis="y", labelsize=8, colors=color)
                ax2.spines["top"].set_visible(False)
            else:
                ax.plot(x, sub[y_cols[1]], marker="s", color=color, linewidth=2,
                        linestyle="--", alpha=0.8, label=y_cols[1])

        ax.set_xticks(x)
        ax.set_xticklabels(RUN_LABELS, rotation=30, ha="right", fontsize=9)
        ax.set_title(SHORT[ds], fontsize=11, fontweight="bold", color=color)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if bounded_01:
            ax.set_ylim(0.0, 1.05)
            ax.axhline(0.5, color="gray", linestyle=":", alpha=0.4, linewidth=0.8)
        _fmt_axis(ax, y_cols[0])

    axes[0].set_ylabel(ylabel)
    # Legend (solid/dashed) — put on first axis
    if len(y_cols) == 2:
        handles = [
            plt.Line2D([0], [0], color="black", linewidth=2, label=y_cols[0]),
            plt.Line2D([0], [0], color="black", linewidth=2, linestyle="--",
                       label=y_cols[1]),
        ]
        axes[0].legend(handles=handles, loc="lower left", fontsize=8, frameon=False)

    suptitle = title
    if add_engreitz_caveat:
        suptitle += f"\n({ENGREITZ_CAVEAT})"
    fig.suptitle(suptitle, fontsize=12, fontweight="bold", y=1.02)
    plt.tight_layout()

    png = OUT / f"{out_stem}.png"
    plt.savefig(png, dpi=200, bbox_inches="tight")
    print(f"  PNG: {png.name}")
    return fig


def plot_per_lane(df, y_col, title, ylabel, out_stem, add_engreitz_caveat=False):
    """Per-lane lines, one panel per dataset. Thin lines per batch; dataset color."""
    datasets = [d for d in DATASET_ORDER if d in df["dataset"].values]
    n = len(datasets)
    fig, axes = plt.subplots(1, n, figsize=(3.0 * n, 4.0), sharey=True)
    if n == 1:
        axes = [axes]

    x = np.arange(len(RUN_LABELS))
    for ax, ds in zip(axes, datasets):
        color = DATASET_COLORS.get(ds, "#444444")
        ds_df = df[(df["dataset"] == ds) & (df["batch"] != "all")]
        lanes = sorted(ds_df["batch"].unique())
        for lane in lanes:
            sub = (ds_df[ds_df["batch"] == lane]
                   .set_index("run_label").reindex(RUN_LABELS))
            ax.plot(x, sub[y_col], marker="o", color=color, linewidth=1.2,
                    alpha=0.55, markersize=4)
        # dataset-level "all" in bold
        all_df = df[(df["dataset"] == ds) & (df["batch"] == "all")]
        if len(all_df):
            sub = all_df.set_index("run_label").reindex(RUN_LABELS)
            ax.plot(x, sub[y_col], marker="o", color=color, linewidth=2.5,
                    label="all")

        ax.set_xticks(x)
        ax.set_xticklabels(RUN_LABELS, rotation=30, ha="right", fontsize=9)
        ax.set_title(f"{SHORT[ds]} (n={len(lanes)} lanes)", fontsize=10,
                     fontweight="bold", color=color)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        _fmt_axis(ax, y_col)

    axes[0].set_ylabel(ylabel)
    suptitle = title
    if add_engreitz_caveat:
        suptitle += f"\n({ENGREITZ_CAVEAT})"
    fig.suptitle(suptitle, fontsize=12, fontweight="bold", y=1.02)
    plt.tight_layout()
    png = OUT / f"{out_stem}.png"
    plt.savefig(png, dpi=200, bbox_inches="tight")
    print(f"  PNG: {png.name}")
    return fig


# Build: guides_per_cell_mean needs guide_all merged with metrics_summary slots
# Merge mean guides/cell into metrics frame
if len(guide_all):
    gpc = guide_all[["dataset", "run_label", "guides_per_cell_mean"]]
    metrics = metrics.merge(gpc, on=["dataset", "run_label"], how="left")

# ---- Build figures ----
figs = []

print("Building figures...")
figs.append(plot_across_params(
    targets, ["auroc", "auprc"],
    "Knockdown efficiency across filtering parameters",
    "auROC / auPRC (dashed)", "sweep_auroc_auprc",
    bounded_01=True, add_engreitz_caveat=True,
))
figs.append(plot_across_params(
    targets, ["frac_significant", "frac_strong_and_significant"],
    "Hit rate across filtering parameters",
    "Fraction (dashed = strong + sig)", "sweep_hit_rate",
    bounded_01=True, add_engreitz_caveat=True,
))
figs.append(plot_across_params(
    targets, ["median_log2fc"],
    "Median knockdown log2FC (targeting guides)",
    "median log2FC", "sweep_median_log2fc",
    add_engreitz_caveat=True,
))
figs.append(plot_across_params(
    metrics, ["n_cells"],
    "Cell yield across filtering parameters",
    "n cells", "sweep_n_cells",
    add_engreitz_caveat=True,
))
figs.append(plot_across_params(
    metrics, ["umi_median", "genes_median"],
    "RNA depth across filtering parameters",
    "median per cell (dashed = genes)", "sweep_rna_depth",
))
figs.append(plot_across_params(
    metrics, ["mito_median"],
    "Median % mito across filtering parameters",
    "% mito", "sweep_mito",
))
figs.append(plot_across_params(
    metrics, ["frac_cells_with_guide", "guide_umi_median"],
    "Guide assignment across filtering parameters",
    "frac cells w/ guide (dashed = guide UMI median, right)",
    "sweep_guide_assignment", twin_second=True,
))
figs.append(plot_across_params(
    metrics, ["guides_per_cell_mean"],
    "Mean guides assigned per cell",
    "guides / cell (mean)", "sweep_guides_per_cell_mean",
))

# Cross-dataset correlation distribution across params
print("Cross-dataset correlation distributions...")
from itertools import combinations
from scipy.stats import pearsonr

target_results = load_combined("combined_intended_target_results.tsv")
corr_rows = []
if len(target_results):
    tr = target_results[target_results["targeting"] == True].copy()
    for run, label in RUNS:
        sub = tr[tr["run"] == run]
        if not len(sub):
            continue
        pivot = sub.pivot_table(index="guide_id", columns="dataset",
                                values="log2_fc", aggfunc="mean")
        ds_cols = [d for d in pivot.columns if d in DATASET_ORDER]
        for a, b in combinations(ds_cols, 2):
            shared = pivot[[a, b]].dropna()
            if len(shared) < 10:
                continue
            r, _ = pearsonr(shared[a], shared[b])
            corr_rows.append({"run_label": label, "pair": f"{SHORT[a]}↔{SHORT[b]}",
                              "a": a, "b": b, "r": r, "n": len(shared)})

if corr_rows:
    corr_df = pd.DataFrame(corr_rows)
    corr_df.to_csv(OUT / "cross_dataset_log2fc_correlations.tsv",
                   sep="\t", index=False)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5),
                                    gridspec_kw={"width_ratios": [1, 2]})
    # Left: distribution (strip + box) per param
    xs = np.arange(len(RUN_LABELS))
    for i, lbl in enumerate(RUN_LABELS):
        vals = corr_df[corr_df["run_label"] == lbl]["r"].values
        if len(vals):
            ax1.boxplot([vals], positions=[i], widths=0.5,
                        showfliers=False, patch_artist=True,
                        boxprops=dict(facecolor="#e0e0e0", edgecolor="#444"),
                        medianprops=dict(color="black"))
            jitter = np.random.RandomState(0).uniform(-0.12, 0.12, len(vals))
            ax1.scatter(np.full(len(vals), i) + jitter, vals,
                        color="#333", s=24, alpha=0.7, zorder=3)
    ax1.set_xticks(xs)
    ax1.set_xticklabels(RUN_LABELS, rotation=30, ha="right")
    ax1.set_ylabel("Pearson r (pairwise dataset log2FC)")
    ax1.set_title("Correlation distribution", fontsize=11, fontweight="bold")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.set_ylim(0, 1)

    # Right: per-pair lines across params, colored by "dominant" dataset
    pairs = sorted(corr_df["pair"].unique())
    for pair in pairs:
        sub = corr_df[corr_df["pair"] == pair].set_index("run_label").reindex(RUN_LABELS)
        a = corr_df[corr_df["pair"] == pair]["a"].iloc[0]
        b = corr_df[corr_df["pair"] == pair]["b"].iloc[0]
        # blend of two dataset colors
        ca = DATASET_COLORS.get(a, "#444")
        cb = DATASET_COLORS.get(b, "#444")
        ax2.plot(xs, sub["r"], marker="o", linewidth=1.8,
                 color=ca, alpha=0.9, label=pair,
                 markerfacecolor=cb, markeredgecolor=ca, markersize=8)
    ax2.set_xticks(xs)
    ax2.set_xticklabels(RUN_LABELS, rotation=30, ha="right")
    ax2.set_ylabel("Pearson r")
    ax2.set_title("Per-pair correlation across params", fontsize=11, fontweight="bold")
    ax2.legend(fontsize=7, loc="lower left", ncol=2, frameon=False)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.set_ylim(0, 1)

    fig.suptitle("Cross-dataset log2FC correlation across filtering parameters",
                 fontsize=12, fontweight="bold", y=1.02)
    plt.tight_layout()
    png = OUT / "sweep_cross_dataset_correlation.png"
    plt.savefig(png, dpi=200, bbox_inches="tight")
    print(f"  PNG: {png.name}")
    figs.append(fig)

# Per-lane variants (where lane data exists)
if len(lane_metrics):
    print("Per-lane figures...")
    # lane_metrics already has batch != 'all'; we also want the 'all' dataset-level
    # from metrics_summary for the overlay
    ms_all = metrics[["dataset", "run_label", "n_cells", "umi_median",
                      "genes_median", "mito_median", "guide_umi_median",
                      "guides_per_cell_mean"]].copy()
    ms_all["batch"] = "all"
    lane_plus_all = pd.concat([lane_metrics, ms_all], ignore_index=True)

    for col, title, ylab, stem in [
        ("n_cells", "Cell yield per lane", "n cells", "sweep_lane_n_cells"),
        ("umi_median", "Median UMI per lane", "median UMI/cell",
         "sweep_lane_umi_median"),
        ("mito_median", "Median % mito per lane", "% mito",
         "sweep_lane_mito"),
        ("guide_umi_median", "Median guide UMI per lane",
         "median guide UMI/cell", "sweep_lane_guide_umi"),
        ("guides_per_cell_mean", "Mean guides/cell per lane",
         "guides / cell (mean)", "sweep_lane_guides_per_cell_mean"),
    ]:
        figs.append(plot_per_lane(lane_plus_all, col, title, ylab, stem,
                                  add_engreitz_caveat=(col == "n_cells")))

# Multi-page PDF
pdf_path = OUT / "param_sweep_all_metrics.pdf"
with PdfPages(pdf_path) as pdf:
    for fig in figs:
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
print(f"\nPDF: {pdf_path}")
print(f"Wrote {len(figs)} figures.")
