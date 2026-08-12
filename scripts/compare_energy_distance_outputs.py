#!/usr/bin/env python3
"""
compare_energy_distance.py
---------------------------
For each dataset directory, compare three runs side-by-side:
  - main run   (dataset root)
  - adam_run   (dataset_root/adam_run/)
  - chikara_run(dataset_root/chikara_run/)

Each dataset produces its own set of plots saved to <outdir>/<dataset_label>/.

File schema:
  pval_edist_full.csv
    index (unnamed)  : gene/gRNA ID
    cell_count       : int
    type             : str  ("positive control", "targeting", etc.)
    distance_0..19   : float  (bootstrap replicates)
    pval_0..19       : float
    distance_mean    : float
    pval_mean        : float
    pval_mean_log    : float
    distance_mean_log: float

  targeting_outlier_table.csv / non_targeting_outlier_table.csv
    index (unnamed)  : gRNA ID
    pval_outlier     : float

Usage
-----
    python compare_energy_distance.py \\
        --datasets /path/to/data_cardio /path/to/data_de /path/to/data_stem \\
        [--labels cardio de stem] \\
        [--outdir ./comparison_plots] \\
        [--pval-cutoff 0.05]
"""

import argparse
import sys
import warnings
from pathlib import Path
from typing import Optional
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Style ─────────────────────────────────────────────────────────────────────
# One colour per run — consistent everywhere
RUN_COLORS = {
    "base":    "#4DBBD5",   # main run  — blue
    "adam":    "#E64B35",   # adam_run  — red
    "chikara": "#3CB44B",   # chikara_run — green
}
RUN_LABELS = {
    "base":    "Main run",
    "adam":    "adam_run",
    "chikara": "chikara_run",
}
RUN_LS = {
    "base":    "--",
    "adam":    "-",
    "chikara": ":",
}
RUNS = ["base", "adam", "chikara"]

plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.labelsize":    11,
    "axes.titlesize":    12,
    "xtick.labelsize":   9,
    "ytick.labelsize":   9,
    "figure.dpi":        150,
    "savefig.bbox":      "tight",
    "savefig.dpi":       150,
})

DIST_COLS = [f"distance_{i}" for i in range(20)]
PVAL_COLS = [f"pval_{i}"     for i in range(20)]
MEAN_DIST = "distance_mean"
MEAN_PVAL = "pval_mean"


# ── I/O ───────────────────────────────────────────────────────────────────────
def load_csv(path: Path, tag: str) -> Optional[pd.DataFrame]:
    if not path.exists():
        print(f"    [WARN] {tag}: not found -> {path}")
        return None
    df = pd.read_csv(path, index_col=0)
    df.columns = df.columns.str.strip()
    print(f"    [OK]   {tag}: {len(df)} rows")
    return df


def load_dataset(dataset_dir: Path, label: str) -> dict:
    sources = {
        "base":    dataset_dir,
        "adam":    dataset_dir / "adam_run",
        "chikara": dataset_dir / "chikara_run",
    }
    data = {"label": label, "dir": dataset_dir}
    for source, src_dir in sources.items():
        data[source] = {
            "pval_edist":  load_csv(src_dir / "pval_edist_full.csv",
                                    f"{source}/pval_edist_full"),
            "nt_outlier":  load_csv(src_dir / "non_targeting_outlier_table.csv",
                                    f"{source}/non_targeting_outlier"),
            "tgt_outlier": load_csv(src_dir / "targeting_outlier_table.csv",
                                    f"{source}/targeting_outlier"),
        }
    return data


# ── Helpers ───────────────────────────────────────────────────────────────────
def get_col(ds: dict, source: str, col: str) -> Optional[pd.Series]:
    df = ds[source]["pval_edist"]
    if df is None or col not in df.columns:
        return None
    return df[col]


def join_two(ds: dict, src_a: str, src_b: str, col: str) -> Optional[pd.DataFrame]:
    """Inner-join two runs on index for a single column."""
    a = get_col(ds, src_a, col)
    b = get_col(ds, src_b, col)
    if a is None or b is None:
        return None
    return (a.rename(src_a).to_frame()
             .join(b.rename(src_b).to_frame(), how="inner")
             .dropna())


def sig_set(ds: dict, source: str, cutoff: float) -> set:
    df = ds[source]["pval_edist"]
    if df is None or MEAN_PVAL not in df.columns:
        return set()
    return set(df.index[df[MEAN_PVAL] < cutoff])


def savefig(fig, path: Path) -> None:
    fig.savefig(path)
    plt.close(fig)
    print(f"    Saved: {path.name}")


def run_legend_handles():
    return [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=RUN_COLORS[r], markersize=8, label=RUN_LABELS[r])
        for r in RUNS
    ]


# ── Plot 1: pairwise scatter (3 pairs) ───────────────────────────────────────
def plot_scatter(ds: dict, outdir: Path) -> None:
    """
    3x2 grid: one row per run-pair, one column per metric
    (distance_mean, -log10(pval_mean)).
    """
    pairs = list(combinations(RUNS, 2))
    specs = [
        (MEAN_DIST, "distance_mean",     lambda x: x),
        (MEAN_PVAL, "-log10(pval_mean)", lambda x: -np.log10(x.clip(lower=1e-300))),
    ]
    fig, axes = plt.subplots(len(pairs), 2,
                             figsize=(10, 4 * len(pairs)),
                             squeeze=False)
    for row, (src_a, src_b) in enumerate(pairs):
        for col, (metric, mlabel, tfm) in enumerate(specs):
            ax = axes[row][col]
            m = join_two(ds, src_a, src_b, metric)
            if m is None or m.empty:
                ax.set_title("no data"); continue
            x, y = tfm(m[src_a]), tfm(m[src_b])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rho, pv = spearmanr(x, y)
            color = RUN_COLORS[src_a]
            ax.scatter(x, y, s=8, alpha=0.35, color=color, linewidths=0)
            lim = [min(x.min(), y.min()), max(x.max(), y.max())]
            ax.plot(lim, lim, "k--", lw=0.8, alpha=0.5)
            ax.set_xlabel(f"{RUN_LABELS[src_a]} — {mlabel}")
            ax.set_ylabel(f"{RUN_LABELS[src_b]} — {mlabel}")
            rho_str = f"{rho:.3f}" if not np.isnan(rho) else "N/A"
            pv_str  = f"{pv:.2e}"  if not np.isnan(pv)  else "N/A"
            ax.set_title(f"{RUN_LABELS[src_a]} vs {RUN_LABELS[src_b]}\n"
                         f"$\\rho$ = {rho_str},  p = {pv_str},  n = {len(m)}")

    fig.suptitle(f"{ds['label']} — pairwise scatter", fontsize=13)
    plt.tight_layout()
    savefig(fig, outdir / "scatter_pairwise.pdf")


# ── Plot 2: bootstrap variance violins ────────────────────────────────────────
def plot_bootstrap_variance(ds: dict, outdir: Path) -> None:
    """
    Violin of all 20 bootstrap replicate distances and p-values,
    one violin per run (3 violins per panel).
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, rep_cols, ylabel in [
        (axes[0], DIST_COLS, "Bootstrap distance"),
        (axes[1], PVAL_COLS, "Bootstrap p-value"),
    ]:
        parts, colors, tick_labels = [], [], []
        for source in RUNS:
            df = ds[source]["pval_edist"]
            if df is None:
                continue
            present = [c for c in rep_cols if c in df.columns]
            if not present:
                continue
            vals = df[present].values.flatten()
            vals = vals[~np.isnan(vals)]
            parts.append(vals)
            colors.append(RUN_COLORS[source])
            tick_labels.append(RUN_LABELS[source])

        if not parts:
            ax.set_title(f"{ylabel}\n(no data)"); continue

        vp = ax.violinplot(parts, positions=range(len(parts)),
                           showmedians=True, widths=0.6)
        for pc, col in zip(vp["bodies"], colors):
            pc.set_facecolor(col); pc.set_alpha(0.65)
        for key in ("cmedians", "cbars", "cmaxes", "cmins"):
            vp[key].set_color("black"); vp[key].set_linewidth(1)
        ax.set_xticks(range(len(parts)))
        ax.set_xticklabels(tick_labels, fontsize=9)
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)

    fig.suptitle(f"{ds['label']} — bootstrap variability", fontsize=13)
    plt.tight_layout()
    savefig(fig, outdir / "bootstrap_variance_violin.pdf")


# ── Plot 3: overlaid histograms ───────────────────────────────────────────────
def plot_histograms(ds: dict, outdir: Path) -> None:
    """
    Overlaid density histograms of distance_mean and -log10(pval_mean),
    all three runs on the same axes.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    specs = [
        (MEAN_DIST, "distance_mean",     lambda x: x),
        (MEAN_PVAL, "-log10(pval_mean)", lambda x: -np.log10(x.clip(lower=1e-300))),
    ]
    for ax, (col, label, tfm) in zip(axes, specs):
        for source in RUNS:
            df = ds[source]["pval_edist"]
            if df is None or col not in df.columns:
                continue
            vals = tfm(df[col].dropna())
            ax.hist(vals, bins=50, alpha=0.45,
                    color=RUN_COLORS[source],
                    linestyle=RUN_LS[source],
                    edgecolor="none", density=True,
                    label=RUN_LABELS[source])
            ax.axvline(float(np.median(vals)),
                       color=RUN_COLORS[source],
                       linestyle=RUN_LS[source],
                       lw=1.5, alpha=0.85)
        ax.set_xlabel(label); ax.set_ylabel("Density")
        ax.set_title(label); ax.legend(fontsize=9)

    fig.suptitle(f"{ds['label']} — distributions across runs", fontsize=13)
    plt.tight_layout()
    savefig(fig, outdir / "histograms_all_runs.pdf")


# ── Plot 4: concordance — upset-style bar + pairwise scatter ─────────────────
def plot_concordance(ds: dict, cutoff: float, outdir: Path) -> None:
    """
    Left: stacked bar of significant target counts broken down by which
    combination of runs they are significant in (7 possible non-empty subsets
    of {main, adam, chikara}).
    Right: 3-panel pairwise -log10(pval_mean) scatter coloured by concordance.
    """
    sig = {r: sig_set(ds, r, cutoff) for r in RUNS}
    all_targets = sig["base"] | sig["adam"] | sig["chikara"]

    # Assign each significant target to one of 7 categories
    category_counts = {}
    for target in all_targets:
        key = tuple(r for r in RUNS if target in sig[r])
        category_counts[key] = category_counts.get(key, 0) + 1

    # Build display order: all-three first, then pairs, then singles
    cat_order = [
        ("base", "adam", "chikara"),
        ("base", "adam"),
        ("base", "chikara"),
        ("adam", "chikara"),
        ("base",),
        ("adam",),
        ("chikara",),
    ]
    cat_labels = {
        ("base", "adam", "chikara"): "All three",
        ("base", "adam"):            "Main + adam",
        ("base", "chikara"):         "Main + chikara",
        ("adam", "chikara"):         "adam + chikara",
        ("base",):                   "Main only",
        ("adam",):                   "adam only",
        ("chikara",):                "chikara only",
    }
    # Colour: blend of the constituent run colours (use first run's colour)
    cat_colors = {
        ("base", "adam", "chikara"): "#888888",
        ("base", "adam"):            "#8E8FBF",
        ("base", "chikara"):         "#45BC8C",
        ("adam", "chikara"):         "#3F9070",
        ("base",):                   RUN_COLORS["base"],
        ("adam",):                   RUN_COLORS["adam"],
        ("chikara",):                RUN_COLORS["chikara"],
    }

    fig = plt.figure(figsize=(14, 5))
    # Left panel: bar chart
    ax_bar = fig.add_subplot(1, 4, 1)
    y_pos, bar_labels, bar_vals, bar_clrs = [], [], [], []
    for cat in cat_order:
        cnt = category_counts.get(cat, 0)
        if cnt == 0:
            continue
        y_pos.append(len(y_pos))
        bar_labels.append(cat_labels[cat])
        bar_vals.append(cnt)
        bar_clrs.append(cat_colors[cat])

    if bar_vals:
        bars = ax_bar.barh(y_pos, bar_vals, color=bar_clrs, alpha=0.85)
        ax_bar.set_yticks(y_pos)
        ax_bar.set_yticklabels(bar_labels, fontsize=8)
        ax_bar.set_xlabel("# significant targets")
        ax_bar.set_title(f"Hit overlap\n(pval_mean < {cutoff})")
        for bar, val in zip(bars, bar_vals):
            ax_bar.text(bar.get_width() + max(bar_vals) * 0.02,
                        bar.get_y() + bar.get_height() / 2,
                        str(val), va="center", fontsize=8)

    # Right panels: pairwise scatters (3 pairs)
    pairs = list(combinations(RUNS, 2))
    for idx, (src_a, src_b) in enumerate(pairs):
        ax = fig.add_subplot(1, 4, idx + 2)
        m = join_two(ds, src_a, src_b, MEAN_PVAL)
        if m is None or m.empty:
            ax.set_title("no data"); continue

        x = -np.log10(m[src_a].clip(lower=1e-300))
        y = -np.log10(m[src_b].clip(lower=1e-300))
        thresh = -np.log10(cutoff)

        # Colour by concordance: in both, in a only, in b only, neither
        sig_a = sig[src_a]; sig_b = sig[src_b]
        pt_colors = []
        for tgt in m.index:
            ia, ib = tgt in sig_a, tgt in sig_b
            if ia and ib:
                pt_colors.append("#888888")
            elif ia:
                pt_colors.append(RUN_COLORS[src_a])
            elif ib:
                pt_colors.append(RUN_COLORS[src_b])
            else:
                pt_colors.append("#DDDDDD")

        ax.scatter(x, y, s=7, alpha=0.5, c=pt_colors, linewidths=0)
        ax.axvline(thresh, color=RUN_COLORS[src_a], lw=0.8, ls="--", alpha=0.6)
        ax.axhline(thresh, color=RUN_COLORS[src_b], lw=0.8, ls="--", alpha=0.6)
        lim = [0, max(x.max(), y.max()) * 1.05]
        ax.plot(lim, lim, "k--", lw=0.7, alpha=0.4)
        ax.set_xlabel(f"{RUN_LABELS[src_a]}\n−log₁₀(pval_mean)", fontsize=8)
        ax.set_ylabel(f"{RUN_LABELS[src_b]}\n−log₁₀(pval_mean)", fontsize=8)
        ax.set_title(f"{RUN_LABELS[src_a]}\nvs {RUN_LABELS[src_b]}")

        legend_els = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#888888",
                   markersize=6, label="Both sig."),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=RUN_COLORS[src_a],
                   markersize=6, label=f"{RUN_LABELS[src_a]} only"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=RUN_COLORS[src_b],
                   markersize=6, label=f"{RUN_LABELS[src_b]} only"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#DDDDDD",
                   markersize=6, label="Not sig."),
        ]
        ax.legend(handles=legend_els, fontsize=6, loc="upper left")

    fig.suptitle(f"{ds['label']} — hit concordance across runs", fontsize=13)
    plt.tight_layout()
    savefig(fig, outdir / "concordance_all_runs.pdf")


# ── Plot 5: outlier p-value distributions ────────────────────────────────────
def plot_outliers(ds: dict, outdir: Path) -> None:
    """
    Overlaid pval_outlier histograms for all three runs,
    for both non-targeting and targeting outlier tables.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, key, title in [
        (axes[0], "nt_outlier",  "Non-targeting outliers"),
        (axes[1], "tgt_outlier", "Targeting outliers"),
    ]:
        for source in RUNS:
            df = ds[source][key]
            if df is None or "pval_outlier" not in df.columns:
                continue
            vals = df["pval_outlier"].dropna()
            ax.hist(vals, bins=30, alpha=0.5,
                    color=RUN_COLORS[source],
                    linestyle=RUN_LS[source],
                    edgecolor="none", density=True,
                    label=f"{RUN_LABELS[source]} (n={len(vals)})")
        ax.set_xlabel("pval_outlier"); ax.set_ylabel("Density")
        ax.set_title(title); ax.legend(fontsize=8)

    fig.suptitle(f"{ds['label']} — outlier p-value distributions", fontsize=13)
    plt.tight_layout()
    savefig(fig, outdir / "outlier_pval_distributions.pdf")


# ── Summary CSV ───────────────────────────────────────────────────────────────
def save_summary_csv(ds: dict, cutoff: float, outdir: Path) -> None:
    rows = []
    for source in RUNS:
        df = ds[source]["pval_edist"]
        row = {"run": RUN_LABELS[source]}
        if df is not None:
            row["n_targets"] = len(df)
            if MEAN_DIST in df.columns:
                row["distance_mean_median"] = round(df[MEAN_DIST].median(), 4)
                row["distance_mean_mean"]   = round(df[MEAN_DIST].mean(), 4)
            if MEAN_PVAL in df.columns:
                for thr in (0.1, 0.05, 0.01):
                    row[f"n_sig_{thr}"] = int((df[MEAN_PVAL] < thr).sum())
        nt  = ds[source]["nt_outlier"]
        tgt = ds[source]["tgt_outlier"]
        row["n_nt_outliers"]  = len(nt)  if nt  is not None else "N/A"
        row["n_tgt_outliers"] = len(tgt) if tgt is not None else "N/A"
        rows.append(row)

    out = outdir / "summary_stats.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"    Saved: {out.name}")


# ── Main ──────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Three-way comparison of main / adam_run / chikara_run "
                    "energy distance outputs, per dataset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--datasets", nargs="+", required=True,
                   help="One or more dataset root directories")
    p.add_argument("--labels", nargs="+", default=None,
                   help="Short label per dataset (default: directory basename)")
    p.add_argument("--outdir", default="./edist_comparison",
                   help="Root output directory; a subdirectory is created per dataset")
    p.add_argument("--pval-cutoff", type=float, default=0.05,
                   help="pval_mean threshold for 'significant' hits")
    return p.parse_args()


def main():
    args = parse_args()
    dataset_paths = [Path(d) for d in args.datasets]
    labels = args.labels or [p.name for p in dataset_paths]

    if len(labels) != len(dataset_paths):
        sys.exit(f"ERROR: {len(labels)} labels but {len(dataset_paths)} datasets")

    root_outdir = Path(args.outdir)
    root_outdir.mkdir(parents=True, exist_ok=True)

    for path, label in zip(dataset_paths, labels):
        print(f"\n{'='*55}")
        print(f"  Dataset: {label}  ({path})")
        print(f"{'='*55}")

        ds = load_dataset(path, label)
        outdir = root_outdir / label
        outdir.mkdir(parents=True, exist_ok=True)
        print(f"  -> Output: {outdir}\n")

        plot_scatter(ds, outdir)
        plot_bootstrap_variance(ds, outdir)
        plot_histograms(ds, outdir)
        plot_concordance(ds, args.pval_cutoff, outdir)
        plot_outliers(ds, outdir)
        save_summary_csv(ds, args.pval_cutoff, outdir)

    print(f"\n{'='*55}")
    print(f"  All done. Outputs in: {root_outdir}")
    print(f"{'='*55}\n")


if __name__ == "__main__":
    main()