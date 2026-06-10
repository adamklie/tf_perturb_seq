"""Per-dataset metric summary across filter configs.

Aggregates the user-requested metrics for each (dataset × filter config) into a
single wide table, plus a "best AUROC config per dataset" view, plus a panel
of dataset-colored bar plots showing each metric at the best config.

Metrics:
  Knockdown sensitivity:  auROC, auPRC, median log2fc (intended-target evaluation)
  Per guide class:        median observed FC and median log2fc, broken out by
                          {positive_control, tf_targeting, all_targeting,
                           non_targeting}
  scRNA / cell quality:   n_cells, umi_median (per cell), mito_median,
                          frac_cells_with_guide, ratio_observed_to_expected
  Guide assignment:       guides_per_cell_mean, frac_cells_exactly_1_guide

Caveats:
  - "Negative controls" and "non-targeting" are the same class in this data
    (label='targeting', targeting=False, gene_name='non-targeting')
  - Non-targeting guides have no intended-target FC computed → log2fc is N/A
  - Proportion of cells with exactly 2 / 3 / >3 guides assigned requires h5mu
    access — not in this script. See `compute_guide_assignment_breakdown.py`
    (TODO).

Outputs:
  results/cross_tech_comparison/per_dataset_summary/
    summary_wide.tsv       — one row per (dataset, config), all metrics as columns
    summary_best.tsv       — one row per dataset at its best AUROC config
    per_dataset_best.pdf   — multi-panel bar chart of selected metrics, dataset-colored
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
OUT = RESULTS / "per_dataset_summary"
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

SHORT = {
    "Hon_WTC11-benchmark_TF-Perturb-seq": "Hon",
    "Huangfu_WTC11-benchmark_TF-Perturb-seq": "Huangfu",
    "Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3": "Gersbach\nGEM-Xv3",
    "Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2": "Gersbach\nHTv2",
    "Engreitz_WTC11-benchmark_TF-Perturb-seq": "Engreitz",
}

overview = pd.read_csv(BASE / "dataset_overview.tsv", sep="\t")
expected = dict(zip(overview["dataset"], overview["expected_cells"]))

# Optional: per-cell guide-count breakdown from h5mu
breakdown_path = OUT / "guide_assignment_breakdown.tsv"
breakdown = pd.read_csv(breakdown_path, sep="\t") if breakdown_path.exists() else pd.DataFrame()

POSITIVE_CONTROL_GENES = {"CD151", "CD55", "CD81", "NGFRAP1", "TFRC"}


def classify_guides(gene_name, targeting):
    """Vectorized guide-class assignment. 5 classes."""
    is_target = targeting == True
    is_pc = gene_name.isin(POSITIVE_CONTROL_GENES) & is_target
    # Olfactory receptor genes: gene names starting with "OR" followed by a digit
    # (e.g. OR10J3, OR2A25). These are the negative controls in the TFP3 library.
    is_or = gene_name.str.match(r"^OR[0-9]", na=False) & is_target
    return np.where(
        ~is_target, "non_targeting",
        np.where(is_pc, "positive_control",
                 np.where(is_or, "negative_control_OR", "tf_targeting")))


def load_trans_per_guide(dataset, config):
    """Load per-guide trans summary for one (dataset, config). Returns df or None."""
    p = (BASE.parent / dataset / "runs" / config /
         "pipeline_dashboard" / "additional_qc" / "trans" / "trans_per_guide_summary.tsv")
    if not p.exists():
        return None
    df = pd.read_csv(p, sep="\t")
    df["guide_class"] = classify_guides(df["gene_name"], df["targeting"])
    return df


def per_dataset_class_stats(results_df):
    """Median log2fc and median fc per (dataset, guide_class).

    5 guide classes:
      - non_targeting (NT, targeting=False)
      - negative_control_OR (olfactory receptors, real genes that shouldn't matter)
      - positive_control (CD81/CD55/CD151/NGFRAP1/TFRC)
      - tf_targeting (the TFs of interest)
      - all_targeting (positive_control + tf_targeting; OR and NT excluded)
    """
    df = results_df.copy()
    df["guide_class"] = classify_guides(df["gene_name"], df["targeting"])
    classes = [df[df["guide_class"] == cls].assign(guide_class=cls)
               for cls in ["non_targeting", "negative_control_OR", "positive_control", "tf_targeting"]]
    # all_targeting = positive_control + tf_targeting (excludes OR and NT)
    allt = df[df["guide_class"].isin(["positive_control", "tf_targeting"])].copy()
    allt["guide_class"] = "all_targeting"
    classes.append(allt)
    return pd.concat(classes, ignore_index=True)


rows = []
for run, label in RUNS:
    run_dir = RESULTS / run

    # metrics_summary
    f = run_dir / "metrics_summary.tsv"
    if not f.exists():
        print(f"  skip {label}: no metrics_summary")
        continue
    msum = pd.read_csv(f, sep="\t")

    # intended_target_metrics
    f = run_dir / "combined_intended_target_metrics.tsv"
    itm = pd.read_csv(f, sep="\t") if f.exists() else pd.DataFrame()

    # intended_target_results — per-guide log2fc / fc
    f = run_dir / "combined_intended_target_results.tsv"
    itr = pd.read_csv(f, sep="\t") if f.exists() else pd.DataFrame()

    # guide_metrics — for guides_per_cell_mean, frac_cells_exactly_1_guide
    f = run_dir / "combined_guide_metrics.tsv"
    gm_full = pd.read_csv(f, sep="\t") if f.exists() else pd.DataFrame()
    gm = gm_full[gm_full["batch"] == "all"].copy() if not gm_full.empty else gm_full

    # gene_metrics — for umi_median per dataset (using batch=all rows)
    f = run_dir / "combined_gene_metrics.tsv"
    gn_full = pd.read_csv(f, sep="\t") if f.exists() else pd.DataFrame()
    gn = gn_full[gn_full["batch"] == "all"].copy() if not gn_full.empty else gn_full

    for ds in DATASET_ORDER:
        if ds not in msum["dataset"].values:
            continue
        row = {"dataset": ds, "config_label": label, "config": run}

        # n_cells / umi / mito / cell-level
        ms_row = msum[msum["dataset"] == ds].iloc[0]
        row["n_cells"] = ms_row["n_cells"]
        row["umi_median_per_cell"] = ms_row["umi_median"]
        row["mito_median"] = ms_row["mito_median"]
        row["frac_cells_with_guide"] = ms_row["frac_cells_with_guide"]
        row["guides_per_cell_median"] = ms_row["guides_per_cell_median"]

        # vs expected
        exp = expected.get(ds, np.nan)
        row["expected_cells"] = exp
        row["ratio_obs_to_expected"] = ms_row["n_cells"] / exp if exp else np.nan

        # intended-target metrics
        if not itm.empty and ds in itm["dataset"].values:
            itm_row = itm[itm["dataset"] == ds].iloc[0]
            row["auroc"] = itm_row.get("auroc", np.nan)
            row["auprc"] = itm_row.get("auprc", np.nan)
            row["median_log2fc_intended_eval"] = itm_row.get("median_log2fc", np.nan)

        # guide metrics — frac exactly 1 guide / guides per cell mean
        if not gm.empty and ds in gm["dataset"].values:
            gm_row = gm[gm["dataset"] == ds].iloc[0]
            row["guides_per_cell_mean"] = gm_row.get("guides_per_cell_mean", np.nan)
            row["n_cells_with_guide"] = gm_row.get("n_cells_with_guide", np.nan)
            row["n_cells_exactly_1_guide"] = gm_row.get("n_cells_exactly_1_guide", np.nan)
            n_with = gm_row.get("n_cells_with_guide", np.nan)
            n_1 = gm_row.get("n_cells_exactly_1_guide", np.nan)
            row["frac_cells_exactly_1_guide_of_assigned"] = n_1 / n_with if (n_with and n_with > 0) else np.nan
            row["frac_cells_exactly_1_guide_of_all"] = n_1 / row["n_cells"] if row["n_cells"] else np.nan

        # gene metrics — umi mean as well
        if not gn.empty and ds in gn["dataset"].values:
            gn_row = gn[gn["dataset"] == ds].iloc[0]
            row["umi_mean_per_cell"] = gn_row.get("umi_mean", np.nan)

        # per-class median log2fc and median fc
        if not itr.empty and ds in itr["dataset"].values:
            ds_results = itr[itr["dataset"] == ds].copy()
            classed = per_dataset_class_stats(ds_results)
            for cls in ["non_targeting", "negative_control_OR", "positive_control", "tf_targeting", "all_targeting"]:
                sub = classed[classed["guide_class"] == cls]
                row[f"median_log2fc_{cls}"] = sub["log2_fc"].median()
                row[f"median_fc_{cls}"] = sub["fc"].median()
                row[f"n_guides_{cls}"] = len(sub)

        # Per-guide trans summary — gives us the genome-wide median log2fc per guide,
        # which is the right way to summarize FC for non-targeting / negative controls
        # (no intended target, but every gene was tested as trans null)
        trans_pg = load_trans_per_guide(ds, run)
        if trans_pg is not None:
            for cls in ["non_targeting", "negative_control_OR", "positive_control", "tf_targeting"]:
                sub = trans_pg[trans_pg["guide_class"] == cls]
                # Median (across guides in this class) of the per-guide median log2fc
                row[f"trans_median_log2fc_{cls}"] = sub["median_log2fc"].median() if len(sub) else np.nan
                row[f"trans_median_abs_log2fc_{cls}"] = sub["median_log2fc"].abs().median() if len(sub) else np.nan
                row[f"trans_n_guides_{cls}"] = len(sub)
                row[f"trans_frac_significant_{cls}"] = sub["frac_significant"].median() if len(sub) else np.nan
            # All targeting = positive control + tf_targeting (OR and NT excluded)
            sub_all = trans_pg[trans_pg["guide_class"].isin(["positive_control", "tf_targeting"])]
            row["trans_median_log2fc_all_targeting"] = sub_all["median_log2fc"].median() if len(sub_all) else np.nan
            row["trans_n_guides_all_targeting"] = len(sub_all)

        # Per-cell guide-count breakdown from h5mu (if available)
        if not breakdown.empty:
            b_row = breakdown[(breakdown["dataset"] == ds) & (breakdown["config"] == run)]
            if len(b_row) > 0:
                br = b_row.iloc[0]
                for k in ["frac_cells_0_guides", "frac_cells_exactly_1_guide",
                          "frac_cells_exactly_2_guides", "frac_cells_exactly_3_guides",
                          "frac_cells_more_than_3_guides",
                          "n_cells_exactly_2_guides", "n_cells_exactly_3_guides",
                          "n_cells_more_than_3_guides"]:
                    row[k] = br.get(k, np.nan)

        rows.append(row)

wide = pd.DataFrame(rows)
wide.to_csv(OUT / "summary_wide.tsv", sep="\t", index=False)
print(f"Saved: summary_wide.tsv ({len(wide)} rows)")

# Best AUROC per dataset
best = (wide.dropna(subset=["auroc"])
            .sort_values(["dataset", "auroc"], ascending=[True, False])
            .groupby("dataset", as_index=False).head(1))
best = best.set_index("dataset")
best = best.loc[[d for d in DATASET_ORDER if d in best.index]].reset_index()
best.to_csv(OUT / "summary_best.tsv", sep="\t", index=False)
print(f"Saved: summary_best.tsv ({len(best)} rows)")
print()
print("=== Best AUROC config per dataset ===")
print(best[["dataset", "config_label", "auroc", "auprc", "n_cells", "ratio_obs_to_expected",
            "frac_cells_with_guide", "guides_per_cell_mean",
            "median_log2fc_positive_control", "median_log2fc_tf_targeting"]].to_string(index=False))


# Multi-panel bar chart — selected metrics at best AUROC config per dataset
PANEL_METRICS = [
    ("auroc", "auROC", None, 0.5),
    ("auprc", "auPRC", None, None),
    ("n_cells", "Total cells", "k", None),
    ("ratio_obs_to_expected", "Observed / expected cells", None, 1.0),
    ("umi_median_per_cell", "Median UMI / cell", None, None),
    ("frac_cells_with_guide", "Frac cells with guide", None, None),
    ("guides_per_cell_mean", "Mean guides / cell", None, None),
    ("frac_cells_exactly_1_guide_of_assigned", "Frac w/ exactly 1 guide\n(of cells with any guide)", None, None),
    ("median_log2fc_positive_control", "Median log2fc\npositive controls", None, 0.0),
    ("median_log2fc_tf_targeting", "Median log2fc\nTF-targeting", None, 0.0),
    ("median_log2fc_all_targeting", "Median log2fc\nall targeting", None, 0.0),
    ("median_fc_tf_targeting", "Median observed FC\nTF-targeting", None, 1.0),
    ("frac_cells_exactly_1_guide", "Frac cells\nexactly 1 guide", None, None),
    ("frac_cells_exactly_2_guides", "Frac cells\nexactly 2 guides", None, None),
    ("frac_cells_exactly_3_guides", "Frac cells\nexactly 3 guides", None, None),
    ("frac_cells_more_than_3_guides", "Frac cells\n>3 guides", None, None),
    ("trans_median_log2fc_non_targeting", "Trans median log2fc\nnon-targeting (null)", None, 0.0),
    ("trans_median_log2fc_negative_control_OR", "Trans median log2fc\nOR neg controls", None, 0.0),
    ("trans_median_log2fc_positive_control", "Trans median log2fc\npositive controls", None, 0.0),
    ("trans_median_log2fc_tf_targeting", "Trans median log2fc\nTF-targeting", None, 0.0),
]

ds_in = [d for d in DATASET_ORDER if d in best["dataset"].values]
n_metrics = len(PANEL_METRICS)
ncols = 3
nrows = int(np.ceil(n_metrics / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows))
axes = np.atleast_2d(axes).flatten()

for ax, (col, title, fmt, hline) in zip(axes, PANEL_METRICS):
    if col not in best.columns:
        ax.set_visible(False)
        continue
    x = np.arange(len(ds_in))
    vals = []
    for ds in ds_in:
        v = best[best["dataset"] == ds][col]
        vals.append(v.iloc[0] if len(v) > 0 and pd.notna(v.iloc[0]) else np.nan)
    colors = [DATASET_COLORS[d] for d in ds_in]
    bars = ax.bar(x, vals, color=colors, alpha=0.9, edgecolor="white")
    if hline is not None:
        ax.axhline(hline, color="black", linestyle="--", linewidth=1, alpha=0.5)
    # value labels
    for xi, v in zip(x, vals):
        if pd.notna(v):
            txt = f"{v:.0f}" if abs(v) > 100 else f"{v:.3f}" if abs(v) < 1 else f"{v:.2f}"
            ax.text(xi, v if v >= 0 else v - 0.02, txt, ha="center",
                    va="bottom" if v >= 0 else "top", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT.get(d, d) for d in ds_in], fontsize=8, rotation=20, ha="right")
    ax.set_title(title, fontsize=10)
    if fmt == "k":
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x/1000:.0f}K"))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

for ax in axes[n_metrics:]:
    ax.set_visible(False)

# Title with best-config legend
config_labels = ", ".join(f"{SHORT.get(d, d).replace(chr(10), ' ')}={best[best['dataset']==d]['config_label'].iloc[0]}" for d in ds_in)
fig.suptitle(f"Per-dataset metrics at best-AUROC filter config\n({config_labels})",
             fontsize=12, fontweight="bold", y=1.00)
plt.tight_layout()
plt.savefig(OUT / "per_dataset_best.pdf", dpi=300, bbox_inches="tight")
plt.savefig(OUT / "per_dataset_best.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"\nSaved: per_dataset_best.pdf / .png")

# Same metrics but across all configs (per dataset, x = config) for the same set
fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows))
axes = np.atleast_2d(axes).flatten()
config_labels_all = [lbl for _, lbl in RUNS]

for ax, (col, title, fmt, hline) in zip(axes, PANEL_METRICS):
    if col not in wide.columns:
        ax.set_visible(False)
        continue
    for ds in ds_in:
        sub = wide[wide["dataset"] == ds].set_index("config_label").reindex(config_labels_all)
        ax.plot(config_labels_all, sub[col].values, "o-", color=DATASET_COLORS[ds],
                label=SHORT.get(ds, ds).replace("\n", " "), linewidth=1.8, markersize=6)
    if hline is not None:
        ax.axhline(hline, color="black", linestyle="--", linewidth=1, alpha=0.4)
    ax.set_title(title, fontsize=10)
    if fmt == "k":
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x/1000:.0f}K"))
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

axes[0].legend(fontsize=7, loc="best", framealpha=0.85)
for ax in axes[n_metrics:]:
    ax.set_visible(False)

fig.suptitle("Per-dataset metrics across filter configs", fontsize=12, fontweight="bold", y=1.00)
plt.tight_layout()
plt.savefig(OUT / "per_dataset_across_configs.pdf", dpi=300, bbox_inches="tight")
plt.savefig(OUT / "per_dataset_across_configs.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved: per_dataset_across_configs.pdf / .png")

print(f"\nDone. Outputs in: {OUT}")
