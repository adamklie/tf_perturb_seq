"""Per-dataset small multiples: auROC & auPRC across filtering parameters.

Goal: ID the best parameter setting per dataset.
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

BASE = Path("/cellar/users/aklie/projects/tf_perturb_seq/datasets/technology-benchmark_WTC11_TF-Perturb-seq")
RESULTS = BASE / "results" / "cross_tech_comparison"
OUT = RESULTS / "cross_parameter"

RUNS = [
    ("cleanser_extremes_200_mito_15pc", "200 UMI"),
    ("cleanser_500_mito_15pc", "500 UMI"),
    ("cleanser_800_mito_15pc", "800 UMI"),
    ("cleanser_extremes_2000_mito_15pc", "2000 UMI"),
    ("cleanser_knee2_mito_15pc", "Knee2"),
]
SHORT = {
    'Hon_WTC11-benchmark_TF-Perturb-seq': 'Hon',
    'Huangfu_WTC11-benchmark_TF-Perturb-seq': 'Huangfu',
    'Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3': 'Gersbach GEM-Xv3',
    'Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2': 'Gersbach HTv2',
    'Engreitz_WTC11-benchmark_TF-Perturb-seq': 'Engreitz',
}
DATASET_ORDER = list(SHORT.keys())

rows = []
for run, label in RUNS:
    f = RESULTS / run / "combined_intended_target_metrics.tsv"
    if not f.exists():
        continue
    df = pd.read_csv(f, sep="\t")
    df["run"] = run
    df["run_label"] = label
    rows.append(df)
m = pd.concat(rows, ignore_index=True)

# Also pull n_cells for context
summary_rows = []
for run, label in RUNS:
    f = RESULTS / run / "metrics_summary.tsv"
    if f.exists():
        d = pd.read_csv(f, sep="\t")
        d["run"], d["run_label"] = run, label
        summary_rows.append(d)
s = pd.concat(summary_rows, ignore_index=True)

run_labels = [lbl for _, lbl in RUNS]
datasets_present = [d for d in DATASET_ORDER if d in m["dataset"].values]

fig, axes = plt.subplots(1, len(datasets_present),
                         figsize=(3.2 * len(datasets_present), 4.2), sharey=True)
if len(datasets_present) == 1:
    axes = [axes]

for ax, ds in zip(axes, datasets_present):
    sub = m[m["dataset"] == ds].set_index("run_label").reindex(run_labels)
    x = np.arange(len(run_labels))
    ax.plot(x, sub["auroc"], marker="o", color="#1f77b4", label="auROC", linewidth=2)
    ax.plot(x, sub["auprc"], marker="s", color="#d62728", label="auPRC", linewidth=2)

    # Best parameter per metric
    for metric, color, y_off in [("auroc", "#1f77b4", 0.04), ("auprc", "#d62728", -0.06)]:
        vals = sub[metric]
        if vals.notna().any():
            best_idx = vals.idxmax()
            best_x = run_labels.index(best_idx)
            best_y = vals.loc[best_idx]
            ax.scatter([best_x], [best_y], s=200, facecolors="none",
                       edgecolors=color, linewidths=2, zorder=5)
            ax.annotate(f"★ {best_idx}", (best_x, best_y),
                        xytext=(0, 12 if metric == "auroc" else -16),
                        textcoords="offset points", ha="center",
                        fontsize=8, color=color, fontweight="bold")

    # n_cells on secondary axis
    s_sub = s[s["dataset"] == ds].set_index("run_label").reindex(run_labels)
    ax2 = ax.twinx()
    ax2.bar(x, s_sub["n_cells"] / 1000, alpha=0.15, color="gray", zorder=0, width=0.6)
    ax2.set_ylabel("n_cells (K)", fontsize=8, color="gray")
    ax2.tick_params(axis="y", labelsize=7, colors="gray")
    ax2.spines["top"].set_visible(False)

    ax.set_xticks(x)
    ax.set_xticklabels(run_labels, rotation=30, ha="right", fontsize=9)
    ax.set_title(SHORT[ds], fontsize=11, fontweight="bold")
    ax.set_ylim(0.0, 1.05)
    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.4, linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.set_zorder(ax2.get_zorder() + 1)
    ax.patch.set_visible(False)

axes[0].set_ylabel("auROC / auPRC")
axes[0].legend(loc="lower left", fontsize=9, frameon=False)
fig.suptitle("Best filtering parameter per dataset (★ = max)",
             fontsize=13, fontweight="bold", y=1.02)
plt.tight_layout()

out_pdf = OUT / "best_param_per_dataset.pdf"
out_png = OUT / "best_param_per_dataset.png"
plt.savefig(out_pdf, bbox_inches="tight")
plt.savefig(out_png, dpi=200, bbox_inches="tight")
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")

# Tabular best-param summary
best_tbl = []
for ds in datasets_present:
    sub = m[m["dataset"] == ds].set_index("run_label").reindex(run_labels)
    row = {"dataset": SHORT[ds]}
    for metric in ["auroc", "auprc"]:
        vals = sub[metric]
        if vals.notna().any():
            row[f"best_{metric}_param"] = vals.idxmax()
            row[f"best_{metric}"] = round(vals.max(), 4)
            row[f"worst_{metric}"] = round(vals.min(), 4)
            row[f"range_{metric}"] = round(vals.max() - vals.min(), 4)
    best_tbl.append(row)
best_df = pd.DataFrame(best_tbl)
best_df.to_csv(OUT / "best_param_per_dataset.tsv", sep="\t", index=False)
print(best_df.to_string(index=False))
