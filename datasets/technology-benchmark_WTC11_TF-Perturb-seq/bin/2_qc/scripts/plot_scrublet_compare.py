"""Scrublet on/off comparison — doublet removal axis.

Verified from `pipeline_info/params_*.json` on GCS (NOT folder names — those
are misleading; the `Benchmark_scrublet_cleanser_*` folder actually uses sceptre):

  cleanser_800_mito_15pc                → ENABLE_SCRUBLET=False, method=cleanser
  scrublet_off_cleanser_800_mito_15pc   → ENABLE_SCRUBLET=False, method=sceptre
  scrublet_on_sceptre_800_mito_15pc     → ENABLE_SCRUBLET=True,  method=sceptre

Therefore:
  - cleanser_800 vs scrublet_off (gray vs green): cleanser vs sceptre method
    comparison (no scrublet either, only GUIDE_ASSIGNMENT_method differs)
  - scrublet_off vs scrublet_on (green vs orange): clean scrublet on/off
    comparison — both use sceptre, only ENABLE_SCRUBLET differs

Outputs to: results/cross_tech_comparison/scrublet_compare/
"""
from pathlib import Path
import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

PROJECT_ROOT = Path("/cellar/users/aklie/projects/tf_perturb_seq")
BASE = PROJECT_ROOT / "datasets" / "technology-benchmark_WTC11_TF-Perturb-seq"
RESULTS = BASE / "results" / "cross_tech_comparison"
OUT = RESULTS / "scrublet_compare"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.append(str(PROJECT_ROOT / "config"))
from loader import load_colors
cfg = load_colors("technology-benchmark_WTC11_TF-Perturb-seq")
DATASET_COLORS = cfg["dataset_colors"]
DATASET_ORDER = cfg["dataset_order"]

# Three runs (labels reflect actual GUIDE_ASSIGNMENT_method from params.json, not folder names)
RUNS = [
    ("cleanser_800_mito_15pc",                "cleanser, scrublet OFF",      "#bdbdbd"),  # gray = cleanser reference
    ("scrublet_off_cleanser_800_mito_15pc",   "sceptre, scrublet OFF",       "#1b9e77"),
    ("scrublet_on_sceptre_800_mito_15pc",     "sceptre, scrublet ON",        "#d95f02"),
]
RUN_LABELS = [r[1] for r in RUNS]

POSITIVE_CONTROL_GENES = {"CD151", "CD55", "CD81", "NGFRAP1", "TFRC"}

SHORT = {
    "Hon_WTC11-benchmark_TF-Perturb-seq": "Hon",
    "Huangfu_WTC11-benchmark_TF-Perturb-seq": "Huangfu",
    "Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3": "Gersbach\nGEM-Xv3",
    "Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2": "Gersbach\nHTv2",
    "Engreitz_WTC11-benchmark_TF-Perturb-seq": "Engreitz",
}

CONFOUND_NOTE = ("Verified from pipeline_info/params.json: gray = cleanser, "
                 "green/orange = sceptre. Green→orange (scrublet OFF→ON) is a "
                 "clean scrublet test (both sceptre). Gray→green is cleanser vs sceptre.")


def classify_guides(gene_name, targeting):
    is_target = targeting == True
    is_pc = gene_name.isin(POSITIVE_CONTROL_GENES) & is_target
    is_or = gene_name.str.match(r"^OR[0-9]", na=False) & is_target
    return np.where(
        ~is_target, "non_targeting",
        np.where(is_pc, "positive_control",
                 np.where(is_or, "negative_control_OR", "tf_targeting")))


def load_per_run(run, label):
    """Load and join all per-(dataset) metrics for one run."""
    rows = []
    msum = pd.read_csv(RESULTS / run / "metrics_summary.tsv", sep="\t")
    itm_path = RESULTS / run / "combined_intended_target_metrics.tsv"
    itm = pd.read_csv(itm_path, sep="\t") if itm_path.exists() else pd.DataFrame()
    gm_path = RESULTS / run / "combined_guide_metrics.tsv"
    gm_full = pd.read_csv(gm_path, sep="\t") if gm_path.exists() else pd.DataFrame()
    gm = gm_full[gm_full["batch"] == "all"] if not gm_full.empty else gm_full

    for ds in DATASET_ORDER:
        if ds not in msum["dataset"].values:
            continue
        ms_row = msum[msum["dataset"] == ds].iloc[0]
        row = {"dataset": ds, "config": run, "config_label": label}
        row["n_cells"] = ms_row["n_cells"]
        row["umi_median"] = ms_row["umi_median"]
        row["mito_median"] = ms_row["mito_median"]
        row["frac_cells_with_guide"] = ms_row["frac_cells_with_guide"]
        row["guide_umi_median"] = ms_row.get("guide_umi_median")
        row["guides_per_cell_median"] = ms_row.get("guides_per_cell_median")

        if not itm.empty and ds in itm["dataset"].values:
            r = itm[itm["dataset"] == ds].iloc[0]
            row["auroc"] = r.get("auroc")
            row["auprc"] = r.get("auprc")
            row["median_log2fc_intended_eval"] = r.get("median_log2fc")

        if not gm.empty and ds in gm["dataset"].values:
            r = gm[gm["dataset"] == ds].iloc[0]
            row["guides_per_cell_mean"] = r.get("guides_per_cell_mean")
            row["n_cells_with_guide"] = r.get("n_cells_with_guide")
            row["n_cells_exactly_1_guide"] = r.get("n_cells_exactly_1_guide")

        # Per-class trans median log2fc
        trans_path = (BASE.parent / ds / "runs" / run /
                      "pipeline_dashboard" / "additional_qc" / "trans" /
                      "trans_per_guide_summary.tsv")
        if trans_path.exists():
            tdf = pd.read_csv(trans_path, sep="\t")
            tdf["guide_class"] = classify_guides(tdf["gene_name"], tdf["targeting"])
            for cls in ["non_targeting", "negative_control_OR", "positive_control", "tf_targeting"]:
                sub = tdf[tdf["guide_class"] == cls]
                row[f"trans_median_log2fc_{cls}"] = sub["median_log2fc"].median() if len(sub) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


frames = [load_per_run(run, label) for run, label, _ in RUNS]
wide = pd.concat(frames, ignore_index=True)
wide.to_csv(OUT / "scrublet_summary_wide.tsv", sep="\t", index=False)
print(f"Saved: scrublet_summary_wide.tsv ({len(wide)} rows)")

# === Cells filtered analysis ===
ds_in = [d for d in DATASET_ORDER if d in wide["dataset"].values]
piv_cells = wide.pivot(index="dataset", columns="config_label", values="n_cells")
ref_label = RUN_LABELS[0]
off_label = RUN_LABELS[1]
on_label = RUN_LABELS[2]
piv_cells["delta_off_vs_ref"] = piv_cells[off_label] - piv_cells[ref_label]
piv_cells["delta_on_vs_off"] = piv_cells[off_label] - piv_cells[on_label]
piv_cells["pct_dropped_on_vs_off"] = piv_cells["delta_on_vs_off"] / piv_cells[off_label] * 100
piv_cells = piv_cells.reindex([d for d in DATASET_ORDER if d in piv_cells.index])
piv_cells.to_csv(OUT / "cells_filtered.tsv", sep="\t")
print()
print("=== Cells filtered ===")
print(piv_cells.round(0).to_string())


def fmt_k(x, _):
    return f"{x/1000:.0f}K"


# === Plot 1: n_cells across the 3 runs (with reference) ===
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
x = np.arange(len(ds_in))
width = 0.27
for j, (run, label, color) in enumerate(RUNS):
    sub = wide[wide["config_label"] == label].set_index("dataset").reindex(ds_in)
    offset = (j - 1) * width
    ax1.bar(x + offset, sub["n_cells"].values, width, color=color, alpha=0.85,
            edgecolor="white", label=label)
ax1.set_xticks(x)
ax1.set_xticklabels([SHORT.get(d, d) for d in ds_in], fontsize=9)
ax1.set_title("Total cells per run", fontsize=11)
ax1.yaxis.set_major_formatter(FuncFormatter(fmt_k))
ax1.legend(fontsize=8, loc="upper right")
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# pct dropped (scrublet_on vs scrublet_off, dataset-colored)
pct_vals = piv_cells.loc[ds_in, "pct_dropped_on_vs_off"].values
colors = [DATASET_COLORS[d] for d in ds_in]
bars = ax2.bar(np.arange(len(ds_in)), pct_vals, color=colors, alpha=0.9, edgecolor="white")
for xi, v in zip(np.arange(len(ds_in)), pct_vals):
    if pd.notna(v):
        ax2.text(xi, v, f"{v:.1f}%", ha="center", va="bottom" if v >= 0 else "top", fontsize=9)
ax2.set_xticks(np.arange(len(ds_in)))
ax2.set_xticklabels([SHORT.get(d, d) for d in ds_in], fontsize=9)
ax2.set_ylabel("% cells dropped (scrublet ON vs OFF)")
ax2.set_title("Cells removed by scrublet", fontsize=11)
ax2.axhline(0, color="black", linewidth=0.5)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

fig.text(0.5, 0.01, CONFOUND_NOTE, ha="center", fontsize=8, color="darkred")
plt.tight_layout(rect=[0, 0.04, 1, 1])
plt.savefig(OUT / "scrublet_cells_filtered.pdf", dpi=300, bbox_inches="tight")
plt.savefig(OUT / "scrublet_cells_filtered.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"\nSaved: scrublet_cells_filtered.pdf / .png")


# === Plot 2: full metric grid across the 3 runs ===
PANEL_METRICS = [
    ("n_cells", "Total cells", "k", None),
    ("umi_median", "Median UMI / cell", None, None),
    ("guide_umi_median", "Median guide UMI", None, None),
    ("mito_median", "Median % mito", None, None),
    ("frac_cells_with_guide", "Frac cells with guide", None, None),
    ("guides_per_cell_mean", "Mean guides / cell", None, None),
    ("auroc", "auROC", None, 0.5),
    ("auprc", "auPRC", None, None),
    ("median_log2fc_intended_eval", "Median log2fc\n(intended-target eval)", None, 0.0),
    ("trans_median_log2fc_non_targeting", "Trans median log2fc\nnon-targeting (null)", None, 0.0),
    ("trans_median_log2fc_negative_control_OR", "Trans median log2fc\nOR neg controls", None, 0.0),
    ("trans_median_log2fc_positive_control", "Trans median log2fc\npositive controls", None, 0.0),
]

ncols = 4
nrows = int(np.ceil(len(PANEL_METRICS) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.5 * nrows))
axes = np.atleast_2d(axes).flatten()

for ax, (col, title, fmt, hline) in zip(axes, PANEL_METRICS):
    if col not in wide.columns:
        ax.set_visible(False)
        continue
    x = np.arange(len(ds_in))
    width_inner = 0.27
    for j, (run, label, color) in enumerate(RUNS):
        sub = wide[wide["config_label"] == label].set_index("dataset").reindex(ds_in)
        vals = sub[col].values
        offset = (j - 1) * width_inner
        ax.bar(x + offset, vals, width_inner, color=color, alpha=0.85, edgecolor="white",
               label=label if ax is axes[0] else None)
    if hline is not None:
        ax.axhline(hline, color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT.get(d, d) for d in ds_in], fontsize=8, rotation=20, ha="right")
    ax.set_title(title, fontsize=10)
    if fmt == "k":
        ax.yaxis.set_major_formatter(FuncFormatter(fmt_k))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

axes[0].legend(fontsize=7, loc="best", framealpha=0.85)
for ax in axes[len(PANEL_METRICS):]:
    ax.set_visible(False)

fig.text(0.5, 0.01, CONFOUND_NOTE, ha="center", fontsize=8, color="darkred")
fig.suptitle("Scrublet on/off — full metric grid (3 runs)", fontsize=13, fontweight="bold", y=1.00)
plt.tight_layout(rect=[0, 0.02, 1, 0.99])
plt.savefig(OUT / "scrublet_metric_grid.pdf", dpi=300, bbox_inches="tight")
plt.savefig(OUT / "scrublet_metric_grid.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved: scrublet_metric_grid.pdf / .png")


# === Plot 3: delta panel — for key metrics, show ON-OFF (within each dataset) ===
DELTA_METRICS = [
    ("n_cells", "Δ n_cells"),
    ("auroc", "Δ auROC"),
    ("auprc", "Δ auPRC"),
    ("median_log2fc_intended_eval", "Δ median log2fc (intended)"),
    ("frac_cells_with_guide", "Δ frac cells with guide"),
    ("guides_per_cell_mean", "Δ mean guides / cell"),
]
ncols = 3
nrows = 2
fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
axes = axes.flatten()

for ax, (col, title) in zip(axes, DELTA_METRICS):
    if col not in wide.columns:
        ax.set_visible(False)
        continue
    off = wide[wide["config_label"] == off_label].set_index("dataset").reindex(ds_in)[col]
    on = wide[wide["config_label"] == on_label].set_index("dataset").reindex(ds_in)[col]
    delta = (on - off).values
    colors = [DATASET_COLORS[d] for d in ds_in]
    ax.bar(np.arange(len(ds_in)), delta, color=colors, alpha=0.9, edgecolor="white")
    for xi, v in zip(np.arange(len(ds_in)), delta):
        if pd.notna(v):
            txt = f"{v:.0f}" if abs(v) > 100 else f"{v:.3f}"
            ax.text(xi, v, txt, ha="center", va="bottom" if v >= 0 else "top", fontsize=8)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_xticks(np.arange(len(ds_in)))
    ax.set_xticklabels([SHORT.get(d, d) for d in ds_in], fontsize=8, rotation=20, ha="right")
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.text(0.5, 0.01, CONFOUND_NOTE, ha="center", fontsize=8, color="darkred")
fig.suptitle("Scrublet ON minus OFF (clean: both use sceptre)", fontsize=12, fontweight="bold", y=1.00)
plt.tight_layout(rect=[0, 0.03, 1, 0.97])
plt.savefig(OUT / "scrublet_deltas.pdf", dpi=300, bbox_inches="tight")
plt.savefig(OUT / "scrublet_deltas.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved: scrublet_deltas.pdf / .png")

# Update README
readme = OUT / "README.md"
readme.write_text(f"""# Scrublet on/off comparison

**Verified from `pipeline_info/params_*.json` on GCS** (do NOT trust folder names):

| Run (GCS folder) | ENABLE_SCRUBLET | GUIDE_ASSIGNMENT_method |
|---|---|---|
| `Benchmark_cleanser_800_mito_15pc` | False | cleanser |
| `Benchmark_scrublet_cleanser_800_mito_15pc` ⚠ "cleanser" in name but actually sceptre | False | **sceptre** |
| `REAL_SCRUBLETS_RENAME_OTHER_Benchmark_scrublet_sceptre_800_mito_15pc` | True | sceptre |

Therefore:
- **Green → Orange** (`scrublet_off_cleanser_800` → `scrublet_on_sceptre_800`):
  **Clean scrublet on/off comparison** — both use sceptre, only `ENABLE_SCRUBLET`
  differs. Cell drops + metric deltas in this axis are real scrublet effects.
- **Gray → Green** (`cleanser_800` → `scrublet_off_cleanser_800`): cleanser vs
  sceptre method comparison (both no scrublet, only `GUIDE_ASSIGNMENT_method`
  differs).

## Outputs
- `scrublet_summary_wide.tsv` — per-(dataset × run) metrics (n_cells, AUROC,
  AUPRC, frac_cells_with_guide, guides_per_cell_mean, per-class trans log2fc)
- `cells_filtered.tsv` — pivot showing n_cells per run + delta + pct dropped
- `scrublet_cells_filtered.pdf` — n_cells across 3 runs + % dropped per dataset
- `scrublet_metric_grid.pdf` — 12 panels of metrics across the 3 runs
- `scrublet_deltas.pdf` — ON minus OFF for 6 key metrics (clean scrublet test)
""")
print(f"\nDone. Outputs in: {OUT}")
