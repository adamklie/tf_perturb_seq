"""Compute exactly-1 / exactly-2 / exactly-3 / >3 guide assignment breakdown.

For each (dataset, filter config) with a local h5mu, count cells by their
number of assigned guides using `mod['guide'].layers['guide_assignment']`.
Writes a single long TSV that joins onto `per_dataset_summary/summary_wide.tsv`.

Output:
  results/cross_tech_comparison/per_dataset_summary/guide_assignment_breakdown.tsv
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import mudata as md
import warnings
warnings.filterwarnings("ignore")

PROJECT_ROOT = Path("/cellar/users/aklie/projects/tf_perturb_seq")
BASE = PROJECT_ROOT / "datasets" / "technology-benchmark_WTC11_TF-Perturb-seq"
RESULTS = BASE / "results" / "cross_tech_comparison"
OUT = RESULTS / "per_dataset_summary"
OUT.mkdir(parents=True, exist_ok=True)

DATASETS = [
    "Hon_WTC11-benchmark_TF-Perturb-seq",
    "Huangfu_WTC11-benchmark_TF-Perturb-seq",
    "Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3",
    "Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2",
    "Engreitz_WTC11-benchmark_TF-Perturb-seq",
]
CONFIGS = [
    ("cleanser_extremes_200_mito_15pc", "200 UMI"),
    ("cleanser_500_mito_15pc", "500 UMI"),
    ("cleanser_800_mito_15pc", "800 UMI"),
    ("cleanser_extremes_2000_mito_15pc", "2000 UMI"),
    ("cleanser_knee2_mito_15pc", "Knee2"),
]


def count_assignments(h5mu_path):
    mu = md.read_h5mu(str(h5mu_path), backed="r")
    g = mu.mod["guide"]
    if "guide_assignment" not in g.layers:
        return None
    layer = g.layers["guide_assignment"]
    # Layer is (n_cells, n_guides) sparse. Sum across guides per cell.
    if sp.issparse(layer):
        per_cell = np.asarray(layer.sum(axis=1)).flatten()
    else:
        per_cell = np.asarray(layer).sum(axis=1)
    per_cell = per_cell.astype(int)
    n_cells = per_cell.shape[0]
    return {
        "n_cells_total": n_cells,
        "n_cells_0_guides": int((per_cell == 0).sum()),
        "n_cells_exactly_1_guide": int((per_cell == 1).sum()),
        "n_cells_exactly_2_guides": int((per_cell == 2).sum()),
        "n_cells_exactly_3_guides": int((per_cell == 3).sum()),
        "n_cells_more_than_3_guides": int((per_cell > 3).sum()),
        "mean_guides_per_cell": float(per_cell.mean()),
        "median_guides_per_cell": float(np.median(per_cell)),
    }


rows = []
for ds in DATASETS:
    for cfg, label in CONFIGS:
        h5mu = BASE.parent / ds / "runs" / cfg / "pipeline_dashboard" / "inference_mudata.h5mu"
        if not h5mu.exists():
            print(f"  skip {ds} / {cfg} — no h5mu")
            continue
        print(f"  loading {ds} / {label} ...")
        try:
            stats = count_assignments(h5mu)
        except Exception as e:
            print(f"    FAILED: {e}")
            continue
        if stats is None:
            print(f"    no guide_assignment layer")
            continue
        row = {"dataset": ds, "config": cfg, "config_label": label, **stats}
        # Fractions
        n = row["n_cells_total"]
        for k in ["0_guides", "exactly_1_guide", "exactly_2_guides",
                  "exactly_3_guides", "more_than_3_guides"]:
            row[f"frac_cells_{k}"] = row[f"n_cells_{k}"] / n if n > 0 else np.nan
        rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(OUT / "guide_assignment_breakdown.tsv", sep="\t", index=False)
print(f"\nSaved: {OUT / 'guide_assignment_breakdown.tsv'} ({len(df)} rows)")
print()
print(df[["dataset", "config_label", "n_cells_total",
         "n_cells_exactly_1_guide", "n_cells_exactly_2_guides",
         "n_cells_exactly_3_guides", "n_cells_more_than_3_guides"]].to_string(index=False))
