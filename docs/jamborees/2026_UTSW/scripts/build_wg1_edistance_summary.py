"""Build the WG1-B cross-dataset energy-distance summary table.

Reshapes `reference/cross_dataset_edistance_summary.tsv` into a slide-deck-friendly
wide-format table with: identity columns, run scope, effect-size diagnostics, the
calibration-robust significance count (recommended), the pval-based count
(⚠ caveat), and a calibration_state classification.

Production rows that don't yet have an ED run on Synapse get a placeholder row
with `data_state` = "pending" / "blocked" — symmetric with WG1-A's qc_summary.

Re-run whenever a new dataset's ED bundle lands.

Usage:
    python scripts/build_wg1_edistance_summary.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

JAMB = Path(__file__).resolve().parent.parent / "docs/jamborees/2026_UTSW"
ED_SUMMARY = JAMB / "reference/cross_dataset_edistance_summary.tsv"
EXP_METADATA = JAMB / "reference/experimental_metadata_simplified.tsv"
OUTPUT = JAMB / "working_groups/wg1_data_qc/edistance_summary.tsv"

# Production roster — anything not in cross_dataset_edistance_summary.tsv gets a
# placeholder row. State comes from experimental_metadata.pipeline_status; values
# below override per-dataset where ED has its own blocker.
PRODUCTION_DATASETS = [
    "Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq",
    "Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq",
    "Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq",
    "Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq",
    "Engreitz_WTC11-endothelial-cells_TF-Perturb-seq",
]

ED_DATA_STATE = {
    "Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq":
        "pending — ED run retry was cancelled 2026-05-10; awaiting decision on rerun strategy",
    "Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq":
        "blocked — awaiting Gersbach team (Sara) deliverables",
    "Engreitz_WTC11-endothelial-cells_TF-Perturb-seq":
        "blocked — no portal data",
}


def classify_calibration(row: pd.Series) -> str:
    """Anti-conservative if all NCs have pval_mean == 0 (or NaN if no NCs)."""
    n_nc = row.get("n_targets_negative_control")
    n_nc_pval_eq_0 = row.get("n_nc_pval_eq_0")
    if pd.isna(n_nc) or n_nc == 0:
        return "no_NCs_in_run"
    if n_nc_pval_eq_0 == n_nc:
        return "anti-conservative"
    return "healthy"


def main() -> None:
    ed = pd.read_csv(ED_SUMMARY, sep="\t")
    meta = pd.read_csv(EXP_METADATA, sep="\t").rename(columns={"differentiation": "lineage"})

    # Drop HTv2 testbed if present (user excluded testbed from WG outputs)
    ed = ed[~ed["dataset_id"].str.contains("benchmark", case=False, na=False)]

    rows: list[dict] = []
    for ds in PRODUCTION_DATASETS:
        ed_row = ed[ed["dataset_id"] == ds]
        meta_row = meta[meta["dataset_id"] == ds].iloc[0] if (meta["dataset_id"] == ds).any() else {}
        if not ed_row.empty:
            r = ed_row.iloc[0]
            row = {
                "dataset_id": ds,
                "dataset_name": meta_row.get("dataset_name", ""),
                "lineage": meta_row.get("lineage", ""),
                "data_state": "complete",
                "n_targets_total": int(r["n_targets_total"]),
                "n_targeting": int(r["n_targets_targeting"]),
                "n_NCs": int(r["n_targets_negative_control"]),
                "n_pos_ctrl": int(r["n_targets_positive_control"]),
                "distance_mean_median_targeting": float(r["distance_mean_median_targeting"]),
                "distance_mean_median_NC": float(r["distance_mean_median_negative_control"]),
                "NC_distance_max": float(r["nc_distance_mean_max"]),
                "n_sig_distance_gt_NC_max": int(r["n_targeting_above_nc_max"]),
                "frac_sig_distance_gt_NC_max":
                    int(r["n_targeting_above_nc_max"]) / int(r["n_targets_targeting"]),
                "n_sig_pval_lt_0p05": int(r["n_pval_lt_0p05"]),
                "frac_sig_pval_lt_0p05": int(r["n_pval_lt_0p05"]) / int(r["n_targets_total"]),
                "n_NCs_pval_eq_0": int(r["n_nc_pval_eq_0"]),
                "calibration_state": classify_calibration(r),
                "preferred_criterion":
                    "distance > NC max"
                    if classify_calibration(r) == "anti-conservative"
                    else "pval_mean < 0.05",
            }
        else:
            row = {
                "dataset_id": ds,
                "dataset_name": meta_row.get("dataset_name", ""),
                "lineage": meta_row.get("lineage", ""),
                "data_state": ED_DATA_STATE.get(ds, "pending"),
                "n_targets_total": None,
                "n_targeting": None,
                "n_NCs": None,
                "n_pos_ctrl": None,
                "distance_mean_median_targeting": None,
                "distance_mean_median_NC": None,
                "NC_distance_max": None,
                "n_sig_distance_gt_NC_max": None,
                "frac_sig_distance_gt_NC_max": None,
                "n_sig_pval_lt_0p05": None,
                "frac_sig_pval_lt_0p05": None,
                "n_NCs_pval_eq_0": None,
                "calibration_state": None,
                "preferred_criterion": None,
            }
        rows.append(row)

    out = pd.DataFrame(rows)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT, sep="\t", index=False)
    print(f"wrote {OUTPUT} ({len(out)} rows × {len(out.columns)} cols)")


if __name__ == "__main__":
    main()
