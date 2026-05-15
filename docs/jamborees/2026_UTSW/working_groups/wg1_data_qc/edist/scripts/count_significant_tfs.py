"""Count significant TFs in one dataset's energy-distance results.

Input:  wg1_significant_tfs.tsv (per-dataset table at data/<id>/energy_distance/)
Output: 1-row TSV with significance counts under both criteria.

Usage:
    uv run python count_significant_tfs.py \\
        --input  ../../../data/<id>/energy_distance/wg1_significant_tfs.tsv \\
        --output ../results/per_dataset/<id>.tsv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    # Upstream metadata join can emit duplicate target_id rows (one per gene-symbol alias).
    df = df.drop_duplicates(subset="target_id", keep="first")
    dataset_id = args.input.resolve().parents[1].name

    targeting = df[df["type"] == "targeting"]
    nc = df[df["type"] == "negative control"]

    n_nc_sig_pval = int((nc["pval_mean"] < 0.05).sum()) if len(nc) else 0
    frac_nc_sig_pval = n_nc_sig_pval / len(nc) if len(nc) else float("nan")

    row = {
        "dataset_id": dataset_id,
        "n_targeting": len(targeting),
        "n_negative_control": len(nc),
        "n_sig_distance_gt_NC_max": int(targeting["sig_distance_gt_NC_max"].sum()),
        "n_sig_pval_lt_0p05": int(targeting["sig_pval_lt_0p05"].sum()),
        "nc_distance_max": nc["distance_mean"].max() if len(nc) else None,
        "distance_mean_median_targeting": targeting["distance_mean"].median(),
        "distance_mean_median_NC": nc["distance_mean"].median() if len(nc) else None,
        "n_NCs_pval_lt_0p05": n_nc_sig_pval,
        "frac_NCs_pval_lt_0p05": frac_nc_sig_pval,
        # NC false-positive rate >2x the 0.05 nominal level => p-values are anti-conservative.
        "calibration_state": "anti-conservative" if frac_nc_sig_pval > 0.10 else "healthy",
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(args.output, sep="\t", index=False)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
