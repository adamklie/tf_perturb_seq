"""Build a long cross-dataset table — one row per (target × dataset).

Input:  directory containing per-dataset `wg1_significant_tfs.tsv` files.
Output: single long TSV. Same identity + stat columns as the per-dataset
        tables, plus `dataset_id` and `dataset_short`.

Tidy form, ready for seaborn / groupby / pivot.

Usage:
    uv run python build_long_per_target_table.py \\
        --input  ../../../../data \\
        --output ../results/pairwise_distance_scatter/per_target_long.tsv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

SHORT_NAMES = {
    "Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq": "HonCM",
    "Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq": "HuangfuDE",
    "Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq": "HuangfuESC",
    "Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq": "GersbachHep",
    "Engreitz_WTC11-endothelial-cells_TF-Perturb-seq": "EngreitzEndo",
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    frames = []
    for ds, short in SHORT_NAMES.items():
        p = args.input / ds / "energy_distance" / "wg1_significant_tfs.tsv"
        if not p.is_file():
            continue
        df = pd.read_csv(p, sep="\t")
        # Upstream metadata join can emit duplicate target_id rows (one per gene-symbol alias).
        df = df.drop_duplicates(subset="target_id", keep="first")
        df.insert(0, "dataset_id", ds)
        df.insert(1, "dataset_short", short)
        frames.append(df)

    if not frames:
        raise SystemExit(f"no wg1_significant_tfs.tsv files found under {args.input}")

    long = pd.concat(frames, ignore_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {args.output} ({len(long)} rows × {len(long.columns)} cols, "
          f"datasets: {long['dataset_short'].unique().tolist()})")


if __name__ == "__main__":
    main()
