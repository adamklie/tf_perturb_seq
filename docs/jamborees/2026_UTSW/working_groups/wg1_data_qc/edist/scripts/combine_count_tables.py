"""Concatenate per-dataset count TSVs into a single table.

Input:  directory of 1-row TSVs (one per dataset) produced by count_significant_tfs.py
Output: single combined TSV (one row per dataset).

Usage:
    uv run python combine_count_tables.py \\
        --input  ../results/significant_tf_counts/per_dataset \\
        --output ../results/significant_tf_counts/significant_tf_counts.tsv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, type=Path, help="directory of per-dataset TSVs")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    paths = sorted(args.input.glob("*.tsv"))
    if not paths:
        raise SystemExit(f"no TSVs found in {args.input}")

    combined = pd.concat([pd.read_csv(p, sep="\t") for p in paths], ignore_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {args.output} ({len(combined)} rows)")


if __name__ == "__main__":
    main()
