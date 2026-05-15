"""Pairwise Jaccard between outlier-gRNA sets across datasets.

Outlier gRNAs are guides with `pval_outlier < threshold` in the energy-distance
pipeline's `targeting_outlier_table.csv` (lower pval → more outlier-like vs.
its sibling guides).

Input:  --datasets-root pointing at the `datasets/` directory; per-dataset
        targeting_outlier_table.csv paths are hard-coded below since
        run-label sub-paths differ across datasets.
Output: symmetric Jaccard TSV (datasets × datasets), with set sizes on the
        diagonal-adjacent metadata rows.

Usage:
    uv run python compute_guide_outlier_jaccard.py \\
        --datasets-root ../../../../../../datasets \\
        --output ../results/guide_outlier_jaccard.tsv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

DATASET_PATHS = {
    "HonCM": "Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/2026_04_19_no_spacer/energy_distance/targeting_outlier_table.csv",
    "HuangfuDE": "Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/muddy_penguin/energy_distance/targeting_outlier_table.csv",
    "HuangfuESC": "Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/sceptre_v1/energy_distance/targeting_outlier_table.csv",
    "GersbachHep": "Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/sara_synapse_syn74842722/energy_distance/targeting_outlier_table.csv",
}


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return float("nan")
    return len(a & b) / len(a | b)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--datasets-root", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--pval-threshold", type=float, default=0.05,
                    help="pval_outlier < threshold → outlier (default 0.05)")
    args = ap.parse_args()

    outlier_sets: dict[str, set] = {}
    set_sizes: dict[str, int] = {}
    total_guides: dict[str, int] = {}
    for short, rel in DATASET_PATHS.items():
        p = args.datasets_root / rel
        if not p.is_file():
            print(f"  skip (missing) {short}: {p}")
            continue
        df = pd.read_csv(p, index_col=0)
        outliers = set(df.index[df["pval_outlier"] < args.pval_threshold])
        outlier_sets[short] = outliers
        set_sizes[short] = len(outliers)
        total_guides[short] = len(df)
        print(f"  {short}: {len(outliers)} / {len(df)} outlier gRNAs")

    if not outlier_sets:
        raise SystemExit("no targeting_outlier_table.csv files found")

    keys = list(outlier_sets)
    mat = pd.DataFrame(index=keys, columns=keys, dtype=float)
    for a in keys:
        for b in keys:
            mat.loc[a, b] = jaccard(outlier_sets[a], outlier_sets[b])

    # Save with set sizes / total guides as a sidecar block
    out = mat.copy()
    out.insert(0, "n_outlier_grnas", pd.Series(set_sizes))
    out.insert(0, "n_total_grnas", pd.Series(total_guides))
    out.index.name = "dataset"

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
