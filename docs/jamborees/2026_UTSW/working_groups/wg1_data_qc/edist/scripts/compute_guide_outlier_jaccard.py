"""Pairwise Jaccard between outlier-gRNA sets across datasets.

Outlier gRNAs are guides with `pval_outlier < threshold` in the energy-distance
pipeline's `targeting_outlier_table.csv` (lower pval → more outlier-like vs.
its sibling guides).

Input:  --datasets-root pointing at the `datasets/` directory. For each entry in
        DATASETS, globs `<datasets_root>/<dataset_id>/*/energy_distance/targeting_outlier_table.csv`
        — discovers the run-label automatically. If a dataset has multiple runs
        with energy_distance outputs, pass --run-label to disambiguate.
Output: symmetric Jaccard TSV (datasets × datasets), with set sizes on the
        diagonal-adjacent metadata rows.

Usage:
    uv run python compute_guide_outlier_jaccard.py \\
        --datasets-root ../../../../../../datasets \\
        --output ../results/guide_outlier_jaccard/guide_outlier_jaccard.tsv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

# short_name -> dataset_id (folder under datasets/). Run-label is discovered via glob.
DATASETS = {
    "HonCM":       "Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq",
    "HuangfuDE":   "Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq",
    "HuangfuESC":  "Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq",
    "GersbachHep": "Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq",
}


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return float("nan")
    return len(a & b) / len(a | b)


def find_outlier_table(datasets_root: Path, dataset_id: str, run_label: str | None) -> Path | None:
    pattern = f"{run_label}/energy_distance/targeting_outlier_table.csv" if run_label \
              else "*/energy_distance/targeting_outlier_table.csv"
    matches = sorted((datasets_root / dataset_id).glob(pattern))
    if not matches:
        return None
    if len(matches) > 1:
        print(f"  warning: {dataset_id} has {len(matches)} matches; using {matches[0].parent.parent.name}/")
    return matches[0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--datasets-root", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--pval-threshold", type=float, default=0.05,
                    help="pval_outlier < threshold → outlier (default 0.05)")
    ap.add_argument("--run-label", default=None,
                    help="Filter to a specific run subdir (e.g. 'sceptre_v1'); default is glob over all runs.")
    args = ap.parse_args()

    outlier_sets: dict[str, set] = {}
    set_sizes: dict[str, int] = {}
    total_guides: dict[str, int] = {}
    for short, dataset_id in DATASETS.items():
        p = find_outlier_table(args.datasets_root, dataset_id, args.run_label)
        if p is None:
            print(f"  skip (no targeting_outlier_table.csv): {short}")
            continue
        df = pd.read_csv(p, index_col=0)
        outliers = set(df.index[df["pval_outlier"] < args.pval_threshold])
        outlier_sets[short] = outliers
        set_sizes[short] = len(outliers)
        total_guides[short] = len(df)
        print(f"  {short}: {len(outliers)} / {len(df)} outlier gRNAs  ({p.parent.parent.name}/)")

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
