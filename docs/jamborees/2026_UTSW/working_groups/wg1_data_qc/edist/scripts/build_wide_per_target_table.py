"""Build a wide cross-dataset table — one row per target, dataset stats side-by-side.

Input:  directory containing per-dataset `wg1_significant_tfs.tsv` files
        (one subdirectory per dataset, e.g. data/<dataset_id>/energy_distance/wg1_significant_tfs.tsv).
Output: single wide TSV. Identity columns (target_id, gene_symbol, …) on the
        left; for each dataset, suffixed stat columns
        (`{short}_distance_mean`, `{short}_pval_mean`, `{short}_sig_*`, etc.).
        Targets missing from a dataset come through as NaN (outer join).

Usage:
    uv run python build_wide_per_target_table.py \\
        --input  ../../../../data \\
        --output ../results/distance_heatmap/per_target_wide.tsv
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

IDENTITY_COLS = [
    "ensembl_gene_id", "gene_symbol", "hgnc_approved_symbol",
    "jaspar_tf_family", "lambert_2018_dbd", "locus",
]
STAT_COLS = [
    "type", "cell_count", "distance_mean", "pval_mean",
    "sig_distance_gt_NC_max", "sig_pval_lt_0p05",
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, type=Path, help="data/ root with <dataset>/energy_distance/wg1_significant_tfs.tsv")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    frames: dict[str, pd.DataFrame] = {}
    for ds, short in SHORT_NAMES.items():
        p = args.input / ds / "energy_distance" / "wg1_significant_tfs.tsv"
        if not p.is_file():
            continue
        df = pd.read_csv(p, sep="\t")
        # Upstream metadata join can emit duplicate target_id rows (one per gene-symbol alias);
        # stat columns are identical within a group, so dedupe by first.
        df = df.drop_duplicates(subset="target_id", keep="first").set_index("target_id")
        frames[short] = df

    if not frames:
        raise SystemExit(f"no wg1_significant_tfs.tsv files found under {args.input}")

    # Identity columns: combine_first across datasets (prefer the first one with a value)
    identity = pd.concat(
        [df[IDENTITY_COLS] for df in frames.values()]
    ).groupby(level=0).first()

    # Stat columns: rename per-dataset with short_name prefix, then join
    wide = identity.copy()
    for short, df in frames.items():
        renamed = df[STAT_COLS].rename(columns=lambda c: f"{short}_{c}")
        wide = wide.join(renamed, how="outer")

    # Restore target_id as a column
    wide = wide.reset_index().rename(columns={"index": "target_id"})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {args.output} ({len(wide)} rows × {len(wide.columns)} cols, datasets: {list(frames)})")


if __name__ == "__main__":
    main()
