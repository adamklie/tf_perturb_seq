"""Build the WG1-D cross-lineage TF activity table.

For each TF that has data in ≥1 dataset, pivots the per-dataset
`wg1_significant_tfs.tsv` files into a single wide-format TSV — one row per
(target_id), columns per dataset for distance / rank / significance, plus a
cross-dataset summary classification.

Drives WG1's "compare e-distance runs across overlapping datasets" calibration
debug task: any time you want to know whether TFs called significant in one
lineage are the same as those in another, this is the file.

Datasets are auto-discovered by scanning
`docs/jamborees/2026_UTSW/datasets/<dataset>/energy_distance/wg1_significant_tfs.tsv`,
so the table widens automatically as more datasets land.

Usage:
    python scripts/build_wg1_tf_cross_lineage.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

JAMB = Path(__file__).resolve().parent.parent / "docs/jamborees/2026_UTSW"
DATASETS_DIR = JAMB / "datasets"
OUTPUT = JAMB / "working_groups/wg1_data_qc/tf_cross_lineage.tsv"

# Short tag per dataset for the column suffix (avoids unwieldy column names).
SHORT_TAGS = {
    "Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq": "HonCM",
    "Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq": "HuangfuDE",
    "Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq": "HuangfuESC",
    "Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq": "GersbachHep",
    "Engreitz_WTC11-endothelial-cells_TF-Perturb-seq": "EngreitzEndo",
}

# Identity columns kept from the per-dataset table (taken from the first dataset
# that has data for a given target). These never disagree across datasets for
# the same target_id.
IDENTITY_COLS = [
    "target_id",
    "ensembl_gene_id",
    "gene_symbol",
    "hgnc_approved_symbol",
    "jaspar_tf_family",
    "lambert_2018_dbd",
    "type",
    "locus",
]


def discover_datasets() -> list[tuple[str, Path]]:
    """Find every dataset with a per-target ED TSV. Returns sorted list."""
    found = []
    for ds_dir in sorted(DATASETS_DIR.iterdir()):
        if not ds_dir.is_dir():
            continue
        candidate = ds_dir / "energy_distance" / "wg1_significant_tfs.tsv"
        if candidate.is_file():
            found.append((ds_dir.name, candidate))
    return found


def classify(row: pd.Series, sig_cols: list[str]) -> str:
    """Cross-dataset classification based on the per-dataset significance bools."""
    values = [row[c] for c in sig_cols if not pd.isna(row[c])]
    if not values:
        return "no_data"
    if all(values):
        return "convergent_significant"
    if not any(values):
        return "convergent_nonsignificant"
    # Mixed: identify which dataset(s) it's specific to.
    sig_datasets = [c.replace("sig_dist_gt_NC_max_", "") for c, v in zip(sig_cols, values) if v]
    if len(sig_datasets) == 1:
        return f"{sig_datasets[0]}-specific"
    return "discordant_partial"


def main() -> None:
    found = discover_datasets()
    if not found:
        raise SystemExit("no per-dataset wg1_significant_tfs.tsv found under datasets/")

    print(f"Found {len(found)} dataset(s) with per-target ED data:")
    for ds, path in found:
        tag = SHORT_TAGS.get(ds, ds)
        print(f"  - {tag} ← {path.relative_to(JAMB)}")

    # Build a per-target wide table by sequential left-join.
    merged: pd.DataFrame | None = None
    sig_cols: list[str] = []
    for ds, path in found:
        tag = SHORT_TAGS.get(ds, ds)
        df = pd.read_csv(path, sep="\t")
        keep = (
            IDENTITY_COLS
            + ["distance_mean", "distance_rank_targeting", "sig_distance_gt_NC_max"]
        )
        df = df[keep].rename(
            columns={
                "distance_mean": f"distance_mean_{tag}",
                "distance_rank_targeting": f"distance_rank_{tag}",
                "sig_distance_gt_NC_max": f"sig_dist_gt_NC_max_{tag}",
            }
        )
        sig_cols.append(f"sig_dist_gt_NC_max_{tag}")
        if merged is None:
            merged = df
        else:
            # Right-side has only the per-dataset cols (identity already in left)
            right = df.drop(columns=IDENTITY_COLS, errors="ignore")
            right["target_id"] = df["target_id"]
            merged = merged.merge(right, on="target_id", how="outer")

    assert merged is not None
    # Cross-dataset summary columns
    merged["n_datasets_with_data"] = merged[
        [f"distance_mean_{SHORT_TAGS.get(ds, ds)}" for ds, _ in found]
    ].notna().sum(axis=1)
    merged["n_datasets_significant"] = merged[sig_cols].fillna(False).astype(bool).sum(axis=1)
    merged["classification"] = merged.apply(lambda r: classify(r, sig_cols), axis=1)

    # Column order: identity, per-dataset (grouped by dataset), cross summary
    per_ds_cols: list[str] = []
    for ds, _ in found:
        tag = SHORT_TAGS.get(ds, ds)
        per_ds_cols += [
            f"distance_mean_{tag}",
            f"distance_rank_{tag}",
            f"sig_dist_gt_NC_max_{tag}",
        ]
    summary_cols = ["n_datasets_with_data", "n_datasets_significant", "classification"]
    out = merged[IDENTITY_COLS + per_ds_cols + summary_cols]

    # Sort: convergent_significant first, then discordant_partial, then by max distance
    classification_order = pd.CategoricalDtype(
        categories=[
            "convergent_significant",
            "discordant_partial",
        ]
        + [f"{SHORT_TAGS.get(ds, ds)}-specific" for ds, _ in found]
        + ["convergent_nonsignificant", "no_data"],
        ordered=True,
    )
    out["classification"] = out["classification"].astype(classification_order)
    out["_sort_dist"] = out[[f"distance_mean_{SHORT_TAGS.get(ds, ds)}" for ds, _ in found]].max(axis=1)
    out = out.sort_values(["classification", "_sort_dist"], ascending=[True, False]).drop(columns="_sort_dist")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT, sep="\t", index=False)
    print(f"\nwrote {OUTPUT} ({len(out)} rows × {len(out.columns)} cols)")
    print("classification value counts:")
    print(out["classification"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
