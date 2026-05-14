"""Build the WG3-A disease-TF activity scorecard.

Inner-joins TF metadata + the gene-disease association table, then left-joins
each dataset's per-target ED data (auto-discovered from
`datasets/<id>/energy_distance/wg1_significant_tfs.tsv`).

Output: `working_groups/wg3_disease_gwas/examples/disease_tf_activity.tsv` — one
row per TF flagged as disease-relevant (per HPO MONDO + OMIM associations), with
per-dataset distance / rank / significance and a cross-dataset classification.

Re-run after `fetch_hpo_gene_disease.py` (annually-ish) and whenever new
per-dataset ED outputs land.

Usage:
    python working_groups/wg3_disease_gwas/examples/build_wg3_disease_tf_activity.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

JAMB = Path(__file__).resolve().parents[3]
TF_METADATA = JAMB / "reference/tf_metadata.tsv"
GENE_DISEASE = JAMB / "reference/gene_disease_associations.tsv"
DATASETS_DIR = JAMB / "datasets"
OUTPUT = JAMB / "working_groups/wg3_disease_gwas/examples/disease_tf_activity.tsv"

SHORT_TAGS = {
    "Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq": "HonCM",
    "Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq": "HuangfuDE",
    "Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq": "HuangfuESC",
    "Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq": "GersbachHep",
    "Engreitz_WTC11-endothelial-cells_TF-Perturb-seq": "EngreitzEndo",
}


def discover_datasets() -> list[tuple[str, Path]]:
    found = []
    for ds_dir in sorted(DATASETS_DIR.iterdir()):
        if not ds_dir.is_dir():
            continue
        candidate = ds_dir / "energy_distance" / "wg1_significant_tfs.tsv"
        if candidate.is_file():
            found.append((ds_dir.name, candidate))
    return found


def classify_cross_dataset(row: pd.Series, sig_cols: list[str]) -> str:
    values = [row[c] for c in sig_cols if not pd.isna(row[c])]
    if not values:
        return "no_data"
    if all(values):
        return "convergent_significant"
    if not any(values):
        return "convergent_nonsignificant"
    sig_datasets = [c.replace("sig_dist_gt_NC_max_", "") for c, v in zip(sig_cols, values) if v]
    if len(sig_datasets) == 1:
        return f"{sig_datasets[0]}-specific"
    return "discordant_partial"


def main() -> None:
    tf = pd.read_csv(TF_METADATA, sep="\t")[
        ["gene_symbol", "hgnc_approved_symbol", "ensembl_gene_id", "jaspar_tf_family", "lambert_2018_dbd"]
    ]
    gd = pd.read_csv(GENE_DISEASE, sep="\t")[
        ["gene_symbol", "n_disease_associations", "mondo_ids", "omim_ids"]
    ]

    # Inner join: only TFs with ≥1 disease association.
    disease_tfs = tf.merge(gd, on="gene_symbol", how="inner")
    print(f"TFs in library with ≥1 disease association: {len(disease_tfs)} / {len(tf)}")

    # Per-dataset wide-join — same auto-discovery as build_wg1_tf_cross_lineage
    found = discover_datasets()
    print(f"Found {len(found)} dataset(s) with per-target ED data:")
    for ds, p in found:
        print(f"  - {SHORT_TAGS.get(ds, ds)} ← {p.relative_to(JAMB)}")

    sig_cols: list[str] = []
    for ds, path in found:
        tag = SHORT_TAGS.get(ds, ds)
        df = pd.read_csv(path, sep="\t")
        # Aggregate per-target ED to per-gene (in case of multiple promoters per gene):
        # take the strongest distance per ensembl_gene_id.
        per_gene = (
            df[df["type"] == "targeting"]
            .groupby("ensembl_gene_id")
            .agg(
                **{
                    f"distance_mean_{tag}": ("distance_mean", "max"),
                    f"distance_rank_{tag}": ("distance_rank_targeting", "min"),
                    f"sig_dist_gt_NC_max_{tag}": ("sig_distance_gt_NC_max", "max"),
                }
            )
            .reset_index()
        )
        disease_tfs = disease_tfs.merge(per_gene, on="ensembl_gene_id", how="left")
        sig_cols.append(f"sig_dist_gt_NC_max_{tag}")

    # Cross-dataset summary
    dist_cols = [f"distance_mean_{SHORT_TAGS.get(ds, ds)}" for ds, _ in found]
    disease_tfs["n_datasets_with_data"] = disease_tfs[dist_cols].notna().sum(axis=1)
    disease_tfs["n_datasets_significant"] = disease_tfs[sig_cols].fillna(False).astype(bool).sum(axis=1)
    disease_tfs["max_distance_across_datasets"] = disease_tfs[dist_cols].max(axis=1)
    disease_tfs["classification"] = disease_tfs.apply(
        lambda r: classify_cross_dataset(r, sig_cols), axis=1
    )

    # Column order: identity → disease → per-dataset → summary
    identity = ["ensembl_gene_id", "gene_symbol", "hgnc_approved_symbol", "jaspar_tf_family", "lambert_2018_dbd"]
    disease_cols = ["n_disease_associations", "mondo_ids", "omim_ids"]
    per_ds: list[str] = []
    for ds, _ in found:
        tag = SHORT_TAGS.get(ds, ds)
        per_ds += [f"distance_mean_{tag}", f"distance_rank_{tag}", f"sig_dist_gt_NC_max_{tag}"]
    summary = ["n_datasets_with_data", "n_datasets_significant", "max_distance_across_datasets", "classification"]

    out = disease_tfs[identity + disease_cols + per_ds + summary].sort_values(
        ["n_datasets_significant", "max_distance_across_datasets"],
        ascending=[False, False],
        na_position="last",
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT, sep="\t", index=False)
    print(f"\nwrote {OUTPUT} ({len(out)} rows × {len(out.columns)} cols)")
    print("classification value counts:")
    print(out["classification"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
