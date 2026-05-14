"""Build the WG5-A TF family activity scorecard.

Groups TFs by `jaspar_tf_family` (falls back to `lambert_2018_dbd` for TFs not
in JASPAR core, then "unannotated"). For each family with ≥3 members in the
library, computes:
  - n_members (library size)
  - per-dataset: n_members_significant + max_distance
  - cross-dataset: n_sig_in_any_lineage + n_sig_in_all_lineages_with_data
  - n_disease_genes / fraction_disease_genes (via HPO MONDO+OMIM)
  - candidate_for_deepdive flag (≥2 sig members in any lineage AND ≥1 disease gene)

Re-runs auto-discover datasets via
`datasets/<id>/energy_distance/wg1_significant_tfs.tsv` (so the per-dataset
columns auto-widen as new datasets land).

Output: `working_groups/wg5_tf_family_case_studies/examples/family_activity_scorecard.tsv`

Usage:
    python working_groups/wg5_tf_family_case_studies/examples/build_wg5_tf_family_scorecard.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

JAMB = Path(__file__).resolve().parents[3]
TF_METADATA = JAMB / "reference/tf_metadata.tsv"
GENE_DISEASE = JAMB / "reference/gene_disease_associations.tsv"
DATASETS_DIR = JAMB / "datasets"
OUTPUT = JAMB / "working_groups/wg5_tf_family_case_studies/examples/family_activity_scorecard.tsv"

MIN_FAMILY_SIZE = 3

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


def family_label(row: pd.Series) -> str:
    """Group key: JASPAR family if present, else Lambert DBD, else 'unannotated'."""
    j = row.get("jaspar_tf_family")
    if isinstance(j, str) and j.strip():
        return j.strip()
    d = row.get("lambert_2018_dbd")
    if isinstance(d, str) and d.strip() and d.strip() != "Unknown":
        return f"DBD:{d.strip()}"
    return "unannotated"


def main() -> None:
    tf = pd.read_csv(TF_METADATA, sep="\t")
    tf["family"] = tf.apply(family_label, axis=1)

    gd = pd.read_csv(GENE_DISEASE, sep="\t")
    disease_symbols = set(gd["gene_symbol"].unique())
    tf["is_disease_gene"] = tf["gene_symbol"].isin(disease_symbols)

    # Family base stats
    families = tf.groupby("family", as_index=False).agg(
        n_members=("gene_symbol", "nunique"),
        n_disease_genes=("is_disease_gene", "sum"),
    )
    families["fraction_disease_genes"] = families["n_disease_genes"] / families["n_members"]
    families = families[families["n_members"] >= MIN_FAMILY_SIZE].reset_index(drop=True)
    print(f"Families with ≥{MIN_FAMILY_SIZE} members: {len(families)}")

    found = discover_datasets()
    print(f"Found {len(found)} dataset(s) with per-target ED data:")
    for ds, p in found:
        print(f"  - {SHORT_TAGS.get(ds, ds)} ← {p.relative_to(JAMB)}")

    sig_cols: list[str] = []
    for ds, path in found:
        tag = SHORT_TAGS.get(ds, ds)
        per_target = pd.read_csv(path, sep="\t")
        per_target = per_target[per_target["type"] == "targeting"]
        # Aggregate per ensembl_gene_id (in case of multiple promoters per gene)
        per_gene = (
            per_target.groupby("ensembl_gene_id")
            .agg(
                distance_mean_max=("distance_mean", "max"),
                sig_dist_gt_NC_max=("sig_distance_gt_NC_max", "max"),
            )
            .reset_index()
        )
        tf_ds = tf.merge(per_gene, on="ensembl_gene_id", how="left")

        agg = (
            tf_ds.groupby("family", as_index=False)
            .agg(
                **{
                    f"n_members_with_data_{tag}": ("distance_mean_max", lambda s: int(s.notna().sum())),
                    f"n_members_sig_{tag}": (
                        "sig_dist_gt_NC_max",
                        lambda s: int(pd.to_numeric(s, errors="coerce").fillna(0).astype(int).sum()),
                    ),
                    f"max_distance_{tag}": ("distance_mean_max", "max"),
                    f"mean_distance_{tag}": ("distance_mean_max", "mean"),
                }
            )
        )
        families = families.merge(agg, on="family", how="left")
        sig_cols.append(f"n_members_sig_{tag}")

    # Cross-dataset rollups
    families["n_sig_in_any_lineage"] = families[sig_cols].fillna(0).sum(axis=1).astype(int)
    families["n_lineages_with_any_sig_member"] = (
        families[sig_cols].fillna(0).gt(0).sum(axis=1).astype(int)
    )
    # Exclude the "unannotated" catch-all — it's a coverage gap, not a coherent family.
    families["candidate_for_deepdive"] = (
        (families["n_sig_in_any_lineage"] >= 2)
        & (families["n_disease_genes"] >= 1)
        & (families["family"] != "unannotated")
    )

    # Sort: candidates first, then by n_sig_in_any_lineage desc
    families = families.sort_values(
        ["candidate_for_deepdive", "n_sig_in_any_lineage", "n_members"],
        ascending=[False, False, False],
    ).reset_index(drop=True)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    families.to_csv(OUTPUT, sep="\t", index=False)
    print(f"\nwrote {OUTPUT} ({len(families)} rows × {len(families.columns)} cols)")
    print(f"  candidates_for_deepdive: {families['candidate_for_deepdive'].sum()}")
    print("\nTop 15 by n_sig_in_any_lineage:")
    cols_show = ["family", "n_members", "n_disease_genes", "n_sig_in_any_lineage", "candidate_for_deepdive"]
    print(families[cols_show].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
