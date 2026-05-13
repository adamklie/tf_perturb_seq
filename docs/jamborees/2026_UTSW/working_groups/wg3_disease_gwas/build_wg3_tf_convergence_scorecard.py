"""Build the WG3-B TF convergence / divergence scorecard.

Filters WG3-A's `disease_tf_activity.tsv` to TFs with energy-distance data in
≥2 lineages, and adds a quantitative cross-lineage magnitude/ratio analysis.

Classification refinement (orthogonal to WG3-A's sig-based labels):
- `convergent_high`: significant in all lineages with data (high-magnitude
  agreement)
- `convergent_low`: not significant in any lineage with data, AND max distance
  across datasets falls in the bottom 75% (genuinely quiet TFs across lineages)
- `convergent_moderate`: not significant in any lineage but max distance ≥75th
  percentile (probably suppressed by calibration / borderline cases worth
  watching)
- `divergent_<lineage>`: significant in exactly one lineage, not others

Magnitude-ratio columns help WG3 sort the "TFs that flip on/off across
lineages" hit list for deep-dive selection.

Output: `working_groups/wg3_disease_gwas/tf_convergence_scorecard.tsv`

Usage:
    python scripts/build_wg3_tf_convergence_scorecard.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

JAMB = Path(__file__).resolve().parent.parent
SOURCE = JAMB / "working_groups/wg3_disease_gwas/disease_tf_activity.tsv"
OUTPUT = JAMB / "working_groups/wg3_disease_gwas/tf_convergence_scorecard.tsv"


def main() -> None:
    src = pd.read_csv(SOURCE, sep="\t")
    print(f"WG3-A source: {len(src)} rows × {len(src.columns)} cols")

    # Identify per-dataset distance and significance columns
    dist_cols = [c for c in src.columns if c.startswith("distance_mean_")]
    sig_cols = [c for c in src.columns if c.startswith("sig_dist_gt_NC_max_")]
    if not dist_cols:
        raise SystemExit("no distance_mean_* columns in WG3-A source")

    # Subset: data in ≥2 lineages
    df = src[src["n_datasets_with_data"] >= 2].copy().reset_index(drop=True)
    print(f"TFs with data in ≥2 lineages: {len(df)} / {len(src)}")

    # Per-TF magnitude stats across the lineages that have data
    df["min_distance_across_datasets"] = df[dist_cols].min(axis=1)
    df["distance_range"] = df["max_distance_across_datasets"] - df["min_distance_across_datasets"]
    # Fold-change-like ratio: max/min, capped at 1e6 to avoid inf when min ~ 0.
    safe_min = df["min_distance_across_datasets"].replace(0, np.nan)
    df["distance_max_over_min"] = (df["max_distance_across_datasets"] / safe_min).clip(upper=1e6)

    # Quantile threshold for "high magnitude" (within disease-TF cohort): 75th
    high_thresh = float(np.nanpercentile(df["max_distance_across_datasets"], 75))
    print(f"75th-percentile max_distance (high-magnitude cutoff): {high_thresh:.4f}")

    def refined_class(row: pd.Series) -> str:
        sig_values = [bool(row[c]) for c in sig_cols if pd.notna(row[c])]
        n_sig = sum(sig_values)
        n_with_data = int(row["n_datasets_with_data"])
        max_dist = float(row["max_distance_across_datasets"])
        if n_sig == n_with_data and n_sig >= 1:
            return "convergent_high"
        if n_sig == 0:
            if max_dist >= high_thresh:
                return "convergent_moderate"
            return "convergent_low"
        # Mixed: find which lineage is sig
        sig_datasets = [
            c.replace("sig_dist_gt_NC_max_", "")
            for c in sig_cols
            if pd.notna(row[c]) and bool(row[c])
        ]
        if len(sig_datasets) == 1:
            return f"divergent_{sig_datasets[0]}"
        return "divergent_partial"

    df["convergence_class"] = df.apply(refined_class, axis=1)

    # Reorder for slide-deck readability: identity → disease → per-dataset → magnitude → class
    identity = ["ensembl_gene_id", "gene_symbol", "hgnc_approved_symbol", "jaspar_tf_family", "lambert_2018_dbd"]
    disease = ["n_disease_associations", "mondo_ids", "omim_ids"]
    per_ds = [c for c in src.columns if c.startswith(("distance_mean_", "distance_rank_", "sig_dist_gt_NC_max_"))]
    magnitude = [
        "n_datasets_with_data",
        "n_datasets_significant",
        "min_distance_across_datasets",
        "max_distance_across_datasets",
        "distance_range",
        "distance_max_over_min",
    ]
    classification = ["classification", "convergence_class"]
    out = df[identity + disease + per_ds + magnitude + classification]

    # Sort: convergent_high first, then divergent (sorted by max_dist), then convergent_moderate, then convergent_low
    sort_key = {"convergent_high": 0, "divergent_partial": 1}
    # divergent_<lineage> entries sort after convergent_high but before convergent_moderate
    out = out.assign(
        _rank=out["convergence_class"].map(
            lambda c: sort_key.get(c, 2 if c.startswith("divergent_") else 3 if c == "convergent_moderate" else 4)
        )
    )
    out = out.sort_values(["_rank", "max_distance_across_datasets"], ascending=[True, False]).drop(columns="_rank")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT, sep="\t", index=False)
    print(f"\nwrote {OUTPUT} ({len(out)} rows × {len(out.columns)} cols)")
    print("convergence_class value counts:")
    print(out["convergence_class"].value_counts().to_string())


if __name__ == "__main__":
    main()
