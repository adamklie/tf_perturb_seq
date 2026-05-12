"""WG3 examples — disease/GWAS-gene TFs, convergent vs divergent activity.

Run from the jamboree root:

    uv run python working_groups/wg3_disease_gwas/examples.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd  # noqa: E402

from _lib import (  # noqa: E402
    landed_datasets,
    load_disease_tf_activity,
    load_gene_disease_associations,
    load_tf_convergence_scorecard,
    load_tf_gene_edges,
)


# §1 — Disease-flagged TFs, sorted by max distance across lineages ------------
print("\n§1 Top-20 disease-flagged TFs by max distance across lineages")
dta = load_disease_tf_activity()
sig_cols = [c for c in dta.columns if c.startswith("sig_dist_gt_NC_max_")]
dta["max_distance"] = dta[
    [c for c in dta.columns if c.startswith("distance_mean_")]
].max(axis=1)
print(dta.sort_values("max_distance", ascending=False).head(20)[
    ["gene_symbol", "jaspar_tf_family", "n_disease_associations",
     "n_datasets_significant", "max_distance"]
].to_string(index=False))


# §2 — Convergence-class breakdown -------------------------------------------
print("\n§2 Convergence-class counts (refined classification)")
sc = load_tf_convergence_scorecard()
print(sc["convergence_class"].value_counts().to_string())


# §3 — Deep-dive: pick a convergent TF that actually has edges --------------
# `intended_target_name` in the edges file is an ENSG ID, so we filter via the
# `tf_gene_symbol` column. Pick the convergent TF with the most edges across
# datasets (avoids picking a high-distance TF that doesn't survive FDR<0.05).
convergent = sc[sc["convergence_class"].str.startswith("convergent", na=False)]
landed_wg4 = landed_datasets("wg4_tf_gene_edges")
target_tf = None
if not convergent.empty and landed_wg4:
    convergent_symbols = set(convergent["gene_symbol"])
    edge_count = pd.Series(dtype=int)
    for short in landed_wg4:
        edges = load_tf_gene_edges(short)
        sub = edges[edges["tf_gene_symbol"].isin(convergent_symbols)]
        c = sub.groupby("tf_gene_symbol").size()
        edge_count = edge_count.add(c, fill_value=0)
    if not edge_count.empty:
        target_tf = edge_count.sort_values(ascending=False).index[0]

if target_tf:
    print(f"\n§3 Downstream targets for {target_tf} per dataset (FDR<0.05)")
    for short in landed_wg4:
        edges = load_tf_gene_edges(short)
        sub = edges[edges["tf_gene_symbol"] == target_tf]
        n_targets = sub["gene_id"].nunique()
        n_up = (sub["log2_fc"] > 0).sum()
        n_down = (sub["log2_fc"] < 0).sum()
        print(f"   {short:14s}  {n_targets:5d} targets  ({n_up} up, {n_down} down)")
else:
    print("\n§3 No convergent TF with landed edges yet.")


# §4 — Disease-gene lookup ---------------------------------------------------
print("\n§4 Disease associations for the deep-dive TF")
gda = load_gene_disease_associations()
if target_tf:
    rows = gda[gda["gene_symbol"] == target_tf].head(10)
    if rows.empty:
        print(f"   No disease associations for {target_tf}.")
    else:
        print(rows.to_string(index=False))
else:
    print("   Skipped (no §3 deep-dive target).")


# §5 — Cross-lineage discordant TFs ------------------------------------------
print("\n§5 Discordant TFs (sig in some lineages but not others)")
disc = sc[sc["convergence_class"].astype(str).str.contains("divergent", na=False)]
print(f"   {len(disc)} divergent TFs total")
print(disc.sort_values("distance_range", ascending=False).head(10)[
    ["gene_symbol", "convergence_class", "min_distance_across_datasets",
     "max_distance_across_datasets", "distance_range", "distance_max_over_min"]
].to_string(index=False))

print("\n⚠ Calibration caveat: Huangfu DE/ESC p-values are anti-conservative. "
      "All classifications use `distance_mean > NC max` as the proxy. "
      "Most 'divergent_HuangfuDE' calls reflect lower absolute distance in DE, "
      "not a real lineage-specific effect. See issues/edistance-calibration/.")
