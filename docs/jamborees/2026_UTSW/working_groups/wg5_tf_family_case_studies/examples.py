"""WG5 examples — TF-family deep-dives.

Run from the jamboree root:

    uv run python working_groups/wg5_tf_family_case_studies/examples.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd  # noqa: E402

from _lib import (  # noqa: E402
    landed_datasets,
    load_family_activity_scorecard,
    load_jaspar_tf_metadata,
    load_tf_gene_edges,
    load_tf_metadata,
)


# §1 — Candidate-for-deep-dive families --------------------------------------
print("\n§1 Candidate-for-deep-dive families")
fam = load_family_activity_scorecard()
candidates = fam[fam["candidate_for_deepdive"]]
sort_col = "n_sig_in_any_lineage" if "n_sig_in_any_lineage" in candidates.columns else "n_members"
print(f"   {len(candidates)} candidate families")
print(candidates.sort_values(sort_col, ascending=False).head(15)[
    ["family", "n_members", "n_disease_genes", "n_sig_in_any_lineage",
     "n_lineages_with_any_sig_member"]
].to_string(index=False))


# §2 — Pick a family, pull its members ---------------------------------------
# Family names in the scorecard mix Lambert DBDs (prefixed "DBD:") with JASPAR
# classes. Resolve both: "DBD:<x>" → match lambert_2018_dbd; otherwise match
# jaspar_tf_family. Default focus: DBD:C2H2 ZF (largest deep-dive candidate).
pick = None
for cand in ("DBD:C2H2 ZF", candidates.iloc[0]["family"]):
    if cand in fam["family"].values:
        pick = cand
        break

tf_meta = load_tf_metadata()
if pick is None:
    members = pd.DataFrame()
    print(f"\n§2 No family resolved.")
elif pick.startswith("DBD:"):
    target = pick[len("DBD:"):]
    members = tf_meta[tf_meta["lambert_2018_dbd"] == target]
    print(f"\n§2 Members of family '{pick}' (Lambert DBD '{target}'): {len(members)} TFs")
    print(members[["gene_symbol", "ensembl_gene_id", "jaspar_tf_family"]].head(15).to_string(index=False))
else:
    members = tf_meta[tf_meta["jaspar_tf_family"] == pick]
    print(f"\n§2 Members of family '{pick}' (JASPAR class): {len(members)} TFs")
    print(members[["gene_symbol", "ensembl_gene_id", "lambert_2018_dbd"]].head(15).to_string(index=False))


# §3 — Pull edges for those members across lineages --------------------------
# Filter on `tf_gene_symbol` (not `intended_target_name`, which is an ENSG ID).
landed = landed_datasets("wg4_tf_gene_edges")
if not members.empty and landed:
    member_symbols = set(members["gene_symbol"])
    print(f"\n§3 Edges for '{pick}' members in landed datasets {landed}")
    for short in landed:
        edges = load_tf_gene_edges(short)
        sub = edges[edges["tf_gene_symbol"].isin(member_symbols)]
        print(f"   {short}: {len(sub):,} edges from "
              f"{sub['tf_gene_symbol'].nunique()} members "
              f"→ {sub['gene_id'].nunique()} unique target genes")
else:
    print("\n§3 No edges to inspect (no family members or no landed datasets).")


# §4 — Motif lookup against JASPAR -------------------------------------------
# JASPAR's `name` column carries the TF symbol (sometimes a complex name like
# "Arnt::Ahr"); use a case-insensitive match.
try:
    jaspar = load_jaspar_tf_metadata()
    if not members.empty:
        upper_names = jaspar["name"].astype(str).str.upper()
        member_upper = {s.upper() for s in members["gene_symbol"]}
        hits = jaspar[upper_names.isin(member_upper)]
        print(f"\n§4 JASPAR motif matches for '{pick}' members: "
              f"{len(hits)} rows ({hits['name'].nunique()} unique TFs)")
        if not hits.empty:
            print(hits[["matrix_id", "name", "class", "family"]].head(10)
                  .to_string(index=False))
    else:
        print("\n§4 Skipping JASPAR lookup (no members).")
except FileNotFoundError as e:
    print(f"\n§4 JASPAR lookup skipped: {e}")


# §5 — GO-enrichment stub ----------------------------------------------------
print("\n§5 GO enrichment of family target genes")
print("""
    # Aggregate per-family target gene list across lineages, then enrich.
    # Suggested: `gseapy.enrichr(gene_list=..., gene_sets=['GO_Biological_Process_2023'])`
    # `uv add gseapy` if needed. Watch out for multiple-testing across families.
""")
