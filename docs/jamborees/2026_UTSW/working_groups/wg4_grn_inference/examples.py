"""WG4 examples — load per-dataset TF→gene edges, inspect network structure.

Run from the jamboree root:

    uv run python working_groups/wg4_grn_inference/examples.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd  # noqa: E402

from _lib import (  # noqa: E402
    landed_datasets,
    load_network_structure_by_lineage,
    load_tf_gene_edges,
)


# §1 — Load edges per dataset as DataFrames ----------------------------------
landed = landed_datasets("wg4_tf_gene_edges")
print(f"\n§1 Edges landed for: {landed}")
edges = {short: load_tf_gene_edges(short) for short in landed}
for short, df in edges.items():
    print(f"   {short:14s}  {len(df):>8,} edges  "
          f"({df['intended_target_name'].nunique()} TFs × "
          f"{df['gene_id'].nunique()} unique target genes)")


# §2 — Cross-lineage network structure (rolled up) ---------------------------
print("\n§2 Per-lineage network structure")
ns = load_network_structure_by_lineage()
print(ns[
    ["dataset_tag", "n_tfs_with_sig_edges", "n_edges",
     "median_tf_outdegree", "max_tf_outdegree",
     "n_tf_tf_edges", "fraction_tf_tf_edges"]
].to_string(index=False))


# §3 — Build NetworkX graphs + degree distributions --------------------------
# Use `tf_gene_symbol` as the source label (intended_target_name is an ENSG).
try:
    import networkx as nx
    print("\n§3 NetworkX degree summaries")
    for short, df in edges.items():
        G = nx.from_pandas_edgelist(
            df, "tf_gene_symbol", "gene_id",
            edge_attr=["log2_fc", "fdr_bh"], create_using=nx.DiGraph,
        )
        out_deg = pd.Series(dict(G.out_degree())).pipe(lambda s: s[s > 0])
        print(f"\n   {short}")
        print(f"     n_nodes={G.number_of_nodes()}, n_edges={G.number_of_edges()}")
        print("     top 10 TFs by out-degree:")
        print(out_deg.sort_values(ascending=False).head(10).to_string())
except ImportError:
    print("\n§3 networkx not installed; `uv add networkx` to enable.")


# §4 — Cross-lineage edge-set Jaccard (sanity check vs network_structure_*.tsv)
print("\n§4 Pairwise edge-set Jaccards (recomputed)")
edge_sets = {
    short: set(zip(df["tf_gene_symbol"], df["gene_id"]))
    for short, df in edges.items()
}
shorts = sorted(edge_sets)
print("   ", "  ".join(f"{s:>12s}" for s in shorts))
for s1 in shorts:
    row = [f"{s1:14s}"]
    for s2 in shorts:
        a, b = edge_sets[s1], edge_sets[s2]
        j = len(a & b) / max(1, len(a | b))
        row.append(f"{j:>12.4f}")
    print("   ".join(row))


# §5 — Filter edges by TF DBD (example: C2H2 zinc-fingers) ------------------
# `tf_dbd` follows Lambert 2018 naming (e.g. "C2H2 ZF", "Homeodomain", "bHLH").
# `tf_family` follows JASPAR class (e.g. "Factors with multiple dispersed
# zinc fingers"). Pick the column that matches your filter vocabulary.
print("\n§5 Filter edges by TF DBD (example: C2H2 ZF)")
for short, df in edges.items():
    if "tf_dbd" not in df.columns:
        continue
    fam = df[df["tf_dbd"].astype(str).str.contains("C2H2 ZF", na=False)]
    print(f"   {short}: {len(fam):,} edges from "
          f"{fam['tf_gene_symbol'].nunique()} C2H2-ZF-family TFs")


# §6 — Multiome integration placeholder --------------------------------------
print("\n§6 Multiome (E2G + ChromBPNet) integration")
print("    Not staged in this folder. Inputs live in ref/scE2G_links/ and the "
      "per-cell-type ChromBPNet outputs from the Engreitz/Mostafavi pipelines. "
      "Plan: project E2G enhancer→gene links onto these edge sets to flag "
      "candidate direct vs. indirect targets per TF. Track in TODO.md.")
