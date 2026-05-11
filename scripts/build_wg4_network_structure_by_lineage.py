"""Build the WG4-B per-lineage network structure summary.

For each lineage, computes degree-distribution and TF→TF connectivity stats off
the WG4-A edge list. Adds cross-lineage edge/TF overlap columns that auto-widen
as more datasets land.

Datasets are auto-discovered by scanning
`docs/jamborees/2026_UTSW/datasets/<dataset>/crispr_pipeline/wg4_tf_gene_edges_FDR05.tsv`.

Output: `working_groups/wg4_grn_inference/network_structure_by_lineage.tsv`

Usage:
    python scripts/build_wg4_network_structure_by_lineage.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

JAMB = Path(__file__).resolve().parent.parent / "docs/jamborees/2026_UTSW"
DATASETS_DIR = JAMB / "datasets"
OUTPUT = JAMB / "working_groups/wg4_grn_inference/network_structure_by_lineage.tsv"

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
        candidate = ds_dir / "crispr_pipeline" / "wg4_tf_gene_edges_FDR05.tsv"
        if candidate.is_file():
            found.append((ds_dir.name, candidate))
    return found


def per_lineage_stats(edges: pd.DataFrame) -> dict[str, float | int]:
    """Compute core network stats from a single dataset's edge table."""
    perturbed_tfs = set(edges["intended_target_name"].unique())
    out_deg = edges.groupby("intended_target_name").size()
    in_deg = edges.groupby("gene_id").size()
    # TF→TF edges: target gene_id is itself a perturbed TF (ensembl id matches
    # any intended_target_name in the same dataset)
    tf_tf_mask = edges["gene_id"].isin(perturbed_tfs)
    return {
        "n_tfs_with_sig_edges": len(perturbed_tfs),
        "n_target_genes": edges["gene_id"].nunique(),
        "n_edges": len(edges),
        "mean_tf_outdegree": float(out_deg.mean()),
        "median_tf_outdegree": float(out_deg.median()),
        "max_tf_outdegree": int(out_deg.max()),
        "p95_tf_outdegree": float(np.percentile(out_deg.values, 95)),
        "mean_gene_indegree": float(in_deg.mean()),
        "median_gene_indegree": float(in_deg.median()),
        "max_gene_indegree": int(in_deg.max()),
        "n_tf_tf_edges": int(tf_tf_mask.sum()),
        "fraction_tf_tf_edges": float(tf_tf_mask.mean()),
        "median_abs_log2fc": float(edges["log2_fc"].abs().median()),
    }


def main() -> None:
    found = discover_datasets()
    if not found:
        raise SystemExit("no per-dataset wg4_tf_gene_edges_FDR05.tsv found under datasets/")

    print(f"Found {len(found)} dataset(s) with WG4-A edge lists:")
    for ds, path in found:
        tag = SHORT_TAGS.get(ds, ds)
        print(f"  - {tag} ← {path.relative_to(JAMB)}")

    edge_tables: dict[str, pd.DataFrame] = {}
    rows: list[dict] = []
    for ds, path in found:
        tag = SHORT_TAGS.get(ds, ds)
        df = pd.read_csv(path, sep="\t")
        edge_tables[tag] = df
        stats = per_lineage_stats(df)
        rows.append({"dataset_tag": tag, "dataset_id": ds, **stats})

    out = pd.DataFrame(rows)

    # Pair-wise overlap columns (auto-widen). Symmetric: filled both directions.
    tags = [r["dataset_tag"] for r in rows]
    edge_sets = {
        t: set(zip(edge_tables[t]["intended_target_name"], edge_tables[t]["gene_id"]))
        for t in tags
    }
    tf_sets = {t: set(edge_tables[t]["intended_target_name"].unique()) for t in tags}

    for other in tags:
        out[f"n_shared_edges_with_{other}"] = out["dataset_tag"].map(
            lambda t: len(edge_sets[t] & edge_sets[other]) if t != other else np.nan
        )
        out[f"n_shared_tfs_with_{other}"] = out["dataset_tag"].map(
            lambda t: len(tf_sets[t] & tf_sets[other]) if t != other else np.nan
        )
        out[f"jaccard_edges_with_{other}"] = out["dataset_tag"].map(
            lambda t: (
                len(edge_sets[t] & edge_sets[other]) / len(edge_sets[t] | edge_sets[other])
                if t != other and (edge_sets[t] | edge_sets[other])
                else np.nan
            )
        )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT, sep="\t", index=False)
    print(f"\nwrote {OUTPUT} ({len(out)} rows × {len(out.columns)} cols)")
    print(out.to_string(index=False, max_colwidth=20))


if __name__ == "__main__":
    main()
