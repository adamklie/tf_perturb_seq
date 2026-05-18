"""CLI entrypoint for meta-consensus across cNMF runs on separate datasets.

Reads config.yaml, loads per-dataset gene_spectra files, finds common genes,
row-normalizes, runs KMeans for each k in meta_k_list, writes per-k outputs
plus an across-k silhouette curve.

Usage:
    python run_meta_consensus.py --config config.yaml
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List

import pandas as pd
import yaml

# Allow running as a script from any cwd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from meta_consensus import (  # noqa: E402
    DatasetSpec,
    build_shared_vs_specific,
    categorize_clusters,
    intersect_genes,
    load_spectra,
    normalize_rows,
    run_kmeans_consensus,
    stack_programs,
    write_category_counts,
)
from plotting import (  # noqa: E402
    plot_category_breakdown,
    plot_clustergram,
    plot_composition_bar,
    plot_k_selection,
    plot_top_genes_heatmap,
)
import plot_cluster_zoom  # noqa: E402


def _is_ensembl(g: str) -> bool:
    return isinstance(g, str) and g.startswith("ENS")


def _log_gene_overlap(dfs: List[pd.DataFrame], dataset_names: List[str], common: pd.Index) -> None:
    print("\n# Gene overlap report")
    for name, df in zip(dataset_names, dfs):
        n_ens = sum(_is_ensembl(g) for g in df.columns)
        n_sym = len(df.columns) - n_ens
        print(f"  {name}: {len(df.columns)} genes  (Ensembl={n_ens}, symbol-like={n_sym})")
    n_ens_c = sum(_is_ensembl(g) for g in common)
    n_sym_c = len(common) - n_ens_c
    print(f"  COMMON: {len(common)} genes  (Ensembl={n_ens_c}, symbol-like={n_sym_c})")
    if len(common) < 1000:
        print(f"  WARNING: small overlap — check for Ensembl/symbol mismatch across datasets.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="Path to config.yaml")
    args = ap.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    datasets = [DatasetSpec(**d) for d in cfg["datasets"]]
    dataset_names = [d.name for d in datasets]
    input_file_type = cfg.get("input_file_type", "gene_spectra_score")
    normalization = cfg.get("normalization", "l2")
    meta_k_list = cfg["meta_k_list"]
    out_root = cfg["output_dir"]
    kmeans_cfg = cfg.get("kmeans", {})
    n_init = int(kmeans_cfg.get("n_init", 10))
    random_state = int(kmeans_cfg.get("random_state", 1))

    os.makedirs(out_root, exist_ok=True)

    print("# Loading spectra")
    dfs = []
    for spec in datasets:
        df = load_spectra(spec, input_file_type)
        print(f"  {spec.name}: shape={df.shape}  from k={spec.k} dt={spec.dt}")
        dfs.append(df)

    common = intersect_genes(dfs)
    _log_gene_overlap(dfs, dataset_names, common)

    S = stack_programs(dfs, common)
    print(f"\n# Stacked: {S.shape[0]} programs x {S.shape[1]} common genes")

    N = normalize_rows(S, normalization)
    print(f"# Normalization: {normalization}")

    silhouettes: Dict[int, float] = {}

    for k in meta_k_list:
        print(f"\n# Meta-consensus k={k}")
        k_dir = os.path.join(out_root, f"k_{k}")
        os.makedirs(k_dir, exist_ok=True)

        labels, meta_spectra, sil = run_kmeans_consensus(
            N, k, n_init=n_init, random_state=random_state
        )
        silhouettes[k] = sil
        print(f"  silhouette = {sil:.4f}")

        cluster_summary = categorize_clusters(labels, dataset_names)
        shared_vs_specific = build_shared_vs_specific(labels, cluster_summary)

        meta_spectra.to_csv(
            os.path.join(k_dir, f"meta_spectra.k_{k}.median.txt"), sep="\t"
        )
        labels.to_csv(
            os.path.join(k_dir, f"cluster_assignments.k_{k}.txt"), sep="\t", header=True
        )
        with open(os.path.join(k_dir, f"silhouette.k_{k}.txt"), "w") as f:
            f.write(f"{sil:.6f}\n")
        cluster_summary.to_csv(
            os.path.join(k_dir, f"cluster_summary.k_{k}.tsv"), sep="\t", index=False
        )
        shared_vs_specific.to_csv(
            os.path.join(k_dir, f"shared_vs_specific.k_{k}.tsv"), sep="\t", index=False
        )
        write_category_counts(
            cluster_summary, dataset_names,
            os.path.join(k_dir, f"category_counts.k_{k}.txt"),
        )

        plot_clustergram(
            N, labels,
            os.path.join(k_dir, f"clustering.k_{k}.png"),
            title=f"Meta-consensus clustergram (k={k}, sil={sil:.3f})",
        )
        plot_composition_bar(
            cluster_summary, dataset_names,
            os.path.join(k_dir, f"composition.k_{k}.png"),
            title=f"Dataset composition per meta-cluster (k={k})",
        )
        plot_category_breakdown(
            cluster_summary, dataset_names,
            os.path.join(k_dir, f"category_breakdown.k_{k}.png"),
            title=f"Meta-cluster categories (k={k})",
        )
        plot_top_genes_heatmap(
            meta_spectra,
            os.path.join(k_dir, f"top_genes.k_{k}.png"),
            n_top=20,
            title=f"Top genes per meta-program (k={k})",
        )

        zoom_clusters = cfg.get("zoom_clusters") or []
        if zoom_clusters:
            go_dir = os.path.join(os.path.dirname(out_root.rstrip("/")), "Shared_GO")
            go_summary_path = os.path.join(go_dir, f"Shared_GO_summary.k_{k}.tsv")
            go_long_path = os.path.join(go_dir, f"Shared_GO_long.k_{k}.tsv")
            if not (os.path.exists(go_summary_path) and os.path.exists(go_long_path)):
                print(f"  [zoom] skipping: Shared_GO outputs for k={k} not found "
                      f"(run Shared_GO/shared_GO_analysis.py first)")
            else:
                go_summary_df = pd.read_csv(go_summary_path, sep="\t")
                go_long_df = pd.read_csv(go_long_path, sep="\t")
                present = set(labels.values)
                for cid in zoom_clusters:
                    cid = int(cid)
                    if cid not in present:
                        print(f"  [zoom] cluster {cid} absent in k={k}, skipping")
                        continue
                    out_png = os.path.join(k_dir, f"cluster_{cid}_zoom.png")
                    print(f"  [zoom] cluster {cid} -> {out_png}")
                    plot_cluster_zoom.plot(
                        N=N,
                        labels=labels,
                        cluster_id=cid,
                        cluster_summary=cluster_summary,
                        go_summary=go_summary_df,
                        go_long=go_long_df,
                        out_path=out_png,
                        title=f"Meta-consensus k={k} — cluster #{cid} zoom",
                    )

    # k-selection curve across all k
    pd.Series(silhouettes).sort_index().to_csv(
        os.path.join(out_root, "k_selection.txt"), sep="\t", header=False
    )
    plot_k_selection(silhouettes, os.path.join(out_root, "k_selection.png"))
    print(f"\n# Done. Outputs in {out_root}")


if __name__ == "__main__":
    main()
