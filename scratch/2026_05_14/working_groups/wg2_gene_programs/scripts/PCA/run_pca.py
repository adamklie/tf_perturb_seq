"""Standalone PCA visualization for the meta-consensus pipeline.

Reads:
  - the meta_consensus config.yaml (for input spectra paths + normalization),
  - per-k cluster_assignments.k_{k}.txt under Result/meta_consensus/k_{k}/.

Re-builds the same normalized program matrix N as the main pipeline, fits
PCA(n_components=2) on it once, and writes:
  - Result/PCA/pca_programs_by_dataset.png   (global, k-independent)
  - Result/PCA/k_{k}/pca_programs_by_cluster.png  (per k)
"""

from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
META_DIR = os.path.normpath(os.path.join(HERE, "..", "meta_consensus"))
sys.path.insert(0, HERE)
sys.path.insert(0, META_DIR)

import yaml  # noqa: E402

from meta_consensus import (  # noqa: E402
    DatasetSpec,
    intersect_genes,
    load_spectra,
    normalize_rows,
    stack_programs,
)
from pca_plots import (  # noqa: E402
    fit_pca,
    plot_pca_programs_by_cluster,
    plot_pca_programs_by_dataset,
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--config",
        default=os.path.join(META_DIR, "config.yaml"),
        help="Path to meta_consensus config.yaml",
    )
    ap.add_argument(
        "--pca-output-dir",
        default="/oak/stanford/groups/engreitz/Users/ymo/Testing_Scripts/Meta_program/Result/PCA",
        help="Where to write PCA plots",
    )
    ap.add_argument(
        "--annotate-clusters",
        action="store_true",
        help="Overlay meta-cluster IDs at each cluster centroid (default: off)",
    )
    args = ap.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    datasets = [DatasetSpec(**d) for d in cfg["datasets"]]
    input_file_type = cfg.get("input_file_type", "gene_spectra_score")
    normalization = cfg.get("normalization", "l2")
    meta_k_list = cfg["meta_k_list"]
    meta_out_dir = cfg["output_dir"]

    os.makedirs(args.pca_output_dir, exist_ok=True)

    print("# Loading and stacking input programs")
    dfs = []
    for spec in datasets:
        df = load_spectra(spec, input_file_type)
        print(f"  {spec.name}: shape={df.shape}  k={spec.k} dt={spec.dt}")
        dfs.append(df)
    common = intersect_genes(dfs)
    S = stack_programs(dfs, common)
    N = normalize_rows(S, normalization)
    print(f"  Stacked + normalized: {N.shape}  (normalization={normalization})")

    print("\n# Fitting PCA on input programs")
    scores, pca = fit_pca(N, n_components=2)
    ev = pca.explained_variance_ratio_
    print(f"  PC1 var: {ev[0]*100:.1f}%   PC2 var: {ev[1]*100:.1f}%")

    global_out = os.path.join(args.pca_output_dir, "pca_programs_by_dataset.png")
    plot_pca_programs_by_dataset(
        scores, pca, N.index, global_out,
        title=f"PCA of input programs ({N.shape[0]} programs, {N.shape[1]} genes)",
    )
    print(f"  wrote {global_out}")

    print("\n# Per-k cluster overlays")
    for k in meta_k_list:
        labels_path = os.path.join(
            meta_out_dir, f"k_{k}", f"cluster_assignments.k_{k}.txt"
        )
        if not os.path.exists(labels_path):
            print(f"  k={k}: skipped (missing {labels_path})")
            continue
        labels_df = pd.read_csv(labels_path, sep="\t", index_col=0)
        labels = labels_df.iloc[:, 0]
        labels.name = "meta_cluster_id"

        k_out_dir = os.path.join(args.pca_output_dir, f"k_{k}")
        os.makedirs(k_out_dir, exist_ok=True)
        out_png = os.path.join(k_out_dir, "pca_programs_by_cluster.png")
        plot_pca_programs_by_cluster(
            scores, pca, N.index, labels, out_png,
            title=f"PCA of input programs colored by dataset (k={k})",
            annotate_centroids=args.annotate_clusters,
        )
        print(f"  k={k}: wrote {out_png}")

    print(f"\n# Done. Plots in {args.pca_output_dir}")


if __name__ == "__main__":
    main()
