"""
Plot sequencing coverage diagnostics for a subset h5ad file.

Produces a 4-panel figure:
  1. Barcode rank plot (knee plot) - UMI counts ranked descending
  2. Histogram of UMI counts per cell
  3. Histogram of genes detected per cell
  4. Scatter: UMI counts vs genes detected per cell

Usage:
    python plot_coverage.py \
        --h5ad   IGVFFI2923WJHO_subset.h5ad \
        --output coverage_plots.png
"""

import argparse
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad",     required=True,
                        help="Subset h5ad (may be log-normalized; used for cell list and obs metadata)")
    parser.add_argument("--raw-h5ad", default=None,
                        help="Original unfiltered h5ad with raw counts (e.g. IGVFFI2923WJHO.h5ad). "
                             "If provided, UMI/gene counts are pulled from here instead of --h5ad.")
    parser.add_argument("--output",  default="coverage_plots.png")
    parser.add_argument("--groupby", default="sample_description",
                        help="obs column to colour scatter by (default: sample_description)")
    return parser.parse_args()


def main():
    args = parse_args()

    print(f"Loading subset h5ad: {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    if args.raw_h5ad:
        print(f"\nNOTE: --h5ad is assumed to be log-normalized. Loading raw counts from")
        print(f"      {args.raw_h5ad} and subsetting to the {adata.n_obs:,} cells in the subset h5ad.")
        print(f"      QC metrics (UMI counts, genes detected) will reflect raw counts.")
        raw = sc.read_h5ad(args.raw_h5ad)
        keep = raw.obs_names.isin(adata.obs_names)
        raw  = raw[keep].copy()
        # Copy obs metadata (sample_description etc.) from annotated subset onto raw
        raw.obs = raw.obs.join(
            adata.obs[["sample_description", "sample_accession", "multiseq_barcode"]
                      if "sample_description" in adata.obs.columns else []],
            how="left"
        )
        sc.pp.calculate_qc_metrics(raw, inplace=True)
        adata = raw
    else:
        print(f"\nNOTE: No --raw-h5ad provided. Assuming --h5ad contains log-normalized counts.")
        print(f"      UMI counts and gene detection numbers will reflect normalized values,")
        print(f"      not raw counts. Pass --raw-h5ad IGVFFI2923WJHO.h5ad for accurate metrics.")
        sc.pp.calculate_qc_metrics(adata, inplace=True)
    umi    = adata.obs["total_counts"].values
    ngenes = adata.obs["n_genes_by_counts"].values

    # Colour by groupby column if present
    groupby = args.groupby if args.groupby in adata.obs.columns else None
    if groupby:
        groups  = adata.obs[groupby].astype(str)
        labels  = sorted(groups.unique())
        cmap    = plt.get_cmap("tab20", len(labels))
        color_map = {l: cmap(i) for i, l in enumerate(labels)}
        colors  = groups.map(color_map).values
    else:
        colors  = ["steelblue"] * len(umi)
        labels  = None

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Sequencing coverage - subset cells", fontsize=15, fontweight="bold")

    # ------------------------------------------------------------------
    # 1. Barcode rank (knee) plot
    # ------------------------------------------------------------------
    ax = axes[0, 0]
    ranked = np.sort(umi)[::-1]
    ax.plot(np.arange(1, len(ranked) + 1), ranked, color="steelblue", linewidth=1.2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Barcode rank", fontsize=11)
    ax.set_ylabel("UMI counts", fontsize=11)
    ax.set_title("Barcode rank plot (knee plot)", fontsize=12, fontweight="bold")
    ax.axhline(np.median(umi), color="tomato", linestyle="--", linewidth=1,
               label=f"Median: {np.median(umi):,.0f}")
    ax.legend(fontsize=9)

    # ------------------------------------------------------------------
    # 2. Histogram - UMI counts per cell
    # ------------------------------------------------------------------
    ax = axes[0, 1]
    ax.hist(umi, bins=80, color="steelblue", edgecolor="white", linewidth=0.3)
    ax.axvline(np.median(umi), color="tomato", linestyle="--", linewidth=1.2,
               label=f"Median: {np.median(umi):,.0f}")
    ax.set_xlabel("UMI counts per cell", fontsize=11)
    ax.set_ylabel("Number of cells", fontsize=11)
    ax.set_title("UMI counts per cell", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)

    # ------------------------------------------------------------------
    # 3. Histogram - genes detected per cell
    # ------------------------------------------------------------------
    ax = axes[1, 0]
    ax.hist(ngenes, bins=80, color="mediumseagreen", edgecolor="white", linewidth=0.3)
    ax.axvline(np.median(ngenes), color="tomato", linestyle="--", linewidth=1.2,
               label=f"Median: {np.median(ngenes):,.0f}")
    ax.set_xlabel("Genes detected per cell", fontsize=11)
    ax.set_ylabel("Number of cells", fontsize=11)
    ax.set_title("Genes detected per cell", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)

    # ------------------------------------------------------------------
    # 4. Scatter - UMI counts vs genes detected
    # ------------------------------------------------------------------
    ax = axes[1, 1]
    if groupby and labels:
        for label in labels:
            mask = groups == label
            ax.scatter(umi[mask], ngenes[mask], s=2, alpha=0.4,
                       color=color_map[label], label=label, rasterized=True)
        ax.legend(title=groupby, fontsize=7, title_fontsize=8,
                  markerscale=4, bbox_to_anchor=(1.01, 1), loc="upper left")
    else:
        ax.scatter(umi, ngenes, s=2, alpha=0.3, color="steelblue", rasterized=True)
    ax.set_xlabel("UMI counts per cell", fontsize=11)
    ax.set_ylabel("Genes detected per cell", fontsize=11)
    ax.set_title("UMI counts vs genes detected", fontsize=12, fontweight="bold")

    plt.tight_layout()
    plt.savefig(args.output, dpi=150, bbox_inches="tight")
    print(f"Saved: {args.output}")

    # Print overall summary stats
    print(f"\n  UMI counts  - median: {np.median(umi):,.0f}  "
          f"mean: {np.mean(umi):,.0f}  "
          f"min: {umi.min():,.0f}  max: {umi.max():,.0f}")
    print(f"  Genes/cell  - median: {np.median(ngenes):,.0f}  "
          f"mean: {np.mean(ngenes):,.0f}  "
          f"min: {ngenes.min():,.0f}  max: {ngenes.max():,.0f}")

    # Per-replicate summary (reads/UMIs per timepoint and replicate)
    if "sample_description" in adata.obs.columns:
        print("\n  Per-replicate summary (UMI counts):")
        print(f"  {'Sample':<25} {'N cells':>8} {'Total UMIs':>12} "
              f"{'Median UMI/cell':>16} {'Median genes/cell':>18}")
        print("  " + "-" * 82)

        obs = adata.obs[["sample_description", "total_counts", "n_genes_by_counts"]].copy()
        for sample, grp in obs.groupby("sample_description", sort=True):
            n_cells    = len(grp)
            total_umi  = grp["total_counts"].sum()
            med_umi    = grp["total_counts"].median()
            med_genes  = grp["n_genes_by_counts"].median()
            print(f"  {sample:<25} {n_cells:>8,} {total_umi:>12,.0f} "
                  f"{med_umi:>16,.0f} {med_genes:>18,.0f}")

        # Also write to a TSV
        tsv_path = args.output.replace(".png", "_per_replicate.tsv")
        summary = (
            obs.groupby("sample_description")
            .agg(
                n_cells=("total_counts", "count"),
                total_umis=("total_counts", "sum"),
                mean_umis_per_cell=("total_counts", "mean"),
                median_umis_per_cell=("total_counts", "median"),
                median_genes_per_cell=("n_genes_by_counts", "median"),
            )
            .reset_index()
        )
        summary.to_csv(tsv_path, sep="\t", index=False, float_format="%.1f")
        print(f"\n  Per-replicate table saved: {tsv_path}")


if __name__ == "__main__":
    main()