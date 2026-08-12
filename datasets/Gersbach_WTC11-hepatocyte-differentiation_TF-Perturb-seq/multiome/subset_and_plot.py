# -*- coding: utf-8 -*-
"""
Subsets an h5ad file to cells in a MULTI-seq barcode mapping and plots
violin plots of hepatocyte marker gene expression.

Usage:
    python subset_and_plot.py \
        --h5ad    IGVFFI2923WJHO.h5ad \
        --mapping cell_barcode_mapping.tsv \
        --output  IGVFFI2923WJHO_subset.h5ad \
        --plot    hepatocyte_markers.png
"""

import argparse
import sys
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_MARKER_GENES = ["AFP", "ALB", "SERPINA1", "APOE", "TTR", "FGB"]

# Ensembl ID -> symbol for display labels
ENSG_TO_SYMBOL = {
    "ENSG00000081051": "AFP",
    "ENSG00000163631": "ALB",
    "ENSG00000197249": "SERPINA1",
    "ENSG00000130203": "APOE",
    "ENSG00000118271": "TTR",
    "ENSG00000171564": "FGB",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad",    required=True, help="Input h5ad file")
    parser.add_argument("--mapping", required=True, help="cell_barcode_mapping.tsv from assign_multiseq_barcodes.py")
    parser.add_argument("--output",  default="IGVFFI2923WJHO_subset.h5ad", help="Output h5ad path")
    parser.add_argument("--plot",    default="hepatocyte_markers.png", help="Output violin plot path")
    parser.add_argument("--skip-norm", action="store_true",
                        help="Skip normalization (use if h5ad is already normalized)")
    parser.add_argument("--skip-subset", action="store_true",
                        help="Skip subsetting and rewriting the h5ad (just plot)")
    parser.add_argument("--min-genes", type=int, default=None,
                        help="For plotting only: minimum genes detected per cell "
                             "(e.g. 300). Cells below this threshold are excluded "
                             "from the violin plots but NOT from the saved h5ad.")
    parser.add_argument("--genes",   nargs="+", default=None,
                        help="Gene IDs to plot (default: symbol names). "
                             "Pass Ensembl IDs if the matrix uses ENSG IDs.")
    return parser.parse_args()


def main():
    args = parse_args()

    # 1. Load mapping and get valid h5ad barcodes
    print(f"[1/5] Loading mapping: {args.mapping}")
    mapping = pd.read_csv(args.mapping, sep="\t")
    # Drop cells without a translated GEX barcode or that are Doublet/Unassigned
    mapping = mapping.dropna(subset=["h5ad_barcode"])
    mapping = mapping[~mapping["multiseq_barcode"].isin(["Doublet", "Unassigned"])]
    valid_barcodes = set(mapping["h5ad_barcode"])
    print(f"      {len(valid_barcodes):,} valid barcodes in mapping")

    # 2. Load h5ad
    print(f"\n[2/5] Loading h5ad: {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    print(f"      Full matrix: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    if args.skip_subset:
        print(f"\n[3/5] Skipping subset and rewrite (--skip-subset set)")
        print(f"[4/5] Skipping normalization and save (--skip-subset set)")
        print(f"[5a/5] Skipping timepoint split (--skip-subset set)")
    else:
        # 3. Subset to cells in mapping
        print(f"\n[3/5] Subsetting to mapped cells")
        keep = adata.obs_names.isin(valid_barcodes)
        adata = adata[keep].copy()
        print(f"      Subset matrix: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

        n_missing = len(valid_barcodes) - adata.n_obs
        if n_missing > 0:
            print(f"      Note: {n_missing:,} barcodes in mapping not found in h5ad "
                  f"(expected if mapping includes unfiltered barcodes)", file=sys.stderr)

        # Attach sample metadata to obs
        mapping_indexed = mapping.set_index("h5ad_barcode")[
            ["multiseq_barcode", "sample_accession", "sample_description"]
        ]
        cols_to_add = mapping_indexed.columns.tolist()
        adata.obs = adata.obs.drop(columns=[c for c in cols_to_add if c in adata.obs.columns])
        adata.obs = adata.obs.join(mapping_indexed, how="left")

        # 4. Normalize
        if args.skip_norm:
            print(f"\n[4/5] Skipping normalization (--skip-norm set)")
        else:
            print(f"\n[4/5] Normalizing (library size + log1p)")
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)

        # 5. Save subset h5ad
        print(f"      Saving subset h5ad: {args.output}")
        adata.write_h5ad(args.output)

        # 6. Split by timepoint and save per-timepoint h5ad files
        print(f"\n[5a/5] Splitting by timepoint (sample_description)")
        if "sample_description" in adata.obs.columns:
            adata.obs["timepoint"] = (
                adata.obs["sample_description"]
                .str.extract(r"(t\d+)", expand=False)
            )
            for tp, group in adata.obs.groupby("timepoint"):
                tp_adata = adata[group.index].copy()
                tp_path  = args.output.replace(".h5ad", f"_{tp}.h5ad")
                tp_adata.write_h5ad(tp_path)
                print(f"      {tp}: {tp_adata.n_obs:,} cells -> {tp_path}")
        else:
            print("      WARNING: sample_description not in obs, skipping timepoint split.",
                  file=sys.stderr)

    # 7. Plot violin plots — optionally filter to high-quality cells for plotting only
    print(f"\n[5b/5] Plotting hepatocyte markers: {args.plot}")
    adata_plot = adata
    if args.min_genes is not None:
        sc.pp.calculate_qc_metrics(adata_plot, inplace=True)
        before = adata_plot.n_obs
        adata_plot = adata_plot[adata_plot.obs["n_genes_by_counts"] >= args.min_genes].copy()
        after = adata_plot.n_obs
        print(f"      --min-genes {args.min_genes}: {before:,} -> {after:,} cells "
              f"({before - after:,} removed, {after/before*100:.1f}% retained for plotting)")
    marker_genes = args.genes if args.genes else DEFAULT_MARKER_GENES

    # Build unversioned -> versioned lookup (e.g. ENSG00000081051 -> ENSG00000081051.2)
    unversioned_to_full = {v.split(".")[0]: v for v in adata_plot.var_names}

    # Resolve each requested ID to the full versioned var_name
    resolved = {g: unversioned_to_full.get(g, g) for g in marker_genes}
    present  = [resolved[g] for g in marker_genes if resolved[g] in adata_plot.var_names]
    missing  = [g for g in marker_genes if resolved[g] not in adata_plot.var_names]

    if missing:
        print(f"      WARNING: genes not found in matrix and will be skipped: {missing}",
              file=sys.stderr)
    if not present:
        print("      ERROR: none of the marker genes found in matrix.", file=sys.stderr)
        sys.exit(1)

    # Use sample_description as groupby if available, else plot all cells together
    groupby = None
    if "sample_description" in adata_plot.obs.columns and adata_plot.obs["sample_description"].notna().any():
        groupby = "sample_description"
        print(f"      Grouping by: {groupby}")

    n_genes = len(present)
    fig, axes = plt.subplots(n_genes, 1, figsize=(10, 3 * n_genes))
    if n_genes == 1:
        axes = [axes]

    y_label = "log1p(CPM)" if not args.skip_norm else "log-normalized expression"

    for ax, gene in zip(axes, present):
        sc.pl.violin(
            adata_plot,
            keys=gene,
            groupby=groupby,
            ax=ax,
            show=False,
            stripplot=False,
            inner="box",
        )
        # Remove legend — x-axis tick labels already show the group names
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()
        # Title: "AFP (ENSG00000081051)"
        unversioned = gene.split(".")[0]
        symbol = ENSG_TO_SYMBOL.get(unversioned, None)
        title = f"{symbol} ({unversioned})" if symbol else gene
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel(y_label, fontsize=10)
        if groupby:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)

    fig.suptitle("Hepatocyte marker gene expression",
                 fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig(args.plot, dpi=150, bbox_inches="tight")
    print(f"      Saved: {args.plot}")
    print("\nDone.")


if __name__ == "__main__":
    main()