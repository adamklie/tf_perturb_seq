#!/usr/bin/env python

"""
Main QC script for TF Perturb-seq MuData outputs.

- Loads a MuData .h5mu file (expects 'gene' and 'guide' modalities)
- Applies joint cell-level filtering (gene + guide QC)
- Computes & saves pre-filter and post-filter QC plots
  (gene + guide, including knee plots)
- Writes a SINGLE TSV with metrics for pre- and post-filter
- Optionally writes filtered gene and guide .h5ad files
"""

import argparse
import logging
import os
import sys
from typing import Optional

import matplotlib
matplotlib.use("Agg")  # ensure non-interactive backend

import mudata
from anndata import AnnData

import pandas as pd  # noqa: E402

from utils import (  # noqa: E402
    describe_mudata,
    add_sample_prefix_to_barcodes,
    plot_knee_plot,
    plot_gene_qc_histograms,
    plot_cells_per_batch,
    prepare_guide_modality,
    plot_guide_qc,
    apply_filters,
    compute_metrics,
)

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)


# ---------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="QC for TF Perturb-seq MuData outputs")

    parser.add_argument(
        "--input-h5mu",
        required=True,
        help="Path to input MuData .h5mu file",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory for QC outputs",
    )
    parser.add_argument(
        "--sample-name",
        default=None,
        help="Sample name to prefix barcodes with (also stored in metrics TSV)",
    )
    parser.add_argument(
        "--gene-mod-key",
        default="gene",
        help="Modality key for gene expression in MuData (default: 'gene')",
    )
    parser.add_argument(
        "--guide-mod-key",
        default="guide",
        help="Modality key for guide assignments in MuData (default: 'guide')",
    )
    parser.add_argument(
        "--guide-assignment-layer",
        default="guide_assignment",
        help="Layer name in guide modality with guide assignment matrix (>0 = assigned)",
    )

    # Column names (for flexibility)
    parser.add_argument("--umi-col", default="total_gene_umis")
    parser.add_argument("--genes-col", default="num_expressed_genes")
    parser.add_argument("--mito-col", default="percent_mito")
    parser.add_argument("--batch-col", default="batch")

    parser.add_argument(
        "--guide-label-col",
        default="label",
        help="Column in guide.var with guide labels (e.g. target/NTC/pos/neg)",
    )
    parser.add_argument(
        "--guide-umi-col",
        default="total_guide_umis",
        help="Column in guide.obs with total guide UMIs per cell",
    )

    # Gene-level filtering thresholds
    parser.add_argument("--min-total-umis", type=float, default=None)
    parser.add_argument("--max-total-umis", type=float, default=None)
    parser.add_argument("--min-genes", type=float, default=None)
    parser.add_argument("--max-genes", type=float, default=None)
    parser.add_argument("--max-mito-pct", type=float, default=None)

    # Guide-related filtering (per-cell)
    parser.add_argument(
        "--guides-per-cell-col",
        default="n_guides_per_cell",
        help="Column in guide.obs holding guides per cell (computed from guide_assignment).",
    )
    parser.add_argument("--min-guides-per-cell", type=float, default=None)
    parser.add_argument("--max-guides-per-cell", type=float, default=None)

    parser.add_argument(
        "--min-cells-per-guide",
        type=int,
        default=None,
        help="Reserved for future filtering based on n_cells_per_guide; currently ignored.",
    )

    # Plotting options
    parser.add_argument(
        "--plot-gene-by-batch",
        action="store_true",
        help="Plot gene QC histograms split by batch",
    )

    # Output options
    parser.add_argument(
        "--write-filtered-h5ad",
        action="store_true",
        help="Write filtered gene and guide AnnData objects to disk",
    )
    parser.add_argument(
        "--gene-h5ad-name",
        default="gene.filtered.h5ad",
        help="Filename for filtered gene AnnData (inside outdir)",
    )
    parser.add_argument(
        "--guide-h5ad-name",
        default="guide.filtered.h5ad",
        help="Filename for filtered guide AnnData (inside outdir)",
    )

    # Metrics TSV name
    parser.add_argument(
        "--metrics-tsv",
        default="qc_metrics.tsv",
        help="Filename for combined metrics TSV (written inside outdir)",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def extract_modalities(
    data: mudata.MuData,
    gene_mod_key: str,
    guide_mod_key: str,
) -> tuple[AnnData, AnnData]:
    if gene_mod_key not in data.mod:
        raise ValueError(f"Gene modality '{gene_mod_key}' not found in MuData.mod")
    if guide_mod_key not in data.mod:
        raise ValueError(f"Guide modality '{guide_mod_key}' not found in MuData.mod")
    return data.mod[gene_mod_key], data.mod[guide_mod_key]


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # --- Load MuData ---
    logger.info(f"Loading MuData from {args.input_h5mu}")
    data = mudata.read_h5mu(args.input_h5mu)
    describe_mudata(data)

    # Optional: prefix barcodes once on the MuData (propagates to all modalities)
    add_sample_prefix_to_barcodes(data, args.sample_name)

    # Extract modalities
    gene, guide_raw = extract_modalities(data, args.gene_mod_key, args.guide_mod_key)

    # Prepare guide modality (compute n_guides_per_cell & n_cells_per_guide)
    guide = prepare_guide_modality(
        data,
        guide_mod_key=args.guide_mod_key,
        guide_assignment_layer=args.guide_assignment_layer,
    )

    # Sanity check: obs_names alignment
    if not gene.obs_names.equals(guide.obs_names):
        raise ValueError("gene and guide modalities have different obs_names after prefixing.")

    # -----------------------------------------------------------------
    # PRE-FILTER QC
    # -----------------------------------------------------------------
    logger.info("Running pre-filter QC")

    pre_prefix = "pre_filter"

    # Gene: knee plot (from X)
    plot_knee_plot(
        gene,
        outdir=args.outdir,
        prefix=pre_prefix + "_gene",
        use_X=True,
        count_col=None,
        title="Gene knee plot (UMIs per cell)",
    )

    # Gene QC histograms
    plot_gene_qc_histograms(
        gene,
        outdir=args.outdir,
        prefix=pre_prefix,
        umi_col=args.umi_col,
        genes_col=args.genes_col,
        mito_col=args.mito_col,
        by_batch=args.plot_gene_by_batch,
        batch_col=args.batch_col,
    )

    # Cells per batch (same for gene/guide)
    plot_cells_per_batch(
        gene,
        outdir=args.outdir,
        prefix=pre_prefix,
        batch_col=args.batch_col,
    )

    # Guide: knee plot (from per-cell total_guide_umis)
    plot_knee_plot(
        guide,
        outdir=args.outdir,
        prefix=pre_prefix + "_guide",
        use_X=False,
        count_col=args.guide_umi_col,
        title="Guide knee plot (guide UMIs per cell)",
    )

    # Guide QC plots
    plot_guide_qc(
        guide,
        outdir=args.outdir,
        prefix=pre_prefix,
        batch_col=args.batch_col,
        label_col=args.guide_label_col,
        guide_umi_col=args.guide_umi_col,
        guides_per_cell_col=args.guides_per_cell_col,
    )

    # Pre-filter metrics (joint)
    metrics_pre_list = []

    # Global (cross-batch) metrics: batch = "all"
    metrics_pre_all = compute_metrics(
        gene=gene,
        guide=guide,
        filter_stage="pre",
        sample_name=args.sample_name,
        umi_col=args.umi_col,
        genes_col=args.genes_col,
        mito_col=args.mito_col,
        guide_umi_col=args.guide_umi_col,
        guides_per_cell_col=args.guides_per_cell_col,
        batch_label="all",
    )
    metrics_pre_list.append(metrics_pre_all)

    # Per-batch metrics (if batch_col exists)
    if args.batch_col in gene.obs.columns:
        for b in sorted(gene.obs[args.batch_col].unique()):
            mask_b = gene.obs[args.batch_col] == b
            gene_b = gene[mask_b].copy()
            guide_b = guide[mask_b].copy()

            metrics_pre_b = compute_metrics(
                gene=gene_b,
                guide=guide_b,
                filter_stage="pre",
                sample_name=args.sample_name,
                umi_col=args.umi_col,
                genes_col=args.genes_col,
                mito_col=args.mito_col,
                guide_umi_col=args.guide_umi_col,
                guides_per_cell_col=args.guides_per_cell_col,
                batch_label=str(b),
            )
            metrics_pre_list.append(metrics_pre_b)

    metrics_pre = pd.concat(metrics_pre_list, ignore_index=True)

    # -----------------------------------------------------------------
    # FILTERING (joint gene + guide)
    # -----------------------------------------------------------------
    logger.info("Applying joint filters (gene + guide)")

    mask = apply_filters(
        gene=gene,
        guide=guide,
        umi_col=args.umi_col,
        genes_col=args.genes_col,
        mito_col=args.mito_col,
        guides_per_cell_col=args.guides_per_cell_col,
        min_total_umis=args.min_total_umis,
        max_total_umis=args.max_total_umis,
        min_genes=args.min_genes,
        max_genes=args.max_genes,
        max_mito_pct=args.max_mito_pct,
        min_guides_per_cell=args.min_guides_per_cell,
        max_guides_per_cell=args.max_guides_per_cell,
        min_cells_per_guide=args.min_cells_per_guide,
    )

    gene_filt = gene[mask].copy()
    guide_filt = guide[mask].copy()

    # -----------------------------------------------------------------
    # POST-FILTER QC
    # -----------------------------------------------------------------
    logger.info("Running post-filter QC")

    post_prefix = "post_filter"

    # Gene: knee plot (post-filter)
    plot_knee_plot(
        gene_filt,
        outdir=args.outdir,
        prefix=post_prefix + "_gene",
        use_X=True,
        count_col=None,
        title="Gene knee plot (UMIs per cell, filtered)",
    )

    # Gene QC histograms
    plot_gene_qc_histograms(
        gene_filt,
        outdir=args.outdir,
        prefix=post_prefix,
        umi_col=args.umi_col,
        genes_col=args.genes_col,
        mito_col=args.mito_col,
        by_batch=args.plot_gene_by_batch,
        batch_col=args.batch_col,
    )

    # Cells per batch
    plot_cells_per_batch(
        gene_filt,
        outdir=args.outdir,
        prefix=post_prefix,
        batch_col=args.batch_col,
    )

    # Guide: knee plot (post-filter)
    plot_knee_plot(
        guide_filt,
        outdir=args.outdir,
        prefix=post_prefix + "_guide",
        use_X=False,
        count_col=args.guide_umi_col,
        title="Guide knee plot (guide UMIs per cell, filtered)",
    )

    # Guide QC plots
    plot_guide_qc(
        guide_filt,
        outdir=args.outdir,
        prefix=post_prefix,
        batch_col=args.batch_col,
        label_col=args.guide_label_col,
        guide_umi_col=args.guide_umi_col,
        guides_per_cell_col=args.guides_per_cell_col,
    )

    # Post-filter metrics
    metrics_post_list = []

    # Global (cross-batch) metrics: batch = "all"
    metrics_post_all = compute_metrics(
        gene=gene_filt,
        guide=guide_filt,
        filter_stage="post",
        sample_name=args.sample_name,
        umi_col=args.umi_col,
        genes_col=args.genes_col,
        mito_col=args.mito_col,
        guide_umi_col=args.guide_umi_col,
        guides_per_cell_col=args.guides_per_cell_col,
        batch_label="all",
    )
    metrics_post_list.append(metrics_post_all)

    # Per-batch metrics (if batch_col exists)
    if args.batch_col in gene_filt.obs.columns:
        for b in sorted(gene_filt.obs[args.batch_col].unique()):
            mask_b = gene_filt.obs[args.batch_col] == b
            gene_b = gene_filt[mask_b].copy()
            guide_b = guide_filt[mask_b].copy()

            metrics_post_b = compute_metrics(
                gene=gene_b,
                guide=guide_b,
                filter_stage="post",
                sample_name=args.sample_name,
                umi_col=args.umi_col,
                genes_col=args.genes_col,
                mito_col=args.mito_col,
                guide_umi_col=args.guide_umi_col,
                guides_per_cell_col=args.guides_per_cell_col,
                batch_label=str(b),
            )
            metrics_post_list.append(metrics_post_b)

    metrics_post = pd.concat(metrics_post_list, ignore_index=True)

    # -----------------------------------------------------------------
    # Write combined metrics TSV (single file)
    # -----------------------------------------------------------------
    metrics_all = pd.concat([metrics_pre, metrics_post], ignore_index=True)
    metrics_path = os.path.join(args.outdir, args.metrics_tsv)
    metrics_all.to_csv(metrics_path, sep="\t", index=False)
    logger.info(f"Wrote combined metrics TSV to: {metrics_path}")

    # -----------------------------------------------------------------
    # Save filtered AnnData (optional, split gene + guide)
    # -----------------------------------------------------------------
    if args.write_filtered_h5ad:
        gene_out = os.path.join(args.outdir, args.gene_h5ad_name)
        guide_out = os.path.join(args.outdir, args.guide_h5ad_name)

        logger.info(f"Writing filtered gene AnnData to {gene_out}")
        gene_filt.write(gene_out)

        logger.info(f"Writing filtered guide AnnData to {guide_out}")
        guide_filt.write(guide_out)

    logger.info("QC script finished.")


if __name__ == "__main__":
    main()
