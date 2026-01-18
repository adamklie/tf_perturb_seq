#!/usr/bin/env python

"""
Helper functions for TF Perturb-seq QC on MuData outputs.
"""

import logging
import os
from typing import Optional, Union

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns            # noqa: E402

from anndata import AnnData
from mudata import MuData

logger = logging.getLogger(__name__)

sns.set_style("whitegrid")


# ---------------------------------------------------------------------
# Basic utilities
# ---------------------------------------------------------------------
def describe_mudata(data: MuData) -> None:
    """Log basic info about a MuData object."""
    logger.info(f"MuData: {data.n_obs} cells total, modalities={list(data.mod.keys())}")
    for mod_key in data.mod.keys():
        logger.info(
            f"  Modality '{mod_key}': {data.mod[mod_key].n_obs} cells, "
            f"{data.mod[mod_key].n_vars} features"
        )


def add_sample_prefix_to_barcodes(
    ad: Union[AnnData, MuData],
    sample_name: Optional[str],
) -> None:
    """
    Prefix barcodes with sample name in-place.

    Works for both AnnData and MuData (via `.obs_names`).
    """
    if sample_name is None:
        return
    if not len(ad.obs_names):
        return
    logger.info(f"Adding sample name '{sample_name}' to barcodes")
    logger.info(f"  Example before: {ad.obs_names[0]}")
    ad.obs.index = sample_name + "#" + ad.obs_names.astype(str)
    logger.info(f"  Example after:  {ad.obs_names[0]}")


# ---------------------------------------------------------------------
# Generic plots
# ---------------------------------------------------------------------
def plot_cells_per_batch(
    adata: AnnData,
    outdir: str,
    prefix: str,
    batch_col: str = "batch",
) -> None:
    """Barplot of number of cells per batch / measurement set."""
    if batch_col not in adata.obs:
        logger.warning(f"Column '{batch_col}' not in adata.obs; skipping cells-per-batch plot.")
        return

    os.makedirs(outdir, exist_ok=True)

    with sns.plotting_context("talk"):
        fig, ax = plt.subplots(figsize=(6, 5))
        sns.countplot(x=batch_col, data=adata.obs, ax=ax)
        ax.set_title(f"Number of cells per measurement set (N={adata.n_obs})")
        ax.set_xlabel(batch_col)
        plt.xticks(rotation=90)
        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_cells_per_batch.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


def plot_knee_plot(
    adata: AnnData,
    outdir: str,
    prefix: str,
    use_X: bool = True,
    count_col: Optional[str] = None,
    title: str = "Knee Plot of UMI Counts",
) -> None:
    """
    Knee plot of UMI count (x-axis) vs barcode rank (y-axis).

    If use_X=True: sum counts from adata.X.
    Else: use adata.obs[count_col].
    """
    os.makedirs(outdir, exist_ok=True)

    if use_X:
        # adata.X can be sparse or dense
        sums = np.array(adata.X.sum(axis=1)).ravel()
    else:
        if count_col is None or count_col not in adata.obs.columns:
            logger.warning(
                f"Cannot make knee plot: use_X=False and count_col='{count_col}' "
                "not found in adata.obs."
            )
            return
        sums = adata.obs[count_col].to_numpy()

    knee_df = pd.DataFrame(
        {
            "sum": sums,
            "barcodes": adata.obs_names.values,
        }
    )
    knee_df = knee_df.sort_values("sum", ascending=False).reset_index(drop=True)
    knee_df["sum_log"] = np.log1p(knee_df["sum"])

    with sns.plotting_context("talk"):
        plt.figure(figsize=(10, 5))
        plt.plot(
            knee_df["sum_log"],
            knee_df.index,
            marker="o",
            linestyle="-",
            markersize=3,
        )
        plt.ylabel("Barcode number")
        plt.xlabel("Log of UMI Counts")
        plt.title(title)
        fname = os.path.join(outdir, f"{prefix}_knee_plot.png")
        plt.tight_layout()
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")

    logger.info(
        f"{title} - min count: {knee_df['sum'].min():.0f}, max count: {knee_df['sum'].max():.0f}"
    )


# ---------------------------------------------------------------------
# Gene QC plots
# ---------------------------------------------------------------------
def plot_gene_qc_histograms(
    gene: AnnData,
    outdir: str,
    prefix: str,
    umi_col: str = "total_gene_umis",
    genes_col: str = "num_expressed_genes",
    mito_col: str = "percent_mito",
    bins: int = 100,
    by_batch: bool = False,
    batch_col: str = "batch",
) -> None:
    """
    Plot three histograms:
      - umi_col
      - genes_col
      - mito_col

    If by_batch is False (default): simple global histograms with a single median.
    If by_batch is True and batch_col exists: per-batch overlaid histograms,
    density-scaled, with per-batch medians and colored dashed lines.
    """
    os.makedirs(outdir, exist_ok=True)

    missing = [c for c in [umi_col, genes_col, mito_col] if c not in gene.obs.columns]
    if missing:
        logger.warning(f"Missing QC columns in gene.obs: {missing}. Skipping QC hist plots.")
        return

    # If not splitting by batch, fall back to simple histograms
    if (not by_batch) or (batch_col not in gene.obs.columns):
        with sns.plotting_context("talk", font_scale=1.2):
            fig, ax = plt.subplots(3, 1, figsize=(12, 10))

            # 1) UMI counts
            sns.histplot(gene.obs[umi_col], bins=bins, ax=ax[0])
            median_umis = gene.obs[umi_col].median()
            ax[0].axvline(median_umis, linestyle="--")
            ax[0].text(
                0.95,
                0.95,
                f"Median: {median_umis:.0f} UMIs",
                ha="right",
                va="top",
                transform=ax[0].transAxes,
            )
            ax[0].set_title("Gene UMI counts")
            ax[0].set_xlabel(umi_col)

            # 2) Genes detected
            sns.histplot(gene.obs[genes_col], bins=bins, ax=ax[1])
            median_genes = gene.obs[genes_col].median()
            ax[1].axvline(median_genes, linestyle="--")
            ax[1].text(
                0.95,
                0.95,
                f"Median: {median_genes:.0f} genes",
                ha="right",
                va="top",
                transform=ax[1].transAxes,
            )
            ax[1].set_title("Genes detected (> 0 UMIs)")
            ax[1].set_xlabel(genes_col)

            # 3) Mito %
            sns.histplot(gene.obs[mito_col], bins=bins, ax=ax[2])
            median_mito = gene.obs[mito_col].median()
            ax[2].axvline(median_mito, linestyle="--")
            ax[2].text(
                0.95,
                0.95,
                f"Median: {median_mito:.2f} %",
                ha="right",
                va="top",
                transform=ax[2].transAxes,
            )
            ax[2].set_title("Mitochondrial %")
            ax[2].set_xlabel(mito_col)

            plt.tight_layout()
            fname = os.path.join(outdir, f"{prefix}_gene_qc_histograms.png")
            plt.savefig(fname, dpi=200)
            plt.close()
            logger.info(f"Saved {fname}")
        return

    # ------------------------------------------------------------------
    # Split by batch
    # ------------------------------------------------------------------
    batches = gene.obs[batch_col].astype("category")
    batch_order = list(batches.cat.categories)
    batch_colors = sns.color_palette("tab10", n_colors=len(batch_order))
    batch_color_dict = dict(zip(batch_order, batch_colors))

    def _text_y(i: int) -> float:
        # stagger labels vertically
        return 0.95 - 0.08 * i

    with sns.plotting_context("talk", font_scale=1.0):
        fig, ax = plt.subplots(3, 1, figsize=(10, 10))

        # 1) UMI counts
        sns.histplot(
            data=gene.obs,
            x=umi_col,
            bins=bins,
            hue=batch_col,
            palette=batch_color_dict,
            ax=ax[0],
            element="step",
            stat="density",
            common_norm=False,
        )
        for i, batch in enumerate(batch_order):
            median_val = gene.obs.loc[gene.obs[batch_col] == batch, umi_col].median()
            ax[0].axvline(median_val, color=batch_color_dict[batch], linestyle="--")
            ax[0].text(
                0.6,
                _text_y(i),
                f"{batch}: {median_val:.0f} UMIs",
                ha="right",
                va="top",
                transform=ax[0].transAxes,
                color=batch_color_dict[batch],
                fontsize=12,
            )
        ax[0].set_title("Gene UMI counts")
        ax[0].set_xlabel(umi_col)

        # 2) Genes detected
        sns.histplot(
            data=gene.obs,
            x=genes_col,
            bins=bins,
            hue=batch_col,
            palette=batch_color_dict,
            ax=ax[1],
            element="step",
            stat="density",
            common_norm=False,
        )
        for i, batch in enumerate(batch_order):
            median_val = gene.obs.loc[gene.obs[batch_col] == batch, genes_col].median()
            ax[1].axvline(median_val, color=batch_color_dict[batch], linestyle="--")
            ax[1].text(
                0.8,
                _text_y(i),
                f"{batch}: {median_val:.0f} genes",
                ha="right",
                va="top",
                transform=ax[1].transAxes,
                color=batch_color_dict[batch],
                fontsize=12,
            )
        ax[1].set_title("Genes detected (> 0 UMIs)")
        ax[1].set_xlabel(genes_col)
        if ax[1].legend_ is not None:
            ax[1].legend_.remove()

        # 3) Mito %
        sns.histplot(
            data=gene.obs,
            x=mito_col,
            bins=bins,
            hue=batch_col,
            palette=batch_color_dict,
            ax=ax[2],
            element="step",
            stat="density",
            common_norm=False,
        )
        for i, batch in enumerate(batch_order):
            median_val = gene.obs.loc[gene.obs[batch_col] == batch, mito_col].median()
            ax[2].axvline(median_val, color=batch_color_dict[batch], linestyle="--")
            ax[2].text(
                0.6,
                _text_y(i),
                f"{batch}: {median_val:.2f} %",
                ha="right",
                va="top",
                transform=ax[2].transAxes,
                color=batch_color_dict[batch],
                fontsize=12,
            )
        ax[2].set_title("Mitochondrial %")
        ax[2].set_xlabel(mito_col)
        if ax[2].legend_ is not None:
            ax[2].legend_.remove()

        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_gene_qc_histograms_by_batch.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


# ---------------------------------------------------------------------
# Guide modality
# ---------------------------------------------------------------------
def prepare_guide_modality(
    data: MuData,
    guide_mod_key: str = "guide",
    guide_assignment_layer: str = "guide_assignment",
) -> AnnData:
    """
    Extracts the guide modality and computes:
      - guide.obs['n_guides_per_cell']
      - guide.var['n_cells_per_guide']
    from layers['guide_assignment'] (>0 indicates assignment).

    Assumes the modality exists. Returns the guide AnnData (view).
    """
    if guide_mod_key not in data.mod:
        raise ValueError(f"Guide modality '{guide_mod_key}' not found in MuData.mod")

    guide = data.mod[guide_mod_key]

    if guide_assignment_layer not in guide.layers:
        logger.warning(
            f"No layer '{guide_assignment_layer}' in guide modality; "
            "cannot compute n_guides_per_cell / n_cells_per_guide."
        )
        return guide

    A = guide.layers[guide_assignment_layer]
    A_bin = (A > 0)
    # handle sparse matrices with .A
    if hasattr(A_bin, "A"):
        A_bin = A_bin.A

    guide.obs["n_guides_per_cell"] = np.asarray(A_bin.sum(axis=1)).ravel()
    guide.var["n_cells_per_guide"] = np.asarray(A_bin.sum(axis=0)).ravel()

    logger.info(
        "Guide modality: computed n_guides_per_cell (obs) and n_cells_per_guide (var) "
        f"from layer '{guide_assignment_layer}'"
    )

    return guide


def plot_guide_qc(
    guide: AnnData,
    outdir: str,
    prefix: str,
    batch_col: str = "batch",
    label_col: str = "label",
    guide_umi_col: str = "total_guide_umis",
    guides_per_cell_col: str = "n_guides_per_cell",
    cells_per_guide_col: str = "n_cells_per_guide",
    bins_umis: int = 100,
    bins_guides: int = 30,
) -> None:
    """
    Three-panel guide QC:
      1. total_guide_umis by batch
      2. n_guides_per_cell by batch
      3. n_cells_per_guide by label
    """
    os.makedirs(outdir, exist_ok=True)

    missing_obs = [
        c
        for c in [guide_umi_col, guides_per_cell_col, batch_col]
        if c not in guide.obs.columns
    ]
    missing_var = [
        c
        for c in [cells_per_guide_col, label_col]
        if c not in guide.var.columns
    ]
    if missing_obs:
        logger.warning(
            f"Skipping guide QC plots: missing columns in guide.obs: {missing_obs}"
        )
        return
    if missing_var:
        logger.warning(
            f"Skipping guide QC plots: missing columns in guide.var: {missing_var}"
        )
        return

    # Colors
    batches = guide.obs[batch_col].astype("category")
    batch_order = list(batches.cat.categories)
    batch_colors = sns.color_palette("tab10", n_colors=len(batch_order))
    batch_color_dict = dict(zip(batch_order, batch_colors))

    labels = guide.var[label_col].astype("category")
    label_order = list(labels.cat.categories)
    label_colors = sns.color_palette("tab10", n_colors=len(label_order))
    guide_label_color_dict = dict(zip(label_order, label_colors))

    def _text_y(i: int) -> float:
        return 0.95 - 0.08 * i

    with sns.plotting_context("talk", font_scale=1.0):
        fig, ax = plt.subplots(3, 1, figsize=(10, 10))

        # 1) Total guide UMIs per cell, split by batch
        sns.histplot(
            data=guide.obs,
            x=guide_umi_col,
            bins=bins_umis,
            hue=batch_col,
            palette=batch_color_dict,
            ax=ax[0],
            element="step",
            stat="density",
            common_norm=False,
        )
        for i, batch in enumerate(batch_order):
            median_val = guide.obs.loc[guide.obs[batch_col] == batch, guide_umi_col].median()
            ax[0].axvline(median_val, color=batch_color_dict[batch], linestyle="--")
            ax[0].text(
                0.5,
                _text_y(i),
                f"{batch}: {median_val:.0f} UMIs",
                ha="right",
                va="top",
                transform=ax[0].transAxes,
                color=batch_color_dict[batch],
                fontsize=12,
            )
        ax[0].set_title("Total guide UMI counts")
        ax[0].set_xlabel(guide_umi_col)

        # 2) n_guides_per_cell, split by batch
        sns.histplot(
            data=guide.obs,
            x=guides_per_cell_col,
            bins=bins_guides,
            hue=batch_col,
            palette=batch_color_dict,
            ax=ax[1],
            element="step",
            stat="density",
            common_norm=False,
        )
        for i, batch in enumerate(batch_order):
            mean_val = guide.obs.loc[guide.obs[batch_col] == batch, guides_per_cell_col].mean()
            ax[1].axvline(mean_val, color=batch_color_dict[batch], linestyle="--")
            ax[1].text(
                0.5,
                _text_y(i),
                f"{batch}: {mean_val:.2f} guides",
                ha="right",
                va="top",
                transform=ax[1].transAxes,
                color=batch_color_dict[batch],
                fontsize=12,
            )
        ax[1].set_title("Number of guides ASSIGNED per cell")
        ax[1].set_xlabel(guides_per_cell_col)
        if ax[1].legend_ is not None:
            ax[1].legend_.remove()

        # 3) n_cells_per_guide by guide label
        sns.histplot(
            data=guide.var,
            x=cells_per_guide_col,
            bins=bins_guides,
            hue=label_col,
            palette=guide_label_color_dict,
            ax=ax[2],
            element="step",
            stat="density",
            common_norm=False,
        )
        for i, label in enumerate(label_order):
            median_val = guide.var.loc[
                guide.var[label_col] == label, cells_per_guide_col
            ].median()
            ax[2].axvline(median_val, color=guide_label_color_dict[label], linestyle="--")
            ax[2].text(
                0.6,
                _text_y(i),
                f"{label}: {median_val:.0f} cells",
                ha="right",
                va="top",
                transform=ax[2].transAxes,
                color=guide_label_color_dict[label],
                fontsize=12,
            )
        ax[2].set_title("Number of cells per guide assignment")
        ax[2].set_xlabel(cells_per_guide_col)
        if ax[2].legend_ is not None:
            ax[2].legend_.remove()

        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_guide_qc.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


# ---------------------------------------------------------------------
# Filtering (joint gene + guide)
# ---------------------------------------------------------------------
def apply_filters(
    gene: AnnData,
    guide: AnnData,
    umi_col: str = "total_gene_umis",
    genes_col: str = "num_expressed_genes",
    mito_col: str = "percent_mito",
    guides_per_cell_col: str = "n_guides_per_cell",
    min_total_umis: Optional[float] = None,
    max_total_umis: Optional[float] = None,
    min_genes: Optional[float] = None,
    max_genes: Optional[float] = None,
    max_mito_pct: Optional[float] = None,
    min_guides_per_cell: Optional[float] = None,
    max_guides_per_cell: Optional[float] = None,
    min_cells_per_guide: Optional[int] = None,
) -> pd.Series:
    """
    Apply cell-level filtering based on both gene and guide QC metrics.

    Returns a boolean mask (pd.Series indexed by cell) to apply jointly
    to both gene and guide modalities (joint filtering).

    Filters:
      - UMI, genes, mito from gene.obs
      - guides per cell from guide.obs[guides_per_cell_col]
      - min_cells_per_guide is currently not implemented (logs a warning).
    """
    if not gene.obs_names.equals(guide.obs_names):
        raise ValueError("gene.obs_names and guide.obs_names differ; cannot joint-filter safely.")

    mask = pd.Series(True, index=gene.obs_names, dtype=bool)

    # Gene-based filters
    if umi_col in gene.obs.columns:
        if min_total_umis is not None:
            mask &= gene.obs[umi_col] >= min_total_umis
        if max_total_umis is not None:
            mask &= gene.obs[umi_col] <= max_total_umis

    if genes_col in gene.obs.columns:
        if min_genes is not None:
            mask &= gene.obs[genes_col] >= min_genes
        if max_genes is not None:
            mask &= gene.obs[genes_col] <= max_genes

    if mito_col in gene.obs.columns and max_mito_pct is not None:
        mask &= gene.obs[mito_col] <= max_mito_pct

    # Guide-based filters (per-cell)
    if guides_per_cell_col in guide.obs.columns:
        gpc = guide.obs[guides_per_cell_col]
        if min_guides_per_cell is not None:
            mask &= gpc >= min_guides_per_cell
        if max_guides_per_cell is not None:
            mask &= gpc <= max_guides_per_cell

    # Placeholder for future min_cells_per_guide logic
    if min_cells_per_guide is not None:
        logger.warning(
            "min_cells_per_guide filtering is not implemented yet (requires "
            "propagating guide-level counts to cell-level masks). Ignoring this parameter."
        )

    logger.info(f"Cells before filtering: {mask.size}")
    logger.info(f"Cells after filtering:  {mask.sum()}")

    return mask


# ---------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------
def compute_metrics(
    gene: AnnData,
    guide: AnnData,
    filter_stage: str,
    sample_name: Optional[str],
    umi_col: str = "total_gene_umis",
    genes_col: str = "num_expressed_genes",
    mito_col: str = "percent_mito",
    guide_umi_col: str = "total_guide_umis",
    guides_per_cell_col: str = "n_guides_per_cell",
    cells_per_guide_col: str = "n_cells_per_guide",
    batch_label: str = "all",
) -> pd.DataFrame:
    """
    Compute summary metrics for one (gene, guide) pair and one filter_stage.

    Returns a one-row DataFrame with:
      - sample_name
      - filter_stage ('pre' or 'post' or custom)
      - batch (e.g. 'all' or specific batch name)
      - n_cells_gene
      - median_total_gene_umis
      - median_num_expressed_genes
      - median_percent_mito
      - n_cells_guide
      - median_total_guide_umis
      - mean_n_guides_per_cell
      - frac_cells_with_guide_assignment
      - median_n_cells_per_guide
    """
    metrics = {
        "sample_name": sample_name if sample_name is not None else "",
        "filter_stage": filter_stage,
        "batch": batch_label,
        "n_cells_gene": gene.n_obs,
        "n_cells_guide": guide.n_obs,
    }

    # Gene metrics
    if umi_col in gene.obs:
        metrics["median_total_gene_umis"] = float(gene.obs[umi_col].median())
    else:
        metrics["median_total_gene_umis"] = np.nan

    if genes_col in gene.obs:
        metrics["median_num_expressed_genes"] = float(gene.obs[genes_col].median())
    else:
        metrics["median_num_expressed_genes"] = np.nan

    if mito_col in gene.obs:
        metrics["median_percent_mito"] = float(gene.obs[mito_col].median())
    else:
        metrics["median_percent_mito"] = np.nan

    # Guide metrics
    if guide_umi_col in guide.obs:
        metrics["median_total_guide_umis"] = float(guide.obs[guide_umi_col].median())
    else:
        metrics["median_total_guide_umis"] = np.nan

    if guides_per_cell_col in guide.obs:
        gpc = guide.obs[guides_per_cell_col]
        metrics["mean_n_guides_per_cell"] = float(gpc.mean())
        metrics["frac_cells_with_guide_assignment"] = float((gpc > 0).mean())
    else:
        metrics["mean_n_guides_per_cell"] = np.nan
        metrics["frac_cells_with_guide_assignment"] = np.nan

    if cells_per_guide_col in guide.var:
        metrics["median_n_cells_per_guide"] = float(guide.var[cells_per_guide_col].median())
    else:
        metrics["median_n_cells_per_guide"] = np.nan

    return pd.DataFrame([metrics])
