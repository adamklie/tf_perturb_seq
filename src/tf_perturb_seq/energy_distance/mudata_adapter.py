#!/usr/bin/env python3
"""
Convert MuData to energy distance pipeline inputs.

This script extracts and preprocesses data from a MuData file (.h5mu) to create
the input files required by the energy distance pipeline:

1. adata.h5ad - Preprocessed gene expression AnnData with PCA embeddings
2. grna_dict.pkl - Dictionary mapping cell barcodes to assigned guide IDs
3. annotation.csv - Guide metadata with required columns

Preprocessing steps (for gene expression):
  - Normalize to target sum (default: 10,000)
  - Log1p transform
  - Select highly variable genes (default: 5,000)
  - Regress out confounders (total_gene_umis, percent_mito)
  - Scale (max_value=10)
  - PCA (default: 50 components)

Usage:
    python mudata_adapter.py \\
        --input inference_mudata.h5mu \\
        --outdir ./energy_distance_inputs \\
        [--n-top-genes 5000] \\
        [--n-pcs 50] \\
        [--batch-key batch] \\
        [--use-harmony]
"""

import argparse
import logging
import os
import pickle
import sys
from pathlib import Path

import mudata as mu
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)


def preprocess_gene_expression(
    adata: sc.AnnData,
    target_sum: float = 1e4,
    n_top_genes: int = 5000,
    max_value: float = 10,
    regress_keys: list = None,
) -> sc.AnnData:
    """
    Preprocess gene expression data for energy distance analysis.

    Steps:
        1. Normalize total counts
        2. Log1p transform
        3. Select highly variable genes
        4. Store raw counts
        5. Regress out confounders
        6. Scale

    Args:
        adata: AnnData with raw counts in .X
        target_sum: Target sum for normalization
        n_top_genes: Number of highly variable genes to select
        max_value: Maximum value for scaling
        regress_keys: Columns in .obs to regress out (default: total_gene_umis, percent_mito)

    Returns:
        Preprocessed AnnData (subset to HVGs)
    """
    if regress_keys is None:
        regress_keys = []
        if "total_gene_umis" in adata.obs.columns:
            regress_keys.append("total_gene_umis")
        if "percent_mito" in adata.obs.columns:
            regress_keys.append("percent_mito")

    logger.info(f"Starting preprocessing: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Make a copy to avoid modifying original
    adata = adata.copy()

    # Normalize
    logger.info(f"Normalizing to target_sum={target_sum}")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    # Highly variable genes
    logger.info(f"Selecting top {n_top_genes} highly variable genes")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="cell_ranger")

    # Store raw (normalized, log-transformed, all genes)
    adata.raw = adata.copy()

    # Subset to HVGs
    adata = adata[:, adata.var["highly_variable"]].copy()
    logger.info(f"After HVG selection: {adata.shape[1]} genes")

    # Regress out confounders
    if regress_keys:
        # Check which keys exist
        valid_keys = [k for k in regress_keys if k in adata.obs.columns]
        if valid_keys:
            logger.info(f"Regressing out: {valid_keys}")
            sc.pp.regress_out(adata, keys=valid_keys)
        else:
            logger.warning(f"No valid regress keys found. Requested: {regress_keys}")

    # Scale
    logger.info(f"Scaling with max_value={max_value}")
    sc.pp.scale(adata, max_value=max_value)

    return adata


def compute_pca(
    adata: sc.AnnData,
    n_pcs: int = 50,
    batch_key: str = None,
    use_harmony: bool = False,
) -> sc.AnnData:
    """
    Compute PCA and optionally integrate with Harmony.

    Args:
        adata: Preprocessed AnnData
        n_pcs: Number of principal components
        batch_key: Column in .obs for batch correction (required if use_harmony=True)
        use_harmony: Whether to run Harmony batch correction

    Returns:
        AnnData with PCA in .obsm['X_pca']
    """
    logger.info(f"Computing PCA with {n_pcs} components")
    sc.tl.pca(adata, n_comps=n_pcs)

    if use_harmony:
        if batch_key is None:
            raise ValueError("batch_key required for Harmony integration")
        if batch_key not in adata.obs.columns:
            raise ValueError(f"batch_key '{batch_key}' not found in adata.obs")

        logger.info(f"Running Harmony integration on batch_key='{batch_key}'")
        try:
            import scanpy.external as sce
            sce.pp.harmony_integrate(adata, key=batch_key)
            # Harmony stores result in X_pca_harmony, copy to X_pca for pipeline compatibility
            adata.obsm["X_pca_original"] = adata.obsm["X_pca"].copy()
            adata.obsm["X_pca"] = adata.obsm["X_pca_harmony"]
            logger.info("Harmony integration complete, using X_pca_harmony as X_pca")
        except ImportError:
            logger.error("harmonypy not installed. Run: pip install harmonypy")
            raise

    logger.info(f"PCA complete: X_pca shape = {adata.obsm['X_pca'].shape}")
    return adata


def extract_grna_assignment_matrix(
    guide_adata: sc.AnnData,
    assignment_layer: str = "guide_assignment",
) -> pd.DataFrame:
    """
    Extract guide assignment matrix as a DataFrame for the energy distance pipeline.

    The upstream pipeline (util_functions.load_files) expects a pickle containing
    a DataFrame with guides as rows and cells as columns. It transposes this to
    (cells x guides), then builds a {gRNA_name: [cell_names]} dict from non-zero entries.

    Args:
        guide_adata: Guide modality AnnData
        assignment_layer: Layer containing binary assignment matrix

    Returns:
        DataFrame with guides as rows, cells as columns (guides x cells)
    """
    logger.info("Extracting gRNA assignment matrix from guide_assignment layer")

    if assignment_layer not in guide_adata.layers:
        raise ValueError(f"Layer '{assignment_layer}' not found. Available: {list(guide_adata.layers.keys())}")

    assignment_matrix = guide_adata.layers[assignment_layer]
    if issparse(assignment_matrix):
        assignment_matrix = assignment_matrix.toarray()

    # Create cells x guides DataFrame, then transpose to guides x cells
    # (upstream pipeline expects guides x cells and transposes it)
    df = pd.DataFrame(
        assignment_matrix,
        index=guide_adata.obs_names,
        columns=guide_adata.var_names,
    )

    n_cells_with_guide = (df.sum(axis=1) > 0).sum()
    n_assignments = int((df.values > 0).sum())
    logger.info(f"Assignment matrix: {df.shape[0]} cells x {df.shape[1]} guides")
    logger.info(f"Cells with at least one guide: {n_cells_with_guide}/{df.shape[0]}")
    logger.info(f"Total guide assignments: {n_assignments}")

    # Transpose to guides x cells (expected by upstream pipeline)
    return df.T


def extract_annotation(
    guide_adata: sc.AnnData,
    type_col: str = "type",
    target_col: str = "intended_target_name",
) -> pd.DataFrame:
    """
    Extract guide annotation DataFrame for energy distance pipeline.

    Required columns for pipeline:
        - guide_id: Unique guide identifier
        - type: 'targeting' or 'non-targeting'
        - intended_target_name: Target gene/region name

    Args:
        guide_adata: Guide modality AnnData
        type_col: Column in .var containing guide type
        target_col: Column in .var containing target name

    Returns:
        DataFrame with guide annotations
    """
    logger.info("Extracting guide annotation")

    var_df = guide_adata.var.copy()

    # Ensure guide_id column exists
    if "guide_id" not in var_df.columns:
        var_df["guide_id"] = var_df.index

    # Check required columns
    if type_col not in var_df.columns:
        raise ValueError(f"Type column '{type_col}' not found. Available: {list(var_df.columns)}")
    if target_col not in var_df.columns:
        raise ValueError(f"Target column '{target_col}' not found. Available: {list(var_df.columns)}")

    # Standardize type column for pipeline compatibility
    # Pipeline expects: 'targeting' or 'non-targeting'
    def standardize_type(x):
        x_lower = str(x).lower().strip()
        if "non" in x_lower and "target" in x_lower:
            return "non-targeting"
        elif "negative" in x_lower:
            return "non-targeting"
        elif "target" in x_lower or "positive" in x_lower or "control" in x_lower:
            return "targeting"
        else:
            return "targeting"  # Default to targeting for unknown types

    var_df["type_standardized"] = var_df[type_col].apply(standardize_type)

    # Log type distribution
    logger.info(f"Guide type distribution (original '{type_col}'):")
    for t, count in var_df[type_col].value_counts().items():
        logger.info(f"  {t}: {count}")

    logger.info(f"Guide type distribution (standardized):")
    for t, count in var_df["type_standardized"].value_counts().items():
        logger.info(f"  {t}: {count}")

    # Select and rename columns for pipeline
    annotation_df = var_df[["guide_id", target_col, "type_standardized"]].copy()
    annotation_df = annotation_df.rename(columns={
        target_col: "intended_target_name",
        "type_standardized": "type",
    })

    # Add additional columns if available
    optional_cols = ["gene_name", "label", "spacer"]
    for col in optional_cols:
        if col in var_df.columns:
            annotation_df[col] = var_df[col].values

    logger.info(f"Annotation shape: {annotation_df.shape}")
    return annotation_df


def prepare_energy_distance_inputs(
    input_path: str,
    outdir: str,
    n_top_genes: int = 5000,
    n_pcs: int = 50,
    target_sum: float = 1e4,
    max_value: float = 10,
    batch_key: str = None,
    use_harmony: bool = False,
    gene_modality: str = "gene",
    guide_modality: str = "guide",
    assignment_layer: str = "guide_assignment",
    type_col: str = "type",
    target_col: str = "intended_target_name",
    pca_key: str = "X_pca",
    skip_preprocessing: bool = False,
) -> dict:
    """
    Prepare all inputs for energy distance pipeline from MuData.

    Args:
        input_path: Path to MuData file (.h5mu)
        outdir: Output directory for pipeline inputs
        n_top_genes: Number of highly variable genes
        n_pcs: Number of PCA components
        target_sum: Normalization target sum
        max_value: Scaling max value
        batch_key: Batch column for Harmony (optional)
        use_harmony: Whether to run Harmony integration
        gene_modality: Name of gene expression modality
        guide_modality: Name of guide modality
        assignment_layer: Layer with guide assignments
        type_col: Column in guide.var with guide type
        target_col: Column in guide.var with target name
        pca_key: Key for PCA in obsm (to check if already computed)
        skip_preprocessing: If True and PCA exists, skip preprocessing

    Returns:
        Dictionary with paths to created files
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Loading MuData from: {input_path}")
    mdata = mu.read(input_path)

    # Validate modalities
    if gene_modality not in mdata.mod:
        raise ValueError(f"Gene modality '{gene_modality}' not found. Available: {list(mdata.mod.keys())}")
    if guide_modality not in mdata.mod:
        raise ValueError(f"Guide modality '{guide_modality}' not found. Available: {list(mdata.mod.keys())}")

    gene = mdata.mod[gene_modality]
    guide = mdata.mod[guide_modality]

    logger.info(f"Gene modality: {gene.shape[0]} cells x {gene.shape[1]} genes")
    logger.info(f"Guide modality: {guide.shape[0]} cells x {guide.shape[1]} guides")

    # Check if cells match
    if not gene.obs_names.equals(guide.obs_names):
        logger.warning("Cell barcodes don't match between modalities. Using intersection.")
        common_cells = gene.obs_names.intersection(guide.obs_names)
        gene = gene[common_cells, :].copy()
        guide = guide[common_cells, :].copy()
        logger.info(f"After intersection: {len(common_cells)} cells")

    # Process gene expression
    if skip_preprocessing and pca_key in gene.obsm:
        logger.info(f"Skipping preprocessing: PCA already exists in obsm['{pca_key}']")
        adata = gene.copy()
    else:
        adata = preprocess_gene_expression(
            gene,
            target_sum=target_sum,
            n_top_genes=n_top_genes,
            max_value=max_value,
        )
        adata = compute_pca(
            adata,
            n_pcs=n_pcs,
            batch_key=batch_key,
            use_harmony=use_harmony,
        )

    # Extract gRNA assignment matrix (guides x cells DataFrame)
    grna_matrix = extract_grna_assignment_matrix(guide, assignment_layer=assignment_layer)

    # Extract annotation
    annotation_df = extract_annotation(guide, type_col=type_col, target_col=target_col)

    # Save outputs
    output_paths = {}

    # 1. Save AnnData with PCA
    adata_path = outdir / "adata.h5ad"
    logger.info(f"Saving AnnData to: {adata_path}")
    adata.write(adata_path)
    output_paths["adata"] = str(adata_path)

    # 2. Save gRNA assignment matrix (guides x cells DataFrame)
    grna_path = outdir / "grna_dict.pkl"
    logger.info(f"Saving gRNA assignment matrix to: {grna_path}")
    grna_matrix.to_pickle(grna_path)
    output_paths["grna_dict"] = str(grna_path)

    # 3. Save annotation
    annotation_path = outdir / "annotation.csv"
    logger.info(f"Saving annotation to: {annotation_path}")
    annotation_df.to_csv(annotation_path, index=False)
    output_paths["annotation"] = str(annotation_path)

    logger.info("=" * 50)
    logger.info("Energy distance inputs prepared successfully!")
    logger.info(f"  AnnData: {adata_path} ({adata.shape[0]} cells, {adata.shape[1]} genes)")
    logger.info(f"  gRNA matrix: {grna_path} ({grna_matrix.shape[0]} guides x {grna_matrix.shape[1]} cells)")
    logger.info(f"  Annotation: {annotation_path} ({len(annotation_df)} guides)")
    logger.info("=" * 50)

    return output_paths


def main():
    parser = argparse.ArgumentParser(
        description="Convert MuData to energy distance pipeline inputs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to MuData file (.h5mu)",
    )
    parser.add_argument(
        "--outdir", "-o",
        required=True,
        help="Output directory for pipeline inputs",
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=5000,
        help="Number of highly variable genes (default: 5000)",
    )
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=50,
        help="Number of PCA components (default: 50)",
    )
    parser.add_argument(
        "--target-sum",
        type=float,
        default=1e4,
        help="Normalization target sum (default: 10000)",
    )
    parser.add_argument(
        "--max-value",
        type=float,
        default=10,
        help="Scaling max value (default: 10)",
    )
    parser.add_argument(
        "--batch-key",
        type=str,
        default=None,
        help="Batch column in obs for Harmony integration",
    )
    parser.add_argument(
        "--use-harmony",
        action="store_true",
        help="Run Harmony batch correction (requires --batch-key)",
    )
    parser.add_argument(
        "--gene-modality",
        type=str,
        default="gene",
        help="Name of gene expression modality (default: gene)",
    )
    parser.add_argument(
        "--guide-modality",
        type=str,
        default="guide",
        help="Name of guide modality (default: guide)",
    )
    parser.add_argument(
        "--assignment-layer",
        type=str,
        default="guide_assignment",
        help="Layer with guide assignments (default: guide_assignment)",
    )
    parser.add_argument(
        "--type-col",
        type=str,
        default="type",
        help="Column in guide.var with guide type (default: type)",
    )
    parser.add_argument(
        "--target-col",
        type=str,
        default="intended_target_name",
        help="Column in guide.var with target name (default: intended_target_name)",
    )
    parser.add_argument(
        "--skip-preprocessing",
        action="store_true",
        help="Skip preprocessing if PCA already exists",
    )

    args = parser.parse_args()

    if args.use_harmony and args.batch_key is None:
        parser.error("--use-harmony requires --batch-key")

    prepare_energy_distance_inputs(
        input_path=args.input,
        outdir=args.outdir,
        n_top_genes=args.n_top_genes,
        n_pcs=args.n_pcs,
        target_sum=args.target_sum,
        max_value=args.max_value,
        batch_key=args.batch_key,
        use_harmony=args.use_harmony,
        gene_modality=args.gene_modality,
        guide_modality=args.guide_modality,
        assignment_layer=args.assignment_layer,
        type_col=args.type_col,
        target_col=args.target_col,
        skip_preprocessing=args.skip_preprocessing,
    )


if __name__ == "__main__":
    main()
