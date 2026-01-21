#!/usr/bin/env python3
"""
Export MuData (.h5mu) components to common compressed formats.

Exports:
- Inference results (cis/trans per element/guide) → TSV.gz
- Cell metadata → TSV.gz
- Gene metadata → TSV.gz
- Guide metadata → TSV.gz
- Gene expression matrix → MTX.gz
- Guide assignment matrix → MTX.gz
- HTO counts matrix → MTX.gz
- Guide assignment table → TSV.gz (long format)
"""

import argparse
import gzip
from pathlib import Path

import mudata as md
import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import csr_matrix, issparse


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    """Write a DataFrame to a gzip-compressed TSV file."""
    with gzip.open(path, "wt") as f:
        df.to_csv(f, sep="\t", index=True)
    print(f"  Wrote: {path} ({df.shape[0]:,} rows × {df.shape[1]:,} cols)")


def write_mtx_gz(
    matrix,
    barcodes: np.ndarray,
    features: np.ndarray,
    outdir: Path,
    feature_names: np.ndarray | None = None,
) -> None:
    """
    Write a matrix in 10x-style MTX format with gzipped components.

    Creates:
    - matrix.mtx.gz (the sparse matrix)
    - barcodes.tsv.gz (cell barcodes)
    - features.tsv.gz (feature IDs and optionally names)
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # Ensure sparse CSR format
    if not issparse(matrix):
        matrix = csr_matrix(matrix)

    # Write matrix.mtx.gz
    mtx_path = outdir / "matrix.mtx"
    mmwrite(str(mtx_path), matrix.T)  # 10x format is features × cells
    with open(mtx_path, "rb") as f_in:
        with gzip.open(str(mtx_path) + ".gz", "wb") as f_out:
            f_out.writelines(f_in)
    mtx_path.unlink()  # Remove uncompressed file

    # Write barcodes.tsv.gz
    with gzip.open(outdir / "barcodes.tsv.gz", "wt") as f:
        f.write("\n".join(barcodes) + "\n")

    # Write features.tsv.gz
    with gzip.open(outdir / "features.tsv.gz", "wt") as f:
        if feature_names is not None:
            for fid, fname in zip(features, feature_names):
                f.write(f"{fid}\t{fname}\n")
        else:
            f.write("\n".join(features) + "\n")

    nnz = matrix.nnz if issparse(matrix) else np.count_nonzero(matrix)
    print(f"  Wrote: {outdir}/ ({matrix.shape[0]:,} × {matrix.shape[1]:,}, {nnz:,} nonzero)")


def export_inference_results(mdata: md.MuData, outdir: Path) -> None:
    """Export inference results from mdata.uns to compressed TSV files."""
    print("\nExporting inference results...")
    results_dir = outdir / "inference_results"
    results_dir.mkdir(parents=True, exist_ok=True)

    result_keys = [
        "cis_per_element_results",
        "cis_per_guide_results",
        "trans_per_element_results",
        "trans_per_guide_results",
    ]

    for key in result_keys:
        if key in mdata.uns:
            df = mdata.uns[key]
            if isinstance(df, pd.DataFrame):
                write_tsv_gz(df, results_dir / f"{key}.tsv.gz")
        else:
            print(f"  Skipped: {key} (not found in mdata.uns)")


def export_cell_metadata(mdata: md.MuData, outdir: Path) -> None:
    """Export merged cell metadata to compressed TSV."""
    print("\nExporting cell metadata...")
    metadata_dir = outdir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    # Collect obs from each modality, prefixing columns to avoid collision
    dfs = []
    for mod_name, adata in mdata.mod.items():
        df = adata.obs.copy()
        df.columns = [f"{mod_name}:{col}" for col in df.columns]
        dfs.append(df)

    # Merge all modality obs (they share the same index)
    merged_obs = pd.concat(dfs, axis=1)

    # Add any MuData-level obs columns that aren't prefixed
    for col in mdata.obs.columns:
        if ":" not in col and col not in merged_obs.columns:
            merged_obs[col] = mdata.obs[col]

    # Remove duplicate columns (keep first occurrence)
    merged_obs = merged_obs.loc[:, ~merged_obs.columns.duplicated()]

    write_tsv_gz(merged_obs, metadata_dir / "cell_metadata.tsv.gz")


def export_gene_metadata(mdata: md.MuData, outdir: Path) -> None:
    """Export gene annotations to compressed TSV."""
    print("\nExporting gene metadata...")
    metadata_dir = outdir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    if "gene" in mdata.mod:
        gene_var = mdata.mod["gene"].var.copy()
        write_tsv_gz(gene_var, metadata_dir / "gene_metadata.tsv.gz")
    else:
        print("  Skipped: gene modality not found")


def export_guide_metadata(mdata: md.MuData, outdir: Path) -> None:
    """Export guide annotations to compressed TSV."""
    print("\nExporting guide metadata...")
    metadata_dir = outdir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    if "guide" in mdata.mod:
        guide_var = mdata.mod["guide"].var.copy()
        write_tsv_gz(guide_var, metadata_dir / "guide_metadata.tsv.gz")
    else:
        print("  Skipped: guide modality not found")


def export_gene_matrix(mdata: md.MuData, outdir: Path) -> None:
    """Export gene expression matrix in MTX format."""
    print("\nExporting gene expression matrix...")

    if "gene" not in mdata.mod:
        print("  Skipped: gene modality not found")
        return

    gene = mdata.mod["gene"]
    mtx_dir = outdir / "matrices" / "gene_counts"

    # Get feature names (symbols) if available
    feature_names = None
    if "symbol" in gene.var.columns:
        feature_names = gene.var["symbol"].astype(str).values

    write_mtx_gz(
        matrix=gene.X,
        barcodes=gene.obs_names.values,
        features=gene.var_names.values,
        outdir=mtx_dir,
        feature_names=feature_names,
    )


def export_guide_assignment_matrix(mdata: md.MuData, outdir: Path) -> None:
    """Export guide assignment matrix in MTX format."""
    print("\nExporting guide assignment matrix...")

    if "guide" not in mdata.mod:
        print("  Skipped: guide modality not found")
        return

    guide = mdata.mod["guide"]

    if "guide_assignment" not in guide.layers:
        print("  Skipped: guide_assignment layer not found")
        return

    mtx_dir = outdir / "matrices" / "guide_assignment"

    # Get guide IDs if available
    feature_names = None
    if "guide_id" in guide.var.columns:
        feature_names = guide.var["guide_id"].astype(str).values

    write_mtx_gz(
        matrix=guide.layers["guide_assignment"],
        barcodes=guide.obs_names.values,
        features=guide.var_names.values,
        outdir=mtx_dir,
        feature_names=feature_names,
    )


def export_hto_matrix(mdata: md.MuData, outdir: Path) -> None:
    """Export HTO counts matrix in MTX format."""
    print("\nExporting HTO counts matrix...")

    if "hashing" not in mdata.mod:
        print("  Skipped: hashing modality not found")
        return

    hashing = mdata.mod["hashing"]
    mtx_dir = outdir / "matrices" / "hto_counts"

    write_mtx_gz(
        matrix=hashing.X,
        barcodes=hashing.obs_names.values,
        features=hashing.var_names.values,
        outdir=mtx_dir,
    )


def export_guide_assignment_table(mdata: md.MuData, outdir: Path) -> None:
    """
    Export guide assignment as a long-format TSV table.

    Format: cell_barcode, guide_id, guide_target, assigned
    Only includes cells with at least one guide assigned.
    """
    print("\nExporting guide assignment table (long format)...")

    if "guide" not in mdata.mod:
        print("  Skipped: guide modality not found")
        return

    guide = mdata.mod["guide"]

    if "guide_assignment" not in guide.layers:
        print("  Skipped: guide_assignment layer not found")
        return

    tables_dir = outdir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    # Get assignment matrix
    assignment = guide.layers["guide_assignment"]
    if issparse(assignment):
        assignment = assignment.tocsr()

    # Get guide metadata
    guide_ids = guide.var["guide_id"].astype(str).values if "guide_id" in guide.var.columns else guide.var_names.values
    guide_targets = (
        guide.var["intended_target_name"].astype(str).values
        if "intended_target_name" in guide.var.columns
        else np.full(len(guide_ids), "")
    )
    cell_barcodes = guide.obs_names.values

    # Build long-format table (only assigned cells)
    rows = []
    if issparse(assignment):
        # Efficient iteration for sparse matrix
        coo = assignment.tocoo()
        for cell_idx, guide_idx in zip(coo.row, coo.col):
            if coo.data[coo.row == cell_idx][coo.col[coo.row == cell_idx] == guide_idx][0] > 0:
                rows.append(
                    {
                        "cell_barcode": cell_barcodes[cell_idx],
                        "guide_id": guide_ids[guide_idx],
                        "intended_target_name": guide_targets[guide_idx],
                    }
                )
        # Simpler approach for sparse
        rows = []
        cx = assignment.tocoo()
        for i, j, v in zip(cx.row, cx.col, cx.data):
            if v > 0:
                rows.append(
                    {
                        "cell_barcode": cell_barcodes[i],
                        "guide_id": guide_ids[j],
                        "intended_target_name": guide_targets[j],
                    }
                )
    else:
        # Dense matrix
        cell_indices, guide_indices = np.where(assignment > 0)
        for cell_idx, guide_idx in zip(cell_indices, guide_indices):
            rows.append(
                {
                    "cell_barcode": cell_barcodes[cell_idx],
                    "guide_id": guide_ids[guide_idx],
                    "intended_target_name": guide_targets[guide_idx],
                }
            )

    df = pd.DataFrame(rows)
    df = df.set_index("cell_barcode")
    write_tsv_gz(df, tables_dir / "guide_assignment.tsv.gz")


def format_mudata(
    mudata_path: str,
    outdir: str,
    export_inference: bool = True,
    export_metadata: bool = True,
    export_matrices: bool = True,
    export_tables: bool = True,
) -> None:
    """
    Export MuData components to common compressed formats.

    Parameters
    ----------
    mudata_path
        Path to input .h5mu file.
    outdir
        Output directory for exported files.
    export_inference
        Export inference results (cis/trans per element/guide).
    export_metadata
        Export cell, gene, and guide metadata.
    export_matrices
        Export count matrices in MTX format.
    export_tables
        Export guide assignment table in long format.
    """
    print(f"Loading MuData from: {mudata_path}")
    mdata = md.read_h5mu(mudata_path)
    print(f"  {mdata.n_obs:,} cells × {mdata.n_vars:,} features")
    print(f"  Modalities: {list(mdata.mod.keys())}")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"\nOutput directory: {outdir}")

    if export_inference:
        export_inference_results(mdata, outdir)

    if export_metadata:
        export_cell_metadata(mdata, outdir)
        export_gene_metadata(mdata, outdir)
        export_guide_metadata(mdata, outdir)

    if export_matrices:
        export_gene_matrix(mdata, outdir)
        export_guide_assignment_matrix(mdata, outdir)
        export_hto_matrix(mdata, outdir)

    if export_tables:
        export_guide_assignment_table(mdata, outdir)

    print("\nDone!")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Export MuData (.h5mu) components to common compressed formats.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Output structure:
  <outdir>/
  ├── inference_results/
  │   ├── cis_per_element_results.tsv.gz
  │   ├── cis_per_guide_results.tsv.gz
  │   ├── trans_per_element_results.tsv.gz
  │   └── trans_per_guide_results.tsv.gz
  ├── metadata/
  │   ├── cell_metadata.tsv.gz
  │   ├── gene_metadata.tsv.gz
  │   └── guide_metadata.tsv.gz
  ├── matrices/
  │   ├── gene_counts/
  │   │   ├── matrix.mtx.gz
  │   │   ├── barcodes.tsv.gz
  │   │   └── features.tsv.gz
  │   ├── guide_assignment/
  │   │   └── ...
  │   └── hto_counts/
  │       └── ...
  └── tables/
      └── guide_assignment.tsv.gz
        """,
    )
    p.add_argument("--mudata", required=True, help="Path to input .h5mu file (MuData).")
    p.add_argument("--outdir", required=True, help="Output directory for exported files.")
    p.add_argument(
        "--no-inference",
        action="store_true",
        help="Skip exporting inference results.",
    )
    p.add_argument(
        "--no-metadata",
        action="store_true",
        help="Skip exporting metadata tables.",
    )
    p.add_argument(
        "--no-matrices",
        action="store_true",
        help="Skip exporting count matrices.",
    )
    p.add_argument(
        "--no-tables",
        action="store_true",
        help="Skip exporting guide assignment table.",
    )
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    format_mudata(
        mudata_path=args.mudata,
        outdir=args.outdir,
        export_inference=not args.no_inference,
        export_metadata=not args.no_metadata,
        export_matrices=not args.no_matrices,
        export_tables=not args.no_tables,
    )
