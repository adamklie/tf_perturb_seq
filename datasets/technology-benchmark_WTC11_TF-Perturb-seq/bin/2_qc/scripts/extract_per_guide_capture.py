#!/usr/bin/env python
"""
Extract per-guide capture metrics from MuData files.

For each guide, calculates:
- n_cells: Number of cells the guide is detected in (UMI > 0)
- n_cells_assigned: Number of cells where this guide is assigned (if assignment available)
- total_umi: Total UMI counts across all cells
- mean_umi: Mean UMI per cell (among cells with guide)
- median_umi: Median UMI per cell (among cells with guide)

Usage:
    python extract_per_guide_capture.py <mudata_path> <output_path> <dataset_name>

Example:
    python extract_per_guide_capture.py \
        /path/to/inference_mudata.h5mu \
        /path/to/output/per_guide_capture.tsv \
        Hon_WTC11-benchmark_TF-Perturb-seq
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import mudata as mu
from scipy import sparse


def extract_per_guide_capture(mudata_path: str, output_path: str, dataset_name: str):
    """Extract per-guide capture metrics from MuData."""

    print(f"Loading MuData from: {mudata_path}")
    mdata = mu.read(mudata_path)

    # Get guide modality (try common names)
    guide_mod = None
    for mod_name in ['guide', 'gdo', 'sgRNA', 'crispr']:
        if mod_name in mdata.mod:
            guide_mod = mdata.mod[mod_name]
            print(f"Found guide modality: '{mod_name}'")
            break

    if guide_mod is None:
        print(f"Available modalities: {list(mdata.mod.keys())}")
        raise ValueError("Could not find guide modality")

    # Get the count matrix (guides x cells after transpose, or cells x guides)
    X = guide_mod.X
    if sparse.issparse(X):
        X = X.toarray()

    # Ensure we have cells x guides orientation
    n_cells, n_guides = X.shape
    print(f"Matrix shape: {n_cells} cells x {n_guides} guides")

    # Get guide names
    guide_names = guide_mod.var_names.tolist()

    # Calculate per-guide metrics
    results = []
    for i, guide_id in enumerate(guide_names):
        guide_counts = X[:, i]

        # Cells with any UMI for this guide
        cells_with_guide = guide_counts > 0
        n_cells_detected = cells_with_guide.sum()

        # UMI statistics (among cells with guide)
        if n_cells_detected > 0:
            umi_values = guide_counts[cells_with_guide]
            total_umi = umi_values.sum()
            mean_umi = umi_values.mean()
            median_umi = np.median(umi_values)
            std_umi = umi_values.std()
            max_umi = umi_values.max()
        else:
            total_umi = 0
            mean_umi = np.nan
            median_umi = np.nan
            std_umi = np.nan
            max_umi = 0

        results.append({
            'guide_id': guide_id,
            'n_cells_detected': int(n_cells_detected),
            'total_umi': float(total_umi),
            'mean_umi': float(mean_umi) if not np.isnan(mean_umi) else np.nan,
            'median_umi': float(median_umi) if not np.isnan(median_umi) else np.nan,
            'std_umi': float(std_umi) if not np.isnan(std_umi) else np.nan,
            'max_umi': float(max_umi),
            'frac_cells_detected': n_cells_detected / n_cells,
        })

    # Create DataFrame
    df = pd.DataFrame(results)
    df['dataset'] = dataset_name
    df['total_cells'] = n_cells

    # Try to get guide metadata from var
    if 'gene_name' in guide_mod.var.columns:
        gene_map = guide_mod.var['gene_name'].to_dict()
        df['gene_name'] = df['guide_id'].map(gene_map)

    if 'label' in guide_mod.var.columns:
        label_map = guide_mod.var['label'].to_dict()
        df['label'] = df['guide_id'].map(label_map)

    # Check for assignment information in obs
    assignment_col = None
    for col in ['assigned_guide', 'guide_assignment', 'perturbation']:
        if col in guide_mod.obs.columns:
            assignment_col = col
            break

    if assignment_col:
        print(f"Found assignment column: '{assignment_col}'")
        # Count cells assigned to each guide
        assignment_counts = guide_mod.obs[assignment_col].value_counts()
        df['n_cells_assigned'] = df['guide_id'].map(assignment_counts).fillna(0).astype(int)
    else:
        df['n_cells_assigned'] = np.nan

    # Reorder columns
    col_order = ['guide_id', 'gene_name', 'label', 'n_cells_detected', 'n_cells_assigned',
                 'frac_cells_detected', 'total_umi', 'mean_umi', 'median_umi', 'std_umi',
                 'max_umi', 'total_cells', 'dataset']
    col_order = [c for c in col_order if c in df.columns]
    df = df[col_order]

    # Save
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)
    print(f"Saved {len(df)} guide metrics to: {output_path}")

    # Print summary
    print(f"\nSummary:")
    print(f"  Total guides: {len(df)}")
    print(f"  Guides detected in >0 cells: {(df['n_cells_detected'] > 0).sum()}")
    print(f"  Median cells per guide: {df['n_cells_detected'].median():.0f}")
    print(f"  Median UMI per guide: {df['total_umi'].median():.0f}")

    return df


def main():
    parser = argparse.ArgumentParser(description='Extract per-guide capture metrics from MuData')
    parser.add_argument('mudata_path', help='Path to MuData file')
    parser.add_argument('output_path', help='Output TSV path')
    parser.add_argument('dataset_name', help='Dataset name for labeling')

    args = parser.parse_args()

    extract_per_guide_capture(args.mudata_path, args.output_path, args.dataset_name)


if __name__ == '__main__':
    main()
