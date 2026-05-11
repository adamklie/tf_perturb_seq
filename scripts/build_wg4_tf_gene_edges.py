"""Build the WG4-A per-dataset TF→gene network edge list (FDR-filtered).

Reads the CRISPR pipeline's `perturbo_trans_per_element_output.tsv.gz`, computes
per-TF Benjamini-Hochberg FDR (treating each TF's genome-wide tests as a family),
filters to FDR < threshold (default 0.05), and writes a slim edge list with
TF identity (gene_symbol + family) joined in.

Output: `datasets/<dataset>/crispr_pipeline/wg4_tf_gene_edges_FDR<thr>.tsv`

Per-TF BH is the lenient (and standard) approach for genome-wide perturb-seq
trans tests — vs. one big BH across all (TF × gene) tests, which is harsh and
biases against TFs with strong polygenic effects.

Usage:
    python scripts/build_wg4_tf_gene_edges.py \\
        --source-dir <DIR with perturbo_trans_per_element_output.tsv.gz> \\
        --output-dir <DIR for wg4_*.tsv> \\
        --tf-metadata-path docs/jamborees/2026_UTSW/reference/tf_metadata.tsv \\
        --fdr-threshold 0.05
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment (monotonic step-up)."""
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p, kind="mergesort")
    ranked = np.empty(n, dtype=int)
    ranked[order] = np.arange(1, n + 1)
    adj = p * n / ranked
    # Enforce monotonicity from the largest p-value downward
    sorted_adj = adj[order]
    for i in range(n - 2, -1, -1):
        sorted_adj[i] = min(sorted_adj[i], sorted_adj[i + 1])
    sorted_adj = np.minimum(sorted_adj, 1.0)
    out = np.empty(n)
    out[order] = sorted_adj
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source-dir", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--tf-metadata-path", required=True, type=Path)
    ap.add_argument("--fdr-threshold", type=float, default=0.05)
    args = ap.parse_args()

    src = args.source_dir / "perturbo_trans_per_element_output.tsv.gz"
    if not src.is_file():
        sys.exit(f"missing: {src}")

    print(f"[read] {src}")
    df = pd.read_csv(src, sep="\t", compression="gzip")
    print(f"  loaded {len(df):,} rows × {len(df.columns)} cols")

    expected = {"gene_id", "intended_target_name", "log2_fc", "p_value"}
    missing = expected - set(df.columns)
    if missing:
        sys.exit(f"missing columns in perturbo_trans: {missing}")

    # Per-TF BH adjustment
    print("[bh] computing FDR per TF (genome-wide tests per intended_target_name as a family)…")
    df["fdr_bh"] = (
        df.groupby("intended_target_name")["p_value"]
        .transform(lambda s: pd.Series(bh_adjust(s.values), index=s.index))
    )

    # Filter
    sig = df[df["fdr_bh"] < args.fdr_threshold].copy()
    print(f"[filter] {len(sig):,} edges at FDR < {args.fdr_threshold}")
    print(f"  unique TFs with ≥1 sig edge: {sig['intended_target_name'].nunique()}")
    print(f"  median edges per significant TF: {sig.groupby('intended_target_name').size().median():.0f}")

    # Join TF identity
    tf = pd.read_csv(args.tf_metadata_path, sep="\t")[
        ["ensembl_gene_id", "gene_symbol", "jaspar_tf_family", "lambert_2018_dbd"]
    ].rename(
        columns={
            "ensembl_gene_id": "intended_target_name",
            "gene_symbol": "tf_gene_symbol",
            "jaspar_tf_family": "tf_family",
            "lambert_2018_dbd": "tf_dbd",
        }
    )
    out = sig.merge(tf, on="intended_target_name", how="left")

    # Column order
    out = out[
        [
            "intended_target_name",
            "tf_gene_symbol",
            "tf_family",
            "tf_dbd",
            "gene_id",
            "log2_fc",
            "log2_fc_std",
            "p_value",
            "fdr_bh",
            "intended_target_chr",
            "intended_target_start",
            "intended_target_end",
        ]
    ].sort_values(["intended_target_name", "fdr_bh"])

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_name = f"wg4_tf_gene_edges_FDR{args.fdr_threshold:g}".replace(".", "p") + ".tsv"
    # Cleaner: use a fixed name. FDR threshold goes in a sidecar / row column.
    out_name = f"wg4_tf_gene_edges_FDR{int(args.fdr_threshold*100):02d}.tsv"
    out_path = args.output_dir / out_name
    out.to_csv(out_path, sep="\t", index=False)
    print(f"[write] {out_path} ({len(out):,} rows × {len(out.columns)} cols)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
