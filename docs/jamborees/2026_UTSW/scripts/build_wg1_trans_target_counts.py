"""Build the WG1-E per-perturbation trans-target count summary.

Groups the WG4-A FDR-filtered edge list by perturbation (intended_target_name)
and reports n_significant_trans_targets per TF, plus the top-5 strongest target
genes (by |log2_fc|) for quick scanning.

Output: `datasets/<dataset>/crispr_pipeline/wg1_trans_target_counts.tsv`

Usage:
    python scripts/build_wg1_trans_target_counts.py \\
        --edges-path datasets/<id>/crispr_pipeline/wg4_tf_gene_edges_FDR05.tsv \\
        --output-dir datasets/<id>/crispr_pipeline/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--edges-path", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--top-n", type=int, default=5)
    args = ap.parse_args()

    if not args.edges_path.is_file():
        sys.exit(f"missing: {args.edges_path}")
    df = pd.read_csv(args.edges_path, sep="\t")

    def top_targets(group: pd.DataFrame) -> str:
        top = group.assign(abs_fc=group["log2_fc"].abs()).sort_values("abs_fc", ascending=False).head(args.top_n)
        return ";".join(top["gene_id"].astype(str))

    summary = (
        df.groupby("intended_target_name")
        .agg(
            tf_gene_symbol=("tf_gene_symbol", "first"),
            tf_family=("tf_family", "first"),
            n_sig_trans_targets=("gene_id", "size"),
            median_abs_log2fc=("log2_fc", lambda s: float(s.abs().median())),
            median_log2fc_signed=("log2_fc", "median"),
            min_fdr=("fdr_bh", "min"),
            top_targets=("gene_id", lambda s: ";".join(
                df.loc[s.index].assign(abs_fc=df.loc[s.index, "log2_fc"].abs())
                .sort_values("abs_fc", ascending=False)
                .head(args.top_n)["gene_id"].astype(str)
            )),
        )
        .reset_index()
        .sort_values("n_sig_trans_targets", ascending=False)
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.output_dir / "wg1_trans_target_counts.tsv"
    summary.to_csv(out_path, sep="\t", index=False)
    print(f"wrote {out_path} ({len(summary):,} rows × {len(summary.columns)} cols)")
    print(f"  TFs with >100 sig trans targets: {(summary['n_sig_trans_targets'] > 100).sum()}")
    print(f"  TFs with >10 sig trans targets: {(summary['n_sig_trans_targets'] > 10).sum()}")
    print(f"  median n_sig_trans_targets: {summary['n_sig_trans_targets'].median():.0f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
