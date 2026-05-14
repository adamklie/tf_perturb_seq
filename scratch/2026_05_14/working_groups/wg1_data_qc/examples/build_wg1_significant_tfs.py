"""Build the per-dataset WG1 significant-TFs table.

For a single dataset's energy-distance run, joins:
  - pval_edist_full.csv (per-target row: distance + permutation pvals)
  - targeting_outlier_table.csv (per-guide outlier pvals; aggregated to per-target "any outlier" flag)
  - tf_metadata.tsv (Ensembl ID + gene_symbol + tf_family from JASPAR)

…and writes a single per-TF TSV with both significance criteria attached, ranked
by `distance_mean` (descending). Sized as a slide-friendly view of the run.

Reads pval_edist_full.csv + targeting_outlier_table.csv from --source-dir (typical
HPC path: `/cellar/.../results/energy_distance/<run_label>/`).

Output path is computed from --output-dir (typical:
`docs/jamborees/2026_UTSW/datasets/<dataset>/energy_distance/`).

Usage:
    python working_groups/wg1_data_qc/examples/build_wg1_significant_tfs.py \\
        --source-dir /cellar/.../Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/results/energy_distance/muddy_penguin \\
        --output-dir /cellar/.../docs/jamborees/2026_UTSW/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/energy_distance \\
        --tf-metadata-path docs/jamborees/2026_UTSW/reference/tf_metadata.tsv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


def parse_target_id(target_id: str) -> tuple[str | None, str | None]:
    """Index format: `<ENSG>|<chr>:<start>-<end>`. Extract (ensembl_gene_id, locus)."""
    if "|" not in target_id:
        return None, target_id
    ensg, locus = target_id.split("|", 1)
    if not re.match(r"^ENSG\d+", ensg):
        return None, target_id
    return ensg, locus


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source-dir", required=True, type=Path, help="ED run output dir (contains pval_edist_full.csv etc.)")
    ap.add_argument("--output-dir", required=True, type=Path, help="Where to write wg1_significant_tfs.tsv")
    ap.add_argument("--tf-metadata-path", required=True, type=Path, help="Path to reference/tf_metadata.tsv")
    args = ap.parse_args()

    pval_path = args.source_dir / "pval_edist_full.csv"
    outl_path = args.source_dir / "targeting_outlier_table.csv"
    tf_path = args.tf_metadata_path

    if not pval_path.is_file():
        sys.exit(f"missing: {pval_path}")
    if not tf_path.is_file():
        sys.exit(f"missing: {tf_path}")

    # Load + clean pval_edist_full
    pv = pd.read_csv(pval_path, index_col=0)
    pv.index.name = "target_id"
    pv = pv.reset_index()
    pv[["ensembl_gene_id", "locus"]] = pv["target_id"].apply(
        lambda s: pd.Series(parse_target_id(str(s)))
    )

    # Load TF metadata + join on ensembl_gene_id
    tf = pd.read_csv(tf_path, sep="\t")[
        ["ensembl_gene_id", "gene_symbol", "hgnc_approved_symbol", "jaspar_tf_family", "lambert_2018_dbd"]
    ]
    df = pv.merge(tf, on="ensembl_gene_id", how="left")

    # Outlier flag: per-target "any outlier guide" boolean.
    # targeting_outlier_table.csv is per-gRNA. We don't have a direct gRNA→target
    # map at this layer of the pipeline, so flag at the dataset level only:
    # was_outlier_any = True if ANY targeting guide in the run was flagged as
    # an outlier (low pval_outlier means more outlier-like; threshold p<0.05).
    if outl_path.is_file():
        outl = pd.read_csv(outl_path, index_col=0)
        n_outliers = (outl["pval_outlier"] < 0.05).sum()
        n_grnas = len(outl)
        # Per-row annotation: how many outlier guides in the run (constant across rows)
        df["n_outlier_grnas_in_run"] = int(n_outliers)
        df["n_grnas_in_run"] = int(n_grnas)
    else:
        df["n_outlier_grnas_in_run"] = pd.NA
        df["n_grnas_in_run"] = pd.NA

    # NC max for calibration-robust significance flag
    nc_max = (
        df.loc[df["type"] == "negative control", "distance_mean"].max()
        if (df["type"] == "negative control").any()
        else float("nan")
    )

    df["sig_distance_gt_NC_max"] = (
        df["distance_mean"] > nc_max if pd.notna(nc_max) else pd.NA
    )
    df["sig_pval_lt_0p05"] = df["pval_mean"] < 0.05

    # Rank within targeting subset
    targeting = df[df["type"] == "targeting"].copy()
    targeting["distance_rank_targeting"] = targeting["distance_mean"].rank(ascending=False, method="min").astype(int)
    df = df.merge(
        targeting[["target_id", "distance_rank_targeting"]],
        on="target_id",
        how="left",
    )

    # Curated output columns
    out_cols = [
        "target_id",
        "ensembl_gene_id",
        "gene_symbol",
        "hgnc_approved_symbol",
        "jaspar_tf_family",
        "lambert_2018_dbd",
        "type",
        "cell_count",
        "distance_mean",
        "distance_rank_targeting",
        "sig_distance_gt_NC_max",
        "pval_mean",
        "sig_pval_lt_0p05",
        "n_outlier_grnas_in_run",
        "n_grnas_in_run",
        "locus",
    ]
    out = df[out_cols].sort_values("distance_mean", ascending=False, na_position="last")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.output_dir / "wg1_significant_tfs.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    n_sig_dist = int(out["sig_distance_gt_NC_max"].fillna(False).sum())
    n_sig_pval = int(out["sig_pval_lt_0p05"].fillna(False).sum())
    print(f"wrote {out_path} ({len(out)} rows × {len(out.columns)} cols)")
    print(f"  sig_distance_gt_NC_max: {n_sig_dist} | sig_pval_lt_0p05: {n_sig_pval}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
