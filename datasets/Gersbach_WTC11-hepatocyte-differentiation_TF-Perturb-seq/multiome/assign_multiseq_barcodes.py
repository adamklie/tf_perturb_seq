# -*- coding: utf-8 -*-
"""
Assigns each cell barcode to its highest-count multiseq barcode across one or
more runs, translates ATAC barcodes -> GEX barcodes using the 10x Multiome
whitelist, and joins sample metadata.

Single-run usage:
    python assign_multiseq_barcodes.py \
        --runs   IGVFFI1271TKWE.tsv:IGVFSM2174QALS \
        --sample_map  IGVFFI6468JUDR.tsv \
        --atac_wl     IGVFFI7587TJLC.txt.gz \
        --gex_wl      IGVFFI8751YQRY.txt.gz \
        --output      cell_barcode_mapping.tsv

Multi-run usage:
    python assign_multiseq_barcodes.py \
        --runs   IGVFFI1271TKWE.tsv:IGVFSM2174QALS \
                 IGVFFI0278DRTS.tsv:IGVFSM0825LHRL \
                 IGVFFI4981LXWG.tsv:IGVFSM1624SZQN \
        --sample_map  IGVFFI6468JUDR.tsv \
        --atac_wl     IGVFFI7587TJLC.txt.gz \
        --gex_wl      IGVFFI8751YQRY.txt.gz \
        --output      cell_barcode_mapping_all.tsv

Each --runs entry is <matrix_tsv>:<igvfsm_id>. Per-run mapping files are
written alongside the combined output as <output>_<igvfsm_id>.tsv.

Arguments:
    --runs        One or more <matrix_tsv>:<igvfsm_id> pairs.
    --sample_map  Multiseq barcode -> sample metadata TSV:
                      barcode | sample accession | sample description
    --atac_wl     10x Multiome ATAC barcode whitelist (IGVFFI7587TJLC).
    --gex_wl      10x Multiome GEX barcode whitelist (IGVFFI8751YQRY).
    --output      Combined output TSV path.
    --method      Assignment method: max (default), ratio, clr, zscore.
    --ratio_min   Minimum top/second ratio (default: 3.0, method=ratio only).
    --clr_min     Minimum CLR score (default: 1.0, method=clr only).
    --zscore_min  Minimum z-score (default: 2.0, method=zscore only).
"""

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_matrix(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.dropna(how="all")
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    return df


def load_whitelist(path: str) -> list:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        return [line.strip() for line in f if line.strip()]


def load_sample_map(path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep="\t")
        if df.shape[1] < 3:
            raise ValueError("Too few columns")
    except Exception:
        df = pd.read_csv(path, sep=r"\s+", engine="python")

    df.columns = [c.strip() for c in df.columns]
    col_map = {}
    for col in df.columns:
        lc = col.lower()
        if lc == "barcode":
            col_map[col] = "multiseq_barcode"
        elif "accession" in lc:
            col_map[col] = "sample_accession"
        elif "description" in lc:
            col_map[col] = "sample_description"

    if len(col_map) < 3:
        raise ValueError(
            f"Could not identify required columns in {path}.\n"
            f"Found: {list(df.columns)}"
        )
    return df.rename(columns=col_map)[
        ["multiseq_barcode", "sample_accession", "sample_description"]
    ]


# ---------------------------------------------------------------------------
# ATAC -> GEX translation
# ---------------------------------------------------------------------------

def build_atac_to_gex(atac_wl_path: str, gex_wl_path: str) -> dict:
    print(f"  Loading ATAC whitelist: {atac_wl_path}")
    atac_wl = load_whitelist(atac_wl_path)
    print(f"  Loading GEX  whitelist: {gex_wl_path}")
    gex_wl  = load_whitelist(gex_wl_path)
    if len(atac_wl) != len(gex_wl):
        raise ValueError(
            f"Whitelist length mismatch: ATAC={len(atac_wl):,}, GEX={len(gex_wl):,}."
        )
    print(f"  {len(atac_wl):,} barcode pairs loaded")
    return dict(zip(atac_wl, gex_wl))


# ---------------------------------------------------------------------------
# Barcode helpers
# ---------------------------------------------------------------------------

def strip_col_prefix(col: str) -> str:
    return col.rsplit("_", 1)[-1] if "_" in col else col


def extract_core_atac(cell_bc: str) -> str:
    core = cell_bc.rsplit("_", 1)[-1]
    core = core.rsplit("-", 1)[0]
    return core


def build_h5ad_name(gex_bc: str, igvfsm_id: str) -> str:
    return f"{gex_bc}_{igvfsm_id}"


# ---------------------------------------------------------------------------
# Assignment methods
# ---------------------------------------------------------------------------

def assign_max(matrix: pd.DataFrame) -> pd.Series:
    return matrix.idxmax(axis=1).apply(strip_col_prefix)


def assign_ratio(matrix: pd.DataFrame, ratio_min: float) -> pd.Series:
    sorted_vals = np.sort(matrix.values, axis=1)
    top    = sorted_vals[:, -1].astype(float)
    second = sorted_vals[:, -2].astype(float)
    ratio  = np.where(second > 0, top / second, np.inf)
    result = matrix.idxmax(axis=1).apply(strip_col_prefix)
    result[top == 0]          = "Unassigned"
    result[ratio < ratio_min] = "Doublet"
    return result


def clr_transform(matrix: pd.DataFrame) -> pd.DataFrame:
    X      = matrix.values.astype(float) + 1.0
    log_X  = np.log(X)
    return pd.DataFrame(log_X - log_X.mean(axis=1, keepdims=True),
                        index=matrix.index, columns=matrix.columns)


def assign_clr(matrix: pd.DataFrame, clr_min: float) -> pd.Series:
    clr    = clr_transform(matrix)
    result = clr.idxmax(axis=1).apply(strip_col_prefix)
    result[clr.max(axis=1) < clr_min] = "Unassigned"
    return result


def assign_zscore(matrix: pd.DataFrame, zscore_min: float) -> pd.Series:
    X   = matrix.values.astype(float)
    mu  = X.mean(axis=1, keepdims=True)
    std = X.std(axis=1, keepdims=True)
    std[std == 0] = 1.0
    Z      = pd.DataFrame((X - mu) / std, index=matrix.index, columns=matrix.columns)
    result = Z.idxmax(axis=1).apply(strip_col_prefix)
    result[Z.max(axis=1) < zscore_min] = "Unassigned"
    return result


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def print_diagnostics(matrix: pd.DataFrame, assigned: pd.Series, run_label: str) -> None:
    sorted_vals = np.sort(matrix.values, axis=1)
    top    = sorted_vals[:, -1].astype(float)
    second = sorted_vals[:, -2].astype(float)
    ratio  = np.where(second > 0, top / second, np.inf)

    clr     = clr_transform(matrix)
    clr_top = clr.max(axis=1).values

    X     = matrix.values.astype(float)
    mu    = X.mean(axis=1, keepdims=True)
    std   = X.std(axis=1, keepdims=True)
    std[std == 0] = 1.0
    z_top = ((X - mu) / std).max(axis=1)

    print(f"\n{'='*60}")
    print(f"DIAGNOSTIC: {run_label}")
    print(f"{'='*60}")

    finite_ratio = ratio[np.isfinite(ratio)]
    print(f"\n  Top/second ratio (n={len(finite_ratio):,} cells with second > 0)")
    for pct, val in zip([5, 25, 50, 75, 95],
                        np.percentile(finite_ratio, [5, 25, 50, 75, 95])):
        print(f"    p{pct:02d}: {val:6.2f}")
    for thresh in [2, 3, 5, 10]:
        print(f"    >= {thresh:2d}x : {(finite_ratio >= thresh).mean()*100:5.1f}% of cells")

    print(f"\n  CLR score of top barcode")
    for pct, val in zip([5, 25, 50, 75, 95],
                        np.percentile(clr_top, [5, 25, 50, 75, 95])):
        print(f"    p{pct:02d}: {val:6.2f}")

    print(f"\n  Z-score of top barcode")
    for pct, val in zip([5, 25, 50, 75, 95],
                        np.percentile(z_top, [5, 25, 50, 75, 95])):
        print(f"    p{pct:02d}: {val:6.2f}")

    print(f"\n  Cells with zero counts: {(top == 0).sum():,}")
    print(f"  Method: {assigned.name}")
    vc = assigned.value_counts()
    for label in ["Doublet", "Unassigned"]:
        if label in vc:
            print(f"    {label}: {vc[label]:,} ({vc[label]/len(assigned)*100:.1f}%)")
    print(f"    Assigned: {(~assigned.isin(['Doublet','Unassigned'])).sum():,}")
    print(f"{'='*60}\n")


# ---------------------------------------------------------------------------
# Per-run processing
# ---------------------------------------------------------------------------

def process_run(matrix_path: str, igvfsm_id: str, atac_to_gex: dict,
                sample_map: pd.DataFrame, method: str,
                ratio_min: float, clr_min: float, zscore_min: float,
                out_path: str) -> pd.DataFrame:

    print(f"\n  Matrix: {matrix_path}  |  IGVFSM: {igvfsm_id}")
    matrix = load_matrix(matrix_path)
    print(f"  {matrix.shape[0]:,} cells x {matrix.shape[1]:,} multiseq barcodes")

    if method == "max":
        assigned = assign_max(matrix)
    elif method == "ratio":
        assigned = assign_ratio(matrix, ratio_min)
    elif method == "clr":
        assigned = assign_clr(matrix, clr_min)
    elif method == "zscore":
        assigned = assign_zscore(matrix, zscore_min)
    assigned.name = method

    print_diagnostics(matrix, assigned, run_label=f"{igvfsm_id} ({method})")

    result = assigned.reset_index()
    result.columns = ["cellBC", "multiseq_barcode"]
    result["igvfsm_id"]    = igvfsm_id
    result["atac_barcode"] = result["cellBC"].apply(extract_core_atac)
    result["gex_barcode"]  = result["atac_barcode"].map(atac_to_gex)
    result["h5ad_barcode"] = result["gex_barcode"].apply(
        lambda g: build_h5ad_name(g, igvfsm_id) if pd.notna(g) else None
    )

    n_total       = len(result)
    n_untranslated = result["gex_barcode"].isna().sum()
    print(f"  Translated: {n_total - n_untranslated:,}/{n_total:,} "
          f"({(n_total - n_untranslated)/n_total*100:.1f}%)")
    if n_untranslated:
        print(f"  WARNING: {n_untranslated:,} ATAC barcodes not in whitelist",
              file=sys.stderr)

    result = result.merge(sample_map, on="multiseq_barcode", how="left")

    n_unmapped = result["sample_accession"].isna().sum()
    n_special  = result["multiseq_barcode"].isin(["Doublet", "Unassigned"]).sum()
    if n_unmapped - n_special > 0:
        print(f"  WARNING: {n_unmapped - n_special:,} cells had unrecognised "
              f"multiseq barcode", file=sys.stderr)

    result = result[[
        "cellBC", "igvfsm_id", "atac_barcode", "gex_barcode", "h5ad_barcode",
        "multiseq_barcode", "sample_accession", "sample_description"
    ]]
    result.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved per-run mapping: {out_path}")
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--runs", required=True, nargs="+",
                        metavar="MATRIX:IGVFSM_ID",
                        help="One or more <matrix_tsv>:<igvfsm_id> pairs")
    parser.add_argument("--sample_map", required=True)
    parser.add_argument("--atac_wl",    required=True)
    parser.add_argument("--gex_wl",     required=True)
    parser.add_argument("--output",     default="cell_barcode_mapping_all.tsv")
    parser.add_argument("--method",     default="max",
                        choices=["max", "ratio", "clr", "zscore"])
    parser.add_argument("--ratio_min",  type=float, default=3.0)
    parser.add_argument("--clr_min",    type=float, default=1.0)
    parser.add_argument("--zscore_min", type=float, default=2.0)
    return parser.parse_args()


def main():
    args = parse_args()

    # Parse run specs
    runs = []
    for spec in args.runs:
        parts = spec.rsplit(":", 1)
        if len(parts) != 2:
            print(f"ERROR: --runs entries must be <matrix_tsv>:<igvfsm_id>, got: {spec}",
                  file=sys.stderr)
            sys.exit(1)
        runs.append((parts[0], parts[1]))

    print(f"[1/4] Building ATAC->GEX translation table")
    atac_to_gex = build_atac_to_gex(args.atac_wl, args.gex_wl)

    print(f"\n[2/4] Loading sample map: {args.sample_map}")
    sample_map = load_sample_map(args.sample_map)
    print(f"      {len(sample_map):,} multiseq barcodes")

    print(f"\n[3/4] Processing {len(runs)} run(s)")
    all_results = []
    out_stem = args.output.replace(".tsv", "")
    for matrix_path, igvfsm_id in runs:
        per_run_path = f"{out_stem}_{igvfsm_id}.tsv"
        df = process_run(
            matrix_path, igvfsm_id, atac_to_gex, sample_map,
            args.method, args.ratio_min, args.clr_min, args.zscore_min,
            per_run_path
        )
        all_results.append(df)

    print(f"\n[4/4] Writing combined mapping: {args.output}")
    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(args.output, sep="\t", index=False)
    print(f"      {len(combined):,} total cells across {len(runs)} run(s)")
    print(f"\n  Run summary:")
    for igvfsm_id, grp in combined.groupby("igvfsm_id"):
        n_assigned = (~grp["multiseq_barcode"].isin(["Doublet", "Unassigned"])).sum()
        print(f"    {igvfsm_id}: {len(grp):,} cells, {n_assigned:,} assigned")
    print("\nDone.")


if __name__ == "__main__":
    main()