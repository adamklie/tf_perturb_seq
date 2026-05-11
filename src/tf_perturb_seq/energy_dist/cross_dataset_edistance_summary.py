"""Build a cross-dataset summary TSV of energy-distance results.

Reads `synapse_paths.tsv` to find each dataset's energy-distance Synapse folder,
downloads `pval_edist_full.csv`, and aggregates per-dataset summary stats:
  - n targets (total + by type)
  - distance_mean: median (overall, targeting, NC, positive control)
  - n targets with distance_mean > NC max  (calibration-robust significance proxy)
  - n targets with pval_mean=0  (raw — affected by the known Huangfu calibration issue)
  - n targets significant at pval_mean < 0.05

The "distance > NC max" column is the recommended significance proxy until the
p-value calibration concern (see issues/edistance-calibration.md) is fixed —
because for the Huangfu runs all NCs end up at pval_mean=0, the raw p-value
threshold is uninformative.

Designed to be run from anywhere. Downloads only `pval_edist_full.csv` per
dataset (~1 MB each).

Authentication:
  - Synapse: SYNAPSE_AUTH_TOKEN env var. Run via `zsh -ic '...'` so the token
    from ~/.zshrc is loaded.

Usage:
    python cross_dataset_edistance_summary.py
    python cross_dataset_edistance_summary.py --output reference/cross_dataset_edistance_summary.tsv
    python cross_dataset_edistance_summary.py --cache-dir /tmp/jamboree_ed_cache
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Optional

import pandas as pd
import synapseclient

JAMB = Path(__file__).resolve().parent.parent / "docs/jamborees/2026_UTSW"
SYNAPSE_PATHS = JAMB / "synapse_paths.tsv"
DEFAULT_OUTPUT = JAMB / "reference/cross_dataset_edistance_summary.tsv"
DEFAULT_CACHE = Path("/tmp/jamboree_edistance_cache")


def find_child(syn, parent_id: str, name: str) -> Optional[str]:
    for c in syn.getChildren(parent_id):
        if c["name"] == name:
            return c["id"]
    return None


def download_pval_edist(syn, ed_parent_id: str, dataset: str, cache_dir: Path) -> Optional[Path]:
    dest = cache_dir / dataset
    dest.mkdir(parents=True, exist_ok=True)
    file_id = find_child(syn, ed_parent_id, "pval_edist_full.csv")
    if file_id is None:
        print(f"  warn: pval_edist_full.csv not found under {ed_parent_id}")
        return None
    try:
        f = syn.get(file_id, downloadLocation=str(dest), ifcollision="overwrite.local")
        return Path(f.path)
    except Exception as e:
        print(f"  warn: download failed for {file_id}: {e!r}")
        return None


def extract_summary(dataset: str, df: pd.DataFrame) -> dict:
    row: dict = {"dataset_id": dataset}

    # Total targets + per-type counts
    row["n_targets_total"] = len(df)
    if "type" in df.columns:
        for t in ("targeting", "negative control", "positive control"):
            row[f"n_targets_{t.replace(' ', '_')}"] = int((df["type"] == t).sum())

    # distance_mean stats
    if "distance_mean" in df.columns:
        row["distance_mean_median_all"] = float(df["distance_mean"].median())
        if "type" in df.columns:
            for t in ("targeting", "negative control", "positive control"):
                sub = df[df["type"] == t]["distance_mean"]
                if not sub.empty:
                    row[f"distance_mean_median_{t.replace(' ', '_')}"] = float(sub.median())
        # Calibration-robust significance proxy: distance > NC max
        if "type" in df.columns:
            nc = df[df["type"] == "negative control"]["distance_mean"]
            if not nc.empty:
                nc_max = float(nc.max())
                row["nc_distance_mean_max"] = nc_max
                row["n_targeting_above_nc_max"] = int(
                    ((df["type"] == "targeting") & (df["distance_mean"] > nc_max)).sum()
                )

    # pval_mean stats (calibration concern — see issue 1)
    if "pval_mean" in df.columns:
        row["n_pval_eq_0"] = int((df["pval_mean"] == 0).sum())
        row["n_pval_lt_0p05"] = int((df["pval_mean"] < 0.05).sum())
        if "type" in df.columns:
            row["n_nc_pval_eq_0"] = int(
                ((df["type"] == "negative control") & (df["pval_mean"] == 0)).sum()
            )

    # Cell-count diagnostics
    if "cell_count" in df.columns:
        row["cell_count_median"] = float(df["cell_count"].median())
        row["cell_count_min"] = int(df["cell_count"].min())

    return row


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    ap.add_argument("--datasets", nargs="*", help="Limit to specific dataset_id(s)")
    ap.add_argument("--synapse-paths", type=Path, default=SYNAPSE_PATHS)
    args = ap.parse_args()

    if not args.synapse_paths.is_file():
        sys.exit(f"synapse_paths.tsv not found: {args.synapse_paths}")
    sp = pd.read_csv(args.synapse_paths, sep="\t", dtype=str).fillna("")
    if "energy_distance" not in sp.columns:
        sys.exit("synapse_paths.tsv has no energy_distance column")

    eligible = sp[sp["energy_distance"].str.startswith("syn", na=False)]
    if args.datasets:
        eligible = eligible[eligible["dataset_id"].isin(args.datasets)]
    if eligible.empty:
        sys.exit("no eligible datasets (none with energy_distance = synXXX)")

    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        sys.exit("SYNAPSE_AUTH_TOKEN env var not set. Run via `zsh -ic` so ~/.zshrc loads.")
    syn = synapseclient.Synapse()
    syn.login(authToken=token, silent=True)

    rows: list[dict] = []
    for _, ds_row in eligible.iterrows():
        dataset = ds_row["dataset_id"]
        ed_id = ds_row["energy_distance"]
        print(f"\n=== {dataset}  ({ed_id}) ===")
        local = download_pval_edist(syn, ed_id, dataset, args.cache_dir)
        if local is None or not local.is_file():
            print("  warn: skipping")
            continue
        try:
            df = pd.read_csv(local, index_col=0)
        except Exception as e:
            print(f"  warn: parse failed: {e!r}")
            continue
        rows.append(extract_summary(dataset, df))

    if not rows:
        sys.exit("no rows produced")

    out_df = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"\nwrote {args.output} ({len(out_df)} rows × {len(out_df.columns)} cols)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
