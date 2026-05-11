"""Validate the outputs of an energy distance pipeline run.

Runs four layers of checks against a single dataset's e-distance output dir
(e.g., the per-run folder under datasets/<id>/results/energy_distance/<run>/):

  Layer 1 — file presence and non-emptiness (the deliverable bundle per
            schemas/energy_distance.json + the intermediate files we expect)
  Layer 2 — schema check: columns of the three CSVs match the verified shape
  Layer 3 — value-range sanity: pval_mean in [0,1], cell_count>0, no all-NaN
            rows, expected `type` values
  Layer 4 — cross-reference vs the Gersbach HTv2 verified Synapse example
            (syn74381167 / pval_edist_full = syn74381183) — schema only, not
            values; can be skipped with --no-synapse-ref

Each layer prints PASS/FAIL/WARN per check. Final summary is exit-code 0 if
all REQUIRED checks pass, 1 otherwise. WARNs are non-fatal.

Usage:
    python validate_edistance_outputs.py --source-dir /path/to/.../muddy_penguin
    python validate_edistance_outputs.py --source-dir <DIR> --no-synapse-ref
    python validate_edistance_outputs.py --source-dir <DIR> --download-ref-to /tmp

The Layer 4 cross-reference downloads the HTv2 pval_edist_full.csv from
Synapse on first run (cached at /tmp/htv2_pval_edist_full.csv by default).
Requires SYNAPSE_AUTH_TOKEN in env; if missing, Layer 4 is skipped with WARN.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd

# ----------------- expected file inventory + schemas -----------------

REQUIRED_DELIVERABLES = [
    "pval_edist_full.csv",
    "targeting_outlier_table.csv",
    "non_targeting_outlier_table.csv",
    "config1_2.json",
    "config3.json",
    "image/gRNA_stat.pdf",
    "image/e-dist_distribution.pdf",
    "image/e-dist_cutoff_value.pdf",
    "image/e-dist_cutoff_value_NEG_CONTROL.pdf",
]
REQUIRED_INTERMEDIATES = [
    "preprocessed.h5ad",
    "gRNA_dict.pickle",
    "pca_dataframe.pickle",
    "annotation_table.csv",
]

# pval_edist_full expected columns (per energy_distance.json verified_columns
# and the pipeline's permutation_test config: num_of_bg=20, so distance_0..19,
# pval_0..19, plus mean + log columns)
PVAL_EDIST_REQUIRED_COLS = (
    {"cell_count", "type", "distance_mean", "pval_mean", "pval_mean_log", "distance_mean_log"}
    | {f"distance_{i}" for i in range(20)}
    | {f"pval_{i}" for i in range(20)}
)

OUTLIER_TABLE_REQUIRED_COLS = {"pval_outlier"}

ALLOWED_TYPES = {"targeting", "target", "positive control", "negative control", "non-targeting"}

HTV2_SYNAPSE_ID_PVAL = "syn74381183"  # pval_edist_full.csv from HTv2 verified run

# ----------------- helpers -----------------

ANSI = {
    "PASS": "\033[32mPASS\033[0m",
    "FAIL": "\033[31mFAIL\033[0m",
    "WARN": "\033[33mWARN\033[0m",
    "INFO": "\033[36mINFO\033[0m",
}

class Result:
    def __init__(self):
        self.fails: list[str] = []
        self.warns: list[str] = []

    def report(self, level: str, msg: str) -> None:
        print(f"{ANSI.get(level, level):>4}  {msg}")
        if level == "FAIL":
            self.fails.append(msg)
        elif level == "WARN":
            self.warns.append(msg)

    def summary(self) -> int:
        print()
        if self.fails:
            print(f"{ANSI['FAIL']}  {len(self.fails)} required check(s) failed:")
            for m in self.fails:
                print(f"  - {m}")
        if self.warns:
            print(f"{ANSI['WARN']}  {len(self.warns)} non-fatal warning(s):")
            for m in self.warns:
                print(f"  - {m}")
        if not self.fails and not self.warns:
            print(f"{ANSI['PASS']}  all checks passed")
        return 1 if self.fails else 0


# ----------------- layer 1: presence + non-empty -----------------

def layer1_presence(src: Path, r: Result) -> None:
    print("\n=== Layer 1: file presence + non-emptiness ===")
    for rel in REQUIRED_DELIVERABLES:
        p = src / rel
        if not p.exists():
            r.report("FAIL", f"missing required deliverable: {rel}")
        elif p.stat().st_size == 0:
            r.report("FAIL", f"empty deliverable: {rel}")
        else:
            r.report("PASS", f"{rel} ({_human_size(p)})")
    for rel in REQUIRED_INTERMEDIATES:
        p = src / rel
        if not p.exists():
            r.report("WARN", f"missing intermediate (not strictly required): {rel}")
        elif p.stat().st_size == 0:
            r.report("WARN", f"empty intermediate: {rel}")
        else:
            r.report("PASS", f"{rel} ({_human_size(p)})")


def _human_size(p: Path) -> str:
    n = p.stat().st_size
    for unit in ("B", "K", "M", "G"):
        if n < 1024:
            return f"{n:.0f}{unit}"
        n /= 1024
    return f"{n:.0f}T"


# ----------------- layer 2: schema -----------------

def layer2_schema(src: Path, r: Result) -> dict:
    print("\n=== Layer 2: CSV schema ===")
    dfs: dict[str, pd.DataFrame] = {}
    for name, required in (
        ("pval_edist_full.csv", PVAL_EDIST_REQUIRED_COLS),
        ("targeting_outlier_table.csv", OUTLIER_TABLE_REQUIRED_COLS),
        ("non_targeting_outlier_table.csv", OUTLIER_TABLE_REQUIRED_COLS),
    ):
        p = src / name
        if not p.exists():
            r.report("FAIL", f"{name}: missing (cannot check schema)")
            continue
        try:
            df = pd.read_csv(p, index_col=0)
        except Exception as e:
            r.report("FAIL", f"{name}: failed to parse ({e!r})")
            continue
        dfs[name] = df
        cols = set(df.columns)
        missing = required - cols
        extra = cols - required
        if missing:
            r.report("FAIL", f"{name}: missing required columns {sorted(missing)}")
        else:
            r.report("PASS", f"{name}: all required columns present ({len(df)} rows × {len(df.columns)} cols)")
        if extra:
            r.report("INFO", f"{name}: extra columns {sorted(extra)} (ok)")
    return dfs


# ----------------- layer 3: value-range sanity -----------------

def layer3_values(dfs: dict, r: Result) -> None:
    print("\n=== Layer 3: value-range sanity ===")

    # pval_edist_full
    df = dfs.get("pval_edist_full.csv")
    if df is not None:
        if df.empty:
            r.report("FAIL", "pval_edist_full.csv: empty")
        else:
            r.report("PASS", f"pval_edist_full.csv: {len(df)} target rows")
        # cell_count > 0
        bad = (df["cell_count"] <= 0).sum() if "cell_count" in df.columns else None
        if bad is None:
            r.report("FAIL", "pval_edist_full.csv: no cell_count column")
        elif bad > 0:
            r.report("FAIL", f"pval_edist_full.csv: {bad} rows with cell_count<=0")
        else:
            r.report("PASS", f"pval_edist_full.csv: cell_count>0 for all {len(df)} rows (range {df['cell_count'].min()}-{df['cell_count'].max()})")
        # pval_mean in [0,1]
        if "pval_mean" in df.columns:
            bad = ((df["pval_mean"] < 0) | (df["pval_mean"] > 1)).sum()
            nans = df["pval_mean"].isna().sum()
            if bad > 0:
                r.report("FAIL", f"pval_edist_full.csv: {bad} pval_mean values outside [0,1]")
            else:
                r.report("PASS", f"pval_edist_full.csv: pval_mean in [0,1] (median {df['pval_mean'].median():.4g}, min {df['pval_mean'].min():.4g})")
            if nans > 0:
                r.report("WARN", f"pval_edist_full.csv: {nans} NaN pval_mean rows")
        # distance_mean > 0
        if "distance_mean" in df.columns:
            bad = (df["distance_mean"] <= 0).sum()
            if bad > 0:
                r.report("WARN", f"pval_edist_full.csv: {bad} distance_mean<=0 rows (might be OK for trivial-effect targets)")
            else:
                r.report("PASS", f"pval_edist_full.csv: distance_mean>0 for all rows (median {df['distance_mean'].median():.4g})")
        # type values
        if "type" in df.columns:
            tvals = set(df["type"].dropna().unique())
            unknown = tvals - ALLOWED_TYPES
            if unknown:
                r.report("WARN", f"pval_edist_full.csv: unexpected type values {sorted(unknown)}")
            else:
                r.report("PASS", f"pval_edist_full.csv: all type values known {sorted(tvals)}")
        # all-NaN rows
        nan_rows = df.isna().all(axis=1).sum()
        if nan_rows > 0:
            r.report("FAIL", f"pval_edist_full.csv: {nan_rows} fully-NaN rows")

    # outlier tables
    for name in ("targeting_outlier_table.csv", "non_targeting_outlier_table.csv"):
        df = dfs.get(name)
        if df is None:
            continue
        if df.empty:
            r.report("FAIL", f"{name}: empty")
            continue
        if "pval_outlier" in df.columns:
            bad = ((df["pval_outlier"] < 0) | (df["pval_outlier"] > 1)).sum()
            if bad > 0:
                r.report("FAIL", f"{name}: {bad} pval_outlier values outside [0,1]")
            else:
                r.report("PASS", f"{name}: {len(df)} gRNAs, pval_outlier in [0,1] (median {df['pval_outlier'].median():.4g})")


# ----------------- layer 4: synapse cross-reference -----------------

def layer4_synapse_ref(dfs: dict, r: Result, ref_path: Path) -> None:
    print("\n=== Layer 4: cross-reference vs HTv2 verified run on Synapse ===")
    df = dfs.get("pval_edist_full.csv")
    if df is None:
        r.report("WARN", "no local pval_edist_full to compare")
        return
    if not ref_path.exists():
        token = os.environ.get("SYNAPSE_AUTH_TOKEN")
        if not token:
            r.report("WARN", "SYNAPSE_AUTH_TOKEN not set; skipping schema cross-reference")
            return
        try:
            import synapseclient
        except ImportError:
            r.report("WARN", "synapseclient not installed; skipping schema cross-reference")
            return
        try:
            print(f"  downloading reference {HTV2_SYNAPSE_ID_PVAL} -> {ref_path.parent} ...")
            syn = synapseclient.Synapse(silent=True)
            syn.login(authToken=token)
            e = syn.get(HTV2_SYNAPSE_ID_PVAL, downloadLocation=str(ref_path.parent))
            ref_path.parent.mkdir(parents=True, exist_ok=True)
            os.replace(e.path, ref_path)
        except Exception as e:
            r.report("WARN", f"Synapse download failed: {e!r}")
            return
    try:
        ref = pd.read_csv(ref_path, index_col=0)
    except Exception as e:
        r.report("WARN", f"failed to read reference {ref_path}: {e!r}")
        return
    our = set(df.columns)
    their = set(ref.columns)
    if our == their:
        r.report("PASS", f"pval_edist_full.csv columns identical to HTv2 reference ({len(our)} cols)")
    else:
        only_ours = our - their
        only_theirs = their - our
        if only_theirs:
            r.report("FAIL", f"pval_edist_full.csv missing HTv2 columns: {sorted(only_theirs)}")
        if only_ours:
            r.report("INFO", f"pval_edist_full.csv has extra columns vs HTv2: {sorted(only_ours)} (ok)")


# ----------------- main -----------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source-dir", required=True, type=Path,
                    help="Path to e-distance run output dir (the one with pval_edist_full.csv)")
    ap.add_argument("--no-synapse-ref", action="store_true",
                    help="Skip Layer 4 (cross-reference vs HTv2 verified run on Synapse)")
    ap.add_argument("--ref-cache", default="/tmp/htv2_pval_edist_full.csv", type=Path,
                    help="Where to cache the downloaded HTv2 reference (default: /tmp/htv2_pval_edist_full.csv)")
    args = ap.parse_args()

    src = args.source_dir.resolve()
    if not src.is_dir():
        print(f"ERROR: --source-dir not a directory: {src}", file=sys.stderr)
        return 2

    print(f"Validating: {src}")
    r = Result()
    layer1_presence(src, r)
    dfs = layer2_schema(src, r)
    layer3_values(dfs, r)
    if not args.no_synapse_ref:
        layer4_synapse_ref(dfs, r, args.ref_cache)
    return r.summary()


if __name__ == "__main__":
    sys.exit(main())
