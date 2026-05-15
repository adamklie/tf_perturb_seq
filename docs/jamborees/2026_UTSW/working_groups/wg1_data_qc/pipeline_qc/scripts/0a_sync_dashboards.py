"""Sync each production run's pipeline_dashboard/dashboard.html to the local path
listed in manifests/production_manifest.tsv (dashboard_local column).

Source types:
  - gcs     -> gsutil cp
  - synapse -> synapse get (via synapseclient)

Idempotent: skips files that already exist locally (override with --force).
Currently the manifest points all dashboards at scratch/2026_05_14/upstream_qc/data/
so a one-time bootstrap there is enough; this script lets you re-pull if needed.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MANIFEST = ROOT / "manifests" / "production_manifest.tsv"


def sync_gcs(uri: str, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["gsutil", "cp", uri, str(dst)], check=True)


def sync_synapse(syn_id: str, dst: Path) -> None:
    import synapseclient

    dst.parent.mkdir(parents=True, exist_ok=True)
    syn = synapseclient.Synapse(silent=True)
    syn.login(authToken=os.environ["SYNAPSE_AUTH_TOKEN"])
    ent = syn.get(syn_id, downloadLocation=str(dst.parent), ifcollision="overwrite.local")
    downloaded = Path(ent.path)
    if downloaded != dst:
        shutil.move(str(downloaded), str(dst))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true", help="Re-download even if local file exists")
    ap.add_argument("--only", help="Only sync this short_name (e.g. 'Huangfu DE')")
    args = ap.parse_args()

    df = pd.read_csv(MANIFEST, sep="\t")
    if args.only:
        df = df[df["short_name"] == args.only]
        if df.empty:
            print(f"No row matching short_name={args.only!r}", file=sys.stderr)
            return 1

    failures = []
    for _, row in df.iterrows():
        dst = Path(row["dashboard_local"])
        if dst.exists() and not args.force:
            print(f"[skip] {row['short_name']}: {dst} already exists")
            continue
        print(f"[sync] {row['short_name']}  ({row['dashboard_source_type']})  -> {dst}")
        try:
            if row["dashboard_source_type"] == "gcs":
                sync_gcs(row["dashboard_source_uri"], dst)
            elif row["dashboard_source_type"] == "synapse":
                sync_synapse(row["dashboard_source_uri"], dst)
            else:
                raise ValueError(f"unknown dashboard_source_type: {row['dashboard_source_type']!r}")
        except Exception as e:
            print(f"  FAILED: {e}", file=sys.stderr)
            failures.append(row["short_name"])

    if failures:
        print(f"\n{len(failures)} sync failures: {failures}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
