"""Snapshot IGVF portal state for the TF Perturb-seq Project (production datasets).

Mirrors the saved query from the user's portal URL:
    https://data.igvf.org/search/?type=<type>&collections=TF+Perturb-seq+Project&collections!=Benchmark

Each run writes:
  - `portal_snapshots/<utc_iso>/<type>.json`            : raw API response (full @graph)
  - `portal_snapshots/<utc_iso>/<type>.tsv`             : flat TSV of the most useful fields
  - `portal_snapshots/<utc_iso>/manifest.tsv`           : query summary (type, count, url)
  - `portal_snapshots/latest -> <utc_iso>`              : symlink to the newest snapshot

Run periodically (cron / manual) to build a history of how the portal state changes.
A separate `diff_snapshots.py` can later compare any two snapshots to highlight new /
removed / changed accessions.

Usage:
    python query_igvf_portal.py
    python query_igvf_portal.py --types MeasurementSet,AnalysisSet
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys
import time
from pathlib import Path
from urllib.parse import urlencode

import pandas as pd
import requests

JAMB = Path(__file__).resolve().parent.parent
SNAPSHOT_DIR = JAMB / "portal_snapshots"

API_BASE = "https://api.data.igvf.org/search/"

# Default query: anything in the TF Perturb-seq Project that's NOT in Benchmark.
DEFAULT_FILTERS: list[tuple[str, str]] = [
    ("collections", "TF Perturb-seq Project"),
    ("collections!", "Benchmark"),
]

# Types to snapshot. Order matters only for output ordering.
DEFAULT_TYPES = [
    "MeasurementSet",
    "AuxiliarySet",
    "AnalysisSet",
    "ConstructLibrarySet",
    "PredictionSet",
]

# Per-type field projections for the TSV summary. Keys = column name, values =
# JSONPath-ish dotted path or callable. Missing fields → empty string.
TSV_FIELDS: dict[str, dict[str, str]] = {
    "MeasurementSet": {
        "accession": "accession",
        "lab": "lab.title",
        "samples_summaries": "samples.0.summary",
        "preferred_assay_titles": "preferred_assay_titles.0",
        "control_type": "control_type",
        "auxiliary_sets": "auxiliary_sets",
        "construct_library_sets": "construct_library_sets",
        "files": "files_count",  # we'll inject this manually if available
        "status": "status",
        "creation_timestamp": "creation_timestamp",
    },
    "AnalysisSet": {
        "accession": "accession",
        "lab": "lab.title",
        "samples_summaries": "samples.0.summary",
        "input_file_sets": "input_file_sets",
        "files_count": "files_count",
        "status": "status",
        "creation_timestamp": "creation_timestamp",
    },
    "AuxiliarySet": {
        "accession": "accession",
        "lab": "lab.title",
        "auxiliary_type": "auxiliary_type",
        "measurement_sets": "measurement_sets",
        "status": "status",
        "creation_timestamp": "creation_timestamp",
    },
    "ConstructLibrarySet": {
        "accession": "accession",
        "lab": "lab.title",
        "scope": "scope",
        "selection_criteria": "selection_criteria",
        "associated_phenotypes": "associated_phenotypes",
        "status": "status",
        "creation_timestamp": "creation_timestamp",
    },
    "PredictionSet": {
        "accession": "accession",
        "lab": "lab.title",
        "input_file_sets": "input_file_sets",
        "status": "status",
        "creation_timestamp": "creation_timestamp",
    },
}


def _dig(obj, path: str):
    """Walk a dotted path through nested dicts/lists; return '' on any miss."""
    if obj is None:
        return ""
    cur = obj
    for part in path.split("."):
        if part.isdigit():
            try:
                cur = cur[int(part)]
            except (IndexError, TypeError):
                return ""
        else:
            if isinstance(cur, dict):
                cur = cur.get(part)
            else:
                return ""
        if cur is None:
            return ""
    if isinstance(cur, list):
        # Render lists of ID-bearing dicts as a `;`-joined accession list when possible.
        out = []
        for v in cur:
            if isinstance(v, dict):
                out.append(v.get("accession") or v.get("@id") or v.get("name") or json.dumps(v))
            else:
                out.append(str(v))
        return ";".join(out)
    if isinstance(cur, dict):
        return cur.get("accession") or cur.get("@id") or json.dumps(cur)
    return cur


def fetch_type(typ: str, filters: list[tuple[str, str]]) -> tuple[str, dict]:
    """Fetch all records of a given type. Uses limit=all to bypass pagination."""
    params = [("type", typ)] + filters + [("limit", "all"), ("format", "json")]
    url = API_BASE + "?" + urlencode(params)
    headers = {"Accept": "application/json"}
    r = requests.get(url, headers=headers, timeout=120)
    r.raise_for_status()
    return url, r.json()


def to_tsv_rows(typ: str, payload: dict) -> pd.DataFrame:
    fields = TSV_FIELDS.get(typ, {"accession": "accession", "status": "status"})
    rows = []
    for entry in payload.get("@graph", []):
        rows.append({col: _dig(entry, path) for col, path in fields.items()})
    return pd.DataFrame(rows, columns=list(fields.keys()))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--types", default=",".join(DEFAULT_TYPES), help="Comma-separated types to query")
    ap.add_argument("--out-root", default=str(SNAPSHOT_DIR), type=Path)
    args = ap.parse_args()

    types = [t.strip() for t in args.types.split(",") if t.strip()]
    out_root = Path(args.out_root)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H-%M-%SZ")
    snap_dir = out_root / stamp
    snap_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    for typ in types:
        try:
            url, payload = fetch_type(typ, DEFAULT_FILTERS)
        except requests.HTTPError as e:
            print(f"[{typ}] HTTP error: {e}; skipping")
            continue
        graph = payload.get("@graph", [])
        total = payload.get("total", len(graph))
        (snap_dir / f"{typ}.json").write_text(json.dumps(payload, indent=2))
        df = to_tsv_rows(typ, payload)
        df.to_csv(snap_dir / f"{typ}.tsv", sep="\t", index=False)
        print(f"[{typ}] {len(graph)} of {total} records → {snap_dir / f'{typ}.tsv'}")
        manifest_rows.append({"type": typ, "n_records": len(graph), "total": total, "url": url})
        time.sleep(0.5)  # polite

    pd.DataFrame(manifest_rows).to_csv(snap_dir / "manifest.tsv", sep="\t", index=False)

    # Update `latest` symlink
    latest = out_root / "latest"
    if latest.exists() or latest.is_symlink():
        latest.unlink()
    latest.symlink_to(stamp)
    print(f"\nSnapshot dir: {snap_dir}")
    print(f"Latest symlink → {latest}")


if __name__ == "__main__":
    main()
