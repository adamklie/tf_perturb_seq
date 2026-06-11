#!/usr/bin/env python3
"""Collect crispr_validator analysis_summary.json files into one summary TSV.

DRAFT — for review. Walks every
  datasets/<ds>/setup/seqspec_validation/results/**/analysis_summary.json
(or a provided manifest's analysis_root dirs) and flattens each summary's
`comparison_rows` into a long TSV plus a per-(dataset,modality) pass/fail rollup.

Usage:
    python collect_results.py --repo-root /carter/users/aklie/projects/tf_perturb_seq \
        --output scripts/seqspec_validation/runs/2026_06_10/summary.tsv
    # or scope to a manifest:
    python collect_results.py --manifest .../manifest.tsv --output summary.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

# worst-first ordering, mirrors FLAG_ORDER in seqspec_parser.py
FLAG_ORDER = {
    "perfect_match": 0,
    "close_enough": 1,
    "very_distant": 2,
    "missing_prediction": 3,
    "missing_seqspec": 4,
}
PASS_FLAGS = {"perfect_match", "close_enough"}

LONG_COLS = [
    "dataset", "group_label", "modality", "region", "flag",
    "max_distance_bp", "note",
]
ROLLUP_COLS = ["dataset", "modality", "n_regions", "worst_flag", "pass"]


def dataset_from_path(p: Path) -> str:
    """…/datasets/<ds>/setup/seqspec_validation/results/<group>/analysis_summary.json"""
    parts = p.parts
    if "datasets" in parts:
        return parts[parts.index("datasets") + 1]
    return "<unknown>"


def find_summaries(repo_root: Path | None, manifest: Path | None) -> list[Path]:
    found: list[Path] = []
    if manifest:
        with open(manifest, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                root = Path(row["analysis_root"])
                found.extend(sorted(root.rglob("analysis_summary.json")))
    elif repo_root:
        glob = "datasets/*/setup/seqspec_validation/results/**/analysis_summary.json"
        found.extend(sorted(repo_root.glob(glob)))
    return found


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root", type=Path)
    ap.add_argument("--manifest", type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    if not (args.repo_root or args.manifest):
        ap.error("provide --repo-root or --manifest")

    summaries = find_summaries(args.repo_root, args.manifest)
    if not summaries:
        print("No analysis_summary.json found.", file=sys.stderr)
        return 1

    long_rows: list[dict] = []
    rollup: dict[tuple[str, str], dict] = {}

    for sp in summaries:
        ds = dataset_from_path(sp)
        try:
            data = json.loads(sp.read_text())
        except Exception as exc:  # noqa: BLE001
            print(f"WARN: cannot parse {sp}: {exc}", file=sys.stderr)
            continue
        group_label = data.get("group_label", sp.parent.name)
        for cr in data.get("comparison_rows", []):
            modality = cr.get("modality", "")
            flag = cr.get("flag", "")
            long_rows.append({
                "dataset": ds,
                "group_label": group_label,
                "modality": modality,
                "region": cr.get("region", ""),
                "flag": flag,
                "max_distance_bp": cr.get("max_distance_bp", ""),
                "note": cr.get("note", ""),
            })
            key = (ds, modality)
            agg = rollup.setdefault(key, {"n": 0, "worst": "perfect_match"})
            agg["n"] += 1
            if FLAG_ORDER.get(flag, 99) > FLAG_ORDER.get(agg["worst"], -1):
                agg["worst"] = flag

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=LONG_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(long_rows)

    rollup_path = args.output.with_name(args.output.stem + "_rollup.tsv")
    with open(rollup_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ROLLUP_COLS, delimiter="\t")
        w.writeheader()
        for (ds, modality), agg in sorted(rollup.items()):
            w.writerow({
                "dataset": ds,
                "modality": modality,
                "n_regions": agg["n"],
                "worst_flag": agg["worst"],
                "pass": agg["worst"] in PASS_FLAGS,
            })

    n_pass = sum(1 for a in rollup.values() if a["worst"] in PASS_FLAGS)
    print(f"Parsed {len(summaries)} summary file(s) -> {len(long_rows)} comparison rows.")
    print(f"Per-(dataset,modality) units: {len(rollup)}  pass: {n_pass}  review: {len(rollup) - n_pass}")
    print(f"  long TSV:   {args.output}")
    print(f"  rollup TSV: {rollup_path}")
    print("NOTE: a gRNA 'missing_seqspec' may be the region-type labeling artifact "
          "(region_type: cdna mislabeled rna), not a real mismatch — inspect before failing.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
