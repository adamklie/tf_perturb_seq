#!/usr/bin/env python3
"""Build tiny per-dataset seqspec-validation samplesheets + prestage assets.

DRAFT — for review. Designed to run on the NRNB login node before submitting the
SLURM array (it needs gsutil + portal access to prestage onlist/guide_design files,
and it bakes target-absolute paths into the emitted CSVs).

For each validatable dataset (one with seqspec YAMLs under setup/seqspec/), this:
  1. picks a source samplesheet,
  2. selects ONE representative complete lane (a measurement_sets value carrying both
     a scRNA and a gRNA row; falls back to first scRNA + first gRNA row),
  3. prestages the gRNA row's barcode_onlist + guide_design to local files
     (workaround for the validator's csv.gz/tsv.gz guide-design 404 bug),
  4. emits a 2-row CSV (scRNA + gRNA) with seqspec rewritten to the dataset's
     absolute seqspec path and onlist/guide_design rewritten to the prestaged paths,
  5. appends a row to a manifest TSV consumed by submit_seqspec_validation.sh.

Validation against fastqs / running the validator is NOT done here.

Usage (on NRNB):
    python build_validation_samplesheets.py \
        --repo-root /carter/users/aklie/projects/tf_perturb_seq \
        --target nrnb \
        --manifest /carter/users/aklie/projects/tf_perturb_seq/scripts/seqspec_validation/runs/2026_06_10/manifest.tsv \
        [--datasets Hon_WTC11-benchmark_TF-Perturb-seq ...] \
        [--keypair /path/igvf_key.json]
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import urllib.request
from base64 import b64encode
from json import loads
from pathlib import Path

PORTAL = "https://api.data.igvf.org"

# Canonical seqspec column = 10-column samplesheet schema used across all datasets.
COLUMNS = [
    "R1_path", "R2_path", "file_modality", "measurement_sets", "sequencing_run",
    "lane", "seqspec", "barcode_onlist", "guide_design", "barcode_hashtag_map",
]

MODALITY_CANON = {
    "scrna sequencing": "scRNA", "scrna": "scRNA",
    "grna sequencing": "gRNA", "grna": "gRNA",
    "cell hashing barcode sequencing": "hash", "hash": "hash",
}

# Per-target absolute repo root used when baking the seqspec / asset paths into the CSV.
TARGET_ROOTS = {
    "nrnb": "/carter/users/aklie/projects/tf_perturb_seq",
    "local": "/Users/adamklie/Desktop/tfp3/tf_perturb_seq",
}

SEQSPEC_BY_MODALITY = {"scRNA": "rna_seqspec.yml", "gRNA": "guide_seqspec.yml"}


def canon_modality(raw: str) -> str:
    return MODALITY_CANON.get(str(raw or "").strip().lower(), str(raw or "").strip())


def is_igvf_accession(v: str) -> bool:
    t = str(v or "").strip().upper()
    return t.startswith("IGVFFI") and t[6:].isalnum()


def is_gcs(v: str) -> bool:
    return str(v or "").strip().lower().startswith("gs://")


def auth_header(keypair: str | None) -> str | None:
    key = secret = None
    if keypair:
        kp = loads(Path(keypair).read_text())
        key, secret = kp["key"], kp["secret"]
    else:
        key, secret = os.getenv("IGVF_API_KEY"), os.getenv("IGVF_SECRET_KEY")
    if not (key and secret):
        return None
    return "Basic " + b64encode(f"{key}:{secret}".encode()).decode()


def portal_get(path: str, auth: str | None) -> dict:
    req = urllib.request.Request(f"{PORTAL}{path}")
    if auth:
        req.add_header("Authorization", auth)
    with urllib.request.urlopen(req, timeout=60) as r:  # noqa: S310 (trusted portal)
        return loads(r.read().decode())


def prestage_asset(value: str, dest_dir: Path, auth: str | None, target_root: str) -> str:
    """Resolve an onlist/guide_design cell to a *target-absolute* local path, downloading
    it if needed. Returns the absolute path string to bake into the CSV."""
    value = str(value or "").strip()
    if not value:
        return ""
    dest_dir.mkdir(parents=True, exist_ok=True)

    if os.path.isabs(value) and os.path.exists(value):
        return value  # already a local asset

    if is_igvf_accession(value):
        obj = portal_get(f"/tabular-files/{value}/@@object?format=json", auth)
        href = obj.get("href")  # e.g. /tabular-files/IGVFFI..../@@download/IGVFFI....csv.gz
        if not href:
            raise RuntimeError(f"No href for {value}; file_format={obj.get('file_format')}")
        ext = "".join(Path(href).suffixes)  # preserves real .csv.gz vs .tsv.gz
        dest = dest_dir / f"{value}{ext}"
        if not dest.exists():
            req = urllib.request.Request(f"{PORTAL}{href}")
            if auth:
                req.add_header("Authorization", auth)
            with urllib.request.urlopen(req, timeout=300) as r, open(dest, "wb") as fh:  # noqa: S310
                fh.write(r.read())
        return _retarget(dest, dest_dir, target_root)

    if is_gcs(value):
        name = value.rstrip("/").split("/")[-1]
        dest = dest_dir / name
        if not dest.exists():
            subprocess.run(["gsutil", "cp", value, str(dest)], check=True)
        return _retarget(dest, dest_dir, target_root)

    # Unknown (relative path?) — leave as-is and warn.
    print(f"  WARN: cannot prestage asset {value!r}; leaving as-is", file=sys.stderr)
    return value


def _retarget(dest: Path, assets_dir: Path, target_root: str) -> str:
    """Express a freshly-downloaded asset path under the target's absolute repo root.

    dest is written under <real-repo>/datasets/.../assets/<name>; we only know the path
    relative to assets_dir's repo root by string. Here dest is already absolute on the
    machine running this script; when --target matches that machine, return it directly."""
    return str(dest.resolve())


def select_lane(rows: list[dict]) -> tuple[dict, dict, str, str]:
    """Return (scrna_row, grna_row, measurement_set, note)."""
    for r in rows:
        r["file_modality"] = canon_modality(r.get("file_modality", ""))

    by_ms: dict[str, set[str]] = {}
    for r in rows:
        by_ms.setdefault(r.get("measurement_sets", ""), set()).add(r["file_modality"])

    complete = sorted(ms for ms, mods in by_ms.items()
                      if {"scRNA", "gRNA"} <= mods and ms)
    if complete:
        ms = complete[0]
        sc = next(r for r in rows if r["file_modality"] == "scRNA" and r["measurement_sets"] == ms)
        gr = next(r for r in rows if r["file_modality"] == "gRNA" and r["measurement_sets"] == ms)
        return sc, gr, ms, "paired_by_measurement_sets"

    # fallback: first scRNA + first gRNA, possibly different samples
    sc = next((r for r in rows if r["file_modality"] == "scRNA"), None)
    gr = next((r for r in rows if r["file_modality"] == "gRNA"), None)
    if not (sc and gr):
        raise RuntimeError("dataset lacks both scRNA and gRNA rows")
    return sc, gr, sc.get("measurement_sets", ""), "FALLBACK_unpaired_first_rows"


def pick_source_samplesheet(ss_dir: Path, override: str | None) -> Path:
    if override:
        p = ss_dir / override if not os.path.isabs(override) else Path(override)
        if not p.exists():
            raise FileNotFoundError(p)
        return p
    cands = sorted(ss_dir.glob("sample_metadata*.csv"))
    local = [c for c in cands if "_localseqspec" in c.name]
    if local:
        return sorted(local)[-1]
    plain = [c for c in cands if "_patched" not in c.name and c.name != "sample_metadata.csv"]
    if plain:
        return sorted(plain)[-1]
    base = ss_dir / "sample_metadata.csv"
    if base.exists():
        return base
    raise FileNotFoundError(f"no sample_metadata*.csv in {ss_dir}")


def build_for_dataset(ds: str, repo_root: Path, target_root: str, auth: str | None,
                      source_override: str | None) -> dict | None:
    ds_dir = repo_root / "datasets" / ds
    seqspec_dir = ds_dir / "setup" / "seqspec"
    if not (seqspec_dir / "rna_seqspec.yml").exists() or not (seqspec_dir / "guide_seqspec.yml").exists():
        print(f"SKIP {ds}: missing rna/guide seqspec under {seqspec_dir}", file=sys.stderr)
        return None

    src = pick_source_samplesheet(ds_dir / "setup" / "samplesheets", source_override)
    with open(src, newline="") as fh:
        rows = list(csv.DictReader(fh))

    sc, gr, ms, note = select_lane(rows)
    print(f"{ds}: source={src.name} ms={ms} ({note})")

    out_dir = ds_dir / "setup" / "seqspec_validation"
    assets_dir = out_dir / "assets"

    # absolute seqspec paths expressed for the target machine
    rel_seqspec = lambda name: f"{target_root}/datasets/{ds}/setup/seqspec/{name}"

    onlist = prestage_asset(gr.get("barcode_onlist", ""), assets_dir, auth, target_root)
    guide = prestage_asset(gr.get("guide_design", ""), assets_dir, auth, target_root)

    def make_row(src_row: dict, modality: str) -> dict:
        return {
            "R1_path": src_row.get("R1_path", ""),
            "R2_path": src_row.get("R2_path", ""),
            "file_modality": modality,
            "measurement_sets": src_row.get("measurement_sets", ""),
            "sequencing_run": src_row.get("sequencing_run", ""),
            "lane": src_row.get("lane", ""),
            "seqspec": rel_seqspec(SEQSPEC_BY_MODALITY[modality]),
            "barcode_onlist": onlist,
            "guide_design": guide if modality == "gRNA" else "",
            "barcode_hashtag_map": "",
        }

    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"{ds}__validation_lane.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerow(make_row(sc, "scRNA"))
        w.writerow(make_row(gr, "gRNA"))

    csv_target = f"{target_root}/datasets/{ds}/setup/seqspec_validation/{csv_path.name}"
    analysis_root = f"{target_root}/datasets/{ds}/setup/seqspec_validation/results"
    downloads_dir = f"{target_root}/datasets/{ds}/setup/seqspec_validation/downloads"
    return {
        "dataset": ds,
        "validation_csv": csv_target,
        "analysis_root": analysis_root,
        "downloads_dir": downloads_dir,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root", required=True, type=Path)
    ap.add_argument("--target", choices=sorted(TARGET_ROOTS), default="nrnb")
    ap.add_argument("--manifest", required=True, type=Path,
                    help="Output manifest TSV consumed by submit_seqspec_validation.sh")
    ap.add_argument("--datasets", nargs="*", help="Subset; default = all validatable datasets")
    ap.add_argument("--source-samplesheet", help="Override source sheet name (applies to all)")
    ap.add_argument("--keypair", help="IGVF keypair JSON {key,secret}")
    args = ap.parse_args()

    target_root = TARGET_ROOTS[args.target]
    auth = auth_header(args.keypair)
    if auth is None:
        print("WARN: no IGVF credentials; accession prestaging will fail", file=sys.stderr)

    ds_root = args.repo_root / "datasets"
    if args.datasets:
        datasets = args.datasets
    else:
        datasets = sorted(
            d.name for d in ds_root.iterdir()
            if (d / "setup" / "seqspec" / "rna_seqspec.yml").exists()
            and (d / "setup" / "seqspec" / "guide_seqspec.yml").exists()
        )

    manifest_rows = []
    for ds in datasets:
        try:
            row = build_for_dataset(ds, args.repo_root, target_root, auth, args.source_samplesheet)
        except Exception as exc:  # noqa: BLE001 — draft; surface and continue
            print(f"ERROR {ds}: {exc}", file=sys.stderr)
            continue
        if row:
            manifest_rows.append(row)

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    with open(args.manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["dataset", "validation_csv", "analysis_root", "downloads_dir"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(manifest_rows)

    print(f"\nWrote manifest with {len(manifest_rows)} dataset(s): {args.manifest}")
    print(f"Submit with: sbatch --array=1-{len(manifest_rows)} "
          f"{target_root}/scripts/seqspec_validation/submit_seqspec_validation.sh "
          f"{args.manifest} {args.repo_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
