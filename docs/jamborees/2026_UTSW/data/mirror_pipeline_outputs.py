"""Mirror a dataset's CRISPR pipeline output bundle from GCS to Synapse.

For a given dataset, this script:
  1. Reads `gcs_output_path` from `reference/experimental_metadata.tsv`.
  2. Downloads the three terminal pipeline directories (`pipeline_dashboard/`,
     `pipeline_info/`, `pipeline_outputs/`) from GCS to a temp dir.
  3. Uploads each to Synapse under `2026_UTSW/datasets/<dataset_id>/crispr_pipeline/`.
  4. Records resulting Synapse folder IDs in `synapse_paths.tsv` (column
     `crispr_pipeline`).
  5. Cleans up the temp dir.

Authentication:
  - GCS: `gcloud auth` set to the IGVF service account
    (`adamklie@igvf-pertub-seq-pipeline.iam.gserviceaccount.com`).
  - Synapse: `SYNAPSE_AUTH_TOKEN` env var. Run via `zsh -ic '...'` so the token from
    ~/.zshrc is loaded.

Usage:
    python mirror_pipeline_outputs.py --dataset <id>                  # default workdir: /tmp/jamboree_pipeline_<id>
    python mirror_pipeline_outputs.py --dataset <id> --workdir <DIR>  # use this if /tmp is small
    python mirror_pipeline_outputs.py --dataset <id> --dry-run        # list only, don't transfer
    python mirror_pipeline_outputs.py --dataset <id> --skip-download  # use existing temp dir

DISK SPACE: bundles are ~60 GB each. The default `/tmp/jamboree_pipeline_<id>` works on
the laptop only when the bundle fits. **On the HPC, `/tmp` is just 20 GB — pass
`--workdir /cellar/users/aklie/scratch/jamboree_pipeline_staging/<id>` to use the
roomy `/cellar` filesystem.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import synapseclient
from synapseclient import File, Folder

JAMB = Path(__file__).resolve().parent.parent
EXP_METADATA = JAMB / "reference/experimental_metadata.tsv"
SYNAPSE_PATHS = JAMB / "synapse_paths.tsv"
SYNAPSE_PARENT = "syn64423137"

BUNDLE_DIRS = ["pipeline_dashboard", "pipeline_info", "pipeline_outputs"]


def lookup_gcs_path(dataset_id: str) -> str:
    df = pd.read_csv(EXP_METADATA, sep="\t")
    row = df[df["dataset_id"] == dataset_id]
    if row.empty:
        sys.exit(f"dataset_id not found in {EXP_METADATA}: {dataset_id}")
    gcs = row["gcs_output_path"].iloc[0]
    if not isinstance(gcs, str) or not gcs.startswith("gs://"):
        sys.exit(f"gcs_output_path not set for {dataset_id}: {gcs!r}")
    return gcs.rstrip("/")


def gcs_du(path: str) -> int:
    out = subprocess.run(
        ["gcloud", "storage", "du", "-s", path],
        capture_output=True, text=True, check=False,
    ).stdout
    for line in out.splitlines():
        parts = line.strip().split()
        if parts and parts[0].isdigit():
            return int(parts[0])
    return -1


def gcs_download(src: str, dst: Path) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    print(f"  downloading {src} -> {dst}")
    subprocess.check_call(["gcloud", "storage", "rsync", "-r", src, str(dst)])


def synapse_get_or_create_folder(syn: synapseclient.Synapse, name: str, parent_id: str) -> str:
    for child in syn.getChildren(parent_id):
        if child["type"].endswith(".Folder") and child["name"] == name:
            return child["id"]
    folder = syn.store(Folder(name, parent=parent_id))
    return folder.id


def synapse_ensure_path(syn: synapseclient.Synapse, parent_id: str, parts: list[str]) -> str:
    cur = parent_id
    for p in parts:
        cur = synapse_get_or_create_folder(syn, p, cur)
    return cur


def synapse_upload_dir(
    syn: synapseclient.Synapse, local_dir: Path, parent_id: str, dry_run: bool = False
) -> str:
    """Upload local_dir's contents into a Synapse folder named local_dir.name under parent_id."""
    target_id = synapse_ensure_path(syn, parent_id, [local_dir.name])
    print(f"  Synapse target folder for {local_dir.name}/ -> {target_id}")
    if dry_run:
        print("  (dry-run; skipping uploads)")
        return target_id
    for path in sorted(local_dir.rglob("*")):
        if path.is_dir():
            continue
        rel = path.relative_to(local_dir)
        if path.stat().st_size == 0:
            # Synapse rejects 0-byte files with `400: File size must be at least one byte`.
            print(f"    skipped (empty) {rel}")
            continue
        sub_parts = list(rel.parent.parts)  # subdirectory parts
        sub_parent = synapse_ensure_path(syn, target_id, sub_parts) if sub_parts else target_id
        try:
            f = syn.store(File(str(path), parent=sub_parent))
            print(f"    uploaded {f.id}  {rel}")
        except Exception as e:
            print(f"    ERROR {rel}: {e!r}")
    return target_id


def update_synapse_paths(dataset_id: str, syn_folder_id: str) -> None:
    df = pd.read_csv(SYNAPSE_PATHS, sep="\t", dtype=str, keep_default_na=False)
    if "crispr_pipeline" not in df.columns:
        df["crispr_pipeline"] = ""
    mask = df["dataset_id"] == dataset_id
    if not mask.any():
        sys.exit(f"dataset_id not in {SYNAPSE_PATHS}: {dataset_id}")
    df.loc[mask, "crispr_pipeline"] = syn_folder_id
    df.to_csv(SYNAPSE_PATHS, sep="\t", index=False)
    print(f"  recorded {syn_folder_id} in {SYNAPSE_PATHS}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, help="dataset_id from experimental_metadata.tsv")
    ap.add_argument("--workdir", default=None, help="Local staging dir (default: /tmp/jamboree_pipeline_<dataset>)")
    ap.add_argument("--keep-workdir", action="store_true", help="Don't delete the temp dir after upload")
    ap.add_argument("--skip-download", action="store_true", help="Assume workdir is already populated")
    ap.add_argument("--dry-run", action="store_true", help="List sources + sizes, do not download or upload")
    args = ap.parse_args()

    gcs_root = lookup_gcs_path(args.dataset)
    print(f"Dataset: {args.dataset}")
    print(f"GCS root: {gcs_root}")

    print("\n--- size check ---")
    for d in BUNDLE_DIRS:
        size = gcs_du(f"{gcs_root}/{d}/")
        print(f"  {d}/: {size / 1e9:.2f} GB" if size > 0 else f"  {d}/: (empty / not found)")

    if args.dry_run:
        print("\n(dry-run; exiting)")
        return

    workdir = Path(args.workdir or f"/tmp/jamboree_pipeline_{args.dataset}")
    workdir.mkdir(parents=True, exist_ok=True)
    print(f"\nStaging dir: {workdir}")

    if not args.skip_download:
        print("\n--- downloading from GCS ---")
        for d in BUNDLE_DIRS:
            gcs_download(f"{gcs_root}/{d}/", workdir / d)
    else:
        print("--skip-download set; using existing files in workdir")

    print("\n--- uploading to Synapse ---")
    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        sys.exit("SYNAPSE_AUTH_TOKEN not set in env (run via `zsh -ic`)")
    syn = synapseclient.Synapse(silent=True)
    syn.login(authToken=token)

    crispr_parent = synapse_ensure_path(
        syn, SYNAPSE_PARENT, ["2026_UTSW", "datasets", args.dataset, "crispr_pipeline"]
    )
    print(f"Per-dataset crispr_pipeline folder on Synapse: {crispr_parent}")

    for d in BUNDLE_DIRS:
        local = workdir / d
        if not local.exists():
            print(f"  WARN: {d}/ missing in workdir, skipping")
            continue
        synapse_upload_dir(syn, local, crispr_parent)

    update_synapse_paths(args.dataset, crispr_parent)

    if not args.keep_workdir:
        print(f"\nCleaning up {workdir}")
        shutil.rmtree(workdir)
    else:
        print(f"\nKeeping {workdir} (--keep-workdir set)")


if __name__ == "__main__":
    main()
