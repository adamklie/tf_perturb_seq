"""Mirror a CRISPR pipeline output bundle from a local HPC directory to Synapse.

Companion to `mirror_pipeline_outputs.py` (GCS source). Use this variant when the
canonical run lives on the HPC rather than GCS — currently the case for the
Gersbach HTv2 benchmark reference run.

Uploads the contents of `--source-dir` directly under
`2026_UTSW/datasets/<dataset_id>/crispr_pipeline/` on Synapse (no run-label
subfolder — describe the run in a README inside the source dir or alongside).

What gets uploaded:
  Canonical layout (whichever subset is present):
    - pipeline_dashboard/   (~1 GB; dashboard.html, inference_mudata.h5mu, additional_qc/, evaluation_output/, figures/, ...)
    - pipeline_info/        (params JSON + software-versions YAML)
    - pipeline_outputs/     (~150 MB; perturbo cis/trans per-element/per-guide TSVs)
  Optional auxiliary dirs (included with --include-aux if present):
    - calibration/          (DEG-calibrated TSVs from CRISPR pipeline post-processing)
    - tf/                   (TF benchmark outputs)
  Optional README:
    - --readme-path <FILE>  (uploaded as README.md alongside the dirs)

Skipped:
  - anndata/                (preprocessing intermediates, ~4 GB on HTv2)

Authentication:
  - Synapse: `SYNAPSE_AUTH_TOKEN` env var. Run via `zsh -ic '...'` so the token
    from `~/.zshrc` is loaded.

Usage:
    python mirror_pipeline_outputs_hpc.py \\
        --dataset Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2 \\
        --source-dir /cellar/.../HTv2/runs/cleanser_800_mito_15pc \\
        --readme-path docs/jamborees/2026_UTSW/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/crispr_pipeline/README.md \\
        --include-aux

    python mirror_pipeline_outputs_hpc.py --dataset <id> --source-dir <DIR> --dry-run
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import synapseclient
from synapseclient import File, Folder

SYNAPSE_PARENT = "syn64423137"

CANONICAL_DIRS = ["pipeline_dashboard", "pipeline_info", "pipeline_outputs"]
OPTIONAL_AUX_DIRS = ["calibration", "tf"]


def synapse_get_or_create_folder(syn, name, parent_id):
    for child in syn.getChildren(parent_id):
        if child["type"].endswith(".Folder") and child["name"] == name:
            return child["id"]
    return syn.store(Folder(name, parent=parent_id)).id


def ensure_folder_path(syn, parent_id, parts):
    cur = parent_id
    for p in parts:
        cur = synapse_get_or_create_folder(syn, p, cur)
    return cur


def upload_file(syn, path, parent_id, dry_run):
    if path.stat().st_size == 0:
        print(f"    skipped (empty) {path.name}")
        return
    if dry_run:
        print(f"    [dry-run] would upload {path.name} ({path.stat().st_size:,} B)")
        return
    syn.store(File(str(path), parent=parent_id))
    print(f"    uploaded {path.name}")


def upload_dir_recursive(syn, source_dir, parent_id, dry_run):
    """Upload every file under source_dir, mirroring the directory tree in Synapse."""
    for entry in sorted(source_dir.iterdir()):
        if entry.is_file():
            upload_file(syn, entry, parent_id, dry_run)
        elif entry.is_dir():
            sub_parent = (
                ensure_folder_path(syn, parent_id, [entry.name])
                if not dry_run
                else parent_id
            )
            if dry_run:
                print(f"  [dry-run] would create folder {entry.name}/ under {parent_id}")
            else:
                print(f"  folder: {entry.name}/")
            upload_dir_recursive(syn, entry, sub_parent, dry_run)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", required=True, help="Dataset ID (used as the Synapse subfolder name under 2026_UTSW/datasets/)")
    ap.add_argument("--source-dir", required=True, type=Path, help="HPC path to the run directory containing pipeline_dashboard/ etc.")
    ap.add_argument("--readme-path", type=Path, help="Optional path to a README file to upload as README.md alongside the dirs.")
    ap.add_argument("--include-aux", action="store_true", help="Also upload calibration/ and tf/ subdirs if present.")
    ap.add_argument("--dry-run", action="store_true", help="List what would be uploaded, don't transfer.")
    args = ap.parse_args()

    if not args.source_dir.is_dir():
        sys.exit(f"--source-dir does not exist or is not a directory: {args.source_dir}")

    dirs_to_upload = []
    for name in CANONICAL_DIRS:
        sub = args.source_dir / name
        if sub.is_dir():
            dirs_to_upload.append(sub)
        else:
            print(f"  note: {name}/ NOT PRESENT in source — skipping")

    if args.include_aux:
        for name in OPTIONAL_AUX_DIRS:
            sub = args.source_dir / name
            if sub.is_dir():
                dirs_to_upload.append(sub)

    if not dirs_to_upload:
        sys.exit("No canonical pipeline subdirs found under --source-dir. Aborting.")

    syn = synapseclient.Synapse()
    if not args.dry_run:
        token = os.environ.get("SYNAPSE_AUTH_TOKEN")
        if not token:
            sys.exit("SYNAPSE_AUTH_TOKEN env var not set. Run via zsh -ic so ~/.zshrc loads.")
        syn.login(authToken=token, silent=True)

    target_parts = ["2026_UTSW", "datasets", args.dataset, "crispr_pipeline"]
    if args.dry_run:
        print(f"[dry-run] target: syn:{SYNAPSE_PARENT}/{'/'.join(target_parts)}")
        target_parent = SYNAPSE_PARENT
    else:
        target_parent = ensure_folder_path(syn, SYNAPSE_PARENT, target_parts)
        print(f"target Synapse folder: {target_parent}")

    if args.readme_path:
        if not args.readme_path.is_file():
            sys.exit(f"--readme-path does not exist: {args.readme_path}")
        print("=== README.md ===")
        upload_file(syn, args.readme_path, target_parent, args.dry_run)

    for sub in dirs_to_upload:
        print(f"=== {sub.name}/ ===")
        if args.dry_run:
            sub_parent = target_parent
        else:
            sub_parent = ensure_folder_path(syn, target_parent, [sub.name])
        upload_dir_recursive(syn, sub, sub_parent, args.dry_run)

    print()
    print(f"DONE. Synapse folder for run: {target_parent}")
    print(f"  https://www.synapse.org/Synapse:{target_parent}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
