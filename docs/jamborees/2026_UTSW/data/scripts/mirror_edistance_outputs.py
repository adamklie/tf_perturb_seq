"""Mirror a dataset's energy distance pipeline outputs from the HPC to Synapse.

For a given dataset, this script:
  1. Walks `--source-dir` (an e-distance run output directory on the HPC).
  2. Uploads a whitelisted set of deliverable artifacts to Synapse under
     `2026_UTSW/datasets/<dataset_id>/energy_distance/`.
  3. Records the resulting Synapse folder ID in `synapse_paths.tsv`'s
     `energy_distance` column when the local TSV is reachable.

Designed to be run **on the HPC** where the e-distance outputs live, so the
upload reads from local disk rather than going GCS-roundtrip. Only deliverables
are uploaded — intermediates / source clones / the source MuData are skipped:

  Skipped (not deliverables):
    - inference_mudata.h5mu          (already mirrored under crispr_pipeline/)
    - preprocessed.h5ad              (intermediate; reproducible from MuData)
    - gRNA_dict.pickle               (intermediate)
    - pca_dataframe.pickle           (intermediate)
    - annotation_table.csv           (intermediate; built from MuData)
    - preprocess_mudata_local.py     (auto-generated wrapper script)
    - energy_dist_pipeline/          (cloned source repo)

  Uploaded (deliverables):
    - pval_edist_full.csv
    - targeting_outlier_table.csv
    - non_targeting_outlier_table.csv
    - target_by_target_matrix.csv     (only if step 3 ran)
    - edist_embedding_info.csv        (only if step 3 ran)
    - discordance_gRNA.csv            (if present)
    - config1_2.json
    - config3.json
    - image/                          (whole directory of diagnostic PDFs)
    - logs/                           (slurm .out / .err — useful for debugging)

Authentication:
  - Synapse: SYNAPSE_AUTH_TOKEN env var. Run via `zsh -ic '...'` so the token
    from ~/.zshrc is loaded.

Usage:
    python mirror_edistance_outputs.py \\
        --dataset Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq \\
        --source-dir <your-hpc-scratch>/results/energy_distance/muddy_penguin

    python mirror_edistance_outputs.py --dataset <id> --source-dir <DIR> --dry-run
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd
import synapseclient
from synapseclient import File, Folder

JAMB = Path(__file__).resolve().parent.parent
SYNAPSE_PATHS = JAMB / "synapse_paths.tsv"
SYNAPSE_PARENT = "syn64423137"

DELIVERABLE_FILES = {
    "pval_edist_full.csv",
    "targeting_outlier_table.csv",
    "non_targeting_outlier_table.csv",
    "target_by_target_matrix.csv",
    "edist_embedding_info.csv",
    "discordance_gRNA.csv",
    "config1_2.json",
    "config3.json",
}
DELIVERABLE_DIRS = {"image", "logs"}


def collect_uploads(source_dir: Path) -> tuple[list[Path], list[Path]]:
    """Return (top-level files, top-level dirs) to upload from source_dir."""
    files = sorted(p for p in source_dir.iterdir() if p.is_file() and p.name in DELIVERABLE_FILES)
    dirs = sorted(p for p in source_dir.iterdir() if p.is_dir() and p.name in DELIVERABLE_DIRS)
    return files, dirs


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


def synapse_upload_file(syn: synapseclient.Synapse, path: Path, parent_id: str) -> None:
    if path.stat().st_size == 0:
        # Synapse rejects 0-byte files with `400: File size must be at least one byte`.
        print(f"    skipped (empty) {path.name}")
        return
    try:
        f = syn.store(File(str(path), parent=parent_id))
        print(f"    uploaded {f.id}  {path.name}")
    except Exception as e:
        print(f"    ERROR {path.name}: {e!r}")


def synapse_upload_dir(syn: synapseclient.Synapse, local_dir: Path, parent_id: str) -> None:
    """Upload local_dir's contents into a Synapse folder named local_dir.name under parent_id."""
    target_id = synapse_ensure_path(syn, parent_id, [local_dir.name])
    print(f"  Synapse target folder for {local_dir.name}/ -> {target_id}")
    for path in sorted(local_dir.rglob("*")):
        if path.is_dir():
            continue
        rel = path.relative_to(local_dir)
        sub_parts = list(rel.parent.parts)
        sub_parent = synapse_ensure_path(syn, target_id, sub_parts) if sub_parts else target_id
        synapse_upload_file(syn, path, sub_parent)


def update_synapse_paths(dataset_id: str, syn_folder_id: str) -> None:
    if not SYNAPSE_PATHS.exists():
        print(f"  NOTE: {SYNAPSE_PATHS} not reachable from this host — record {syn_folder_id} manually")
        return
    df = pd.read_csv(SYNAPSE_PATHS, sep="\t", dtype=str, keep_default_na=False)
    if "energy_distance" not in df.columns:
        df["energy_distance"] = ""
    mask = df["dataset_id"] == dataset_id
    if not mask.any():
        sys.exit(f"dataset_id not in {SYNAPSE_PATHS}: {dataset_id}")
    df.loc[mask, "energy_distance"] = syn_folder_id
    df.to_csv(SYNAPSE_PATHS, sep="\t", index=False)
    print(f"  recorded {syn_folder_id} in {SYNAPSE_PATHS}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, help="dataset_id (matches synapse_paths.tsv)")
    ap.add_argument("--source-dir", required=True, type=Path,
                    help="HPC e-distance run output dir (e.g., .../results/energy_distance/muddy_penguin)")
    ap.add_argument("--dry-run", action="store_true", help="List planned uploads, don't transfer")
    args = ap.parse_args()

    src = args.source_dir.resolve()
    if not src.is_dir():
        sys.exit(f"source-dir not a directory: {src}")

    files, dirs = collect_uploads(src)
    print(f"Dataset: {args.dataset}")
    print(f"Source:  {src}")
    print(f"\n--- planned uploads ---")
    for p in files:
        size_mb = p.stat().st_size / 1e6
        print(f"  file {p.name}  ({size_mb:.2f} MB)")
    for d in dirs:
        n = sum(1 for _ in d.rglob("*") if _.is_file())
        print(f"  dir  {d.name}/  ({n} file(s))")

    if not files and not dirs:
        sys.exit("Nothing to upload — no deliverable files/dirs found in source-dir.")

    if args.dry_run:
        print("\n(dry-run; exiting)")
        return

    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        sys.exit("SYNAPSE_AUTH_TOKEN not set in env (run via `zsh -ic`)")
    syn = synapseclient.Synapse(silent=True)
    syn.login(authToken=token)

    edist_parent = synapse_ensure_path(
        syn, SYNAPSE_PARENT, ["2026_UTSW", "datasets", args.dataset, "energy_distance"]
    )
    print(f"\nPer-dataset energy_distance folder on Synapse: {edist_parent}")

    print("\n--- uploading ---")
    for p in files:
        synapse_upload_file(syn, p, edist_parent)
    for d in dirs:
        synapse_upload_dir(syn, d, edist_parent)

    update_synapse_paths(args.dataset, edist_parent)
    print(f"\n[done] energy_distance folder: {edist_parent}")


if __name__ == "__main__":
    main()
