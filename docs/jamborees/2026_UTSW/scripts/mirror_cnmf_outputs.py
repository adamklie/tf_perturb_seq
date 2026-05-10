"""Mirror a dataset's cNMF run outputs from the HPC to Synapse.

Implements the curation rule documented in `schemas/cnmf.json` -> `bundle_inclusion_rule`:

  (1) Selected-k full data for downstream analysis (WG1 + WG2):
        - all loading variants at the selected k
        - cell usages at the selected k
        - integrated MuData at the selected k
        - full Eval/<sel>_<dt>/ TXT bundle
        - selected-k Plot/Program_*/ + Plot/Perturb_gene_*/
        - Annotation/<sel>_<dt>.xlsx
        - Interpretation/Summary_table/<sel>_<dt>/

  (2) Sweep-as-provenance (revisit-the-k-decision):
        - <run>.gene_spectra_score.k_<X>.dt_<dt>.txt for every k
        - <run>.clustering.k_<X>.dt_<Y>.png for every (k, dt)
        - Eval/<k>_<dt>/ for every (k, dt) (full TXT bundles)
        - Plot/k_selection_*/ folder

  (3) Run-level reproducibility:
        - <run>.k_selection.png
        - <run>.k_selection_stats.df.npz
        - <run>.overdispersed_genes.txt
        - config_*.yml
        - logs/
        - README.txt (selected k + rationale)

Excluded (intermediates / duplicates / sweep-only bulk):
        - cnmf_tmp/, Inference/
        - prog_data/, loading/ (duplicate of selected-k loadings)
        - per-k gene_spectra_tpm / spectra.consensus / starcat_spectra / usages.consensus
        - per-k adata/cNMF_<k>_<dt>.h5mu (the dominant ~70 GB sweep bulk)
        - Evaluation/ legacy folder

Designed to be run **on the HPC** where the cNMF run output lives. After running,
rsync the updated `synapse_paths.tsv` back to the laptop.

Authentication:
  - Synapse: SYNAPSE_AUTH_TOKEN env var. Run via `zsh -ic '...'` so the token
    from ~/.zshrc is loaded.

Usage:
    python mirror_cnmf_outputs.py \\
        --dataset Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq \\
        --source-dir /cellar/users/aklie/.../PerturbNMF/Result/050926_HuangfuDE_20iter_5KHVG_torch_halsvar_batch \\
        --selected-k 50

    python mirror_cnmf_outputs.py --dataset <id> --source-dir <DIR> --selected-k <K> --dry-run
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd
import synapseclient
from synapseclient import File, Folder

JAMB = Path(__file__).resolve().parent.parent
SYNAPSE_PATHS = JAMB / "synapse_paths.tsv"
SYNAPSE_PARENT = "syn64423137"

# Selected-k flat files (one per k, but only the chosen k goes here)
SELECTED_K_FLAT_TEMPLATES = [
    "{run}.gene_spectra_tpm.k_{k}.dt_{dt}.txt",
    "{run}.spectra.k_{k}.dt_{dt}.consensus.txt",
    "{run}.starcat_spectra.k_{k}.dt_{dt}.txt",
    "{run}.usages.k_{k}.dt_{dt}.consensus.txt",
    "{run}.clustering.k_{k}.dt_{dt}.png",
]
# All-k sweep flat files (one per k for transparency)
SWEEP_FLAT_PATTERN = re.compile(r"^.*\.gene_spectra_score\.k_\d+\.dt_\d+_\d+\.txt$")
SWEEP_CLUSTERING_PATTERN = re.compile(r"^.*\.clustering\.k_\d+\.dt_\d+_\d+\.png$")
# Run-level
RUN_LEVEL_FILES = [
    "{run}.k_selection.png",
    "{run}.k_selection_stats.df.npz",
    "{run}.overdispersed_genes.txt",
]
RUN_LEVEL_GLOB = ["config_*.yml", "README.txt"]
# Selected-k folders
SELECTED_K_FOLDER_TEMPLATES = [
    "Annotation/{k}_{dt}.xlsx",
    "Interpretation/Summary_table/{k}_{dt}",
    "Plot/Program_{k}_{dt}",
    "Plot/Perturb_gene_{k}_{dt}",
    "adata/cNMF_{k}_{dt}.h5mu",
]
# All-k subdirs of Eval/
EVAL_DIR_NAME = "Eval"
# Top-level folders we always include
ALWAYS_INCLUDE_DIRS = ["logs"]
# Whitelisted under Plot/ — only k_selection_*/ from Plot is generic (rest are selected-k)
PLOT_KSELECTION_PATTERN = re.compile(r"^k_selection.*$")


def fmt_dt(dt: float) -> str:
    return f"{int(dt)}_{int(round((dt - int(dt)) * 10))}"


def collect_uploads(src: Path, run: str, sel_k: int, dt: float) -> list[tuple[Path, str]]:
    """Return list of (local_path, synapse_relative_path) tuples to upload."""
    sel_dt = fmt_dt(dt)
    uploads: list[tuple[Path, str]] = []

    # Selected-k flat files
    for tpl in SELECTED_K_FLAT_TEMPLATES:
        f = src / tpl.format(run=run, k=sel_k, dt=sel_dt)
        if f.is_file():
            uploads.append((f, f.name))

    # All-k sweep flat files (gene_spectra_score for every k)
    for f in sorted(src.iterdir()):
        if f.is_file() and SWEEP_FLAT_PATTERN.match(f.name):
            uploads.append((f, f.name))
        elif f.is_file() and SWEEP_CLUSTERING_PATTERN.match(f.name):
            # All-k clustering pngs (sweep-as-provenance)
            uploads.append((f, f.name))

    # Run-level flat files
    for tpl in RUN_LEVEL_FILES:
        f = src / tpl.format(run=run)
        if f.is_file():
            uploads.append((f, f.name))
    for pattern in RUN_LEVEL_GLOB:
        for f in sorted(src.glob(pattern)):
            if f.is_file():
                uploads.append((f, f.name))

    # Selected-k folders / files
    for tpl in SELECTED_K_FOLDER_TEMPLATES:
        target = src / tpl.format(k=sel_k, dt=sel_dt)
        if target.is_file():
            uploads.append((target, str(target.relative_to(src))))
        elif target.is_dir():
            for sub in target.rglob("*"):
                if sub.is_file():
                    uploads.append((sub, str(sub.relative_to(src))))

    # All-k Eval/<k>_<dt>/ folders
    eval_dir = src / EVAL_DIR_NAME
    if eval_dir.is_dir():
        for k_dt_dir in sorted(eval_dir.iterdir()):
            if not k_dt_dir.is_dir():
                continue
            if not re.match(r"^\d+_\d+_\d+$", k_dt_dir.name) and not re.match(r"^\d+_\d+$", k_dt_dir.name):
                continue
            for sub in sorted(k_dt_dir.iterdir()):
                if sub.is_file():
                    uploads.append((sub, str(sub.relative_to(src))))

    # Plot/k_selection_*/
    plot_dir = src / "Plot"
    if plot_dir.is_dir():
        for sub in sorted(plot_dir.iterdir()):
            if sub.is_dir() and PLOT_KSELECTION_PATTERN.match(sub.name):
                for f in sub.rglob("*"):
                    if f.is_file():
                        uploads.append((f, str(f.relative_to(src))))

    # logs/
    for d in ALWAYS_INCLUDE_DIRS:
        d_path = src / d
        if d_path.is_dir():
            for f in d_path.rglob("*"):
                if f.is_file():
                    uploads.append((f, str(f.relative_to(src))))

    # Dedup (preserve order)
    seen = set()
    deduped = []
    for src_p, rel in uploads:
        if str(src_p) in seen:
            continue
        seen.add(str(src_p))
        deduped.append((src_p, rel))
    return deduped


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


def update_synapse_paths_tsv(dataset: str, syn_id: str) -> None:
    if not SYNAPSE_PATHS.is_file():
        print(f"  warning: {SYNAPSE_PATHS} not found; skipping log")
        return
    df = pd.read_csv(SYNAPSE_PATHS, sep="\t", dtype=str).fillna("")
    if "cnmf" not in df.columns:
        print(f"  warning: 'cnmf' column not in {SYNAPSE_PATHS}; skipping log")
        return
    mask = df["dataset_id"] == dataset
    if not mask.any():
        print(f"  warning: dataset_id '{dataset}' not in {SYNAPSE_PATHS}; skipping log")
        return
    df.loc[mask, "cnmf"] = syn_id
    df.to_csv(SYNAPSE_PATHS, sep="\t", index=False)
    print(f"  logged: {SYNAPSE_PATHS} cnmf <- {syn_id}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", required=True, help="Dataset ID (Synapse subfolder name under 2026_UTSW/datasets/)")
    ap.add_argument("--source-dir", required=True, type=Path, help="HPC path to the cNMF run dir (PerturbNMF/Result/<run_name>/)")
    ap.add_argument("--selected-k", required=True, type=int, help="Group-selected k for downstream analysis")
    ap.add_argument("--density-threshold", type=float, default=2.0, help="Density threshold for selected-k bundle (default: 2.0)")
    ap.add_argument("--dry-run", action="store_true", help="List uploads without transferring")
    args = ap.parse_args()

    if not args.source_dir.is_dir():
        sys.exit(f"--source-dir does not exist or is not a directory: {args.source_dir}")
    run_name = args.source_dir.name

    uploads = collect_uploads(args.source_dir, run_name, args.selected_k, args.density_threshold)
    if not uploads:
        sys.exit(f"No files matched the curation rule under {args.source_dir}")

    total_bytes = sum(p.stat().st_size for p, _ in uploads)
    print(f"run_name: {run_name}")
    print(f"selected_k: {args.selected_k} dt: {args.density_threshold}")
    print(f"files to upload: {len(uploads)} ({total_bytes / 1024 / 1024 / 1024:.2f} GB)")

    if args.dry_run:
        for p, rel in uploads[:30]:
            print(f"  [dry-run] {rel} ({p.stat().st_size:,} B)")
        if len(uploads) > 30:
            print(f"  [dry-run] ... and {len(uploads) - 30} more")
        return 0

    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        sys.exit("SYNAPSE_AUTH_TOKEN env var not set. Run via zsh -ic so ~/.zshrc loads.")
    syn = synapseclient.Synapse()
    syn.login(authToken=token, silent=True)

    target_parts = ["2026_UTSW", "datasets", args.dataset, "cnmf", run_name]
    target_root = ensure_folder_path(syn, SYNAPSE_PARENT, target_parts)
    print(f"target Synapse folder: {target_root}")

    folder_cache: dict[str, str] = {"": target_root}

    def folder_for(rel: str) -> str:
        if rel in folder_cache:
            return folder_cache[rel]
        parent = "/".join(rel.split("/")[:-1])
        parent_id = folder_for(parent)
        new_id = synapse_get_or_create_folder(syn, rel.split("/")[-1], parent_id)
        folder_cache[rel] = new_id
        return new_id

    for src_p, rel in uploads:
        if src_p.stat().st_size == 0:
            print(f"  skipped (empty) {rel}")
            continue
        parent_rel = "/".join(rel.split("/")[:-1])
        parent_id = folder_for(parent_rel)
        syn.store(File(str(src_p), parent=parent_id))
        print(f"  uploaded {rel}")

    update_synapse_paths_tsv(args.dataset, target_root)
    print()
    print(f"DONE. Synapse folder: https://www.synapse.org/Synapse:{target_root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
