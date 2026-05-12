"""Per-dataset audit: local filesystem vs GitHub-tracked vs Synapse mirror.

Usage:
    uv run python scripts/audit_dataset.py <dataset_id> [--synapse]

Without --synapse: skips the Synapse query (fast, offline).
With --synapse:    walks the dataset's Synapse mirror via SYNAPSE_AUTH_TOKEN.

Emits a markdown report on stdout that distinguishes:
  - scripts / configs / docs  (KEEP on GitHub)
  - pipeline results          (REMOVE from GitHub; keep local + Synapse)
  - samplesheets              (REMOVE from GitHub; canonical lives in GCS pipeline_config)
  - bulk binaries (.h5ad/.h5mu, .pickle, .tsv.gz, large .csv)
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SYNAPSE_PATHS_TSV = REPO / "docs" / "jamborees" / "2026_UTSW" / "synapse_paths.tsv"

SCRIPT_EXTS = {"sh", "py", "ipynb", "R", "config", "yaml", "yml", "toml", "json"}
DOC_EXTS = {"md", "txt"}


def classify(path: str, size: int) -> str:
    """Bucket a tracked file path."""
    fname = os.path.basename(path)
    ext = fname.rsplit(".", 1)[-1].lower() if "." in fname else ""

    # Definite results
    if "/calibration/" in path and ext in {"tsv", "csv", "npy", "pickle", "pkl", "pdf", "png"}:
        return "result-calibration"
    if "/pipeline_qc/" in path:
        return "result-pipeline_qc"
    if "/energy_distance/" in path and ext in {"tsv", "csv", "pickle", "pkl", "pdf", "png"}:
        return "result-energy_distance"
    if "/energy_distance/logs/" in path or "/energy_distance/image/" in path:
        return "result-energy_distance"
    if "/cnmf/" in path and ("/Result/" in path or "/Data/" in path):
        return "result-cnmf"
    if "/PerturbNMF/" in path and ("/Result/" in path or "/Data/" in path):
        return "result-cnmf"
    if fname.startswith("perturbo_") and "_per_" in fname and (ext.startswith("tsv") or ext == "gz"):
        return "result-perturbo"
    if ext in {"err", "out"} and "/logs/" in path:
        return "log"
    if ext in {"err", "out"}:
        return "log"
    if ext in {"h5ad", "h5mu", "pickle", "pkl"}:
        return "binary"
    if ext == "gz" and (".tsv" in fname or ".csv" in fname):
        return "binary"

    # Samplesheets / sample_metadata
    if fname.startswith("sample_metadata") and ext == "csv":
        return "samplesheet"
    if "/samplesheets/" in path:
        return "samplesheet"

    # Scripts and docs
    if ext in SCRIPT_EXTS:
        return "script"
    if ext in DOC_EXTS or fname == "README":
        return "doc"

    # Plots
    if ext in {"pdf", "png", "svg", "jpg", "jpeg"}:
        return "plot"

    return "other"


def fmt_size(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}TB"


def git_tracked(dataset_id: str) -> list[tuple[str, int]]:
    """Return list of (path, size_bytes) tracked under datasets/<dataset_id>/."""
    out = subprocess.check_output(
        [os.environ.get("GIT_BIN", "/cm/shared/apps/git/2.33.1/bin/git"), "-C", str(REPO), "ls-tree", "-r", "--long", "HEAD", f"datasets/{dataset_id}/"],
        stderr=subprocess.DEVNULL,
    ).decode()
    files = []
    for line in out.splitlines():
        parts = line.split(maxsplit=4)
        if len(parts) < 5:
            continue
        try:
            sz = int(parts[3])
        except ValueError:
            continue
        files.append((parts[4], sz))
    return files


def local_tree(dataset_id: str, max_depth: int = 3) -> list[tuple[int, str, int]]:
    """Walk the local dataset and return (depth, rel_path, n_files) per directory."""
    base = REPO / "datasets" / dataset_id
    if not base.exists():
        return []
    rows = []
    for root, dirs, files in os.walk(base, followlinks=False):
        dirs.sort()
        rel = os.path.relpath(root, base)
        depth = 0 if rel == "." else rel.count(os.sep) + 1
        if depth > max_depth:
            dirs.clear()
            continue
        # skip ignored heavy dirs from display (but still count their files)
        rows.append((depth, "." if rel == "." else rel, len(files)))
    return rows


def synapse_walk(synapse_id: str) -> list[tuple[int, str, str]]:
    """Walk a Synapse folder and return (depth, name, id) — flat list."""
    try:
        import synapseclient
        import synapseutils
    except ImportError:
        return []
    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        return []
    syn = synapseclient.Synapse(silent=True)
    syn.login(authToken=token, silent=True)
    rows = []
    for dirpath, dirnames, filenames in synapseutils.walk(syn, synapse_id):
        # dirpath is (name, syn_id) per synapseutils.walk
        depth = dirpath[0].count("/")
        rows.append((depth, dirpath[0], dirpath[1]))
        for fname, fid in filenames:
            rows.append((depth + 1, dirpath[0] + "/" + fname, fid))
    return rows


GCS_ROOT = "gs://igvf-pertub-seq-pipeline-data"


def gcs_pipeline_config(dataset_id: str) -> list[str]:
    """List the top-level files at GCS pipeline_config for this dataset (config + samplesheets)."""
    try:
        out = subprocess.check_output(
            ["gsutil", "ls", f"{GCS_ROOT}/{dataset_id}/"],
            stderr=subprocess.DEVNULL,
            timeout=30,
        ).decode()
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
        return []
    return [line.strip() for line in out.splitlines() if line.strip()]


def lookup_synapse_root(dataset_id: str) -> dict[str, str]:
    """Read synapse_paths.tsv and return analysis_col → syn_id mapping for this dataset."""
    if not SYNAPSE_PATHS_TSV.exists():
        return {}
    rows = SYNAPSE_PATHS_TSV.read_text().splitlines()
    header = rows[0].split("\t")
    for row in rows[1:]:
        cells = row.split("\t")
        if not cells or cells[0] != dataset_id:
            continue
        return {h: v for h, v in zip(header, cells) if v and v != "-"}
    return {}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset_id")
    parser.add_argument("--synapse", action="store_true", help="Query Synapse (slow)")
    parser.add_argument("--top", type=int, default=15, help="Top-N largest per category")
    args = parser.parse_args()

    ds = args.dataset_id
    print(f"# Audit: `{ds}`\n")

    # --- Local ---
    print("## Local structure (HPC)\n")
    print("```")
    for depth, rel, n in local_tree(ds, max_depth=3):
        indent = "  " * depth
        print(f"{indent}{rel}/  ({n} files)")
    print("```\n")

    # --- GitHub ---
    print("## On GitHub (`origin/main`)\n")
    tracked = git_tracked(ds)
    by_cat: dict[str, list[tuple[str, int]]] = {}
    total_size = 0
    for path, sz in tracked:
        cat = classify(path, sz)
        by_cat.setdefault(cat, []).append((path, sz))
        total_size += sz
    print(f"**{len(tracked)} files, {fmt_size(total_size)} total**\n")

    order = [
        ("result-calibration", "REMOVE"),
        ("result-pipeline_qc", "REMOVE"),
        ("result-energy_distance", "REMOVE"),
        ("result-cnmf", "REMOVE"),
        ("result-perturbo", "REMOVE"),
        ("binary", "REMOVE (data)"),
        ("log", "REMOVE (slurm log)"),
        ("samplesheet", "REMOVE? (canonical = GCS pipeline_config)"),
        ("plot", "REVIEW"),
        ("other", "REVIEW"),
        ("script", "KEEP"),
        ("doc", "KEEP"),
    ]
    for cat, verdict in order:
        items = by_cat.get(cat, [])
        if not items:
            continue
        cat_size = sum(s for _, s in items)
        print(f"### {cat}  ({len(items)} files, {fmt_size(cat_size)}) → **{verdict}**\n")
        for path, sz in sorted(items, key=lambda x: -x[1])[: args.top]:
            print(f"  - `{fmt_size(sz):>8}`  {path}")
        if len(items) > args.top:
            print(f"  - …({len(items) - args.top} more)")
        print()

    # --- GCS pipeline_config ---
    print("## GCS pipeline_config\n")
    gcs_files = gcs_pipeline_config(ds)
    if not gcs_files:
        print(f"_No GCS dir at `{GCS_ROOT}/{ds}/`._\n")
    else:
        print(f"`{GCS_ROOT}/{ds}/`:\n")
        for line in gcs_files:
            print(f"  - `{line.replace(GCS_ROOT + '/' + ds + '/', '')}`")
        print()

    # --- Synapse ---
    syn_map = lookup_synapse_root(ds)
    print("## Synapse mirror\n")
    if not syn_map:
        print("_No entry in `synapse_paths.tsv`._\n")
    else:
        print("From `docs/jamborees/2026_UTSW/synapse_paths.tsv`:\n")
        for col, sid in syn_map.items():
            if col in ("dataset_id", "dataset_name"):
                continue
            print(f"  - **{col}**: [{sid}](https://www.synapse.org/Synapse:{sid})")
        print()
        if args.synapse:
            print("\n### Synapse tree\n")
            for col, sid in syn_map.items():
                if col in ("dataset_id", "dataset_name"):
                    continue
                if not sid.startswith("syn"):
                    continue
                print(f"\n**{col} ({sid}):**")
                print("```")
                rows = synapse_walk(sid)
                for depth, name, fid in rows[:200]:
                    indent = "  " * depth
                    print(f"{indent}{name}  ({fid})")
                if len(rows) > 200:
                    print(f"… ({len(rows) - 200} more)")
                print("```")


if __name__ == "__main__":
    main()
