#!/usr/bin/env python3
"""
Build the 3-column samples.tsv that drives `scripts/qc_array.sh`.

Scans datasets/<DS>/<RUN_LABEL>/crispr_pipeline/pipeline_outputs/inference_mudata.h5mu
and emits one row per existing file (or --include-missing to surface gaps).

Output schema (header mandatory):
    input <TAB> outdir <TAB> run_name

Usage:
    python3 build_samples_tsv.py \\
        --repo-root /cellar/users/aklie/projects/tf_perturb_seq \\
        --output samples.tsv

    # Restrict to specific datasets / runs:
    python3 build_samples_tsv.py --repo-root <...> --output sw.tsv \\
        --datasets Hon_WTC11-benchmark_TF-Perturb-seq \\
        --runs cleanser_800 cleanser_500

    # Audit what's missing locally:
    python3 build_samples_tsv.py --repo-root <...> --output audit.tsv --include-missing
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


MUDATA_REL = "crispr_pipeline/pipeline_outputs/inference_mudata.h5mu"
DEFAULT_QC_SUBDIR = "qc"


def short_name(dataset: str) -> str:
    """Strip the `_TF-Perturb-seq` suffix from a dataset dir name."""
    suffix = "_TF-Perturb-seq"
    if dataset.endswith(suffix):
        return dataset[: -len(suffix)]
    return dataset


def discover(
    datasets_dir: Path,
    dataset_filter: set[str] | None,
    run_filter: set[str] | None,
) -> list[tuple[Path, Path, str, bool]]:
    """
    Return list of (input_path, outdir_path, run_name, exists) tuples,
    one per (dataset × run) candidate.
    """
    rows: list[tuple[Path, Path, str, bool]] = []

    if not datasets_dir.is_dir():
        raise SystemExit(f"ERROR: datasets dir not found: {datasets_dir}")

    for ds_dir in sorted(datasets_dir.iterdir()):
        if not ds_dir.is_dir():
            continue
        ds_name = ds_dir.name
        if dataset_filter and ds_name not in dataset_filter:
            continue

        # Each subdir of the dataset that isn't `setup` or `bin` is a candidate run dir.
        for run_dir in sorted(ds_dir.iterdir()):
            if not run_dir.is_dir():
                continue
            run_name = run_dir.name
            if run_name in {"setup", "bin", "logs"} or run_name.startswith("."):
                continue
            if run_filter and run_name not in run_filter:
                continue

            mudata_path = run_dir / MUDATA_REL
            outdir = run_dir / DEFAULT_QC_SUBDIR
            tsv_run_name = f"{short_name(ds_name)}_{run_name}"
            rows.append((mudata_path, outdir, tsv_run_name, mudata_path.is_file()))

    return rows


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--repo-root", required=True, help="Repo root containing datasets/")
    ap.add_argument("--output", required=True, help="Output TSV path")
    ap.add_argument(
        "--datasets",
        nargs="+",
        help="Restrict to specific dataset dir names (basename match)",
    )
    ap.add_argument(
        "--runs",
        nargs="+",
        help="Restrict to specific run label dir names (basename match)",
    )
    ap.add_argument(
        "--include-missing",
        action="store_true",
        help="Emit rows even when inference_mudata.h5mu is missing (for audit)",
    )
    ap.add_argument(
        "--qc-subdir",
        default=DEFAULT_QC_SUBDIR,
        help=f"Outdir basename within each run dir (default: {DEFAULT_QC_SUBDIR})",
    )
    args = ap.parse_args()

    repo_root = Path(args.repo_root).resolve()
    datasets_dir = repo_root / "datasets"
    output_path = Path(args.output)

    dataset_filter = set(args.datasets) if args.datasets else None
    run_filter = set(args.runs) if args.runs else None

    rows = discover(datasets_dir, dataset_filter, run_filter)
    if args.qc_subdir != DEFAULT_QC_SUBDIR:
        rows = [(i, o.parent / args.qc_subdir, r, e) for (i, o, r, e) in rows]

    output_path.parent.mkdir(parents=True, exist_ok=True)

    n_total = len(rows)
    n_existing = sum(1 for *_, e in rows if e)
    emitted = []

    with output_path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["input", "outdir", "run_name"])
        for input_path, outdir, run_name, exists in rows:
            if not exists and not args.include_missing:
                continue
            writer.writerow([str(input_path), str(outdir), run_name])
            emitted.append((input_path, exists))

    print(f"Repo root:     {repo_root}")
    print(f"Datasets dir:  {datasets_dir}")
    if dataset_filter:
        print(f"Datasets:      {sorted(dataset_filter)}")
    if run_filter:
        print(f"Runs:          {sorted(run_filter)}")
    print()
    print(f"Candidates scanned:  {n_total}")
    print(f"  with mudata:       {n_existing}")
    print(f"  missing mudata:    {n_total - n_existing}")
    print(f"Rows written:        {len(emitted)}")
    print(f"Output:              {output_path}")
    print()

    if not args.include_missing and n_total - n_existing > 0:
        print("Skipped (missing inference_mudata.h5mu):", file=sys.stderr)
        for input_path, outdir, run_name, exists in rows:
            if not exists:
                print(f"  {input_path}", file=sys.stderr)
        print(file=sys.stderr)
        print("Re-run with --include-missing to surface these in the TSV.", file=sys.stderr)

    if emitted:
        print("Next:")
        print(f"  N=$(($(wc -l < {output_path}) - 1))")
        print(f"  sbatch --array=1-${{N}} {repo_root}/scripts/qc_array.sh {output_path} {repo_root}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
