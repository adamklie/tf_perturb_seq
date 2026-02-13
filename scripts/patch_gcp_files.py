#!/usr/bin/env python3
"""
Decompress .gz files on GCS and produce a patched sample metadata CSV.

The pipeline requires uncompressed files for certain columns (seqspec, barcode_onlist,
guide_design, barcode_hashtag_map). This script:

1. Reads the input samplesheet CSV
2. Identifies all unique .gz GCS paths in the specified columns
3. Decompresses each .gz file to a patch/ directory on GCS
4. Writes a patched CSV with updated paths

Usage:
    # Dry run (show what would be done):
    python patch_gcp_files.py --input sample_metadata_gcp.csv --output patched.csv --dry-run

    # Actually decompress and patch:
    python patch_gcp_files.py --input sample_metadata_gcp.csv --output patched.csv

    # Custom patch directory:
    python patch_gcp_files.py --input sample_metadata_gcp.csv --output patched.csv \
        --patch-dir gs://bucket/dataset/date/patch
"""

import argparse
import csv
import os
import subprocess
import sys


DEFAULT_COLUMNS = [
    "seqspec",
    "barcode_onlist",
    "guide_design",
    "barcode_hashtag_map",
]


def derive_patch_dir(gcs_paths: list[str]) -> str:
    """
    Derive the patch directory from GCS paths.

    Paths follow: gs://bucket/dataset/date/...
    Returns: gs://bucket/dataset/date/patch
    """
    prefixes = set()
    for path in gcs_paths:
        # gs://bucket/dataset/date/rest...
        parts = path.replace("gs://", "").split("/")
        if len(parts) >= 3:
            # bucket / dataset / date
            prefixes.add(f"gs://{parts[0]}/{parts[1]}/{parts[2]}")

    if len(prefixes) != 1:
        raise ValueError(
            f"Cannot auto-derive patch directory: found {len(prefixes)} distinct base paths: {prefixes}. "
            "Use --patch-dir to specify explicitly."
        )

    return f"{prefixes.pop()}/patch"


def check_gcs_exists(gcs_path: str) -> bool:
    """Check if a GCS path exists."""
    try:
        result = subprocess.run(
            ["gsutil", "stat", gcs_path],
            capture_output=True,
            text=True,
            timeout=30,
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, Exception):
        return False


def decompress_gcs_file(src: str, dest: str) -> bool:
    """Decompress a .gz file on GCS via streaming."""
    try:
        # gsutil cat src.gz | gunzip | gsutil cp - dest
        cat_proc = subprocess.Popen(
            ["gsutil", "cat", src],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        gunzip_proc = subprocess.Popen(
            ["gunzip"],
            stdin=cat_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        cat_proc.stdout.close()

        cp_proc = subprocess.Popen(
            ["gsutil", "cp", "-", dest],
            stdin=gunzip_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        gunzip_proc.stdout.close()

        _, cp_err = cp_proc.communicate(timeout=300)
        cat_proc.wait()
        gunzip_proc.wait()

        return cp_proc.returncode == 0
    except Exception as e:
        print(f"    ERROR: {e}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Decompress .gz files on GCS and produce a patched sample metadata CSV",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input sample metadata CSV file with GCS paths",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output patched CSV file path",
    )
    parser.add_argument(
        "--columns",
        nargs="+",
        default=DEFAULT_COLUMNS,
        help=f"Columns to patch (default: {DEFAULT_COLUMNS})",
    )
    parser.add_argument(
        "--patch-dir",
        help="GCS patch directory (auto-derived if not specified)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without executing",
    )

    args = parser.parse_args()

    # Read input CSV
    print(f"Reading: {args.input}")
    with open(args.input, "r", newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    print(f"Found {len(rows)} rows")
    print(f"Checking columns: {args.columns}")

    # Collect unique .gz paths that need decompression
    gz_paths = set()
    all_gcs_paths = []
    for row in rows:
        for col in args.columns:
            if col not in row:
                continue
            value = row[col].strip()
            if not value or not value.startswith("gs://"):
                continue
            all_gcs_paths.append(value)
            if value.endswith(".gz"):
                gz_paths.add(value)

    if not gz_paths:
        print("\nNo .gz files found — nothing to patch.")
        print(f"Writing unmodified CSV to: {args.output}")
        if not args.dry_run:
            with open(args.output, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows)
        sys.exit(0)

    print(f"\nFound {len(gz_paths)} unique .gz files to decompress")

    # Determine patch directory
    if args.patch_dir:
        patch_dir = args.patch_dir.rstrip("/")
    else:
        patch_dir = derive_patch_dir(all_gcs_paths)

    print(f"Patch directory: {patch_dir}")

    # Build mapping: original .gz path -> patched path
    path_map = {}
    for gz_path in sorted(gz_paths):
        filename = os.path.basename(gz_path)
        if filename.endswith(".gz"):
            filename = filename[:-3]  # strip .gz
        patched_path = f"{patch_dir}/{filename}"
        path_map[gz_path] = patched_path

    # Decompress files
    print(f"\n{'=' * 60}")
    if args.dry_run:
        print("DRY RUN — no files will be decompressed")
    print(f"{'=' * 60}\n")

    skipped = 0
    decompressed = 0
    failed = 0

    for i, (src, dest) in enumerate(sorted(path_map.items()), 1):
        src_name = os.path.basename(src)
        dest_name = os.path.basename(dest)
        print(f"  [{i}/{len(path_map)}] {src_name} -> {dest_name}", end="", flush=True)

        if args.dry_run:
            print(" (dry run)")
            continue

        # Check if destination already exists
        if check_gcs_exists(dest):
            print(" SKIP (exists)")
            skipped += 1
            continue

        # Decompress
        print(" ...", end=" ", flush=True)
        if decompress_gcs_file(src, dest):
            print("OK")
            decompressed += 1
        else:
            print("FAILED")
            failed += 1

    # Summary
    print(f"\n{'=' * 60}")
    print("DECOMPRESSION SUMMARY")
    print(f"{'=' * 60}")
    print(f"Total .gz files: {len(path_map)}")
    if not args.dry_run:
        print(f"Decompressed: {decompressed}")
        print(f"Skipped (exists): {skipped}")
        print(f"Failed: {failed}")

    if failed > 0:
        print(f"\nERROR: {failed} files failed to decompress — aborting CSV generation.")
        sys.exit(1)

    # Write patched CSV
    print(f"\nWriting patched CSV: {args.output}")

    patched_rows = []
    for row in rows:
        patched_row = dict(row)
        for col in args.columns:
            if col not in patched_row:
                continue
            value = patched_row[col].strip()
            if value in path_map:
                patched_row[col] = path_map[value]
        patched_rows.append(patched_row)

    if not args.dry_run:
        with open(args.output, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(patched_rows)
        print("Done!")
    else:
        print("(dry run — no file written)")

    print(f"\n{'=' * 60}")
    print("PATH MAPPINGS")
    print(f"{'=' * 60}")
    for src, dest in sorted(path_map.items()):
        print(f"  {os.path.basename(src)} -> {dest}")


if __name__ == "__main__":
    main()
