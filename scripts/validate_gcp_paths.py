#!/usr/bin/env python3
"""
Validate that all GCS paths in a sample metadata file exist in GCP.

This script reads a sample metadata CSV file, extracts all GCS paths,
and verifies each one exists in the bucket.

Usage:
    python validate_gcp_paths.py --input sample_metadata_gcp.CSV
    python validate_gcp_paths.py --input sample_metadata_gcp.CSV --columns R1_path R2_path
"""

import argparse
import csv
import subprocess
import sys
from collections import defaultdict


# Columns that may contain file references
DEFAULT_FILE_COLUMNS = [
    "R1_path",
    "R2_path",
    "seqspec",
    "barcode_onlist",
    "guide_design",
    "barcode_hashtag_map",
]


def is_gcs_path(value: str) -> bool:
    """Check if a value is a GCS path."""
    return value.strip().startswith("gs://")


def collect_gcs_paths(input_file: str, file_columns: list[str]) -> dict[str, list[str]]:
    """
    Read the input CSV and collect all unique GCS paths.

    Returns:
        Dictionary mapping GCS path to list of (row_num, column) tuples where it appears
    """
    gcs_paths = defaultdict(list)

    with open(input_file, "r", newline="") as f:
        reader = csv.DictReader(f)

        for row_num, row in enumerate(reader, start=2):  # Start at 2 (header is row 1)
            for col in file_columns:
                if col not in row:
                    continue
                value = row[col].strip()
                if not value:
                    continue

                if is_gcs_path(value):
                    gcs_paths[value].append((row_num, col))

    return gcs_paths


def check_gcs_path_exists(gcs_path: str) -> bool:
    """Check if a GCS path exists using gsutil."""
    try:
        result = subprocess.run(
            ["gsutil", "stat", gcs_path],
            capture_output=True,
            text=True,
            timeout=30,
        )
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        return False
    except Exception:
        return False


def check_gcs_paths_batch(gcs_paths: list[str]) -> dict[str, bool]:
    """
    Check multiple GCS paths efficiently using gsutil ls.

    Returns dict mapping path to exists boolean.
    """
    results = {}

    # Group paths by bucket/prefix for more efficient checking
    for path in gcs_paths:
        results[path] = check_gcs_path_exists(path)

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Validate GCS paths in a sample metadata file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input sample metadata CSV file with GCS paths",
    )
    parser.add_argument(
        "--columns",
        nargs="+",
        default=DEFAULT_FILE_COLUMNS,
        help=f"Columns to check (default: {DEFAULT_FILE_COLUMNS})",
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Only show errors, not progress",
    )

    args = parser.parse_args()

    # Collect GCS paths
    if not args.quiet:
        print(f"Reading: {args.input}")
        print(f"Checking columns: {args.columns}")

    gcs_paths = collect_gcs_paths(args.input, args.columns)

    if not gcs_paths:
        print("No GCS paths found in the specified columns.")
        sys.exit(0)

    if not args.quiet:
        print(f"\nFound {len(gcs_paths)} unique GCS paths to validate\n")

    # Check each path
    missing = []
    found = []

    for i, (path, locations) in enumerate(gcs_paths.items(), 1):
        if not args.quiet:
            print(f"  [{i}/{len(gcs_paths)}] Checking {path.split('/')[-1]}...", end=" ", flush=True)

        exists = check_gcs_path_exists(path)

        if exists:
            found.append(path)
            if not args.quiet:
                print("OK")
        else:
            missing.append((path, locations))
            if not args.quiet:
                print("MISSING")

    # Summary
    print(f"\n{'=' * 60}")
    print("VALIDATION SUMMARY")
    print(f"{'=' * 60}")
    print(f"Total paths checked: {len(gcs_paths)}")
    print(f"Found: {len(found)}")
    print(f"Missing: {len(missing)}")

    if missing:
        print(f"\n{'=' * 60}")
        print("MISSING FILES")
        print(f"{'=' * 60}")
        for path, locations in missing:
            print(f"\n  {path}")
            print(f"    Referenced in:")
            for row_num, col in locations:
                print(f"      - Row {row_num}, column '{col}'")

        print(f"\n{'=' * 60}")
        print("VALIDATION FAILED")
        print(f"{'=' * 60}")
        sys.exit(1)
    else:
        print(f"\n{'=' * 60}")
        print("VALIDATION PASSED - All files exist")
        print(f"{'=' * 60}")
        sys.exit(0)


if __name__ == "__main__":
    main()
