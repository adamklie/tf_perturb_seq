#!/usr/bin/env python3
"""
Upload files from IGVF portal and local paths to GCP.

This script reads a sample metadata CSV file, identifies file references
(both IGVF accessions and local paths), and uploads them to GCP.
IGVF accessions are transferred directly from S3 to GCS using gcloud transfer jobs.
Local files are uploaded using gsutil.

Usage:
    python upload_to_gcp.py \
        --input sample_metadata.csv \
        --output sample_metadata_gcp.csv \
        --gcs-bucket igvf-pertub-seq-pipeline-data \
        --gcs-prefix WTC11_CM_TF_PerturbSeq/2025_01_18/
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from urllib.parse import urlparse

import requests


# IGVF Portal API base URL
IGVF_API_BASE = "https://api.data.igvf.org"

# Authentication - can be set via environment variables or command line
IGVF_API_KEY = os.environ.get("IGVF_API_KEY", "")
IGVF_SECRET_KEY = os.environ.get("IGVF_SECRET_KEY", "")

# Columns that may contain file references
DEFAULT_FILE_COLUMNS = [
    "R1_path",
    "R2_path",
    "seqspec",
    "barcode_onlist",
    "guide_design",
    "barcode_hashtag_map",
]

# Pattern to identify IGVF file accessions
IGVF_FILE_PATTERN = re.compile(r"^IGVFFI[A-Z0-9]+$")


@dataclass
class FileInfo:
    """Information about a file to be uploaded."""
    original_ref: str  # Original reference (accession or path)
    source_type: str   # "igvf" or "local"
    s3_uri: Optional[str] = None  # S3 URI for IGVF files
    s3_bucket: Optional[str] = None  # S3 bucket name
    s3_prefix: Optional[str] = None  # S3 prefix (path within bucket)
    local_path: Optional[str] = None  # Local path for local files
    filename: Optional[str] = None  # Filename
    gcs_path: Optional[str] = None  # Destination GCS path


def is_igvf_accession(value: str) -> bool:
    """Check if a value is an IGVF file accession."""
    return bool(IGVF_FILE_PATTERN.match(value.strip()))


def is_local_path(value: str) -> bool:
    """Check if a value is a local file path."""
    value = value.strip()
    return value.startswith("/") or value.startswith("./") or value.startswith("../")


def get_igvf_file_info(accession: str, api_key: str = "", secret_key: str = "") -> dict:
    """Query IGVF portal API for file information."""
    # Use provided credentials or fall back to module-level/env vars
    auth = None
    if api_key and secret_key:
        auth = (api_key, secret_key)
    elif IGVF_API_KEY and IGVF_SECRET_KEY:
        auth = (IGVF_API_KEY, IGVF_SECRET_KEY)

    # Try different file type endpoints
    file_types = ["sequence-files", "tabular-files", "reference-files", "configuration-files"]

    for file_type in file_types:
        url = f"{IGVF_API_BASE}/{file_type}/{accession}/"
        try:
            response = requests.get(url, auth=auth, timeout=30)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.HTTPError:
            continue

    raise ValueError(f"Could not find file info for accession: {accession}")


def parse_s3_uri(s3_uri: str) -> tuple[str, str]:
    """Parse S3 URI into bucket and prefix."""
    parsed = urlparse(s3_uri)
    bucket = parsed.netloc
    prefix = parsed.path.lstrip("/")
    return bucket, prefix


def collect_file_references(
    input_file: str,
    file_columns: list[str],
) -> tuple[dict[str, FileInfo], list[dict]]:
    """
    Read the input CSV and collect all unique file references.

    Returns:
        - Dictionary mapping original reference to FileInfo
        - List of rows from the input file
    """
    file_refs = {}
    rows = []

    with open(input_file, "r", newline="") as f:
        reader = csv.DictReader(f)  # CSV format (comma-separated)
        fieldnames = reader.fieldnames

        for row in reader:
            rows.append(row)
            for col in file_columns:
                if col not in row:
                    continue
                value = row[col].strip()
                if not value:
                    continue

                if value in file_refs:
                    continue

                if is_igvf_accession(value):
                    file_refs[value] = FileInfo(
                        original_ref=value,
                        source_type="igvf",
                    )
                elif is_local_path(value):
                    file_refs[value] = FileInfo(
                        original_ref=value,
                        source_type="local",
                        local_path=value,
                        filename=os.path.basename(value),
                    )

    return file_refs, rows


def resolve_igvf_files(file_refs: dict[str, FileInfo]) -> None:
    """Query IGVF portal to resolve S3 URIs for IGVF accessions."""
    igvf_refs = [ref for ref, info in file_refs.items() if info.source_type == "igvf"]

    print(f"\nResolving {len(igvf_refs)} IGVF accessions...")

    for i, ref in enumerate(igvf_refs, 1):
        print(f"  [{i}/{len(igvf_refs)}] Querying {ref}...", end=" ", flush=True)
        try:
            file_info = get_igvf_file_info(ref)
            s3_uri = file_info.get("s3_uri")

            if not s3_uri:
                print("WARNING: No S3 URI found")
                continue

            bucket, prefix = parse_s3_uri(s3_uri)
            filename = os.path.basename(prefix)

            file_refs[ref].s3_uri = s3_uri
            file_refs[ref].s3_bucket = bucket
            file_refs[ref].s3_prefix = prefix
            file_refs[ref].filename = filename

            print(f"OK ({bucket})")
        except Exception as e:
            print(f"ERROR: {e}")


def compute_gcs_paths(
    file_refs: dict[str, FileInfo],
    gcs_bucket: str,
    gcs_prefix: str,
) -> None:
    """Compute GCS destination paths for all files.

    For IGVF files, the full S3 prefix is preserved (this matches how gcloud
    transfer jobs copy files). For local files, only the filename is used.
    """
    # Ensure prefix ends with /
    if gcs_prefix and not gcs_prefix.endswith("/"):
        gcs_prefix += "/"

    for ref, info in file_refs.items():
        if info.source_type == "igvf" and info.s3_prefix:
            # IGVF files: include full S3 prefix (matches gcloud transfer behavior)
            info.gcs_path = f"gs://{gcs_bucket}/{gcs_prefix}{info.s3_prefix}"
        elif info.filename:
            # Local files: just use filename
            info.gcs_path = f"gs://{gcs_bucket}/{gcs_prefix}{info.filename}"


def print_summary(file_refs: dict[str, FileInfo]) -> None:
    """Print a summary of files to be uploaded."""
    igvf_files = [f for f in file_refs.values() if f.source_type == "igvf" and f.s3_uri]
    local_files = [f for f in file_refs.values() if f.source_type == "local"]
    failed_igvf = [f for f in file_refs.values() if f.source_type == "igvf" and not f.s3_uri]

    # Group IGVF files by S3 bucket
    by_bucket = defaultdict(list)
    for f in igvf_files:
        by_bucket[f.s3_bucket].append(f)

    print("\n" + "=" * 60)
    print("UPLOAD SUMMARY")
    print("=" * 60)

    print(f"\nIGVF Portal Files (S3 → GCS transfer): {len(igvf_files)}")
    for bucket, files in by_bucket.items():
        print(f"  From s3://{bucket}/: {len(files)} files")
        for f in files[:3]:  # Show first 3
            print(f"    - {f.original_ref} → {f.gcs_path}")
        if len(files) > 3:
            print(f"    ... and {len(files) - 3} more")

    print(f"\nLocal Files (gsutil upload): {len(local_files)}")
    for f in local_files[:3]:
        print(f"  - {f.local_path} → {f.gcs_path}")
    if len(local_files) > 3:
        print(f"  ... and {len(local_files) - 3} more")

    if failed_igvf:
        print(f"\nWARNING: {len(failed_igvf)} IGVF accessions could not be resolved:")
        for f in failed_igvf:
            print(f"  - {f.original_ref}")

    print("\n" + "=" * 60)


def confirm_upload() -> bool:
    """Ask user to confirm the upload."""
    while True:
        response = input("\nProceed with upload? [y/N]: ").strip().lower()
        if response in ("y", "yes"):
            return True
        if response in ("n", "no", ""):
            return False
        print("Please enter 'y' or 'n'")


def create_role_arn_file(output_dir: str) -> str:
    """Create the role ARN JSON file for S3 access."""
    role_arn_path = os.path.join(output_dir, "role-arn.json")
    role_arn_data = {
        "roleArn": "arn:aws:iam::407227577691:role/S3toPertubSeqGoogleCloudTransfer"
    }
    with open(role_arn_path, "w") as f:
        json.dump(role_arn_data, f, indent=2)
    return role_arn_path


def upload_igvf_files(
    file_refs: dict[str, FileInfo],
    gcs_bucket: str,
    gcs_prefix: str,
    project: str,
    role_arn_file: str,
    dry_run: bool = False,
) -> dict[str, bool]:
    """
    Upload IGVF files from S3 to GCS using gcloud transfer jobs.

    Returns dict mapping accession to success status.
    """
    igvf_files = [f for f in file_refs.values() if f.source_type == "igvf" and f.s3_uri]

    if not igvf_files:
        return {}

    # Group by S3 bucket for batch transfers
    by_bucket = defaultdict(list)
    for f in igvf_files:
        by_bucket[f.s3_bucket].append(f)

    results = {}
    gcs_dest = f"gs://{gcs_bucket}/{gcs_prefix}"

    for bucket, files in by_bucket.items():
        s3_src = f"s3://{bucket}/"

        # Create include prefixes for all files from this bucket
        include_prefixes = [f.s3_prefix for f in files]

        print(f"\nCreating transfer job from s3://{bucket}/ ({len(files)} files)...")

        # Create a unique job name
        job_name = f"igvf-upload-{bucket}-{int(time.time())}"

        cmd = [
            "gcloud", "transfer", "jobs", "create",
            s3_src, gcs_dest,
            f"--name=transferJobs/{job_name}",
            f"--description=IGVF Portal upload - {len(files)} files",
            f"--source-creds-file={role_arn_file}",
            f"--include-prefixes={','.join(include_prefixes)}",
            "--overwrite-when=never",
        ]

        if dry_run:
            print(f"  [DRY RUN] Would run: {' '.join(cmd)}")
            for f in files:
                results[f.original_ref] = True
        else:
            try:
                print(f"  Running transfer job...")
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                print(f"  Transfer job created successfully")
                for f in files:
                    results[f.original_ref] = True
            except subprocess.CalledProcessError as e:
                print(f"  ERROR: {e.stderr}")
                for f in files:
                    results[f.original_ref] = False

    return results


def upload_local_files(
    file_refs: dict[str, FileInfo],
    gcs_bucket: str,
    gcs_prefix: str,
    dry_run: bool = False,
) -> dict[str, bool]:
    """
    Upload local files to GCS using gsutil.

    Returns dict mapping local path to success status.
    """
    local_files = [f for f in file_refs.values() if f.source_type == "local"]

    if not local_files:
        return {}

    results = {}
    gcs_dest = f"gs://{gcs_bucket}/{gcs_prefix}"

    print(f"\nUploading {len(local_files)} local files to {gcs_dest}...")

    for i, f in enumerate(local_files, 1):
        print(f"  [{i}/{len(local_files)}] {f.filename}...", end=" ", flush=True)

        if not os.path.exists(f.local_path):
            print("ERROR: File not found")
            results[f.original_ref] = False
            continue

        cmd = ["gsutil", "cp", f.local_path, f.gcs_path]

        if dry_run:
            print(f"[DRY RUN] Would run: {' '.join(cmd)}")
            results[f.original_ref] = True
        else:
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                print("OK")
                results[f.original_ref] = True
            except subprocess.CalledProcessError as e:
                print(f"ERROR: {e.stderr}")
                results[f.original_ref] = False

    return results


def update_metadata_file(
    input_file: str,
    output_file: str,
    rows: list[dict],
    file_refs: dict[str, FileInfo],
    file_columns: list[str],
    upload_results: dict[str, bool],
) -> None:
    """Update the metadata file with GCS paths."""
    # Get fieldnames from input file
    with open(input_file, "r") as f:
        reader = csv.DictReader(f)  # CSV format
        fieldnames = reader.fieldnames

    updated_rows = []
    for row in rows:
        new_row = row.copy()
        for col in file_columns:
            if col not in row:
                continue
            value = row[col].strip()
            if not value:
                continue

            if value in file_refs and file_refs[value].gcs_path:
                # Only update if upload was successful
                if upload_results.get(value, False):
                    new_row[col] = file_refs[value].gcs_path

        updated_rows.append(new_row)

    # Write output file (CSV format)
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(updated_rows)

    print(f"\nUpdated metadata written to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Upload files from IGVF portal and local paths to GCP",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input sample metadata CSV file",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output sample metadata CSV file with updated GCS paths",
    )
    parser.add_argument(
        "--gcs-bucket",
        required=True,
        help="GCS bucket name (without gs:// prefix)",
    )
    parser.add_argument(
        "--gcs-prefix",
        required=True,
        help="GCS prefix/folder path within the bucket",
    )
    parser.add_argument(
        "--project",
        default="igvf-pertub-seq-pipeline",
        help="GCP project ID (default: igvf-pertub-seq-pipeline)",
    )
    parser.add_argument(
        "--columns",
        nargs="+",
        default=DEFAULT_FILE_COLUMNS,
        help=f"Columns to process (default: {DEFAULT_FILE_COLUMNS})",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without actually uploading",
    )
    parser.add_argument(
        "--no-confirm",
        action="store_true",
        help="Skip confirmation prompt",
    )

    args = parser.parse_args()

    # Validate input file
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)

    # Set GCP project
    print(f"Setting GCP project to: {args.project}")
    subprocess.run(
        ["gcloud", "config", "set", "project", args.project],
        check=True,
        capture_output=True,
    )

    # Collect file references
    print(f"\nReading input file: {args.input}")
    print(f"Processing columns: {args.columns}")
    file_refs, rows = collect_file_references(args.input, args.columns)

    print(f"\nFound {len(file_refs)} unique file references:")
    igvf_count = sum(1 for f in file_refs.values() if f.source_type == "igvf")
    local_count = sum(1 for f in file_refs.values() if f.source_type == "local")
    print(f"  - IGVF accessions: {igvf_count}")
    print(f"  - Local paths: {local_count}")

    # Resolve IGVF files
    resolve_igvf_files(file_refs)

    # Compute GCS paths
    compute_gcs_paths(file_refs, args.gcs_bucket, args.gcs_prefix)

    # Print summary
    print_summary(file_refs)

    # Confirm
    if not args.no_confirm and not args.dry_run:
        if not confirm_upload():
            print("Upload cancelled.")
            sys.exit(0)

    # Create role ARN file in temp directory
    with tempfile.TemporaryDirectory() as tmpdir:
        role_arn_file = create_role_arn_file(tmpdir)

        # Upload files
        upload_results = {}

        # Upload IGVF files
        igvf_results = upload_igvf_files(
            file_refs,
            args.gcs_bucket,
            args.gcs_prefix,
            args.project,
            role_arn_file,
            dry_run=args.dry_run,
        )
        upload_results.update(igvf_results)

        # Upload local files
        local_results = upload_local_files(
            file_refs,
            args.gcs_bucket,
            args.gcs_prefix,
            dry_run=args.dry_run,
        )
        upload_results.update(local_results)

    # Update metadata file
    update_metadata_file(
        args.input,
        args.output,
        rows,
        file_refs,
        args.columns,
        upload_results,
    )

    # Summary
    success_count = sum(1 for v in upload_results.values() if v)
    fail_count = sum(1 for v in upload_results.values() if not v)

    print(f"\n{'=' * 60}")
    print("UPLOAD COMPLETE")
    print(f"{'=' * 60}")
    print(f"Successful: {success_count}")
    print(f"Failed: {fail_count}")

    if fail_count > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
