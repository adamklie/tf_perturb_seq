#!/usr/bin/env python3
"""
Scaffold setup/scripts/ for a new TFP3 dataset.

Generates 1_generate_per_sample_metadata.sh, 2_upload_to_gcp.sh, and
3_patch_gcp_files.sh under datasets/<name>/setup/scripts/ with the
DATASET_NAME and ACCESSION substituted in. Also creates empty
setup/configs/ and setup/samplesheets/ directories.

Idempotent: existing files are not overwritten unless --force is passed.

Usage:
    python3 scaffold_setup_scripts.py \\
        --dataset-name Hon_WTC11-newcondition_TF-Perturb-seq \\
        --accession IGVFDS12345678
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


SCRIPT_1 = r"""#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# CONFIGURATION - Only change BASE_DIR for your environment
# =============================================================================

BASE_DIR=__BASE_DIR__
DATASET_NAME=__DATASET_NAME__
ACCESSION=__ACCESSION__

# =============================================================================
# DERIVED PATHS (no need to change)
# =============================================================================

DATASET_DIR=${BASE_DIR}/datasets/${DATASET_NAME}
SCRIPT=${BASE_DIR}/src/tf_perturb_seq/portal/generate_per_sample.py

# Uncomment + edit if the portal accession is missing linked seqspecs:
# SEQSPEC_DIR=${DATASET_DIR}/bin/1_CRISPR_pipeline/seqspec/yaml_files
# SEQSPEC_FLAGS=(
#   --hash_seqspec  ${SEQSPEC_DIR}/hash_seqspec.yml
#   --rna_seqspec   ${SEQSPEC_DIR}/rna_seqspec.yml
#   --sgrna_seqspec ${SEQSPEC_DIR}/guide_seqspec.yml
# )

# =============================================================================
# RUN
# =============================================================================

echo "=========================================="
echo "Generate Per-Sample Metadata"
echo "=========================================="
echo "Dataset:    ${DATASET_NAME}"
echo "Accession:  ${ACCESSION}"
echo "Output:     ${DATASET_DIR}/setup/samplesheets/sample_metadata.csv"
echo ""

python3 ${SCRIPT} \
  --accession ${ACCESSION} \
  --output ${DATASET_DIR}/setup/samplesheets/sample_metadata.csv \
  "${SEQSPEC_FLAGS[@]:-}"

echo ""
echo "Done! Output: ${DATASET_DIR}/setup/samplesheets/sample_metadata.csv"
"""

SCRIPT_2 = r"""#!/usr/bin/env bash
#
# Upload files from IGVF portal and local paths to GCP.
#
# Usage:
#   ./2_upload_to_gcp.sh              # Normal run
#   DRY_RUN=true ./2_upload_to_gcp.sh # Dry run (recommended first)
#
set -euo pipefail

# =============================================================================
# CONFIGURATION - Only change BASE_DIR for your environment
# =============================================================================

BASE_DIR=__BASE_DIR__
DATASET_NAME=__DATASET_NAME__

PROJECT=igvf-pertub-seq-pipeline
GCS_BUCKET=igvf-pertub-seq-pipeline-data

# =============================================================================
# DERIVED PATHS (no need to change)
# =============================================================================

DATASET_DIR=${BASE_DIR}/datasets/${DATASET_NAME}
UPLOAD_SCRIPT=${BASE_DIR}/src/tf_perturb_seq/gcp/upload_to_gcp.py

INPUT_FILE=${DATASET_DIR}/setup/samplesheets/sample_metadata.csv
OUTPUT_FILE=${DATASET_DIR}/setup/samplesheets/sample_metadata_gcp_$(date +%Y_%m_%d).csv
GCS_PREFIX=${DATASET_NAME}/$(date +%Y_%m_%d)/

# Default columns: R1_path,R2_path,seqspec,barcode_onlist,guide_design,barcode_hashtag_map
# Uncomment to restrict:
# COLUMNS="R1_path R2_path"

# =============================================================================
# EXECUTION
# =============================================================================

echo "=========================================="
echo "IGVF to GCP Upload"
echo "=========================================="
echo "  Input:      ${INPUT_FILE}"
echo "  Output:     ${OUTPUT_FILE}"
echo "  GCS Dest:   gs://${GCS_BUCKET}/${GCS_PREFIX}"
echo "  Project:    ${PROJECT}"
echo ""

command -v gcloud >/dev/null || { echo "gcloud not found"; exit 1; }
command -v gsutil >/dev/null || { echo "gsutil not found"; exit 1; }
python3 -c "import requests" 2>/dev/null || { echo "Install: pip install requests"; exit 1; }

if ! gcloud auth print-access-token >/dev/null 2>&1; then
    gcloud auth login
fi
gcloud config set project "${PROJECT}"

CMD=(
    python3 "${UPLOAD_SCRIPT}"
    --input "${INPUT_FILE}"
    --output "${OUTPUT_FILE}"
    --gcs-bucket "${GCS_BUCKET}"
    --gcs-prefix "${GCS_PREFIX}"
    --project "${PROJECT}"
)
[[ -n "${COLUMNS:-}" ]] && CMD+=(--columns ${COLUMNS})
[[ "${DRY_RUN:-false}" == "true" ]] && CMD+=(--dry-run) && echo "DRY RUN"

"${CMD[@]}"

echo ""
echo "Next: review ${OUTPUT_FILE}, then run 3_patch_gcp_files.sh"
"""

SCRIPT_3 = r"""#!/usr/bin/env bash
#
# Decompress .tsv.gz on GCS (seqspec, barcode_onlist, guide_design, barcode_hashtag_map)
# and produce a patched samplesheet for the CRISPR pipeline.
#
# Usage:
#   ./3_patch_gcp_files.sh
#   DRY_RUN=true ./3_patch_gcp_files.sh
#
# IMPORTANT: update INPUT_FILE / OUTPUT_FILE dates before running.
#
set -euo pipefail

# =============================================================================
# CONFIGURATION - Only change BASE_DIR for your environment
# =============================================================================

BASE_DIR=__BASE_DIR__
DATASET_NAME=__DATASET_NAME__

# =============================================================================
# DERIVED PATHS (UPDATE DATES BELOW each run)
# =============================================================================

DATASET_DIR=${BASE_DIR}/datasets/${DATASET_NAME}
PATCH_SCRIPT=${BASE_DIR}/src/tf_perturb_seq/gcp/patch_gcp_files.py

# UPDATE THESE TO MATCH THE STEP 2 OUTPUT DATE:
INPUT_FILE=${DATASET_DIR}/setup/samplesheets/sample_metadata_gcp_YYYY_MM_DD.csv
OUTPUT_FILE=${DATASET_DIR}/setup/samplesheets/sample_metadata_gcp_YYYY_MM_DD_patched.csv

# =============================================================================
# RUN
# =============================================================================

echo "=========================================="
echo "Patch GCP Files (decompress .gz)"
echo "=========================================="
echo "Dataset:    ${DATASET_NAME}"
echo "Input:      ${INPUT_FILE}"
echo "Output:     ${OUTPUT_FILE}"
echo ""

CMD=(
    python3 "${PATCH_SCRIPT}"
    --input "${INPUT_FILE}"
    --output "${OUTPUT_FILE}"
)
[[ "${DRY_RUN:-false}" == "true" ]] && CMD+=(--dry-run) && echo "DRY RUN"

"${CMD[@]}"

echo ""
echo "Done! Patched samplesheet: ${OUTPUT_FILE}"
echo "Next: edit 4_run_CRISPR_pipeline.sh to consume this file."
"""


DATASET_NAME_RE = re.compile(r"^[A-Z][A-Za-z0-9]*_[A-Za-z0-9-]+_TF-Perturb-seq(_[A-Za-z0-9-]+)?$")
ACCESSION_RE = re.compile(r"^IGVFDS[A-Z0-9]+$")


def find_repo_root(start: Path) -> Path:
    cur = start.resolve()
    for parent in [cur, *cur.parents]:
        if (parent / "datasets").is_dir() and (parent / "src" / "tf_perturb_seq").is_dir():
            return parent
    raise RuntimeError(
        f"Could not locate repo root (looking for datasets/ + src/tf_perturb_seq/) above {start}"
    )


def render(template: str, subs: dict[str, str]) -> str:
    out = template
    for k, v in subs.items():
        out = out.replace(f"__{k}__", v)
    return out


def write_if_absent(path: Path, content: str, force: bool) -> str:
    if path.exists() and not force:
        return "skipped (exists)"
    path.write_text(content)
    path.chmod(0o755)
    return "wrote"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset-name", required=True, help="e.g. Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq")
    ap.add_argument("--accession", required=True, help="IGVF analysis set accession (IGVFDS...)")
    ap.add_argument("--repo-root", help="Repo root (default: auto-detect from CWD)")
    ap.add_argument("--base-dir", help="BASE_DIR to write into the scripts (default: same as --repo-root)")
    ap.add_argument("--force", action="store_true", help="Overwrite existing scripts")
    args = ap.parse_args()

    if not DATASET_NAME_RE.match(args.dataset_name):
        print(
            f"WARN: dataset name '{args.dataset_name}' does not match the canonical pattern "
            "<Lab>_<CellLine>-<cond>_TF-Perturb-seq[_<chem>] - continuing anyway.",
            file=sys.stderr,
        )
    if not ACCESSION_RE.match(args.accession):
        print(f"ERROR: accession '{args.accession}' is not IGVFDS<ALNUM>", file=sys.stderr)
        return 2

    repo_root = Path(args.repo_root) if args.repo_root else find_repo_root(Path.cwd())
    base_dir = args.base_dir or str(repo_root)

    dataset_dir = repo_root / "datasets" / args.dataset_name
    scripts_dir = dataset_dir / "setup" / "scripts"
    configs_dir = dataset_dir / "setup" / "configs"
    samplesheets_dir = dataset_dir / "setup" / "samplesheets"

    for d in (scripts_dir, configs_dir, samplesheets_dir):
        d.mkdir(parents=True, exist_ok=True)

    subs = {
        "BASE_DIR": base_dir,
        "DATASET_NAME": args.dataset_name,
        "ACCESSION": args.accession,
    }

    files = [
        (scripts_dir / "1_generate_per_sample_metadata.sh", SCRIPT_1),
        (scripts_dir / "2_upload_to_gcp.sh", SCRIPT_2),
        (scripts_dir / "3_patch_gcp_files.sh", SCRIPT_3),
    ]

    print(f"Repo root:    {repo_root}")
    print(f"Dataset dir:  {dataset_dir}")
    print(f"BASE_DIR:     {base_dir}")
    print(f"Accession:    {args.accession}")
    print()

    for path, tmpl in files:
        action = write_if_absent(path, render(tmpl, subs), args.force)
        print(f"  {action:20s} {path.relative_to(repo_root)}")

    print()
    print("Next steps:")
    print(f"  1. Export IGVF_API_KEY / IGVF_SECRET_KEY")
    print(f"  2. bash {scripts_dir.relative_to(repo_root)}/1_generate_per_sample_metadata.sh")
    print(f"  3. bash {scripts_dir.relative_to(repo_root)}/2_upload_to_gcp.sh   (gcloud auth login first if needed)")
    print(f"  4. Update INPUT_FILE/OUTPUT_FILE dates, then  bash {scripts_dir.relative_to(repo_root)}/3_patch_gcp_files.sh")
    return 0


if __name__ == "__main__":
    sys.exit(main())
