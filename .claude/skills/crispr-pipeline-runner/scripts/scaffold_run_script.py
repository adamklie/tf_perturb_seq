#!/usr/bin/env python3
"""
Scaffold a CRISPR Pipeline run: writes (or refreshes) the per-run Nextflow
.config and the 4_run_CRISPR_pipeline.sh driver under
datasets/<DATASET>/setup/{configs,scripts}/ for a new RUN_LABEL.

Behavior:
  - configs/<DATASET>_<RUN_LABEL>.config:
      - If --base-config is given, copies that as a starting point.
      - Otherwise emits a TFP3 default skeleton.
      - Never overwrites unless --force.
  - scripts/4_run_CRISPR_pipeline.sh:
      - If absent, writes the canonical template (newer pattern).
      - If present, prints a per-knob diff summary (does NOT overwrite without --force).

Usage:
    python3 scaffold_run_script.py \\
        --dataset-name Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq \\
        --run-label seqspec_v4 \\
        --data-date 2026_04_15

    # Copy + adjust an existing config:
    python3 scaffold_run_script.py \\
        --dataset-name Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq \\
        --run-label seqspec_v4 \\
        --data-date 2026_04_15 \\
        --base-config datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/setup/configs/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq_seqspec_v3.config
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path


DRIVER = r"""#!/usr/bin/env bash
set -euo pipefail

# Pin Nextflow version (Nextflow downloads this at startup)
export NXF_VER=__NXF_VER__

# GCS auth for the pipeline (service account key for igvf-pertub-seq-pipeline)
# Update this path to your local service-account JSON.
export GOOGLE_APPLICATION_CREDENTIALS=__GAC_PATH__

# Optional: Tower monitoring (https://tower.nf)
# export TOWER_ACCESS_TOKEN=...

# =============================================================================
# CONFIGURATION
# =============================================================================

DATASET_NAME=__DATASET_NAME__
BASE_DIR=__BASE_DIR__/datasets/${DATASET_NAME}

# Data date (GCS subfolder under the dataset; matches sample_metadata_gcp_<DATA_DATE>.csv)
DATA_DATE=__DATA_DATE__

# Sample metadata with GCS paths (Stage 1 output)
SAMPLE_METADATA=$BASE_DIR/setup/samplesheets/sample_metadata_gcp_${DATA_DATE}_patched.csv

# CRISPR Pipeline path (clone of https://github.com/IGVF/CRISPR_Pipeline)
PIPELINE_PATH=__PIPELINE_PATH__

# Run label (Nextflow run name) — bump this per run
RUN_LABEL=__RUN_LABEL__

# Dataset-specific config
CONFIG=$BASE_DIR/setup/configs/${DATASET_NAME}_${RUN_LABEL}.config

# Output directory on GCS
OUTDIR=gs://igvf-pertub-seq-pipeline-data/${DATASET_NAME}/${DATA_DATE}/outs/${RUN_LABEL}

# Log file with dataset name, run label, and timestamp
LOG_FILE=$BASE_DIR/logs/${DATASET_NAME}_${RUN_LABEL}_$(date +%Y%m%d_%H%M%S).log

# Run in background? Set to true to run with nohup
RUN_IN_BACKGROUND=${RUN_IN_BACKGROUND:-false}

# =============================================================================
# RUN PIPELINE
# =============================================================================

mkdir -p $BASE_DIR/logs
cd $PIPELINE_PATH

NF_CMD="nextflow run main.nf \
    -profile google \
    -c $CONFIG \
    --input $SAMPLE_METADATA \
    --outdir $OUTDIR \
    -resume \
    -with-tower"

echo "============================================="
echo " CRISPR pipeline launch"
echo "============================================="
echo " Dataset:           $DATASET_NAME"
echo " Sample metadata:   $SAMPLE_METADATA"
echo " Config:            $CONFIG"
echo " RUN_LABEL:         $RUN_LABEL"
echo " OUTDIR (GCS):      $OUTDIR"
echo " Pipeline path:     $PIPELINE_PATH"
echo "============================================="

if [ "$RUN_IN_BACKGROUND" = true ]; then
    echo "Running CRISPR pipeline in background..."
    echo "Log file: $LOG_FILE"
    echo "Monitor with: tail -f $LOG_FILE"
    nohup bash -c "$NF_CMD" > $LOG_FILE 2>&1 &
    echo "Pipeline started with PID: $!"
else
    echo "Running CRISPR pipeline in foreground..."
    echo "To run in background, use: RUN_IN_BACKGROUND=true $0"
    $NF_CMD
fi
"""


CONFIG_SKELETON = r"""// Per-run Nextflow config for __DATASET_NAME__ / __RUN_LABEL__
// Loaded via:  nextflow run main.nf -profile google -c <this file> ...
//
// Edit the params{} block below to match the chemistry/protocol.
// See .claude/skills/crispr-pipeline-runner/references/01-config-spec.md for full param docs.

params {
    input = null

    // ---- Chemistry / protocol -------------------------------------------------
    ENABLE_DATA_HASHING = false       // true for HTO datasets
    ENABLE_SCRUBLET = false
    use_igvf_reference = true
    is_10x3v3 = false                 // true for 10x 3' v3 chemistry
    reverse_complement_guides = true
    spacer_tag = ""                   // e.g. "TAGCTCTTAAAC" for Hon cardio

    DUAL_GUIDE = false
    REFERENCE_transcriptome = 'human'
    REFERENCE_gtf_download_path = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz'
    REFERENCE_gtf_local_path = '/path/to/gencode_gtf.gtf.gz'

    // ---- QC ------------------------------------------------------------------
    QC_min_genes_per_cell = 500
    QC_min_cells_per_gene = 0.05
    QC_pct_mito = 15
    QC_barcode_filter = 'knee'        // 'knee' | 'cleanser_500' | 'cleanser_800' | etc.

    // ---- Guide assignment ----------------------------------------------------
    Multiplicity_of_infection = 'high'

    GUIDE_ASSIGNMENT_method = 'sceptre'        // 'sceptre' | 'cleanser'
    GUIDE_ASSIGNMENT_capture_method = 'CROP-seq'
    GUIDE_ASSIGNMENT_cleanser_probability_threshold = 1
    GUIDE_ASSIGNMENT_SCEPTRE_probability_threshold = 'default'
    GUIDE_ASSIGNMENT_SCEPTRE_n_em_rep = 'default'

    // ---- Inference -----------------------------------------------------------
    INFERENCE_method = 'default'
    INFERENCE_target_guide_pairing_strategy = 'default'
    INFERENCE_predefined_pairs_to_test = "path/to/file.csv"
    INFERENCE_max_target_distance_bp = 1000000
    INFERENCE_SCEPTRE_side = 'both'
    INFERENCE_SCEPTRE_grna_integration_strategy = 'union'
    INFERENCE_SCEPTRE_resampling_approximation = 'skew_normal'
    INFERENCE_SCEPTRE_control_group = 'default'
    INFERENCE_SCEPTRE_resampling_mechanism = 'default'
    INFERENCE_SCEPTRE_formula_object = 'default'

    NETWORK_custom_central_nodes = 'undefined'
    NETWORK_central_nodes_num = 1

    // ---- Dashboard ---------------------------------------------------------
    css = "assets/css"
    js = "assets/js"
    svg = "assets/svg"

    // ---- Resource caps -----------------------------------------------------
    max_cpus = 128
    max_memory = 256.GB

    // ---- Containers --------------------------------------------------------
    containers {
        base     = 'sjiang9/conda-docker:0.3'
        cleanser = 'ghcr.io/gersbachlab-bioinformatics/cleanser:1.2.1'
        sceptre  = 'sjiang9/sceptre-igvf:0.1'
        perturbo = 'ghcr.io/pinellolab/perturbo:sha-f3dc8ca'
        aria2    = 'biasofpriene/aria2c'
    }

    // ---- Google Cloud ------------------------------------------------------
    google_bucket  = 'gs://igvf-pertub-seq-pipeline-data'
    google_project = 'igvf-pertub-seq-pipeline'
    google_region  = 'us-central1'

    // ---- Boilerplate -------------------------------------------------------
    outdir = "./pipeline_outputs"
    publish_dir_mode = 'copy'
    email = null
    email_on_fail = null
    plaintext_email = false
    monochrome_logs = false
    hook_url = null
    version = false
    pipelines_testdata_base_path = 'https://raw.githubusercontent.com/nf-core/test-datasets/'

    config_profile_name = null
    config_profile_description = null
    custom_config_version = 'master'
    custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}"
    config_profile_contact = null
    config_profile_url = null
}

// For profiles{}, process{}, singularity{}, tower{}, report{} -
// copy from a sibling dataset's .config or upstream pipeline defaults.
// Those rarely change per-dataset.
"""


DATASET_NAME_RE = re.compile(r"^[A-Z][A-Za-z0-9]*_[A-Za-z0-9-]+_TF-Perturb-seq(_[A-Za-z0-9-]+)?$")
RUN_LABEL_RE = re.compile(r"^[A-Za-z0-9_]+$")
DATE_RE = re.compile(r"^\d{4}_\d{2}_\d{2}$")


def find_repo_root(start: Path) -> Path:
    cur = start.resolve()
    for parent in [cur, *cur.parents]:
        if (parent / "datasets").is_dir() and (parent / "src" / "tf_perturb_seq").is_dir():
            return parent
    raise RuntimeError(f"Could not locate repo root above {start}")


def render(template: str, subs: dict[str, str]) -> str:
    out = template
    for k, v in subs.items():
        out = out.replace(f"__{k}__", v)
    return out


def diff_driver(existing: Path, expected: dict[str, str]) -> list[str]:
    """Cheap line-level scan: find lines for the known knobs and report current vs expected."""
    text = existing.read_text()
    notes = []
    patterns = {
        "DATA_DATE": r"^DATA_DATE=(.+)$",
        "RUN_LABEL": r"^RUN_LABEL=(.+)$",
    }
    for knob, pat in patterns.items():
        m = re.search(pat, text, re.MULTILINE)
        if m:
            current = m.group(1).strip()
            wanted = expected.get(knob, "")
            if wanted and current != wanted:
                notes.append(f"  - {knob}:  current={current!r}  →  expected={wanted!r}")
    return notes


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset-name", required=True)
    ap.add_argument("--run-label", required=True, help="e.g. seqspec_v3, cleanser_800")
    ap.add_argument("--data-date", required=True, help="YYYY_MM_DD")
    ap.add_argument("--base-config", help="Existing .config to copy + edit instead of using the default skeleton")
    ap.add_argument("--repo-root", help="Repo root (default: auto-detect)")
    ap.add_argument("--base-dir", help="BASE_DIR to write into the driver (default: --repo-root)")
    ap.add_argument("--pipeline-path", default="/Users/adamklie/Desktop/tfp3/CRISPR_Pipeline",
                    help="Local CRISPR_Pipeline checkout (default: ~/Desktop/tfp3/CRISPR_Pipeline)")
    ap.add_argument("--nxf-ver", default="26.03.2-edge")
    ap.add_argument("--gac-path", default="/Users/adamklie/Desktop/tfp3/igvf-pertub-seq-pipeline.json",
                    help="GOOGLE_APPLICATION_CREDENTIALS path")
    ap.add_argument("--force", action="store_true", help="Overwrite existing files")
    args = ap.parse_args()

    if not DATASET_NAME_RE.match(args.dataset_name):
        print(f"WARN: dataset name '{args.dataset_name}' is not canonical.", file=sys.stderr)
    if not RUN_LABEL_RE.match(args.run_label):
        print(f"ERROR: run-label '{args.run_label}' must be [A-Za-z0-9_]+", file=sys.stderr)
        return 2
    if not DATE_RE.match(args.data_date):
        print(f"ERROR: data-date '{args.data_date}' is not YYYY_MM_DD", file=sys.stderr)
        return 2

    repo_root = Path(args.repo_root) if args.repo_root else find_repo_root(Path.cwd())
    base_dir = args.base_dir or str(repo_root)

    dataset_dir = repo_root / "datasets" / args.dataset_name
    if not dataset_dir.is_dir():
        print(f"ERROR: dataset dir not found: {dataset_dir}", file=sys.stderr)
        print(f"Run igvf-portal-staging/scripts/scaffold_setup_scripts.py first.", file=sys.stderr)
        return 2

    configs_dir = dataset_dir / "setup" / "configs"
    scripts_dir = dataset_dir / "setup" / "scripts"
    configs_dir.mkdir(parents=True, exist_ok=True)
    scripts_dir.mkdir(parents=True, exist_ok=True)

    config_path = configs_dir / f"{args.dataset_name}_{args.run_label}.config"
    driver_path = scripts_dir / "4_run_CRISPR_pipeline.sh"

    print(f"Repo root:    {repo_root}")
    print(f"Dataset dir:  {dataset_dir}")
    print(f"RUN_LABEL:    {args.run_label}")
    print(f"DATA_DATE:    {args.data_date}")
    print()

    # 1. Config
    if config_path.exists() and not args.force:
        print(f"  skipped (exists)    {config_path.relative_to(repo_root)}")
    else:
        if args.base_config:
            src = Path(args.base_config)
            if not src.is_file():
                print(f"ERROR: --base-config not found: {src}", file=sys.stderr)
                return 2
            shutil.copyfile(src, config_path)
            print(f"  copied from base    {config_path.relative_to(repo_root)}")
            print(f"                        (source: {src})")
        else:
            config_path.write_text(render(CONFIG_SKELETON, {
                "DATASET_NAME": args.dataset_name,
                "RUN_LABEL": args.run_label,
            }))
            print(f"  wrote skeleton      {config_path.relative_to(repo_root)}")
        print()
        print(f"  REVIEW: edit {config_path.relative_to(repo_root)} to match chemistry/protocol")
        print(f"          (see .claude/skills/crispr-pipeline-runner/references/01-config-spec.md)")

    # 2. Driver
    driver_subs = {
        "NXF_VER": args.nxf_ver,
        "GAC_PATH": args.gac_path,
        "DATASET_NAME": args.dataset_name,
        "BASE_DIR": base_dir,
        "DATA_DATE": args.data_date,
        "RUN_LABEL": args.run_label,
        "PIPELINE_PATH": args.pipeline_path,
    }
    if driver_path.exists() and not args.force:
        notes = diff_driver(driver_path, {"DATA_DATE": args.data_date, "RUN_LABEL": args.run_label})
        print(f"  exists (kept)       {driver_path.relative_to(repo_root)}")
        if notes:
            print(f"  EDIT before run:")
            for n in notes:
                print(n)
    else:
        driver_path.write_text(render(DRIVER, driver_subs))
        driver_path.chmod(0o755)
        print(f"  wrote driver        {driver_path.relative_to(repo_root)}")

    print()
    print("Next steps:")
    print(f"  1. Edit {config_path.relative_to(repo_root)} to match the dataset's chemistry.")
    print(f"  2. Verify NXF_VER + GOOGLE_APPLICATION_CREDENTIALS in {driver_path.relative_to(repo_root)}.")
    print(f"  3. RUN_IN_BACKGROUND=true bash {driver_path.relative_to(repo_root)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
