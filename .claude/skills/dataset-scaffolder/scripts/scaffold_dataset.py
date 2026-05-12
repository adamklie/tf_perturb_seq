#!/usr/bin/env python3
"""
Scaffold a new TFP3 dataset under datasets/<DATASET_NAME>/ with the canonical
layout.

Creates:
  datasets/<DATASET_NAME>/
    ├── README.md               (from docs/data/dataset_template/README.md, headers filled)
    ├── dataset_config.yaml     (from template, dataset_name + igvf_accession filled)
    ├── setup/
    │   ├── scripts/            (Stage 1 scripts via igvf-portal-staging scaffolder)
    │   │   ├── 1_generate_per_sample_metadata.sh
    │   │   ├── 2_upload_to_gcp.sh
    │   │   └── 3_patch_gcp_files.sh
    │   ├── configs/            (empty)
    │   └── samplesheets/       (empty)
    └── (run dirs added later, per CRISPR pipeline run)

Idempotent: existing files are skipped unless --force is passed.

Usage:
    python3 scaffold_dataset.py \\
        --dataset-name Hon_WTC11-newcondition_TF-Perturb-seq \\
        --accession IGVFDS12345678 \\
        --lab Hon \\
        --cell-line WTC11 \\
        --condition newcondition
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path


DATASET_NAME_RE = re.compile(r"^[A-Z][A-Za-z0-9]*_[A-Za-z0-9-]+_TF-Perturb-seq(_[A-Za-z0-9-]+)?$")
ACCESSION_RE = re.compile(r"^IGVFDS[A-Z0-9]+$")


def find_repo_root(start: Path) -> Path:
    cur = start.resolve()
    for parent in [cur, *cur.parents]:
        if (parent / "datasets").is_dir() and (parent / "src" / "tf_perturb_seq").is_dir():
            return parent
    raise RuntimeError(f"Could not locate repo root above {start}")


def write_if_absent(path: Path, content: str, force: bool) -> str:
    if path.exists() and not force:
        return "skipped (exists)"
    path.write_text(content)
    return "wrote"


def render_template(template: Path, subs: dict[str, str]) -> str:
    text = template.read_text()
    for k, v in subs.items():
        text = text.replace("{" + k + "}", v)
    return text


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--dataset-name", required=True)
    ap.add_argument("--accession", required=True, help="IGVFDS<ALNUM> analysis set accession")
    ap.add_argument("--lab", help="Lab name (default: parsed from dataset-name)")
    ap.add_argument("--cell-line", help="Cell line (default: parsed)")
    ap.add_argument("--condition", help="Condition slot (default: parsed)")
    ap.add_argument("--repo-root", help="Repo root (default: auto-detect)")
    ap.add_argument("--base-dir", help="BASE_DIR for Stage 1 scripts (default: --repo-root)")
    ap.add_argument("--force", action="store_true", help="Overwrite existing files")
    ap.add_argument("--skip-setup-scripts", action="store_true",
                    help="Don't invoke igvf-portal-staging scaffolder")
    args = ap.parse_args()

    # Validate
    if not DATASET_NAME_RE.match(args.dataset_name):
        print(f"WARN: dataset name '{args.dataset_name}' is not canonical pattern", file=sys.stderr)
    if not ACCESSION_RE.match(args.accession):
        print(f"ERROR: accession '{args.accession}' is not IGVFDS<ALNUM>", file=sys.stderr)
        return 2

    # Parse identity fields from dataset-name if not given
    parts = args.dataset_name.split("_", 1)
    if len(parts) >= 2:
        default_lab = parts[0]
        rest = parts[1]
        if rest.endswith("_TF-Perturb-seq") or "_TF-Perturb-seq_" in rest:
            cellcond = rest.split("_TF-Perturb-seq", 1)[0]
            if "-" in cellcond:
                default_cell, default_cond = cellcond.split("-", 1)
            else:
                default_cell, default_cond = cellcond, ""
        else:
            default_cell, default_cond = "", ""
    else:
        default_lab, default_cell, default_cond = "", "", ""

    lab = args.lab or default_lab
    cell_line = args.cell_line or default_cell
    condition = args.condition or default_cond

    repo_root = Path(args.repo_root) if args.repo_root else find_repo_root(Path.cwd())
    base_dir = args.base_dir or str(repo_root)
    template_dir = repo_root / "docs" / "data" / "dataset_template"

    if not template_dir.is_dir():
        print(f"ERROR: template dir not found: {template_dir}", file=sys.stderr)
        return 2

    dataset_dir = repo_root / "datasets" / args.dataset_name

    print(f"Repo root:    {repo_root}")
    print(f"Dataset dir:  {dataset_dir}")
    print(f"Identity:     lab={lab!r}  cell_line={cell_line!r}  condition={condition!r}")
    print(f"Accession:    {args.accession}")
    print()

    # Create directories
    dataset_dir.mkdir(parents=True, exist_ok=True)
    for sub in ("setup/scripts", "setup/configs", "setup/samplesheets"):
        (dataset_dir / sub).mkdir(parents=True, exist_ok=True)

    # README.md (lightweight string substitution; template has {KEY} placeholders)
    readme_src = template_dir / "README.md"
    readme_dst = dataset_dir / "README.md"
    subs = {
        "DATASET_NAME": args.dataset_name,
        "LAB_NAME": lab,
        "CELL_LINE": cell_line,
        "DIFFERENTIATION_STATE": condition,
        "ACCESSION": args.accession,
        "SYNAPSE_ID or N/A": "N/A",
    }
    readme_content = render_template(readme_src, subs)
    print(f"  {write_if_absent(readme_dst, readme_content, args.force):20s} {readme_dst.relative_to(repo_root)}")

    # dataset_config.yaml
    cfg_src = template_dir / "dataset_config.yaml"
    cfg_dst = dataset_dir / "dataset_config.yaml"
    cfg_content = cfg_src.read_text()
    cfg_content = re.sub(r'(dataset_name:\s*)"[^"]*"', f'\\1"{args.dataset_name}"', cfg_content)
    cfg_content = re.sub(r'(igvf_accession:\s*)"[^"]*"', f'\\1"{args.accession}"', cfg_content)
    # TFP3 uses 15% mito, not template's 20%
    cfg_content = re.sub(r'^qc_pct_mito:\s*20\s*$', 'qc_pct_mito: 15', cfg_content, flags=re.MULTILINE)
    print(f"  {write_if_absent(cfg_dst, cfg_content, args.force):20s} {cfg_dst.relative_to(repo_root)}")

    # Stage 1 scripts via sibling scaffolder
    if not args.skip_setup_scripts:
        print()
        print("Delegating to igvf-portal-staging scaffolder for setup/scripts/ ...")
        sibling = repo_root / ".claude" / "skills" / "igvf-portal-staging" / "scripts" / "scaffold_setup_scripts.py"
        if not sibling.is_file():
            print(f"  WARN: {sibling.relative_to(repo_root)} not found — Stage 1 scripts not written", file=sys.stderr)
        else:
            cmd = [
                sys.executable, str(sibling),
                "--dataset-name", args.dataset_name,
                "--accession", args.accession,
                "--repo-root", str(repo_root),
                "--base-dir", base_dir,
            ]
            if args.force:
                cmd.append("--force")
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                print(f"  WARN: sibling scaffolder failed: {e}", file=sys.stderr)

    # Empty placeholders to make directories visible in git (only if --force; otherwise idempotent skips)
    placeholders = [
        dataset_dir / "setup" / "configs" / ".gitkeep",
        dataset_dir / "setup" / "samplesheets" / ".gitkeep",
    ]
    for p in placeholders:
        if not p.exists():
            p.write_text("")
            print(f"  wrote (placeholder)  {p.relative_to(repo_root)}")

    print()
    print("Next steps:")
    print(f"  1. Review/edit {dataset_dir.relative_to(repo_root)}/README.md")
    print(f"  2. Review/edit {dataset_dir.relative_to(repo_root)}/dataset_config.yaml")
    print(f"  3. Set IGVF_API_KEY / IGVF_SECRET_KEY (~/.bashrc or .env)")
    print(f"  4. Use the igvf-portal-staging skill to run Stages 1.1-1.3")
    print(f"  5. After Stage 1.3 produces the patched samplesheet, use crispr-pipeline-runner")
    print(f"     to scaffold the per-run 4_run_CRISPR_pipeline.sh + .config")
    print(f"  6. Add a row to docs/data/DATA.md inventory")
    return 0


if __name__ == "__main__":
    sys.exit(main())
