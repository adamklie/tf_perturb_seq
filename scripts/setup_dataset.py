#!/usr/bin/env python3
"""
Setup a new dataset directory structure for the TF Perturb-seq pipeline.

This script creates the standard directory structure and copies template files
for a new dataset, pre-filling configuration values where possible.

Usage:
    python scripts/setup_dataset.py --name "Lab_CellLine-differentiation_TF-Perturb-seq" --accession IGVFDSXXXXXX

Example:
    python scripts/setup_dataset.py \
        --name "Hon_H9-neuron-differentiation_TF-Perturb-seq" \
        --accession IGVFDS1234ABCD \
        --lab "Hon" \
        --cell-line "H9" \
        --differentiation "neuron"
"""

import argparse
import os
import shutil
from pathlib import Path


def get_repo_root():
    """Get the repository root directory."""
    script_dir = Path(__file__).resolve().parent
    return script_dir.parent


def setup_dataset(
    name: str,
    accession: str,
    lab: str = None,
    cell_line: str = None,
    differentiation: str = None,
    base_dir: str = None,
    force: bool = False,
):
    """Create dataset directory structure from template."""
    repo_root = get_repo_root()
    template_dir = repo_root / "templates" / "dataset_template"
    datasets_dir = repo_root / "datasets"
    target_dir = datasets_dir / name

    if not template_dir.exists():
        raise FileNotFoundError(f"Template directory not found: {template_dir}")

    if target_dir.exists() and not force:
        raise FileExistsError(
            f"Dataset directory already exists: {target_dir}\n"
            "Use --force to overwrite existing files."
        )

    # Create directory structure
    print(f"Creating dataset: {name}")
    print(f"Target directory: {target_dir}")

    # Directories to create
    dirs_to_create = [
        "bin/1_CRISPR_pipeline/seqspec/yaml_files",
        "bin/2_qc",
        "bin/3_gene_program_discovery",
        "logs",
        "results",
    ]

    for dir_path in dirs_to_create:
        (target_dir / dir_path).mkdir(parents=True, exist_ok=True)
        print(f"  Created: {dir_path}/")

    # Determine base_dir for scripts
    if base_dir is None:
        base_dir = str(repo_root)

    # Template substitutions
    substitutions = {
        "{DATASET_NAME}": name,
        "{ACCESSION}": accession,
        "{LAB_NAME}": lab or _extract_lab_from_name(name),
        "{CELL_LINE}": cell_line or _extract_cell_line_from_name(name),
        "{DIFFERENTIATION_STATE}": differentiation or _extract_diff_from_name(name),
        "TEMPLATE_DATASET": name,
        "IGVFDSXXXXXXXX": accession,
        "/Users/adamklie/Desktop/projects/tf_perturb_seq": base_dir,
    }

    # Copy and process template files
    template_files = [
        ("README.md", "README.md"),
        ("dataset_config.yaml", "dataset_config.yaml"),
        ("bin/1_CRISPR_pipeline/1_generate_per_sample_metadata.sh", "bin/1_CRISPR_pipeline/1_generate_per_sample_metadata.sh"),
        ("bin/1_CRISPR_pipeline/2_upload_to_gcp.sh", "bin/1_CRISPR_pipeline/2_upload_to_gcp.sh"),
        ("bin/1_CRISPR_pipeline/3_run_CRISPR_pipeline.sh", "bin/1_CRISPR_pipeline/3_run_CRISPR_pipeline.sh"),
    ]

    for src_rel, dst_rel in template_files:
        src_path = template_dir / src_rel
        dst_path = target_dir / dst_rel

        if src_path.exists():
            content = src_path.read_text()

            # Apply substitutions
            for old, new in substitutions.items():
                content = content.replace(old, new)

            dst_path.write_text(content)
            print(f"  Created: {dst_rel}")

            # Make shell scripts executable
            if dst_rel.endswith(".sh"):
                os.chmod(dst_path, 0o755)

    # Create placeholder seqspec files
    seqspec_dir = target_dir / "bin/1_CRISPR_pipeline/seqspec/yaml_files"
    for seqspec_name in ["rna_seqspec.yml", "guide_seqspec.yml", "hash_seqspec.yml"]:
        seqspec_path = seqspec_dir / seqspec_name
        if not seqspec_path.exists():
            seqspec_path.write_text(f"# {seqspec_name}\n# TODO: Add seqspec configuration\n")
            print(f"  Created placeholder: bin/1_CRISPR_pipeline/seqspec/yaml_files/{seqspec_name}")

    print()
    print("Dataset setup complete!")
    print()
    print("Next steps:")
    print(f"  1. cd datasets/{name}")
    print("  2. Add seqspec YAML files to bin/1_CRISPR_pipeline/seqspec/yaml_files/")
    print("  3. Review and update dataset_config.yaml")
    print("  4. Run: ./bin/1_CRISPR_pipeline/1_generate_per_sample_metadata.sh")
    print()
    print(f"See docs/CRISPR_PIPELINE.md for detailed instructions.")


def _extract_lab_from_name(name: str) -> str:
    """Extract lab name from dataset name (first part before underscore)."""
    parts = name.split("_")
    return parts[0] if parts else "Unknown"


def _extract_cell_line_from_name(name: str) -> str:
    """Extract cell line from dataset name."""
    parts = name.split("_")
    if len(parts) >= 2:
        # Second part usually contains cell line
        cell_part = parts[1].split("-")[0]
        return cell_part
    return "Unknown"


def _extract_diff_from_name(name: str) -> str:
    """Extract differentiation state from dataset name."""
    if "benchmark" in name.lower():
        return "benchmark (undifferentiated)"
    parts = name.split("_")
    if len(parts) >= 2:
        # Look for differentiation info in second part
        for part in parts[1].split("-"):
            if "differentiation" in part.lower():
                continue
            if part.lower() not in ["tf", "perturb", "seq"]:
                return part
    return "Unknown"


def main():
    parser = argparse.ArgumentParser(
        description="Setup a new dataset directory for TF Perturb-seq pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python scripts/setup_dataset.py --name "Hon_H9-neuron_TF-Perturb-seq" --accession IGVFDS1234ABCD

    # With metadata
    python scripts/setup_dataset.py \\
        --name "Engreitz_WTC11-endothelial_TF-Perturb-seq" \\
        --accession IGVFDS5678EFGH \\
        --lab "Engreitz" \\
        --cell-line "WTC11" \\
        --differentiation "endothelial"
        """,
    )
    parser.add_argument(
        "--name",
        required=True,
        help="Dataset name (should follow convention: Lab_CellLine-state_TF-Perturb-seq)",
    )
    parser.add_argument(
        "--accession",
        required=True,
        help="IGVF Analysis Set accession (e.g., IGVFDS4761PYUO)",
    )
    parser.add_argument(
        "--lab",
        help="Lab name (auto-detected from name if not provided)",
    )
    parser.add_argument(
        "--cell-line",
        help="Cell line (auto-detected from name if not provided)",
    )
    parser.add_argument(
        "--differentiation",
        help="Differentiation state (auto-detected from name if not provided)",
    )
    parser.add_argument(
        "--base-dir",
        help="Base directory for tf_perturb_seq repo (default: auto-detect)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing files in target directory",
    )

    args = parser.parse_args()

    try:
        setup_dataset(
            name=args.name,
            accession=args.accession,
            lab=args.lab,
            cell_line=args.cell_line,
            differentiation=args.differentiation,
            base_dir=args.base_dir,
            force=args.force,
        )
    except (FileNotFoundError, FileExistsError) as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
