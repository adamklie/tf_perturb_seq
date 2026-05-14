"""Generate experimental metadata table(s) for the 2026 UTSW jamboree.

The bulk of these fields are hand-curated from per-dataset documentation in
`docs/DATA.md`, each dataset's `*.md` / `README.md`, and the IGVF portal. A few
columns (chemistry / spacer_tag / guide_assignment / etc.) are extracted directly
from each dataset's pipeline config file.

Outputs:
  - /tmp/experimental_metadata.tsv                              (comprehensive — for Synapse)
  - <jamboree>/reference/experimental_metadata_simplified.tsv   (simplified — kept locally)

Marker conventions in the curated table:
  - "?"  : value not yet known / needs to be filled in by the team
  - "-"  : not applicable / does not exist
  - blank : same meaning as "?" but reserved for future automated fills
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

# Script lives at docs/jamborees/2026_UTSW/reference/scripts/<script>.py
# parents: [0]=scripts [1]=reference [2]=2026_UTSW [3]=jamborees [4]=docs [5]=repo root
REPO_ROOT = Path(__file__).resolve().parents[5]
JAMB = REPO_ROOT / "docs/jamborees/2026_UTSW"
DATASETS_DIR = REPO_ROOT / "datasets"

OUT_FULL = JAMB / "reference/experimental_metadata.tsv"
OUT_SIMPLE = JAMB / "reference/experimental_metadata_simplified.tsv"

# Config-file param names that we want to surface in the metadata.
CONFIG_KEYS = [
    "ENABLE_DATA_HASHING",
    "is_10x3v3",
    "reverse_complement_guides",
    "spacer_tag",
    "GUIDE_ASSIGNMENT_method",
    "GUIDE_ASSIGNMENT_capture_method",
    "Multiplicity_of_infection",
    "DUAL_GUIDE",
]

PARAM_RE = re.compile(r"^\s*([A-Za-z0-9_]+)\s*=\s*(.+?)\s*$")


# Hand-curated metadata derived from DATA.md, per-dataset .md files, and IGVF portal.
# When in doubt, a "?" placeholder is used so the team can fill in.
DATASETS: list[dict] = [
    {
        "dataset_id": "Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq",
        "dataset_name": "Hon WTC11 Cardiomyocyte",
        "lab": "Hon",
        "cell_line": "WTC11",
        "differentiation": "Cardiomyocyte (12-day)",
        "igvf_analysis_set": "IGVFDS6332VCTO",
        "igvf_construct_library_set": "IGVFDS3299AXST",
        "igvf_guide_file": "IGVFFI8270UPKB",
        "guide_pools": "ABCDF",
        "perturbation_method": "CRISPRi",
        "assay_technology": "10x 5'-Perturb (Sigma backbone, HT-like chemistry)",
        "multiplexing_method": "HTO",
        "n_measurement_sets": 28,
        "pipeline_config": "Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq_initial_run.config",
        "pipeline_run_label": "initial_run",
        "gcs_output_path": "gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/2026_04_15/outs/initial_run",
        "pipeline_status": "Pipeline troubleshooting (Weizhou) — seqspec hash modality issues",
        "notes": "Production yaml's HTO library_spec lacks a clean tag-region position. seqspec_v2 run with stripped i7/i5 was abandoned. See dataset README for full context.",
    },
    {
        "dataset_id": "Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq",
        "dataset_name": "Huangfu HUES8 Definitive Endoderm",
        "lab": "Huangfu",
        "cell_line": "HUES8",
        "differentiation": "Definitive endoderm",
        "igvf_analysis_set": "IGVFDS9951KTRR",
        "igvf_construct_library_set": "IGVFDS3299AXST",
        "igvf_guide_file": "IGVFFI8270UPKB",
        "guide_pools": "ABCD",
        "perturbation_method": "CRISPRi",
        "assay_technology": "10x 3' v3",
        "multiplexing_method": "none",
        "n_measurement_sets": 8,
        "pipeline_config": "Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq_muddy_penguin.config",
        "pipeline_run_label": "muddy_penguin",
        "gcs_output_path": "gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/2026_04_09/outs/muddy_penguin",
        "pipeline_status": "Pipeline running (muddy_penguin run with 13bp spacer_tag GAGTACATGGGGG)",
        "notes": "Two runs compared spacer_tag length: entertaining_hamster=GAGTACATGGGG (12bp), muddy_penguin=GAGTACATGGGGG (13bp). Production data may carry an extra leading G vs the benchmark.",
    },
    {
        "dataset_id": "Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq",
        "dataset_name": "Huangfu HUES8 Embryonic Stem Cell",
        "lab": "Huangfu",
        "cell_line": "HUES8",
        "differentiation": "Embryonic stem cell (undifferentiated)",
        "igvf_analysis_set": "IGVFDS1216AEWT",
        "igvf_construct_library_set": "IGVFDS3299AXST",
        "igvf_guide_file": "IGVFFI8270UPKB",
        "guide_pools": "ABCD",
        "perturbation_method": "CRISPRi",
        "assay_technology": "10x 3' v3",
        "multiplexing_method": "none",
        "n_measurement_sets": 8,
        "pipeline_config": "Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq_2026_04_13.config",
        "pipeline_run_label": "sceptre_v1",
        "gcs_output_path": "gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/2026_04_13/outs/sceptre_v1",
        "pipeline_status": "Setting up pipeline scripts; portal had blockers (analysis set + seqspecs) audited 2026-03-25",
        "notes": "Shares construct library set + guide file with the definitive-endoderm dataset.",
    },
    {
        "dataset_id": "Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq",
        "dataset_name": "Gersbach WTC11 Hepatocyte",
        "lab": "Gersbach",
        "cell_line": "WTC11",
        "differentiation": "Hepatocyte (22-day)",
        "igvf_analysis_set": "?",
        "igvf_construct_library_set": "IGVFDS3299AXST",
        "igvf_guide_file": "IGVFFI8270UPKB",
        "guide_pools": "ABCDF",
        "perturbation_method": "CRISPRi",
        "assay_technology": "10x Perturb-seq on NovaSeq X Plus, 25B kit; R1 ~28bp + R2 ~90bp (consistent with 10x 3' v3)",
        "multiplexing_method": "none (no HTO indicators in file metadata)",
        "n_measurement_sets": 47,
        "pipeline_config": "nextflow.config",
        "pipeline_run_label": "?",
        "gcs_output_path": "?",
        "pipeline_status": "Pipeline troubleshooting (Sara). Analysis set is in progress (Gersbach team).",
        "notes": "47 measurement sets, all `in progress` status (not yet released). Sub pool `10XLane1-8_S{1..47}`, 22-day hepatocyte differentiation. Lab-specific pipeline layout — uses cleanser / direct-capture rather than sceptre / CROP-seq.",
    },
    {
        "dataset_id": "Engreitz_WTC11-endothelial-cells_TF-Perturb-seq",
        "dataset_name": "Engreitz WTC11 Endothelial",
        "lab": "Engreitz",
        "cell_line": "WTC11",
        "differentiation": "Endothelial",
        "igvf_analysis_set": "-",
        "igvf_construct_library_set": "-",
        "igvf_guide_file": "-",
        "guide_pools": "?",
        "perturbation_method": "CRISPRi",
        "assay_technology": "CC Perturb-seq (Engreitz lab; same technology as the Engreitz_WTC11-benchmark dataset)",
        "multiplexing_method": "?",
        "n_measurement_sets": 0,
        "pipeline_config": "-",
        "pipeline_run_label": "-",
        "gcs_output_path": "-",
        "pipeline_status": "Not on the IGVF portal yet; deferred — see Engreitz_WTC11-benchmark for the related benchmark dataset and CC Perturb-seq technology details.",
        "notes": "Listed in 2026_05_07_state.png but no IGVF portal entry and no entry in tf_perturb_seq/datasets/ yet. Treat as a placeholder for now.",
    },
]


def parse_config(path: Path) -> dict[str, str]:
    if not path or not path.exists():
        return {}
    out: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            m = PARAM_RE.match(line)
            if not m:
                continue
            k, v = m.group(1), m.group(2).rstrip(",").strip()
            v = v.strip("'\"")
            out[k] = v
    return out


def main() -> None:
    rows = []
    for entry in DATASETS:
        cfg_name = entry.get("pipeline_config")
        cfg_path = DATASETS_DIR / entry["dataset_id"] / cfg_name if cfg_name and cfg_name != "?" else None
        cfg_params = parse_config(cfg_path) if cfg_path else {}

        rec = dict(entry)
        for k in CONFIG_KEYS:
            rec[f"config_{k.lower()}"] = cfg_params.get(k, "?")
        rows.append(rec)

    df = pd.DataFrame(rows)

    column_order = [
        "dataset_id",
        "dataset_name",
        "lab",
        "cell_line",
        "differentiation",
        "igvf_analysis_set",
        "igvf_construct_library_set",
        "igvf_guide_file",
        "guide_pools",
        "perturbation_method",
        "assay_technology",
        "multiplexing_method",
        "n_measurement_sets",
        "pipeline_run_label",
        "pipeline_config",
        "gcs_output_path",
        "config_enable_data_hashing",
        "config_is_10x3v3",
        "config_reverse_complement_guides",
        "config_spacer_tag",
        "config_guide_assignment_method",
        "config_guide_assignment_capture_method",
        "config_multiplicity_of_infection",
        "config_dual_guide",
        "pipeline_status",
        "notes",
    ]
    df = df[column_order]

    OUT_FULL.parent.mkdir(parents=True, exist_ok=True)
    OUT_SIMPLE.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_FULL, sep="\t", index=False)
    print(f"Wrote {OUT_FULL} ({len(df)} rows × {len(df.columns)} cols)")

    simple = df[
        [
            "dataset_id",
            "dataset_name",
            "lab",
            "cell_line",
            "differentiation",
            "igvf_analysis_set",
            "guide_pools",
            "perturbation_method",
            "assay_technology",
            "multiplexing_method",
            "n_measurement_sets",
            "pipeline_status",
        ]
    ]
    simple.to_csv(OUT_SIMPLE, sep="\t", index=False)
    print(f"Wrote {OUT_SIMPLE} ({len(simple)} rows × {len(simple.columns)} cols)")

    print("\n--- preview (simplified) ---")
    print(simple.to_string(index=False))


if __name__ == "__main__":
    main()
