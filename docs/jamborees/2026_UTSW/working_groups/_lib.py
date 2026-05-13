"""Shared loaders for working-group examples.

All paths resolve relative to the jamboree folder (`docs/jamborees/2026_UTSW/`).
Each loader returns a pandas DataFrame. Synapse-auth helpers are lazy — the
local artifacts are what the examples need; Synapse only matters when you go
fetch the upstream comprehensive bundles.

Run from any examples.py via:

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from _lib import load_tf_cross_lineage, DATASETS, ...
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import pandas as pd

JAMBOREE_ROOT = Path(__file__).resolve().parent.parent

DATASETS: dict[str, str] = {
    "HonCM": "Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq",
    "HuangfuDE": "Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq",
    "HuangfuESC": "Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq",
    "GersbachHep": "Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq",
    "EngreitzEndo": "Engreitz_WTC11-endothelial-cells_TF-Perturb-seq",
}


def _read_tsv(rel: str, **kwargs) -> pd.DataFrame:
    path = JAMBOREE_ROOT / rel
    if not path.exists():
        raise FileNotFoundError(
            f"{rel} not found under {JAMBOREE_ROOT}. "
            "Check the working-groups README for landing status."
        )
    return pd.read_csv(path, sep="\t", **kwargs)


# -------- Reference data (cross-dataset) --------

def load_tf_metadata(simplified: bool = False) -> pd.DataFrame:
    name = "tf_metadata_simplified.tsv" if simplified else "tf_metadata.tsv"
    return _read_tsv(f"reference/{name}")


def load_experimental_metadata(simplified: bool = False) -> pd.DataFrame:
    name = "experimental_metadata_simplified.tsv" if simplified else "experimental_metadata.tsv"
    return _read_tsv(f"reference/{name}")


def load_cross_dataset_pipeline_summary() -> pd.DataFrame:
    return _read_tsv("reference/cross_dataset_pipeline_summary.tsv")


def load_cross_dataset_edistance_summary() -> pd.DataFrame:
    return _read_tsv("reference/cross_dataset_edistance_summary.tsv")


def load_gene_disease_associations() -> pd.DataFrame:
    return _read_tsv("reference/gene_disease_associations.tsv")


def load_jaspar_tf_metadata() -> pd.DataFrame:
    return _read_tsv("reference/jaspar_core_tf_metadata.tsv")


# -------- WG1 roll-ups --------

def load_qc_summary() -> pd.DataFrame:
    return _read_tsv("working_groups/wg1_data_qc/qc_summary.tsv")


def load_edistance_summary_wg1() -> pd.DataFrame:
    return _read_tsv("working_groups/wg1_data_qc/edistance_summary.tsv")


def load_tf_cross_lineage() -> pd.DataFrame:
    return _read_tsv("working_groups/wg1_data_qc/tf_cross_lineage.tsv")


# -------- WG3 roll-ups --------

def load_disease_tf_activity() -> pd.DataFrame:
    return _read_tsv("working_groups/wg3_disease_gwas/disease_tf_activity.tsv")


def load_tf_convergence_scorecard() -> pd.DataFrame:
    return _read_tsv("working_groups/wg3_disease_gwas/tf_convergence_scorecard.tsv")


# -------- WG4 roll-ups --------

def load_network_structure_by_lineage() -> pd.DataFrame:
    return _read_tsv("working_groups/wg4_grn_inference/network_structure_by_lineage.tsv")


# -------- WG5 roll-ups --------

def load_family_activity_scorecard() -> pd.DataFrame:
    return _read_tsv("working_groups/wg5_tf_family_case_studies/family_activity_scorecard.tsv")


# -------- Per-dataset companions --------

def _dataset_dir(dataset: str) -> Path:
    if dataset in DATASETS:
        full = DATASETS[dataset]
    elif dataset in DATASETS.values():
        full = dataset
    else:
        raise KeyError(
            f"Unknown dataset '{dataset}'. Expected one of: "
            + ", ".join(DATASETS.keys())
        )
    return JAMBOREE_ROOT / "data" / full


def load_significant_tfs(dataset: str) -> pd.DataFrame:
    """Per-target energy distance + significance + TF metadata join (WG1-B detail)."""
    path = _dataset_dir(dataset) / "energy_distance" / "wg1_significant_tfs.tsv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path.relative_to(JAMBOREE_ROOT)} not landed yet — "
            "check working_groups/wg1_data_qc/README.md for status."
        )
    return pd.read_csv(path, sep="\t")


def load_trans_target_counts(dataset: str) -> pd.DataFrame:
    """Per-perturbation count of significant trans targets (WG1-E)."""
    path = _dataset_dir(dataset) / "crispr_pipeline" / "wg1_trans_target_counts.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path.relative_to(JAMBOREE_ROOT)} not landed yet.")
    return pd.read_csv(path, sep="\t")


def load_tf_gene_edges(dataset: str) -> pd.DataFrame:
    """Per-dataset TF→gene edges at per-TF BH FDR<0.05 (WG4-A)."""
    path = _dataset_dir(dataset) / "crispr_pipeline" / "wg4_tf_gene_edges_FDR05.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path.relative_to(JAMBOREE_ROOT)} not landed yet.")
    return pd.read_csv(path, sep="\t")


# -------- Synapse (optional, for fetching upstream bundles) --------

def synapse_login():
    """Authenticate to Synapse using SYNAPSE_AUTH_TOKEN env var.

    Only needed when an example needs to pull a comprehensive bundle that
    isn't mirrored locally. Returns a `synapseclient.Synapse` instance.
    """
    try:
        import synapseclient
    except ImportError as e:
        raise ImportError(
            "synapseclient not installed. `uv add synapseclient` or skip this section."
        ) from e
    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        raise RuntimeError(
            "SYNAPSE_AUTH_TOKEN not set. Get one at "
            "https://www.synapse.org/Profile:v/settings (Personal Access Tokens)."
        )
    syn = synapseclient.Synapse()
    syn.login(authToken=token, silent=True)
    return syn


# -------- Misc helpers --------

def landed_datasets(scope: str) -> list[str]:
    """Return short-names with the requested companion file landed locally.

    scope: "wg1_significant_tfs", "wg1_trans_target_counts", "wg4_tf_gene_edges".
    """
    rel = {
        "wg1_significant_tfs": "energy_distance/wg1_significant_tfs.tsv",
        "wg1_trans_target_counts": "crispr_pipeline/wg1_trans_target_counts.tsv",
        "wg4_tf_gene_edges": "crispr_pipeline/wg4_tf_gene_edges_FDR05.tsv",
    }
    if scope not in rel:
        raise KeyError(f"Unknown scope '{scope}'. Pick one of {list(rel)}.")
    return [
        short for short, full in DATASETS.items()
        if (JAMBOREE_ROOT / "data" / full / rel[scope]).exists()
    ]
