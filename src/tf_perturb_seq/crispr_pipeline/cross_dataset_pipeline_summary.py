"""Build a cross-dataset summary TSV of CRISPR pipeline metrics.

Reads `synapse_paths.tsv` to find each production dataset's mirrored
`crispr_pipeline/` Synapse folder, walks down to
`pipeline_dashboard/additional_qc/{gene,guide,intended_target,trans}/*_metrics.tsv`,
and aggregates a small set of per-dataset summary numbers (cell counts, knockdown
efficiency, perturbo significance counts, AUROC/AUPRC). Output is one row per
dataset, written to `reference/cross_dataset_pipeline_summary.tsv`.

Useful for Working Group 1's data-summarization role (Topic 1, Figure 1).

Designed to be run from anywhere (laptop or HPC). Downloads only the
~10 KB metrics TSVs, not the 63 GB bundle.

Authentication:
  - Synapse: SYNAPSE_AUTH_TOKEN env var. Run via `zsh -ic '...'` so the token
    from ~/.zshrc is loaded.

Usage:
    python cross_dataset_pipeline_summary.py
    python cross_dataset_pipeline_summary.py --output reference/cross_dataset_pipeline_summary.tsv
    python cross_dataset_pipeline_summary.py --cache-dir /tmp/jamboree_qc_cache
    python cross_dataset_pipeline_summary.py --datasets Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Optional

import pandas as pd
import synapseclient

JAMB = Path(__file__).resolve().parent.parent / "docs/jamborees/2026_UTSW"
SYNAPSE_PATHS = JAMB / "synapse_paths.tsv"
DEFAULT_OUTPUT = JAMB / "reference/cross_dataset_pipeline_summary.tsv"
DEFAULT_CACHE = Path("/tmp/jamboree_pipeline_qc_cache")

QC_FILES = {
    "gene": "gene_metrics.tsv",
    "guide": "guide_metrics.tsv",
    "intended_target": "intended_target_metrics.tsv",
    "trans": "trans_metrics.tsv",
}


def find_child(syn, parent_id: str, name: str) -> Optional[str]:
    for c in syn.getChildren(parent_id):
        if c["name"] == name:
            return c["id"]
    return None


def find_path(syn, parent_id: str, parts: list[str]) -> Optional[str]:
    cur = parent_id
    for p in parts:
        cur = find_child(syn, cur, p)
        if cur is None:
            return None
    return cur


def download_file(syn, file_id: str, dest_dir: Path) -> Optional[Path]:
    try:
        f = syn.get(file_id, downloadLocation=str(dest_dir), ifcollision="overwrite.local")
        return Path(f.path)
    except Exception as e:
        print(f"  warn: download failed for {file_id}: {e!r}")
        return None


def get_qc_tables(syn, crispr_parent_id: str, dataset: str, cache_dir: Path) -> dict[str, pd.DataFrame]:
    """Walk Synapse: <crispr_parent>/pipeline_dashboard/additional_qc/{kind}/{kind}_metrics.tsv."""
    tables: dict[str, pd.DataFrame] = {}
    dest = cache_dir / dataset
    dest.mkdir(parents=True, exist_ok=True)

    # pipeline_dashboard or dashboard? (Hon CM uses "dashboard", Huangfu DE/ESC use "pipeline_dashboard")
    dashboard_id = find_child(syn, crispr_parent_id, "pipeline_dashboard") or find_child(syn, crispr_parent_id, "dashboard")
    if dashboard_id is None:
        print(f"  warn: no pipeline_dashboard/ or dashboard/ under {crispr_parent_id}")
        return tables

    aq_id = find_child(syn, dashboard_id, "additional_qc")
    if aq_id is None:
        print(f"  warn: no additional_qc/ under {dashboard_id}")
        return tables

    for kind, fname in QC_FILES.items():
        kind_id = find_child(syn, aq_id, kind)
        if kind_id is None:
            print(f"  note: {kind}/ not found under additional_qc")
            continue
        file_id = find_child(syn, kind_id, fname)
        if file_id is None:
            print(f"  note: {fname} not found in {kind}/")
            continue
        local = download_file(syn, file_id, dest)
        if local and local.is_file():
            try:
                tables[kind] = pd.read_csv(local, sep="\t")
            except Exception as e:
                print(f"  warn: parse failed for {fname}: {e!r}")
    return tables


def extract_summary(dataset: str, tables: dict[str, pd.DataFrame]) -> dict:
    row = {"dataset_id": dataset}

    # gene_metrics: take the "all" row (first row where batch == "all")
    gm = tables.get("gene")
    if gm is not None and "batch" in gm.columns:
        all_row = gm[gm["batch"] == "all"]
        if not all_row.empty:
            r = all_row.iloc[0]
            row["n_cells"] = int(r.get("n_cells", 0))
            row["gene_umi_median"] = float(r.get("umi_median", 0))
            row["mito_pct_median"] = float(r.get("mito_median", 0))

    # guide_metrics: "all" row
    gum = tables.get("guide")
    if gum is not None and "batch" in gum.columns:
        all_row = gum[gum["batch"] == "all"]
        if not all_row.empty:
            r = all_row.iloc[0]
            row["guide_umi_median"] = float(r.get("guide_umi_median", 0))
            row["guides_per_cell_mean"] = float(r.get("guides_per_cell_mean", 0))
            row["frac_cells_with_guide"] = float(r.get("frac_cells_with_guide", 0))
            row["n_guides_total"] = int(r.get("n_guides_total", 0))

    # intended_target_metrics: single row
    it = tables.get("intended_target")
    if it is not None and not it.empty:
        r = it.iloc[0]
        row["intended_n_guides_tested"] = int(r.get("n_guides_tested", 0))
        row["intended_n_significant"] = int(r.get("n_significant", 0))
        row["intended_frac_significant"] = float(r.get("frac_significant", 0))
        row["intended_median_log2fc"] = float(r.get("median_log2fc", 0))
        row["intended_auroc"] = float(r.get("auroc", 0))
        row["intended_auprc"] = float(r.get("auprc", 0))

    # trans_metrics: single row
    tr = tables.get("trans")
    if tr is not None and not tr.empty:
        r = tr.iloc[0]
        row["trans_n_targeting_guides"] = int(r.get("n_targeting_guides", 0))
        row["trans_mean_significant_per_guide"] = float(r.get("mean_significant_per_guide_targeting", 0))
        row["trans_total_significant_tests"] = int(r.get("total_significant_tests", 0))
        # AUROC/AUPRC may be empty for trans — coerce
        for k in ("auroc", "auprc"):
            v = r.get(k)
            if isinstance(v, str) and not v.strip():
                v = 0.0
            row[f"trans_{k}"] = float(v) if v == v and v != "" else 0.0  # NaN-safe

    return row


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output TSV path")
    ap.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE, help="Where to cache the downloaded metrics TSVs")
    ap.add_argument("--datasets", nargs="*", help="Limit to specific dataset_id(s); default = all production with a crispr_pipeline Synapse ID")
    ap.add_argument("--synapse-paths", type=Path, default=SYNAPSE_PATHS, help="Path to synapse_paths.tsv")
    args = ap.parse_args()

    if not args.synapse_paths.is_file():
        sys.exit(f"synapse_paths.tsv not found: {args.synapse_paths}")
    sp = pd.read_csv(args.synapse_paths, sep="\t", dtype=str).fillna("")

    if "crispr_pipeline" not in sp.columns:
        sys.exit("synapse_paths.tsv has no crispr_pipeline column")

    # Pick datasets that have a real Synapse ID (skip "-" and empty)
    eligible = sp[sp["crispr_pipeline"].str.startswith("syn", na=False)]
    if args.datasets:
        eligible = eligible[eligible["dataset_id"].isin(args.datasets)]
    if eligible.empty:
        sys.exit("no eligible datasets (none with crispr_pipeline = synXXX)")

    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        sys.exit("SYNAPSE_AUTH_TOKEN env var not set. Run via `zsh -ic` so ~/.zshrc loads.")
    syn = synapseclient.Synapse()
    syn.login(authToken=token, silent=True)

    rows: list[dict] = []
    for _, ds_row in eligible.iterrows():
        dataset = ds_row["dataset_id"]
        crispr_id = ds_row["crispr_pipeline"]
        print(f"\n=== {dataset}  ({crispr_id}) ===")
        tables = get_qc_tables(syn, crispr_id, dataset, args.cache_dir)
        if not tables:
            print(f"  warn: no QC tables found")
            continue
        row = extract_summary(dataset, tables)
        rows.append(row)

    if not rows:
        sys.exit("no rows produced")

    df = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"\nwrote {args.output} ({len(df)} rows × {len(df.columns)} cols)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
