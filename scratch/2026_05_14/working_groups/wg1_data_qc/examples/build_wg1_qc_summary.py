"""Build the WG1-A cross-dataset QC summary table.

Left-joins `reference/cross_dataset_pipeline_summary.tsv` (per-dataset CRISPR
pipeline metrics, currently 3/5 datasets) onto `reference/experimental_metadata_simplified.tsv`
(identity + pipeline status, 5/5 datasets), selects a curated set of QC columns,
and writes the result to `working_groups/wg1_data_qc/examples/qc_summary.tsv`.

Re-run whenever a new dataset gets a CRISPR pipeline mirror on
Synapse. Output rows for datasets not yet in `cross_dataset_pipeline_summary.tsv`
will have empty QC columns (`data_state = pending` / `blocked`).

Usage:
    python working_groups/wg1_data_qc/examples/build_wg1_qc_summary.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

JAMB = Path(__file__).resolve().parents[3]
PIPELINE_SUMMARY = JAMB / "reference/cross_dataset_pipeline_summary.tsv"
EXP_METADATA = JAMB / "reference/experimental_metadata_simplified.tsv"
OUTPUT = JAMB / "working_groups/wg1_data_qc/examples/qc_summary.tsv"

# Output column order (grouped). Renames applied below.
IDENTITY_COLS = [
    "dataset_id", "dataset_name", "lab", "cell_line", "lineage", "n_measurement_sets",
]
QC_COLS = [
    "n_cells", "gene_umi_median", "mito_pct_median",
    "guide_umi_median", "guides_per_cell_mean", "frac_cells_with_guide",
]
KNOCKDOWN_COLS = [
    "intended_n_guides_tested", "intended_n_significant",
    "intended_frac_significant", "intended_auroc",
]
STATE_COL = ["data_state"]
STATUS_COL = ["pipeline_status"]

OUTPUT_COLUMNS = IDENTITY_COLS + STATE_COL + QC_COLS + KNOCKDOWN_COLS + STATUS_COL


def main() -> None:
    pipe = pd.read_csv(PIPELINE_SUMMARY, sep="\t")
    meta = pd.read_csv(EXP_METADATA, sep="\t").rename(columns={"differentiation": "lineage"})

    df = meta.merge(pipe, on="dataset_id", how="left")

    # data_state — derived from whether per-dataset QC metrics are present
    df["data_state"] = df["n_cells"].apply(
        lambda x: "complete" if pd.notna(x) else "pending"
    )
    # Engreitz is a different state — distinguish from "pending"
    df.loc[df["dataset_id"].str.startswith("Engreitz_"), "data_state"] = "blocked"

    out = df[OUTPUT_COLUMNS]
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT, sep="\t", index=False)
    print(f"wrote {OUTPUT} ({len(out)} rows × {len(out.columns)} cols)")


if __name__ == "__main__":
    main()
