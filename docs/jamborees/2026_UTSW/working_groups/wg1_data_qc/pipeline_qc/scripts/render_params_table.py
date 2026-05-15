"""Render the manifest.tsv pipeline-params columns as a visual table (PDF + PNG).

Highlights cells whose value is the outlier across the four datasets.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REPO = ROOT.parents[5]

MANIFEST = ROOT / "manifests" / "manifest.tsv"
COLORS_YAML = REPO / "config" / "colors" / "production_TF-Perturb-seq.yaml"
OUT_DIR = ROOT / "results" / "cross_dataset_metrics"

# Param rows to display (in order) and a pretty label for each.
ROWS = [
    ("operator_platform", "Operator / platform"),
    ("ENABLE_DATA_HASHING", "ENABLE_DATA_HASHING"),
    ("use_igvf_reference", "use_igvf_reference"),
    ("is_10x3v3", "is_10x3v3"),
    ("reverse_complement_guides", "reverse_complement_guides"),
    ("spacer_tag", "spacer_tag"),
    ("QC_min_genes_per_cell", "QC_min_genes_per_cell"),
    ("QC_min_cells_per_gene", "QC_min_cells_per_gene"),
    ("QC_pct_mito", "QC_pct_mito"),
    ("QC_barcode_filter", "QC_barcode_filter"),
    ("Multiplicity_of_infection", "Multiplicity_of_infection"),
    ("GUIDE_ASSIGNMENT_method", "GUIDE_ASSIGNMENT_method"),
    ("GUIDE_ASSIGNMENT_capture_method", "GUIDE_ASSIGNMENT_capture_method"),
    ("INFERENCE_SCEPTRE_control_group", "INFERENCE_SCEPTRE_control_group"),
    ("INFERENCE_SCEPTRE_GENE_CHUNK_SIZE", "INFERENCE_SCEPTRE_GENE_CHUNK_SIZE"),
    ("base_container_short", "Base container"),
    ("compute", "max_cpus / max_memory"),
]


CONTAINER_ALIASES = {
    "ghcr.io/pinellolab/crispr_pipeline/conda-docker:latest": "pinellolab/crispr_pipeline\nconda-docker:latest",
    "ghcr.io/pinellolab/crispr_pipeline/conda-docker:latest (sif)": "pinellolab/crispr_pipeline\nconda-docker:latest (sif)",
    "sjiang9/conda-docker:0.3": "sjiang9/conda-docker:0.3",
}


def load_params() -> pd.DataFrame:
    df = pd.read_csv(MANIFEST, sep="\t", dtype=str).fillna("")
    df["operator_platform"] = df["params_source"].str.extract(r"\(([^)]+)\)").fillna(df["params_source"])
    df["compute"] = df["max_cpus"] + " / " + df["max_memory_GB"] + " GB"
    df["spacer_tag"] = df["spacer_tag"].replace("", '""')
    df["base_container_short"] = df["base_container"].map(lambda v: CONTAINER_ALIASES.get(v, v))
    return df


def lighten(hex_color: str, alpha: float = 0.18) -> tuple:
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255
    return (1 - alpha + alpha * r, 1 - alpha + alpha * g, 1 - alpha + alpha * b)


def is_outlier(values: list[str]) -> list[bool]:
    """A cell is an outlier if its value appears < 2 times across the row."""
    counts = {v: values.count(v) for v in set(values)}
    return [counts[v] < 2 for v in values]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_params()
    with open(COLORS_YAML) as f:
        palette = yaml.safe_load(f)["dataset_colors"]

    datasets = df["dataset"].tolist()
    short_names = df["short_name"].tolist()
    runs = df["run"].tolist()
    n_cols = len(datasets)

    # Build cell text + outlier mask.
    cell_text = []
    outlier_mask = []
    for key, _label in ROWS:
        values = df[key].tolist()
        cell_text.append(values)
        outlier_mask.append(is_outlier(values))

    row_labels = [label for _, label in ROWS]

    fig_h = 0.42 * (len(ROWS) + 3) + 0.6
    fig_w = 18.0
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_axis_off()

    col_headers = [f"{sn}\n({run})" for sn, run in zip(short_names, runs)]

    tbl = ax.table(
        cellText=cell_text,
        rowLabels=row_labels,
        colLabels=col_headers,
        cellLoc="center",
        rowLoc="right",
        loc="center",
        colWidths=[0.21] * n_cols,
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.0, 1.55)

    container_row = next(i for i, (k, _) in enumerate(ROWS) if k == "base_container_short")
    for j in range(n_cols):
        tbl[container_row + 1, j].set_height(tbl[container_row + 1, j].get_height() * 1.6)
    tbl[container_row + 1, -1].set_height(tbl[container_row + 1, 0].get_height())

    header_colors = [palette[d] for d in datasets]
    for j, color in enumerate(header_colors):
        cell = tbl[0, j]
        cell.set_facecolor(color)
        cell.set_text_props(color="white", weight="bold")
        cell.set_height(cell.get_height() * 1.5)

    for i in range(len(ROWS)):
        row = i + 1
        rowlabel_cell = tbl[row, -1]
        rowlabel_cell.set_text_props(weight="bold", family="monospace")
        rowlabel_cell.set_facecolor("#F4F4F4")
        for j in range(n_cols):
            cell = tbl[row, j]
            cell.set_text_props(family="monospace")
            if outlier_mask[i][j]:
                cell.set_facecolor(lighten(header_colors[j], alpha=0.35))
                cell.set_text_props(family="monospace", weight="bold")
            else:
                cell.set_facecolor("#FFFFFF" if i % 2 == 0 else "#FAFAFA")

    plt.suptitle(
        "TFP3 production CRISPR pipeline parameters by dataset",
        fontsize=13,
        weight="bold",
        y=0.995,
    )
    fig.text(
        0.5,
        0.02,
        "Highlighted cells differ from the majority of datasets for that parameter.  Source: "
        "docs/jamborees/2026_UTSW/working_groups/wg1_data_qc/pipeline_qc/manifests/manifest.tsv",
        ha="center",
        fontsize=8,
        style="italic",
        color="#555555",
    )

    plt.subplots_adjust(left=0.16, right=0.99, top=0.94, bottom=0.06)

    for ext in ("pdf", "png"):
        out = OUT_DIR / f"pipeline_params_table.{ext}"
        fig.savefig(out, bbox_inches="tight", dpi=200 if ext == "png" else None)
        print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
