"""Build a folks-ready filtered trans-DE results table per dataset.

Modeled on IGVF's `global differential expression` shape (e.g. IGVFFI5989UAVX):
one row per significant (perturbation, target_gene) trans hit, with the
perturbed TF's identity + the target gene's symbol/location + effect size
+ p-value + FDR.

Reads each dataset's `wg4_tf_gene_edges_FDR05.tsv` (already filtered at
per-TF BH FDR<0.05) and left-joins target-gene annotations from the IGVF
GTF (cached at `reference/gene_annotations.tsv`).

Output: `datasets/<dataset>/crispr_pipeline/trans_de_results.tsv.gz`

Re-runs auto-discover datasets via the input file's presence.

Usage:
    python scripts/build_trans_de_results.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

JAMB = Path(__file__).resolve().parent.parent / "docs/jamborees/2026_UTSW"
DATASETS_DIR = JAMB / "datasets"
GENE_ANNOT = JAMB / "reference/gene_annotations.tsv"

COLUMN_ORDER = [
    "perturbation_target_id",
    "perturbation_gene_symbol",
    "perturbation_tf_family",
    "perturbation_tf_dbd",
    "perturbation_chr",
    "perturbation_start",
    "perturbation_end",
    "target_gene_id",
    "target_gene_symbol",
    "target_gene_type",
    "target_chr",
    "target_start",
    "target_end",
    "target_strand",
    "log2_fc",
    "log2_fc_std",
    "p_value",
    "fdr_bh",
]


def discover_datasets() -> list[tuple[str, Path]]:
    found = []
    for ds_dir in sorted(DATASETS_DIR.iterdir()):
        if not ds_dir.is_dir():
            continue
        candidate = ds_dir / "crispr_pipeline" / "wg4_tf_gene_edges_FDR05.tsv"
        if candidate.is_file():
            found.append((ds_dir.name, candidate))
    return found


def main() -> None:
    if not GENE_ANNOT.is_file():
        raise SystemExit(f"missing {GENE_ANNOT} — run scripts/build_gene_annotations.py first")
    annot = pd.read_csv(GENE_ANNOT, sep="\t")[
        ["gene_id", "gene_symbol", "gene_type", "chr", "start", "end", "strand"]
    ].rename(
        columns={
            "gene_id": "target_gene_id",
            "gene_symbol": "target_gene_symbol",
            "gene_type": "target_gene_type",
            "chr": "target_chr",
            "start": "target_start",
            "end": "target_end",
            "strand": "target_strand",
        }
    )

    found = discover_datasets()
    if not found:
        raise SystemExit("no wg4_tf_gene_edges_FDR05.tsv found under datasets/")

    for ds, path in found:
        edges = pd.read_csv(path, sep="\t").rename(
            columns={
                "intended_target_name": "perturbation_target_id",
                "tf_gene_symbol": "perturbation_gene_symbol",
                "tf_family": "perturbation_tf_family",
                "tf_dbd": "perturbation_tf_dbd",
                "intended_target_chr": "perturbation_chr",
                "intended_target_start": "perturbation_start",
                "intended_target_end": "perturbation_end",
                "gene_id": "target_gene_id",
            }
        )

        # Strip any version suffix on the target gene_id (defensive — they're
        # usually unversioned already in perturbo trans output).
        edges["target_gene_id"] = edges["target_gene_id"].astype(str).str.split(".").str[0]

        merged = edges.merge(annot, on="target_gene_id", how="left")

        # Sort by perturbation, then ascending FDR (most significant first per TF)
        merged = merged.sort_values(
            ["perturbation_gene_symbol", "fdr_bh"], na_position="last"
        )[COLUMN_ORDER]

        out_path = path.parent / "trans_de_results.tsv.gz"
        merged.to_csv(out_path, sep="\t", index=False, compression="gzip")

        n_total = len(merged)
        n_sym = merged["target_gene_symbol"].notna().sum()
        size_kb = out_path.stat().st_size / 1024
        print(
            f"wrote {out_path.relative_to(JAMB)}\n"
            f"  {n_total:,} rows × {len(merged.columns)} cols  ({size_kb:,.1f} KB compressed)\n"
            f"  target gene-symbol annotation: {n_sym:,} / {n_total:,} "
            f"({100 * n_sym / n_total:.1f}%)"
        )


if __name__ == "__main__":
    main()
