"""Build an Element Universe BED from a harmonized guide library.

The DACC `Element Universe` spec is a BED file listing every distinct
genomic element targeted by a guide library (promoters, enhancers, etc).

Columns produced (BED6+ extension):
    chrom                  — element chromosome
    start                  — element start (0-based, half-open BED convention)
    end                    — element end
    name                   — `<genomic_element>:<symbol>` (e.g. `promoter:SOX17`)
    score                  — n guides targeting this element
    strand                 — `.` (strand not always defined for promoters)
    genomic_element        — element type (`promoter`, `enhancer`, ...)
    intended_target_gene   — ENSG of the gene the element regulates
    intended_target_symbol — gene symbol

Skips `non-targeting` guides (no genomic coordinates). Includes positive +
negative controls when they have valid coordinates (drops the synthetic
chrPC/chrNC entries used in our library for positive/negative controls).

Usage:
    PYTHONPATH=src uv run python -m tf_perturb_seq.dacc.build_element_universe \\
        --guides ref/guide_libraries/harmonized/harmonized_guide_file_poolabcdf_ensg.tsv \\
        --out docs/jamborees/2026_UTSW/reference/element_universe.bed
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def build_element_universe(guides_path: Path) -> pd.DataFrame:
    """Return the Element Universe DataFrame from a harmonized guide library TSV."""
    guides = pd.read_csv(guides_path, sep="\t")

    expected = {
        "intended_target_chr",
        "intended_target_start",
        "intended_target_end",
        "intended_target_name",
        "gene_name",
        "type",
        "genomic_element",
    }
    missing = expected - set(guides.columns)
    if missing:
        raise ValueError(
            f"Guide library {guides_path} missing required columns: {missing}."
        )

    df = guides.copy()
    df = df[df["type"] != "non-targeting"]
    df = df[df["intended_target_chr"].astype(str).str.startswith("chr", na=False)]
    df = df[~df["intended_target_chr"].isin(["chrPC", "chrNC"])]
    df = df.dropna(subset=["intended_target_chr", "intended_target_start",
                           "intended_target_end"])

    df["intended_target_start"] = df["intended_target_start"].astype(int)
    df["intended_target_end"] = df["intended_target_end"].astype(int)

    grouped = (
        df.groupby(
            ["intended_target_chr", "intended_target_start", "intended_target_end",
             "genomic_element", "intended_target_name", "gene_name"],
            dropna=False,
        )
        .size()
        .reset_index(name="n_guides")
    )

    grouped["name"] = grouped["genomic_element"].astype(str) + ":" + grouped["gene_name"].astype(str)
    grouped["strand"] = "."

    grouped = grouped.rename(columns={
        "intended_target_chr": "chrom",
        "intended_target_start": "start",
        "intended_target_end": "end",
        "n_guides": "score",
        "intended_target_name": "intended_target_gene",
        "gene_name": "intended_target_symbol",
    })

    return grouped[[
        "chrom", "start", "end", "name", "score", "strand",
        "genomic_element", "intended_target_gene", "intended_target_symbol",
    ]].sort_values(["chrom", "start", "end"]).reset_index(drop=True)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--guides", type=Path, required=True,
                   help="Path to harmonized guide library TSV (with ENSG mapping).")
    p.add_argument("--out", type=Path, required=True,
                   help="Output BED path.")
    args = p.parse_args(argv)

    df = build_element_universe(args.guides)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False, header=False)

    print(f"Wrote {len(df)} elements to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
