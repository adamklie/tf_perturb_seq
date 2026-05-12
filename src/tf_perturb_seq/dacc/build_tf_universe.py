"""Build a TF Universe TSV from a harmonized guide library.

The DACC `TF Universe` spec lists every TF whose promoter is targeted by a
guide library. One row per unique TF.

Columns produced:
    gene          — Ensembl gene ID (GENCODE V43; from intended_target_name col)
    gene_symbol   — HGNC-style gene symbol (from gene_name col)
    genomic_element — element type targeted (typically `promoter` for our library)
    n_guides      — how many guides target this TF (across all elements)

Usage:
    python -m tf_perturb_seq.dacc.build_tf_universe \\
        --guides ref/guide_libraries/harmonized/harmonized_guide_file_poolabcdf_ensg.tsv \\
        --out ref/dacc/tf_universe.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def build_tf_universe(guides_path: Path) -> pd.DataFrame:
    """Return the TF Universe DataFrame from a harmonized guide library TSV."""
    guides = pd.read_csv(guides_path, sep="\t")

    expected = {"intended_target_name", "gene_name", "type", "genomic_element"}
    missing = expected - set(guides.columns)
    if missing:
        raise ValueError(
            f"Guide library {guides_path} missing required columns: {missing}. "
            "Expected a harmonized + ensg'd library (see ref/guide_libraries/harmonized/)."
        )

    targeting = guides[guides["type"] == "targeting"].copy()
    if targeting.empty:
        raise ValueError(f"No rows with type=='targeting' in {guides_path}.")

    grouped = (
        targeting.groupby(["intended_target_name", "gene_name", "genomic_element"])
        .size()
        .reset_index(name="n_guides")
        .rename(columns={
            "intended_target_name": "gene",
            "gene_name": "gene_symbol",
        })
        .sort_values(["gene_symbol", "gene"])
        .reset_index(drop=True)
    )

    return grouped[["gene", "gene_symbol", "genomic_element", "n_guides"]]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--guides", type=Path, required=True,
                   help="Path to harmonized guide library TSV (with ENSG mapping).")
    p.add_argument("--out", type=Path, required=True,
                   help="Output TSV path.")
    args = p.parse_args(argv)

    df = build_tf_universe(args.guides)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    print(f"Wrote {len(df)} TFs to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
