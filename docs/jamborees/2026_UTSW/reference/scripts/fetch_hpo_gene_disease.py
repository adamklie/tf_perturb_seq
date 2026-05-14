"""Fetch + process HPO's `genes_to_disease.txt` into a per-gene disease-association table.

Output: `docs/jamborees/2026_UTSW/reference/gene_disease_associations.tsv` —
one row per (gene_symbol), columns:
  - gene_symbol
  - ncbi_gene_id
  - n_disease_associations  (count of MONDO + OMIM disease IDs)
  - mondo_ids               (semicolon-separated MONDO:xxxxxxx, up to 10)
  - omim_ids                (semicolon-separated OMIM:xxxxxx, up to 10)
  - all_disease_ids         (semicolon-separated, full list)

Source: https://github.com/obophenotype/human-phenotype-ontology/releases —
the `genes_to_disease.txt` file is a stable artifact mapping gene symbols to
disease IDs (MONDO / OMIM / DECIPHER / ORPHA). We keep only MONDO + OMIM
(the user's preferred sources for the jamboree).

Idempotent — re-fetches only if `--force` is passed.

Usage:
    python scripts/fetch_hpo_gene_disease.py
    python scripts/fetch_hpo_gene_disease.py --force
"""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

import pandas as pd

# Script lives at docs/jamborees/2026_UTSW/reference/scripts/<script>.py
# parents: [0]=scripts [1]=reference [2]=2026_UTSW
JAMB = Path(__file__).resolve().parents[2]
RAW_PATH = JAMB / "reference/_cache/hpo_genes_to_disease.txt"
OUTPUT = JAMB / "reference/gene_disease_associations.tsv"

# Stable mirror of the HPO release artifact.
HPO_URL = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_disease.txt"


def fetch_raw(force: bool = False) -> Path:
    if RAW_PATH.is_file() and not force:
        print(f"[cache] using existing {RAW_PATH}")
        return RAW_PATH
    RAW_PATH.parent.mkdir(parents=True, exist_ok=True)
    print(f"[fetch] {HPO_URL} -> {RAW_PATH}")
    urllib.request.urlretrieve(HPO_URL, RAW_PATH)
    print(f"[fetch] done ({RAW_PATH.stat().st_size:,} B)")
    return RAW_PATH


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true", help="Re-fetch even if cache exists")
    args = ap.parse_args()

    raw = fetch_raw(force=args.force)

    # Expected schema (HPO release): ncbi_gene_id, gene_symbol, association_type, disease_id, source
    df = pd.read_csv(raw, sep="\t", comment="#")
    expected = {"gene_symbol", "disease_id"}
    if not expected.issubset(df.columns):
        sys.exit(f"unexpected schema in {raw}: got cols {list(df.columns)}, need at least {expected}")

    # Filter to MONDO + OMIM (Mondo-based source per user choice; OMIM is what HPO
    # uses for most clinical-disease associations and Mondo cross-references it).
    keep_prefix = ("MONDO:", "OMIM:")
    filt = df[df["disease_id"].str.startswith(keep_prefix, na=False)].copy()
    n_total = len(df)
    n_filt = len(filt)
    print(f"[filter] kept {n_filt:,} of {n_total:,} rows (MONDO + OMIM IDs)")

    def cap_list(items: list[str], n: int = 10) -> str:
        return ";".join(items[:n])

    grouped = (
        filt.groupby("gene_symbol")
        .agg(
            ncbi_gene_id=("ncbi_gene_id", "first"),
            n_disease_associations=("disease_id", "nunique"),
            mondo_ids=("disease_id", lambda s: cap_list(sorted({x for x in s if x.startswith("MONDO:")}))),
            omim_ids=("disease_id", lambda s: cap_list(sorted({x for x in s if x.startswith("OMIM:")}))),
            all_disease_ids=("disease_id", lambda s: cap_list(sorted(set(s)), n=1_000_000)),
        )
        .reset_index()
        .sort_values("n_disease_associations", ascending=False)
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    grouped.to_csv(OUTPUT, sep="\t", index=False)
    print(f"[write] {OUTPUT} ({len(grouped):,} unique gene symbols × {len(grouped.columns)} cols)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
