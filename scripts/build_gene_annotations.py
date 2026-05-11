"""Build a flat gene annotation TSV from the IGVF GTF.

One row per gene_id (version-stripped to match perturbo outputs), with the
gene's symbol, chromosome, start, end, strand, and gene_type. Cached at
`reference/gene_annotations.tsv` so downstream scripts don't re-parse the
GTF.

Re-run when the upstream GTF (`reference/IGVFFI9573KOZR.gtf.gz`) changes.

Usage:
    python scripts/build_gene_annotations.py
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd

JAMB = Path(__file__).resolve().parent.parent / "docs/jamborees/2026_UTSW"
GTF = JAMB / "reference/IGVFFI9573KOZR.gtf.gz"
OUTPUT = JAMB / "reference/gene_annotations.tsv"

ATTR_RE = re.compile(r'(\w+) "([^"]+)"')


def main() -> None:
    rows: list[dict] = []
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = dict(ATTR_RE.findall(parts[8]))
            raw_id = attrs.get("gene_id", "")
            gene_id = raw_id.split(".")[0]  # strip version suffix
            rows.append(
                {
                    "gene_id": gene_id,
                    "gene_id_versioned": raw_id,
                    "gene_symbol": attrs.get("gene_name", ""),
                    "gene_type": attrs.get("gene_type", ""),
                    "chr": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                }
            )

    df = pd.DataFrame(rows).drop_duplicates(subset="gene_id", keep="first")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT, sep="\t", index=False)
    print(f"wrote {OUTPUT} ({len(df):,} rows × {len(df.columns)} cols)")
    print(f"  unique gene_ids: {df['gene_id'].nunique():,}")
    print(f"  unique gene_symbols: {df['gene_symbol'].nunique():,}")
    print(f"  gene_types: {df['gene_type'].value_counts().head(8).to_dict()}")


if __name__ == "__main__":
    main()
