"""Build a DACC `Gene Universe` TSV from a cNMF HVG list + the IGVF GTF.

The DACC `Gene Universe` spec requires two columns — `gene` (Ensembl ID,
GENCODE V43) and `gene_symbol` — listing the genes used as input to gene-
program inference (the HVGs cNMF was fit on).

The reference IGVF GTF (`ref/genome/IGVFFI9573KOZR.gtf.gz`) is the canonical
source for symbol→ENSG mapping. The Hon submission to the portal previously
failed validation because 286 HVG symbols couldn't be resolved to ENSGs —
this generator surfaces those misses explicitly so they can be fixed before
submission instead of after.

Usage:
    PYTHONPATH=src uv run python -m tf_perturb_seq.dacc.build_gene_universe \\
        --hvg datasets/<ds>/<run>/cnmf/Result/Inference/Inference.overdispersed_genes.txt \\
        --gtf ref/genome/IGVFFI9573KOZR.gtf.gz \\
        --out docs/jamborees/2026_UTSW/datasets/<ds>/dacc/gene_universe.tsv \\
        --misses-out docs/jamborees/2026_UTSW/datasets/<ds>/dacc/gene_universe.misses.tsv

The HVG file is a plain one-symbol-per-line text file (cNMF default output).
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path

import pandas as pd

_ATTR_RE = re.compile(r'(\w+) "([^"]+)"')


def parse_gtf_symbol_map(gtf_path: Path) -> pd.DataFrame:
    """Return a DataFrame with `gene` (ENSG, no version) and `gene_symbol` for
    every `gene` feature in the GTF.
    """
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    rows: list[tuple[str, str]] = []
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = dict(_ATTR_RE.findall(fields[8]))
            gid = attrs.get("gene_id")
            name = attrs.get("gene_name")
            if gid is None or name is None:
                continue
            # Strip Ensembl version suffix (e.g. "ENSG00000139618.13" → "ENSG00000139618")
            gid = gid.split(".")[0]
            rows.append((gid, name))

    df = pd.DataFrame(rows, columns=["gene", "gene_symbol"]).drop_duplicates()
    if df.empty:
        raise ValueError(f"No 'gene' features parsed from {gtf_path}.")
    return df


def build_gene_universe(hvg_path: Path, gtf_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (universe_df, misses_df).

    universe_df columns: gene, gene_symbol.
    misses_df columns: gene_symbol (unmapped HVG symbols).
    """
    hvg_syms = [s.strip() for s in hvg_path.read_text().splitlines() if s.strip()]
    if not hvg_syms:
        raise ValueError(f"No HVG symbols found in {hvg_path}.")

    gtf_map = parse_gtf_symbol_map(gtf_path)
    # Build two lookups. Symbols can collide across multiple ENSGs — take the
    # first occurrence (deterministic by GTF order) for reproducibility.
    sym_to_gene = (gtf_map.drop_duplicates(subset=["gene_symbol"])
                   .set_index("gene_symbol")["gene"].to_dict())
    gene_to_sym = gtf_map.set_index("gene")["gene_symbol"].to_dict()

    rows: list[tuple[str, str]] = []
    misses: list[str] = []
    seen_inputs: set[str] = set()
    for sym in hvg_syms:
        if sym in seen_inputs:
            continue
        seen_inputs.add(sym)
        if sym in sym_to_gene:
            rows.append((sym_to_gene[sym], sym))
            continue
        # Try ENSG-with-version form
        stripped = re.sub(r"\.\d+$", "", sym)
        if stripped in gene_to_sym:
            rows.append((stripped, gene_to_sym[stripped]))
            continue
        misses.append(sym)

    universe = pd.DataFrame(rows, columns=["gene", "gene_symbol"]).drop_duplicates()
    misses_df = pd.DataFrame({"gene_symbol": misses})

    return universe.reset_index(drop=True), misses_df.reset_index(drop=True)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hvg", type=Path, required=True,
                   help="Path to cNMF overdispersed_genes.txt (one symbol per line).")
    p.add_argument("--gtf", type=Path, required=True,
                   help="Path to IGVF GTF (gzipped is OK).")
    p.add_argument("--out", type=Path, required=True,
                   help="Output gene_universe TSV path.")
    p.add_argument("--misses-out", type=Path, default=None,
                   help="Optional path to write unmapped symbols.")
    args = p.parse_args(argv)

    universe, misses = build_gene_universe(args.hvg, args.gtf)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    universe.to_csv(args.out, sep="\t", index=False)

    if args.misses_out is not None:
        misses.to_csv(args.misses_out, sep="\t", index=False)
        print(f"Wrote {len(misses)} unmapped symbols to {args.misses_out}", file=sys.stderr)

    print(
        f"Wrote {len(universe)} HVGs to {args.out} "
        f"(mapped {len(universe)}/{len(universe) + len(misses)}; "
        f"{len(misses)} unmapped)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
