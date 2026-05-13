"""Generate TF metadata table(s) for the 2026 UTSW jamboree.

For each target gene in the perturb-seq library, build a comprehensive row of TF
annotations by joining:

  - `ref/guide_libraries/target_genes.tsv`                                  -> target list +
                                                              n_promoter rows per gene
                                                              (continuation rows with NaN
                                                              `Gene` are mapped via the
                                                              promoter-ID prefix)
  - `ref/genome/IGVFFI9573KOZR.gtf.gz`                             -> ensembl_gene_id (current)
  - `ref/guide_libraries/harmonized/...poolabcdf_ensg.tsv`  -> alias resolution for old
                                                              symbols + ENSG fallback
  - Lambert et al. 2018 (humantfs.ccbr.utoronto.ca)         -> curated TF flag, DBD,
                                                              assessment (joined by ENSG
                                                              with HGNC-symbol fallback)
  - `jaspar_core_tf_metadata.tsv` (tax_id "9606")           -> class / family / matrix
                                                              IDs / uniprot (joined on
                                                              gene symbol, trying both
                                                              original and approved)

Outputs:
  - /tmp/tf_metadata.tsv                                  (comprehensive — for Synapse)
  - <jamboree>/reference/tf_metadata_simplified.tsv       (simplified — kept locally)
"""

from __future__ import annotations

import gzip
import io
import re
from pathlib import Path

import pandas as pd
import requests

REPO_ROOT = Path("/Users/adamklie/Desktop/tfp3/tf_perturb_seq")
JAMB = REPO_ROOT / "docs/jamborees/2026_UTSW"

TARGET_GENES = REPO_ROOT / "ref/guide_libraries/target_genes.tsv"
GTF = REPO_ROOT / "ref/genome/IGVFFI9573KOZR.gtf.gz"
HARMONIZED_GUIDES = REPO_ROOT / "ref/guide_libraries/harmonized/harmonized_guide_file_poolabcdf_ensg.tsv"
JASPAR = JAMB / "jaspar_core_tf_metadata.tsv"
LAMBERT_URL = "https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv"
LAMBERT_REFERER = "https://humantfs.ccbr.utoronto.ca/download.php"
LAMBERT_CACHE = Path("/tmp/jamboree_cache/lambert_2018.csv")
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
HGNC_CACHE = Path("/tmp/jamboree_cache/hgnc_complete_set.txt")

OUT_FULL = JAMB / "reference/tf_metadata.tsv"
OUT_SIMPLE = JAMB / "reference/tf_metadata_simplified.tsv"

GENE_LINE_RE = re.compile(r'gene_id "([^"]+)".*?gene_name "([^"]+)"')
HGNC_ID_RE = re.compile(r'hgnc_id "([^"]+)"')


# ---------- loaders ----------


def load_target_genes() -> tuple[pd.DataFrame, pd.Series]:
    """Return (deduped target rows keyed by gene_symbol, n-promoter-rows series)."""
    raw = pd.read_csv(TARGET_GENES, sep="\t")
    # Continuation rows leave Gene blank; recover via the promoter-ID prefix.
    inferred = raw["Set A, promoter ID"].str.extract(r"^([^_]+)")[0]
    raw["gene_symbol"] = raw["Gene"].fillna(inferred)
    raw = raw.dropna(subset=["gene_symbol"])
    n_promoters = raw.groupby("gene_symbol").size().rename("n_promoters_in_target_library")
    out = pd.DataFrame({"gene_symbol": sorted(raw["gene_symbol"].unique())})
    out = out.merge(n_promoters, left_on="gene_symbol", right_index=True, how="left")
    out["n_promoters_in_target_library"] = out["n_promoters_in_target_library"].astype("Int64")
    return out, n_promoters


def load_gtf_gene_map() -> pd.DataFrame:
    """Return one row per gene from the GTF: gene_name, ensembl_gene_id, hgnc_id, contig.

    Multiple entries can share a gene_name when the gene appears on alt contigs.
    For symbol-based lookup we prefer the primary chromosome (chr1..chr22, X, Y, M).
    """
    primary = re.compile(r"^chr([0-9]+|[XYM])$")
    rows: list[tuple[str, str, str | None, str]] = []
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            m = GENE_LINE_RE.search(parts[8])
            if not m:
                continue
            gene_id, gene_name = m.group(1), m.group(2)
            hgnc_match = HGNC_ID_RE.search(parts[8])
            hgnc_id = hgnc_match.group(1) if hgnc_match else None
            rows.append((gene_name, gene_id.split(".")[0], hgnc_id, parts[0]))
    df = pd.DataFrame(rows, columns=["gtf_gene_name", "ensembl_gene_id", "hgnc_id", "contig"])
    df["_primary"] = df["contig"].apply(lambda c: bool(primary.match(c)))
    df = df.sort_values(["gtf_gene_name", "_primary"], ascending=[True, False])
    df = df.drop_duplicates("gtf_gene_name").drop(columns=["_primary"])
    return df


def load_hgnc_alias_map() -> dict[str, tuple[str, str | None]]:
    """any_symbol -> (approved_symbol, hgnc_id).

    Includes:
      - approved `symbol` -> itself (so we can bridge to the GTF via hgnc_id when the
        symbol differs between HGNC and the GTF, e.g. SCAND3 (HGNC) vs ZBED9 (GTF)).
      - every value in `alias_symbol` (pipe-separated) -> approved symbol + hgnc_id.
      - every value in `prev_symbol` (pipe-separated) -> approved symbol + hgnc_id.
    `prev_symbol` and the symbol itself take precedence over `alias_symbol` on conflict.
    """
    if not HGNC_CACHE.exists():
        print(f"Fetching HGNC complete set from {HGNC_URL}")
        HGNC_CACHE.parent.mkdir(parents=True, exist_ok=True)
        r = requests.get(HGNC_URL, timeout=180)
        r.raise_for_status()
        HGNC_CACHE.write_bytes(r.content)
    df = pd.read_csv(HGNC_CACHE, sep="\t", dtype=str, low_memory=False)
    df = df[df["status"].str.lower() == "approved"]

    out: dict[str, tuple[str, str | None]] = {}
    # Lowest priority first: alias_symbol, prev_symbol, then symbol itself last.
    for col in ("alias_symbol", "prev_symbol"):
        sub = df[df[col].notna()][["symbol", "hgnc_id", col]]
        for _, r in sub.iterrows():
            for token in str(r[col]).split("|"):
                token = token.strip()
                if token:
                    out[token] = (r["symbol"], r.get("hgnc_id"))
    for _, r in df[["symbol", "hgnc_id"]].iterrows():
        out[r["symbol"]] = (r["symbol"], r.get("hgnc_id"))
    return out


def load_harmonized_alias_map() -> pd.DataFrame:
    """gene_name -> ensembl_gene_id from the harmonized guide file."""
    df = pd.read_csv(HARMONIZED_GUIDES, sep="\t")
    df = df[df["intended_target_name"].astype(str).str.startswith("ENSG")]
    df = df[["gene_name", "intended_target_name"]].drop_duplicates()
    df = df.rename(columns={"gene_name": "harmonized_gene_name", "intended_target_name": "ensembl_gene_id_harmonized"})
    return df


def load_lambert() -> pd.DataFrame:
    if LAMBERT_CACHE.exists():
        print(f"Reading Lambert 2018 from cache {LAMBERT_CACHE}")
        df = pd.read_csv(LAMBERT_CACHE)
    else:
        print(f"Fetching Lambert 2018 from {LAMBERT_URL}")
        LAMBERT_CACHE.parent.mkdir(parents=True, exist_ok=True)
        headers = {"User-Agent": "Mozilla/5.0", "Referer": LAMBERT_REFERER}
        r = requests.get(LAMBERT_URL, headers=headers, timeout=120)
        r.raise_for_status()
        LAMBERT_CACHE.write_bytes(r.content)
        df = pd.read_csv(io.StringIO(r.text))
    cols = {c.strip().lower(): c for c in df.columns}
    keep = pd.DataFrame(
        {
            "lambert_hgnc_symbol": df[cols["hgnc symbol"]],
            "lambert_ensembl_id": df[cols["ensembl id"]],
            "lambert_2018_is_tf": df[cols["is tf?"]],
            "lambert_2018_dbd": df[cols["dbd"]],
            "lambert_2018_tf_assessment": df[cols.get("tf assessment", cols["hgnc symbol"])] if "tf assessment" in cols else pd.NA,
            "lambert_2018_binding_mode": df[cols["binding mode"]] if "binding mode" in cols else pd.NA,
        }
    )
    keep = keep.dropna(subset=["lambert_hgnc_symbol"]).drop_duplicates("lambert_hgnc_symbol")
    keep["in_lambert_2018"] = True
    return keep


def load_jaspar_human() -> pd.DataFrame:
    df = pd.read_csv(JASPAR, sep="\t", dtype={"tax_id": str})
    human = df[df["tax_id"] == "9606"].copy()

    def _join(values: pd.Series) -> str:
        items = sorted({str(v) for v in values if pd.notna(v) and str(v) != ""})
        return ";".join(items)

    grouped = (
        human.groupby("name", as_index=False)
        .agg(
            jaspar_matrix_ids=("matrix_id", lambda s: ";".join(map(str, sorted(s)))),
            jaspar_base_ids=("base_id", _join),
            jaspar_tf_class=("class", _join),
            jaspar_tf_family=("family", _join),
            jaspar_uniprot_ids=("uniprot_ids", _join),
        )
        .rename(columns={"name": "jaspar_name"})
    )
    grouped["in_jaspar_core"] = True
    return grouped


# ---------- resolution ----------


def resolve_ensembl(
    targets: pd.DataFrame,
    gtf_map: pd.DataFrame,
    hgnc_alias_map: dict[str, tuple[str, str | None]],
    harmonized_map: pd.DataFrame,
) -> pd.DataFrame:
    """Resolve each target gene to (ensembl_gene_id, hgnc_approved_symbol, source).

    Resolution priority:
      1. Direct: gene_symbol matches a GTF gene_name.
      2. HGNC alias: gene_symbol -> approved_symbol or hgnc_id (via HGNC complete set),
         then look up approved_symbol or hgnc_id in the GTF.
      3. Harmonized guide file: gene_symbol -> ENSG (lab-curated, last-resort fallback).
    """
    df = targets.copy()
    df["ensembl_gene_id"] = pd.NA
    df["ensembl_gene_id_source"] = pd.NA
    df["hgnc_approved_symbol"] = pd.NA

    by_name = gtf_map.set_index("gtf_gene_name")["ensembl_gene_id"].to_dict()
    by_hgnc = (
        gtf_map.dropna(subset=["hgnc_id"]).set_index("hgnc_id")["ensembl_gene_id"].to_dict()
    )
    ensg_to_name = gtf_map.set_index("ensembl_gene_id")["gtf_gene_name"].to_dict()
    harmonized_lookup = (
        harmonized_map.set_index("harmonized_gene_name")["ensembl_gene_id_harmonized"].to_dict()
    )

    for idx, row in df.iterrows():
        sym = row["gene_symbol"]

        # 1) Direct GTF gene_name match
        if sym in by_name:
            df.at[idx, "ensembl_gene_id"] = by_name[sym]
            df.at[idx, "ensembl_gene_id_source"] = "gtf_direct"
            df.at[idx, "hgnc_approved_symbol"] = sym
            continue

        # 2) HGNC alias resolution
        if sym in hgnc_alias_map:
            approved, hgnc_id = hgnc_alias_map[sym]
            ensg = by_name.get(approved)
            if ensg is None and hgnc_id and hgnc_id in by_hgnc:
                ensg = by_hgnc[hgnc_id]
            if ensg is not None:
                df.at[idx, "ensembl_gene_id"] = ensg
                df.at[idx, "ensembl_gene_id_source"] = "hgnc_alias"
                # Use the GTF's gene_name (which may differ from HGNC's approved if the
                # GTF is older/newer than HGNC) so the approved symbol is consistent
                # with the rest of our data.
                df.at[idx, "hgnc_approved_symbol"] = ensg_to_name.get(ensg, approved)
                continue

        # 3) Harmonized guide file fallback
        if sym in harmonized_lookup:
            ensg = harmonized_lookup[sym]
            df.at[idx, "ensembl_gene_id"] = ensg
            df.at[idx, "ensembl_gene_id_source"] = "harmonized_guides"
            df.at[idx, "hgnc_approved_symbol"] = ensg_to_name.get(ensg, sym)
            continue

    # Approved symbol falls back to input symbol when unresolved
    df["hgnc_approved_symbol"] = df["hgnc_approved_symbol"].fillna(df["gene_symbol"])
    return df


# ---------- main ----------


def main() -> None:
    targets, _ = load_target_genes()
    print(f"Targets: {len(targets)} unique gene symbols")

    gtf_map = load_gtf_gene_map()
    print(f"GTF unique gene names: {len(gtf_map)}")

    hgnc_alias_map = load_hgnc_alias_map()
    print(f"HGNC alias map: {len(hgnc_alias_map)} alias/prev -> approved entries")

    harmonized_map = load_harmonized_alias_map()
    print(f"Harmonized guide map: {len(harmonized_map)} gene_name -> ENSG entries")

    lambert = load_lambert()
    print(f"Lambert 2018: {len(lambert)} unique gene symbols")

    jaspar = load_jaspar_human()
    print(f"JASPAR human: {len(jaspar)} unique gene symbols")

    df = resolve_ensembl(targets, gtf_map, hgnc_alias_map, harmonized_map)

    # ---- Lambert join ----
    # Prefer Ensembl-ID match; fall back to symbol (approved then original).
    lam_by_ensg = lambert.set_index("lambert_ensembl_id")
    lam_by_sym = lambert.set_index("lambert_hgnc_symbol")
    lam_cols = ["lambert_2018_is_tf", "lambert_2018_dbd",
                "lambert_2018_tf_assessment", "lambert_2018_binding_mode"]
    df["in_lambert_2018"] = False
    for c in lam_cols:
        df[c] = pd.NA
    for idx, row in df.iterrows():
        match = None
        ensg = row["ensembl_gene_id"]
        if pd.notna(ensg) and ensg in lam_by_ensg.index:
            match = lam_by_ensg.loc[ensg]
        else:
            for sym in (row["hgnc_approved_symbol"], row["gene_symbol"]):
                if pd.notna(sym) and sym in lam_by_sym.index:
                    match = lam_by_sym.loc[sym]
                    break
        if match is not None:
            df.at[idx, "in_lambert_2018"] = True
            for c in lam_cols:
                df.at[idx, c] = match[c]

    # ---- JASPAR join ----
    # Symbol match: try approved symbol then original.
    jas_by_sym = jaspar.set_index("jaspar_name")
    jas_cols = ["jaspar_matrix_ids", "jaspar_base_ids",
                "jaspar_tf_class", "jaspar_tf_family", "jaspar_uniprot_ids"]
    df["in_jaspar_core"] = False
    for c in jas_cols:
        df[c] = pd.NA
    for idx, row in df.iterrows():
        match = None
        for sym in (row["hgnc_approved_symbol"], row["gene_symbol"]):
            if pd.notna(sym) and sym in jas_by_sym.index:
                match = jas_by_sym.loc[sym]
                break
        if match is not None:
            df.at[idx, "in_jaspar_core"] = True
            for c in jas_cols:
                df.at[idx, c] = match[c]

    column_order = [
        "gene_symbol",
        "hgnc_approved_symbol",
        "n_promoters_in_target_library",
        "ensembl_gene_id",
        "ensembl_gene_id_source",
        "in_lambert_2018",
        "lambert_2018_is_tf",
        "lambert_2018_dbd",
        "lambert_2018_tf_assessment",
        "lambert_2018_binding_mode",
        "in_jaspar_core",
        "jaspar_matrix_ids",
        "jaspar_base_ids",
        "jaspar_tf_class",
        "jaspar_tf_family",
        "jaspar_uniprot_ids",
    ]
    df = df[column_order].sort_values("gene_symbol").reset_index(drop=True)

    OUT_FULL.parent.mkdir(parents=True, exist_ok=True)
    OUT_SIMPLE.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_FULL, sep="\t", index=False)
    print(f"\nWrote {OUT_FULL} ({len(df)} rows)")

    simple = df[
        [
            "gene_symbol",
            "hgnc_approved_symbol",
            "ensembl_gene_id",
            "n_promoters_in_target_library",
            "in_lambert_2018",
            "lambert_2018_dbd",
            "jaspar_tf_class",
            "jaspar_tf_family",
        ]
    ]
    simple.to_csv(OUT_SIMPLE, sep="\t", index=False)
    print(f"Wrote {OUT_SIMPLE} ({len(simple)} rows)")

    print("\n--- coverage ---")
    print(f"Targets total: {len(df)}")
    print(f"With Ensembl ID: {df['ensembl_gene_id'].notna().sum()}")
    print(f"  via GTF direct: {(df['ensembl_gene_id_source']=='gtf_direct').sum()}")
    print(f"  via HGNC alias: {(df['ensembl_gene_id_source']=='hgnc_alias').sum()}")
    print(f"  via harmonized guides: {(df['ensembl_gene_id_source']=='harmonized_guides').sum()}")
    print(f"  unresolved: {df['ensembl_gene_id'].isna().sum()}")
    if df["ensembl_gene_id"].isna().any():
        print(f"  unresolved symbols: {df.loc[df['ensembl_gene_id'].isna(),'gene_symbol'].tolist()}")
    print(f"In Lambert 2018: {df['in_lambert_2018'].sum()}")
    print(f"In JASPAR human: {df['in_jaspar_core'].sum()}")
    print(f"\nTotal promoter rows accounted for: {df['n_promoters_in_target_library'].sum()}")


if __name__ == "__main__":
    main()
