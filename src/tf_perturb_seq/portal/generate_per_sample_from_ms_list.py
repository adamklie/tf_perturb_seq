"""Generate the IGVF pipeline samplesheet from an explicit list of measurement sets.

For datasets without an analysis set (e.g. Gersbach Hep as of 2026-05), we
build the samplesheet by enumerating the measurement sets directly and
passing the guide library file accession explicitly.

Same per-row schema as src/tf_perturb_seq/portal/generate_per_sample.py:
    R1_path, R2_path, file_modality, measurement_sets, sequencing_run,
    lane, seqspec, barcode_onlist, guide_design, barcode_hashtag_map

Limitations:
  - No cell-hashing / auxiliary set handling (auxiliary sets are normally
    pulled via the analysis set's input_file_sets). If the dataset uses
    hashing, extend this script to accept --auxiliary_sets.
  - No multi-guide-library support (we accept a single --guide_design).

Usage:
    python generate_per_sample_from_ms_list.py \\
        --measurement_sets hep_ms_list.txt \\
        --guide_design IGVFFI8270UPKB \\
        --output sample_metadata.csv
"""

import argparse
import csv
import json
import os
import sys
from pathlib import Path

import requests
from requests.auth import HTTPBasicAuth

PORTAL = "https://api.data.igvf.org"

MODALITY_MAP = {
    "scRNA sequencing": "scRNA",
    "gRNA sequencing": "gRNA",
    "cell hashing barcode sequencing": "hash",
}


def get_auth(keypair_path=None):
    if keypair_path:
        with open(keypair_path) as f:
            kp = json.load(f)
        return HTTPBasicAuth(kp["key"], kp["secret"])
    k = os.getenv("IGVF_API_KEY")
    s = os.getenv("IGVF_SECRET_KEY")
    if k and s:
        return HTTPBasicAuth(k, s)
    raise RuntimeError("No IGVF credentials (set IGVF_API_KEY/SECRET_KEY or pass --keypair)")


def fetch(path, auth):
    if not path.startswith("/"):
        path = "/" + path
    r = requests.get(f"{PORTAL}{path}", auth=auth, headers={"Accept": "application/json"}, timeout=30)
    r.raise_for_status()
    return r.json()


def detect_file_modality(f):
    """Pick the pipeline modality for one sequence-file object.

    Tries (in order):
      1. assay_titles on the file (rarely set per-file in practice)
      2. file aliases (Gersbach Hep encodes 'GEX' / 'CRISPR' here)
      3. submitted_file_name path (contains '/GEX_' or '/CRISPR_')
    Returns the pipeline modality ("scRNA" / "gRNA" / "hash") or None.
    """
    for at in f.get("assay_titles", []) or []:
        atl = at.lower()
        if "rna" in atl or "scrna" in atl:
            return "scRNA"
        if "grna" in atl or "guide" in atl or "sgrna" in atl:
            return "gRNA"
        if "hash" in atl:
            return "hash"
    for al in f.get("aliases", []) or []:
        all_lower = al.lower()
        if "-gex-" in all_lower or "gex_" in all_lower or "scrna" in all_lower:
            return "scRNA"
        if "-crispr-" in all_lower or "crispr_" in all_lower or "-grna-" in all_lower or "-sgrna-" in all_lower:
            return "gRNA"
        if "hash" in all_lower:
            return "hash"
    sn = (f.get("submitted_file_name") or "").lower()
    if "/gex_" in sn or "_gex_" in sn or "gex-" in sn:
        return "scRNA"
    if "/crispr_" in sn or "_crispr_" in sn or "crispr-" in sn or "/grna_" in sn or "/sgrna_" in sn:
        return "gRNA"
    if "/hash_" in sn or "_hash_" in sn:
        return "hash"
    return None


def find_paired_auxiliary_set(ms_id, auth):
    """Find the gRNA auxiliary set that references this measurement set, if any."""
    # Direct: aux sets that reference this ms via measurement_sets
    url = (
        f"{PORTAL}/search/?type=AuxiliarySet"
        f"&measurement_sets.accession={ms_id}"
        f"&file_set_type=gRNA+sequencing"
        f"&format=json&limit=10"
        f"&field=accession&field=aliases"
    )
    r = requests.get(url, auth=auth, headers={"Accept": "application/json"}, timeout=30)
    r.raise_for_status()
    items = r.json().get("@graph", [])
    if not items:
        return None
    if len(items) > 1:
        # Multiple — pick the one with matching iHep alias if obvious; else first
        print(f"  [{ms_id}] WARNING: {len(items)} paired aux sets — taking first ({items[0]['accession']})")
    return items[0]["accession"]


def process_fileset_files(fset_id, fset_kind, modality_hint, files, *,
                          measurement_set_acc, guide_design, auth, fallback_seqspecs,
                          barcode_onlist, strand_specificity):
    """Yield rows for all R1/R2 pairs in a fileset (measurement set or aux set).

    modality_hint forces the modality if the fileset is a single-modality aux set
    (e.g., gRNA sequencing). For measurement sets containing mixed modalities, pass
    None and we auto-detect per-file.
    """
    sequence_file_index = {}
    for f_ref in files:
        f = fetch(f_ref + "/@@object?format=json", auth)
        if not f["@id"].startswith("/sequence-files/"):
            continue
        if f.get("status") in ["deleted", "revoked"]:
            continue
        if f.get("content_type") != "reads":
            continue
        if f.get("illumina_read_type") not in ("R1", "R2"):
            continue

        if modality_hint:
            file_modality = modality_hint
        else:
            file_modality = detect_file_modality(f)
            if file_modality is None:
                raise ValueError(
                    f"{fset_id}: could not detect modality for {f['@id']} (aliases={f.get('aliases')}, "
                    f"submitted={f.get('submitted_file_name')})"
                )

        key = (
            f.get("sequencing_run"),
            f.get("lane"),
            f.get("flowcell_id"),
            f.get("index"),
            file_modality,
        )
        sequence_file_index.setdefault(key, {})[f.get("illumina_read_type")] = f

    for key, reads in sequence_file_index.items():
        r1 = reads.get("R1")
        r2 = reads.get("R2")
        if not r1 or not r2:
            continue
        # Seqspec lookup with per-modality fallback
        seqspec_path = ""
        for ss in r1.get("seqspecs", []) or []:
            ss_obj = fetch(ss + "/@@object?format=json", auth)
            if (
                ss_obj.get("upload_status") == "validated"
                and ss_obj.get("status") in ["in progress", "preview", "released"]
            ):
                seqspec_path = ss.split("/")[-2]
                break
        if not seqspec_path:
            seqspec_path = fallback_seqspecs.get(key[4])
            if not seqspec_path:
                raise ValueError(
                    f"{fset_id}: missing seqspec on R1 {r1['@id']} for modality {key[4]} "
                    f"and no fallback provided"
                )

        yield {
            "R1_path": r1["@id"].split("/")[-2],
            "R2_path": r2["@id"].split("/")[-2],
            "file_modality": key[4],
            "measurement_sets": measurement_set_acc,
            "sequencing_run": key[0],
            "lane": key[1],
            "seqspec": seqspec_path,
            "barcode_onlist": barcode_onlist[0].split("/")[-2] if barcode_onlist else "",
            "guide_design": guide_design,
            "barcode_hashtag_map": "",
        }


def rows_for_measurement_set(ms_id, guide_design, auth, fallback_seqspecs, *,
                              barcode_onlist_fallback=None,
                              strand_specificity_fallback=None,
                              aux_map=None):
    """Yield per-(sequencing_run, lane, modality) rows for one measurement set
    AND its paired CRISPR auxiliary set (Hep convention)."""
    ms = fetch(f"/measurement-sets/{ms_id}/@@object?format=json", auth)

    barcode_onlist = ms.get("onlist_files") or []
    onlist_method = ms.get("onlist_method") or ""
    if onlist_method and onlist_method not in ["no combination", "multi"]:
        raise ValueError(f"{ms_id}: unsupported onlist_method={onlist_method}")
    strand_specificity = ms.get("strand_specificity") or ""
    if not strand_specificity:
        if not strand_specificity_fallback:
            raise ValueError(
                f"{ms_id}: missing strand_specificity on portal and no --strand_specificity_fallback"
            )
        strand_specificity = strand_specificity_fallback
        print(f"  [{ms_id}] strand_specificity from fallback: {strand_specificity}")
    if not barcode_onlist:
        if not barcode_onlist_fallback:
            raise ValueError(
                f"{ms_id}: missing onlist_files on portal and no --barcode_onlist_fallback"
            )
        # Format the fallback as a portal-style list so downstream code is uniform
        barcode_onlist = [f"/tabular-files/{barcode_onlist_fallback}/" if not barcode_onlist_fallback.startswith("/") else barcode_onlist_fallback]
        print(f"  [{ms_id}] barcode_onlist from fallback: {barcode_onlist_fallback}")

    # Pass 1: GEX measurement set's files (scRNA modality)
    yield from process_fileset_files(
        ms_id, "measurement_set", None, ms.get("files", []),
        measurement_set_acc=ms_id, guide_design=guide_design, auth=auth,
        fallback_seqspecs=fallback_seqspecs,
        barcode_onlist=barcode_onlist, strand_specificity=strand_specificity,
    )

    # Pass 2: paired CRISPR auxiliary set's files (gRNA modality)
    aux_id = (aux_map or {}).get(ms_id) if aux_map else None
    if not aux_id:
        aux_id = find_paired_auxiliary_set(ms_id, auth)
    if not aux_id:
        print(f"  [{ms_id}] WARNING: no paired CRISPR auxiliary set found")
        return
    aux = fetch(f"/auxiliary-sets/{aux_id}/@@object?format=json", auth)
    print(f"  [{ms_id}] paired aux set: {aux_id} ({aux.get('file_set_type')})")
    yield from process_fileset_files(
        aux_id, "auxiliary_set", "gRNA", aux.get("files", []),
        measurement_set_acc=ms_id, guide_design=guide_design, auth=auth,
        fallback_seqspecs=fallback_seqspecs,
        barcode_onlist=barcode_onlist, strand_specificity=strand_specificity,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--measurement_sets",
        required=True,
        help="Path to a file with one IGVF measurement set accession per line, "
        "OR a comma-separated list of accessions.",
    )
    ap.add_argument(
        "--guide_design",
        required=True,
        help="IGVF guide RNA sequences file accession (e.g. IGVFFI8270UPKB).",
    )
    ap.add_argument("--output", required=True, help="Output samplesheet CSV path")
    ap.add_argument("--keypair", help="Optional IGVF keypair JSON")
    ap.add_argument("--rna_seqspec", help="Fallback seqspec path for scRNA modality")
    ap.add_argument("--sgrna_seqspec", help="Fallback seqspec path for gRNA modality")
    ap.add_argument("--hash_seqspec", help="Fallback seqspec path for hashing modality")
    ap.add_argument(
        "--barcode_onlist",
        help="Fallback barcode onlist accession (e.g. IGVFFI4695IKAL) if portal field is null.",
    )
    ap.add_argument(
        "--strand_specificity",
        help="Fallback strand_specificity if missing on portal (e.g. '5 prime to 3 prime').",
    )
    ap.add_argument(
        "--aux_map",
        help="Optional TSV with explicit gex_measurement_set <TAB> crispr_auxiliary_set rows. "
        "If not provided, each GEX MS is searched against the portal for its paired aux set.",
    )
    args = ap.parse_args()

    # Resolve measurement set list — accept lines like "IGVFDS... # comment" or "IGVFDS...\t# comment"
    if Path(args.measurement_sets).is_file():
        ms_list = []
        for line in Path(args.measurement_sets).read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # take the first token before whitespace or '#'
            acc = line.split("#", 1)[0].split()[0].strip()
            if acc:
                ms_list.append(acc)
    else:
        ms_list = [x.strip() for x in args.measurement_sets.split(",") if x.strip()]
    if not ms_list:
        sys.exit("No measurement sets resolved.")
    print(f"Processing {len(ms_list)} measurement sets...")

    auth = get_auth(args.keypair)
    fallback_seqspecs = {
        "scRNA": args.rna_seqspec,
        "gRNA": args.sgrna_seqspec,
        "hash": args.hash_seqspec,
    }

    # Load aux_map if provided
    aux_map = {}
    if args.aux_map and Path(args.aux_map).is_file():
        with open(args.aux_map) as f:
            for line in f:
                if line.startswith("#") or line.startswith("gex_"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 2 and parts[1] and parts[1] != "MISSING":
                    aux_map[parts[0]] = parts[1]
        print(f"Loaded {len(aux_map)} GEX→CRISPR aux mappings from {args.aux_map}")

    all_rows = []
    for i, ms_id in enumerate(ms_list, 1):
        print(f"[{i}/{len(ms_list)}] {ms_id}")
        try:
            for row in rows_for_measurement_set(
                ms_id, args.guide_design, auth, fallback_seqspecs,
                barcode_onlist_fallback=args.barcode_onlist,
                strand_specificity_fallback=args.strand_specificity,
                aux_map=aux_map,
            ):
                all_rows.append(row)
        except Exception as e:
            print(f"  ERROR on {ms_id}: {e}")
            raise

    if not all_rows:
        sys.exit("No rows produced — every measurement set yielded zero usable read pairs.")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        fieldnames = list(all_rows[0].keys())
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(all_rows)
    print(f"\nWrote {len(all_rows)} rows to {args.output}")
    # Breakdown by modality
    from collections import Counter
    by_mod = Counter(r["file_modality"] for r in all_rows)
    print(f"  Rows by modality: {dict(by_mod)}")
    by_ms = Counter(r["measurement_sets"] for r in all_rows)
    print(f"  Measurement sets covered: {len(by_ms)}")


if __name__ == "__main__":
    main()
