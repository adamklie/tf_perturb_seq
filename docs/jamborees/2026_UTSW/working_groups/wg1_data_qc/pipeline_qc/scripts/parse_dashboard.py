"""Parse each dataset's pipeline_dashboard/dashboard.html and emit three TSVs
under results/upstream_mapping_and_filtering/:

  per_lane_mapping_summary.tsv
    dataset, short_name, modality, measurement_set,
    total_reads, paired_reads_mapped, alignment_pct, detected_barcodes,
    pct_reads_in_onlist   (kallisto-only; NaN where mapper differs)

  dataset_summary_metrics.tsv
    dataset, short_name, n_cells, umi_median, mito_median,
    total_scrna_reads, total_guide_reads, total_hto_reads,
    scrna_pct_aligned, guide_pct_aligned, hto_pct_aligned,
    n_measurement_sets, reads_per_cell

  filtering_funnel.tsv
    dataset, short_name, stage_index, stage, cells, removed

Source paths come from manifests/manifest.tsv (dashboard_local column).
All numeric columns are raw (no "M" / "K" / "%"), suitable for direct plotting.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MANIFEST = ROOT / "manifests" / "manifest.tsv"
RESULTS = ROOT / "results" / "upstream_mapping_and_filtering"

MODALITY_KEYS = {
    "Mapping scRNA": "scRNA",
    "Mapping Guide": "Guide",
    "Mapping Hashing": "Hashing",
}

SUFFIX_FACTOR = {"K": 1e3, "M": 1e6, "G": 1e9, "B": 1e9, "T": 1e12}


def parse_human_number(s: str | None) -> float | None:
    """Parse '641M' / '91.5%' / '19,878,833' / '4900' -> float (no unit). NaN-safe."""
    if s is None:
        return None
    s = s.strip()
    if not s or s.lower() in {"n/a", "na", "none"}:
        return None
    s = s.replace(",", "").rstrip("%")
    if not s:
        return None
    if s[-1] in SUFFIX_FACTOR:
        return float(s[:-1]) * SUFFIX_FACTOR[s[-1]]
    return float(s)


def card_spans(html: str) -> list[tuple[int, int, str]]:
    """(start, end, card_title) for each h3 card, in document order."""
    headers = list(re.finditer(r'<h3 id="card-title-\d+">([^<]+)</h3>', html))
    spans = []
    for i, m in enumerate(headers):
        start = m.start()
        end = headers[i + 1].start() if i + 1 < len(headers) else len(html)
        spans.append((start, end, m.group(1).strip()))
    return spans


def parse_mapping_card(card_html: str) -> dict | None:
    """Top-line value-display paragraph has fields common across all dashboards.
    For kallisto-mapped datasets we also pull pct_reads_in_onlist from the
    paginated table; for STARsolo/alevin-fry datasets that column will be NaN.
    """
    ms = re.search(r"<p>([A-Za-z0-9_-]+)</p>", card_html)
    if not ms:
        return None

    row = {"measurement_set": ms.group(1)}
    for key, label in [
        ("total_reads", r"Total Reads[^<:]*:"),
        ("paired_reads_mapped", r"Paired Reads Mapped:"),
        ("alignment_pct", r"Alignment Percentage:"),
        ("detected_barcodes", r"Total Detected [^<]*Barcodes[^<]*:"),
    ]:
        m = re.search(rf'{label}\s*<span[^>]*>(?:<span[^>]*>)?([\d.,KMBG]+)%?', card_html)
        row[key] = parse_human_number(m.group(1)) if m else None

    m = re.search(
        r"<tr>\s*<td>percentageReadsOnOnlist</td>\s*<td>([^<]+)</td>",
        card_html,
    )
    row["pct_reads_in_onlist"] = parse_human_number(m.group(1)) if m else None
    return row


def parse_mapping(html: str) -> pd.DataFrame:
    cols = [
        "modality",
        "measurement_set",
        "total_reads",
        "paired_reads_mapped",
        "alignment_pct",
        "detected_barcodes",
        "pct_reads_in_onlist",
    ]
    rows = []
    for start, end, title in card_spans(html):
        modality = MODALITY_KEYS.get(title)
        if modality is None:
            continue
        parsed = parse_mapping_card(html[start:end])
        if parsed is None:
            continue
        parsed["modality"] = modality
        rows.append(parsed)
    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(rows)[cols]


def parse_funnel(html: str) -> pd.DataFrame:
    fc = re.search(r'<span class="flowchart">(.+?)</span>\s*<span class="flow-footnote"', html, re.S)
    if not fc:
        return pd.DataFrame(columns=["stage_index", "stage", "cells", "removed"])
    flow = fc.group(1)
    steps = re.findall(
        r'<span class="flow-step">(.+?)(?=<span class="flow-arrow|<span class="flow-footnote|$)',
        flow,
        re.S,
    )
    rows = []
    for i, step in enumerate(steps):
        # Bound stage-title match at </span> so titles with literal '<' survive
        # (e.g. "Mito filter (pct_counts_mt < 20%)").
        title = re.search(r'flow-title">(.+?)</span>', step, re.S)
        cells = re.search(r'flow-count">Cells:\s*([\d,]+)', step)
        removed = re.search(r'flow-removed">Removed:\s*([\d,]+)', step)
        if not (title and cells):  # skip header step with no count
            continue
        rows.append(
            {
                "stage_index": i,
                "stage": title.group(1).strip(),
                "cells": int(cells.group(1).replace(",", "")),
                "removed": int(removed.group(1).replace(",", "")) if removed else None,
            }
        )
    return pd.DataFrame(rows)


def parse_topline_summary(html: str) -> dict[str, float | None]:
    out: dict[str, float | None] = {"n_cells": None, "umi_median": None, "mito_median": None}
    m = re.search(
        r"Cells:\s*<span[^>]*>([\d.,KMB]+)</span>\s*<br>\s*Median UMIs:\s*<span[^>]*>([\d.,]+)</span>",
        html,
    )
    if m:
        out["n_cells"] = parse_human_number(m.group(1))
        out["umi_median"] = parse_human_number(m.group(2))
    m = re.search(r"Median mito%:\s*<span[^>]*><span[^>]*>([\d.]+)</span>", html)
    if m:
        out["mito_median"] = parse_human_number(m.group(1))
    return out


def build_dataset_summary(short: str, dataset: str, mapping: pd.DataFrame, topline: dict) -> dict:
    row = {"dataset": dataset, "short_name": short, **topline}

    def _agg(modality: str, prefix: str) -> None:
        sub = mapping[mapping["modality"] == modality]
        if sub.empty:
            row[f"total_{prefix}_reads"] = None
            row[f"{prefix}_pct_aligned"] = None
            return
        total = sub["total_reads"].sum()
        aligned = sub["paired_reads_mapped"].sum()
        row[f"total_{prefix}_reads"] = total
        row[f"{prefix}_pct_aligned"] = 100 * aligned / total if total else None

    _agg("scRNA", "scrna")
    _agg("Guide", "guide")
    _agg("Hashing", "hto")

    n_ms = mapping[mapping["modality"] == "scRNA"]["measurement_set"].nunique()
    row["n_measurement_sets"] = int(n_ms) if n_ms else None
    row["reads_per_cell"] = (
        row["total_scrna_reads"] / row["n_cells"]
        if row.get("total_scrna_reads") and row.get("n_cells")
        else None
    )
    return row


def main() -> int:
    manifest = pd.read_csv(MANIFEST, sep="\t")
    RESULTS.mkdir(parents=True, exist_ok=True)

    all_mapping, all_funnel, all_summary = [], [], []

    for _, m in manifest.iterrows():
        dash = Path(m["dashboard_local"])
        if not dash.exists():
            print(f"[skip] {m['short_name']}: {dash} not found")
            continue
        print(f"[parse] {m['short_name']}  ({dash.stat().st_size/1e6:.1f} MB)")
        html = dash.read_text()

        mapping = parse_mapping(html)
        mapping.insert(0, "dataset", m["dataset"])
        mapping.insert(1, "short_name", m["short_name"])
        all_mapping.append(mapping)

        funnel = parse_funnel(html)
        funnel.insert(0, "dataset", m["dataset"])
        funnel.insert(1, "short_name", m["short_name"])
        all_funnel.append(funnel)

        topline = parse_topline_summary(html)
        all_summary.append(build_dataset_summary(m["short_name"], m["dataset"], mapping, topline))

        modalities_seen = ",".join(sorted(mapping["modality"].unique())) if not mapping.empty else "(none)"
        print(
            f"        mapping: {len(mapping):3d} rows ({modalities_seen})   "
            f"funnel: {len(funnel)} stages   topline: {topline}"
        )

    mapping_out = pd.concat(all_mapping, ignore_index=True)
    funnel_out = pd.concat(all_funnel, ignore_index=True)
    summary_out = pd.DataFrame(all_summary)

    print()
    for name, df in [
        ("per_lane_mapping_summary.tsv", mapping_out),
        ("filtering_funnel.tsv", funnel_out),
        ("dataset_summary_metrics.tsv", summary_out),
    ]:
        p = RESULTS / name
        df.to_csv(p, sep="\t", index=False)
        print(f"wrote {p.relative_to(ROOT)}  ({len(df)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
