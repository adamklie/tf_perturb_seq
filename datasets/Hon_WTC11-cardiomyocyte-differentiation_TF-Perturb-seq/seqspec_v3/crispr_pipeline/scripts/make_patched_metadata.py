"""Rewrite sample_metadata_gcp_<DATE>.csv → *_patched.csv pointing at patch/."""
import csv, sys, re
from pathlib import Path

DATASET = "Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq"
DATE = "2026_04_15"
BUCKET = "gs://igvf-pertub-seq-pipeline-data"
PATCH_PREFIX = f"{BUCKET}/{DATASET}/{DATE}/patch/"

BASE = Path(__file__).parent
src = BASE / f"sample_metadata_gcp_{DATE}.csv"
dst = BASE / f"sample_metadata_gcp_{DATE}_patched.csv"

BARCODE_ONLIST = PATCH_PREFIX + "IGVFFI9487JPEN.tsv"
GUIDE_DESIGN   = PATCH_PREFIX + "IGVFFI8270UPKB.tsv"     # Hon-style TSV
HASHTAG_MAP    = PATCH_PREFIX + "IGVFFI3028BJFI_with_header.tsv.gz"

def seqspec_patch(url):
    # gs://...<DATE>/<ymd>/<uuid>/IGVFFIxxxx.yaml.gz  → patch/IGVFFIxxxx.yaml
    m = re.search(r"/([^/]+)\.yaml\.gz$", url or "")
    return PATCH_PREFIX + f"{m.group(1)}.yaml" if m else url

rows = list(csv.DictReader(open(src)))
for r in rows:
    r["seqspec"] = seqspec_patch(r["seqspec"])
    if r["barcode_onlist"]:
        r["barcode_onlist"] = BARCODE_ONLIST
    if r["guide_design"]:
        r["guide_design"] = GUIDE_DESIGN
    if r["barcode_hashtag_map"]:
        r["barcode_hashtag_map"] = HASHTAG_MAP

with open(dst, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)

print(f"Wrote {dst} ({len(rows)} rows)")
# sanity
from collections import Counter
for col in ["barcode_onlist","guide_design","seqspec","barcode_hashtag_map"]:
    c = Counter(r[col] for r in rows if r.get(col))
    print(f"  {col}: {len(c)} unique, {sum(c.values())} populated")
