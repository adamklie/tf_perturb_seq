#!/bin/bash
#
# Step 3: Patch Hep CRISPR Pipeline inputs (data-driven from samplesheet).
#
# Inputs (uploaded by Step 2):
#   - All seqspec yaml.gz files referenced by IGVF accession in the samplesheet
#     (~388 portal-validated MS seqspecs + 47 aux seqspecs where present)
#   - barcode_onlist IGVFFI1697XAXI.tsv.gz (10x 5' HT v2 onlist)
#   - guide_design  IGVFFI8270UPKB.csv.gz (TF perturb library)
#
# Outputs (this script writes to <BUCKET>/<DATASET>/<DATE>/patch/):
#   - <acc>.yaml         (decompressed portal seqspecs)
#   - rna_seqspec.yml    (uploaded from local setup/seqspec/ as fallback)
#   - guide_seqspec.yml  (uploaded from local setup/seqspec/ as fallback)
#   - IGVFFI1697XAXI.tsv (decompressed onlist)
#   - IGVFFI8270UPKB.tsv (decompressed guide_design; .csv→.tsv per seqSpecCheck convention)
#
# Also rewrites the samplesheet → sample_metadata_gcp_<DATE>_patched.csv with all paths
# pointing at the /patch/ siblings above.
#
# Usage:
#   DRY_RUN=true bash 3_patch_gcp_files.sh
#   bash 3_patch_gcp_files.sh
#
set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR=/carter/users/aklie/projects/tf_perturb_seq
DATASET=Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq
DATE=2026_06_03   # CHANGE_ME — must match Step 2's GCS_PREFIX date
BUCKET=gs://igvf-pertub-seq-pipeline-data
PATCH_DIR=${BUCKET}/${DATASET}/${DATE}/patch
INPUT_SHEET=${BASE_DIR}/datasets/${DATASET}/setup/samplesheets/sample_metadata_gcp_${DATE}.csv
OUTPUT_SHEET=${BASE_DIR}/datasets/${DATASET}/setup/samplesheets/sample_metadata_gcp_${DATE}_patched.csv

# Local seqspec fallbacks (uploaded as-is, no decompression needed)
RNA_SEQSPEC_LOCAL=${BASE_DIR}/datasets/${DATASET}/setup/seqspec/rna_seqspec.yml
GUIDE_SEQSPEC_LOCAL=${BASE_DIR}/datasets/${DATASET}/setup/seqspec/guide_seqspec.yml
RNA_SEQSPEC_REL=datasets/${DATASET}/setup/seqspec/rna_seqspec.yml
GUIDE_SEQSPEC_REL=datasets/${DATASET}/setup/seqspec/guide_seqspec.yml

# =============================================================================
# PRECHECKS
# =============================================================================
echo "=== Hep Step 3 patch (data-driven) ==="
echo "  Input sheet:  ${INPUT_SHEET}"
echo "  Output sheet: ${OUTPUT_SHEET}"
echo "  Patch dir:    ${PATCH_DIR}"
echo "  Local rna:    ${RNA_SEQSPEC_LOCAL}"
echo "  Local guide:  ${GUIDE_SEQSPEC_LOCAL}"
echo

if [ ! -f "${INPUT_SHEET}" ]; then
    echo "ERROR: Step 2 output samplesheet not found: ${INPUT_SHEET}"
    echo "  Run Step 2 first."
    exit 1
fi
for f in "${RNA_SEQSPEC_LOCAL}" "${GUIDE_SEQSPEC_LOCAL}"; do
    [ -f "${f}" ] || { echo "ERROR: Local seqspec missing: ${f}"; exit 1; }
done

# =============================================================================
# DERIVE WORK LIST FROM SAMPLESHEET
# =============================================================================
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Unique GCS *.yaml.gz seqspec sources (post-Step-2 paths)
awk -F',' 'NR>1 && $7~/\.ya?ml\.gz$/ {print $7}' "${INPUT_SHEET}" | sort -u > ${TMPDIR}/seqspec_gcs.txt
# Unique GCS *.tsv.gz / *.csv.gz sources (onlist + guide)
awk -F',' 'NR>1 && ($8~/\.tsv\.gz$|\.csv\.gz$/) {print $8}' "${INPUT_SHEET}" | sort -u > ${TMPDIR}/onlist_gcs.txt
awk -F',' 'NR>1 && ($9~/\.tsv\.gz$|\.csv\.gz$/) {print $9}' "${INPUT_SHEET}" | sort -u > ${TMPDIR}/guide_gcs.txt

N_SEQSPEC=$(wc -l < ${TMPDIR}/seqspec_gcs.txt)
N_ONLIST=$(wc -l < ${TMPDIR}/onlist_gcs.txt)
N_GUIDE=$(wc -l < ${TMPDIR}/guide_gcs.txt)

echo "=== Work plan ==="
echo "  Portal seqspec yaml.gz to decompress: ${N_SEQSPEC}"
echo "  Onlist tsv.gz to decompress:          ${N_ONLIST}"
echo "  Guide csv.gz to decompress:           ${N_GUIDE}"
echo "  Local yaml uploads:                   2 (rna + guide fallbacks)"
echo

if [[ "${DRY_RUN:-false}" == "true" ]]; then
    echo "=== DRY RUN — would execute the following ==="
    echo "  Decompress ${N_SEQSPEC} portal seqspecs:"
    head -3 ${TMPDIR}/seqspec_gcs.txt | sed 's|^|    |'
    [ "$N_SEQSPEC" -gt 3 ] && echo "    ... and $((N_SEQSPEC-3)) more"
    echo "  Decompress ${N_ONLIST} onlist files:"
    head ${TMPDIR}/onlist_gcs.txt | sed 's|^|    |'
    echo "  Decompress ${N_GUIDE} guide files:"
    head ${TMPDIR}/guide_gcs.txt | sed 's|^|    |'
    echo "  Upload local yamls:"
    echo "    ${RNA_SEQSPEC_LOCAL} -> ${PATCH_DIR}/rna_seqspec.yml"
    echo "    ${GUIDE_SEQSPEC_LOCAL} -> ${PATCH_DIR}/guide_seqspec.yml"
    echo "  Rewrite samplesheet -> ${OUTPUT_SHEET}"
    exit 0
fi

# =============================================================================
# PRECHECK: Step 2 sources present?
# =============================================================================
echo "=== Verifying Step 2 sources exist on GCS (spot check) ==="
for src in $(head -1 ${TMPDIR}/seqspec_gcs.txt) $(head -1 ${TMPDIR}/onlist_gcs.txt) $(head -1 ${TMPDIR}/guide_gcs.txt); do
    if ! gsutil ls "${src}" >/dev/null 2>&1; then
        echo "ERROR: Step 2 source not found: ${src}"
        echo "  Hint: wait for Step 2 transfer to complete. Check status:"
        echo "  gcloud transfer operations list --filter='operationName:transferJobs-igvf-upload-*'"
        exit 1
    fi
done
echo "  Spot check OK."
echo

# =============================================================================
# DECOMPRESS PORTAL SEQSPECS (parallel)
# =============================================================================
decompress_one() {
    local src="$1"
    local dst_base
    dst_base=$(basename "${src}" .gz)
    gsutil -q cat "${src}" | gunzip | gsutil -q cp - "${PATCH_DIR}/${dst_base}"
}
export -f decompress_one
export PATCH_DIR

echo "=== Decompressing ${N_SEQSPEC} portal seqspecs (parallel x 8) ==="
cat ${TMPDIR}/seqspec_gcs.txt | xargs -P 8 -I {} bash -c 'decompress_one "$@"' _ {} 2>&1 | tail -5 || true
echo "  done."

# Onlist (.tsv.gz → .tsv)
echo "=== Decompressing onlist ==="
while read src; do
    dst_base=$(basename "${src}" .gz)
    echo "  ${src} -> ${PATCH_DIR}/${dst_base}"
    gsutil cat "${src}" | gunzip | gsutil cp - "${PATCH_DIR}/${dst_base}"
done < ${TMPDIR}/onlist_gcs.txt

# Guide (.csv.gz → .tsv per pipeline convention)
echo "=== Decompressing guide_design ==="
while read src; do
    dst_base=$(basename "${src}" .gz)
    dst_base=${dst_base%.csv}.tsv
    echo "  ${src} -> ${PATCH_DIR}/${dst_base}"
    gsutil cat "${src}" | gunzip | gsutil cp - "${PATCH_DIR}/${dst_base}"
done < ${TMPDIR}/guide_gcs.txt

# Local yaml fallbacks
echo "=== Uploading local seqspec yamls ==="
gsutil cp "${RNA_SEQSPEC_LOCAL}"   "${PATCH_DIR}/rna_seqspec.yml"
gsutil cp "${GUIDE_SEQSPEC_LOCAL}" "${PATCH_DIR}/guide_seqspec.yml"

# =============================================================================
# REWRITE SAMPLESHEET
# =============================================================================
echo "=== Writing patched samplesheet ${OUTPUT_SHEET} ==="
python3 << PYEOF
import csv, re, os
src = "${INPUT_SHEET}"
dst = "${OUTPUT_SHEET}"
patch = "${PATCH_DIR}"
rna_rel = "${RNA_SEQSPEC_REL}"
guide_rel = "${GUIDE_SEQSPEC_REL}"

def patch_seqspec(v):
    # Portal seqspec .yaml.gz on GCS  →  patch/<acc>.yaml
    m = re.match(r'^gs://.*/(IGVFFI\w+)\.(ya?ml)\.gz$', v)
    if m: return f"{patch}/{m.group(1)}.{m.group(2)}"
    # Local fallback rna_seqspec.yml  →  patch/rna_seqspec.yml
    if v.endswith('rna_seqspec.yml'):   return f"{patch}/rna_seqspec.yml"
    if v.endswith('guide_seqspec.yml'): return f"{patch}/guide_seqspec.yml"
    return v

def patch_tab(v):
    m = re.match(r'^gs://.*/(IGVFFI\w+)\.(tsv|csv)\.gz$', v)
    if m: return f"{patch}/{m.group(1)}.tsv"
    return v

with open(src, newline='') as fin, open(dst, 'w', newline='') as fout:
    r = csv.reader(fin); w = csv.writer(fout)
    hdr = next(r); w.writerow(hdr)
    si = hdr.index('seqspec')
    oi = hdr.index('barcode_onlist')
    gi = hdr.index('guide_design')
    for row in r:
        row[si] = patch_seqspec(row[si])
        row[oi] = patch_tab(row[oi])
        row[gi] = patch_tab(row[gi])
        w.writerow(row)
print(f"  patched samplesheet written: {dst}")
PYEOF

echo
echo "=== Done. Patch dir listing: ==="
gsutil ls "${PATCH_DIR}/" | head -10
echo "  ... ($(gsutil ls ${PATCH_DIR}/ | wc -l) files total)"
echo
echo "Step 4 ready to launch against: ${OUTPUT_SHEET}"
