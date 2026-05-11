#!/bin/bash
#
# Decompress .gz files on GCS to a patch directory
#
# The pipeline requires uncompressed files for:
# - barcode_onlist (.tsv)
# - guide_design (.tsv -- portal serves .csv.gz but content is tab-separated)
# - seqspec (.yaml)  -- 84 unique seqspecs for this dataset, stripped of i7/i5
#                       sample-index reads (see SEQSPEC NOTE below)
# - barcode_hashtag_map (handled separately; see note below)
#
# Usage: ./3_patch_gcp_files.sh
#
# SEQSPEC NOTE — Index-read stripping
# -----------------------------------
# The Hon-supplied seqspecs include an extra `!Read` block for the i7 sample
# index (10bp `primer_id: index7`, sometimes mislabeled as `truseq_read2` with
# max_len=8 in the RNA seqspecs). `parsing_guide_metadata.py` in the pipeline
# feeds every read_id into `seqspec index -t kb`, which then produces a
# chemistry string referencing 3 fastq files per lane. The pipeline only
# stages 2 (R1+R2), so kallisto fails with "Number of files does not match
# number of input files required by technology". We strip those reads here.
set -e

STRIP_PY=$(cat <<'PYEOF'
import sys, re
text = sys.stdin.read()
m = re.search(r'^sequence_spec:\s*\n', text, flags=re.MULTILINE)
if not m:
    sys.stdout.write(text); sys.exit(0)
head = text[:m.end()]; rest = text[m.end():]
lines = rest.split('\n')
blocks = []; cur = None; end_idx = len(lines)
for i, ln in enumerate(lines):
    if ln.startswith('- !Read'):
        if cur is not None: blocks.append(cur)
        cur = [ln]
    elif ln and not ln.startswith(' ') and not ln.startswith('-'):
        if cur is not None: blocks.append(cur); cur = None
        end_idx = i; break
    else:
        if cur is not None: cur.append(ln)
if cur is not None: blocks.append(cur)
epilogue = '\n'.join(lines[end_idx:])
kept = []
for b in blocks:
    body = '\n'.join(b)
    if re.search(r'primer_id:\s*index[57]\b', body): continue
    m_max = re.search(r'max_len:\s*(\d+)', body)
    if m_max and int(m_max.group(1)) < 12: continue
    kept.append(body)
out = head + '\n'.join(kept)
if not out.endswith('\n'): out += '\n'
out += epilogue
sys.stdout.write(out)
PYEOF
)

BUCKET="gs://igvf-pertub-seq-pipeline-data"
DATASET="Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq"
DATE="2026_04_15"
PATCH_DIR="${BUCKET}/${DATASET}/${DATE}/patch"

# =============================================================================
# BARCODE ONLIST & GUIDE DESIGN
# =============================================================================

BARCODE_ONLIST_SRC="${BUCKET}/${DATASET}/${DATE}/2024/12/05/996cf65e-5b5c-494d-b8d3-f1f1e8b0181c/IGVFFI9487JPEN.tsv.gz"
GUIDE_DESIGN_SRC="${BUCKET}/${DATASET}/${DATE}/2024/08/21/2dca8d2e-e59f-425d-b0d4-835ec3fb61f9/IGVFFI8270UPKB.csv.gz"

echo "=== Decompressing barcode_onlist (IGVFFI9487JPEN) ==="
gsutil cat "${BARCODE_ONLIST_SRC}" | gunzip | gsutil cp - "${PATCH_DIR}/IGVFFI9487JPEN.tsv"

echo "=== Decompressing guide_design (IGVFFI8270UPKB) ==="
# NOTE: IGVF stores this file with a .csv extension but the contents are tab-
#       separated. The Hon benchmark pipeline expects .tsv — write both to be safe.
gsutil cat "${GUIDE_DESIGN_SRC}" | gunzip | gsutil cp - "${PATCH_DIR}/IGVFFI8270UPKB.tsv"

# =============================================================================
# SEQSPEC YAML FILES (84 unique; paths listed verbatim from step-2 metadata)
# =============================================================================

SEQSPEC_RELPATHS=(
    "2025/07/21/1532d24f-ecca-4349-a7d7-0080b463b373/IGVFFI9505FIRT.yaml.gz"
    "2025/07/21/1b92e47d-ec27-4723-9c64-52c9bd61fa04/IGVFFI5605QWGQ.yaml.gz"
    "2025/07/21/1fecf456-ee69-4f83-8203-019bcb2c3fa2/IGVFFI6281ZPWU.yaml.gz"
    "2025/07/21/314e935b-f559-45fc-95e8-a8b7961b5e70/IGVFFI7351DFXB.yaml.gz"
    "2025/07/21/3abe1553-6910-49df-a22d-e33c4a028525/IGVFFI9226TQDW.yaml.gz"
    "2025/07/21/400b0efd-2ed5-4030-a4f3-0f7dfb7e5b0f/IGVFFI1313DWZJ.yaml.gz"
    "2025/07/21/4ab35945-52d2-4a28-9a46-3407d5eac8ed/IGVFFI6105GRYO.yaml.gz"
    "2025/07/21/4e37ace5-8d1c-47fd-9a62-b719d439a3e7/IGVFFI4474EBAY.yaml.gz"
    "2025/07/21/53f46a50-bcba-433e-86ce-cefa40c1b3d0/IGVFFI9051ZCVI.yaml.gz"
    "2025/07/21/633a01ca-3dc0-4ae0-8087-99cf72ff0686/IGVFFI9369EFAR.yaml.gz"
    "2025/07/21/66f7fd38-5f86-4708-b36e-c2fe46689131/IGVFFI3145AMWL.yaml.gz"
    "2025/07/21/6d4b667d-af2d-431a-b3ee-e66b86427f61/IGVFFI7979QFLR.yaml.gz"
    "2025/07/21/81f7312d-2067-482e-b326-508bfd3b1703/IGVFFI6264KFSA.yaml.gz"
    "2025/07/21/8b1f3abc-38c6-4211-9095-ec68f7e653ed/IGVFFI4956FUMQ.yaml.gz"
    "2025/07/21/8cbb5176-1a10-41d9-859b-51c5ad946326/IGVFFI9346BRXZ.yaml.gz"
    "2025/07/21/a221c8b2-74e3-46f6-a26e-896006842a9a/IGVFFI3081TOBQ.yaml.gz"
    "2025/07/21/a8a83740-2fee-4449-b84e-cc845c7bcd84/IGVFFI9735GJHE.yaml.gz"
    "2025/07/21/b3d57d66-3de8-4f69-8854-c8275a808b65/IGVFFI1143RPTQ.yaml.gz"
    "2025/07/21/bac6f919-60a9-4bd1-afcd-861af3208709/IGVFFI7117BPUP.yaml.gz"
    "2025/07/21/c369c4d5-06ba-4695-928f-d71307e9522a/IGVFFI4558NWWG.yaml.gz"
    "2025/07/21/c57c7ea9-415a-42e5-9012-109e2cd84f1b/IGVFFI5565FKGG.yaml.gz"
    "2025/07/21/cc2cd454-d964-4e00-bc85-6a4070221236/IGVFFI1941QGQW.yaml.gz"
    "2025/07/21/d2411d89-052e-419e-b1a8-499db3a0cbda/IGVFFI0699UWYY.yaml.gz"
    "2025/07/21/dc1a827f-f276-4c45-a3c5-407c3dea1acc/IGVFFI9894IGJG.yaml.gz"
    "2025/07/21/e72ffb9f-21b9-482c-ab98-c257e0311a92/IGVFFI3803RMXF.yaml.gz"
    "2025/07/21/f39e788e-2b83-4466-8b08-2d0e8ca36ff6/IGVFFI2360MPMH.yaml.gz"
    "2025/07/21/f5dc09b4-03d2-439a-802f-1803c220a5b6/IGVFFI9154XFLC.yaml.gz"
    "2025/07/21/ffc6975c-4aac-485f-8e4e-bcf27b7de389/IGVFFI3585YHTP.yaml.gz"
    "2025/10/13/14c20ebd-bdc5-4b11-a572-991356088230/IGVFFI5781ASKG.yaml.gz"
    "2025/10/13/1521d4f0-5cd2-42b7-a405-eb0f62ea32e9/IGVFFI3387BACM.yaml.gz"
    "2025/10/13/18452234-ec4a-489f-a5b8-5d2bd8dfb533/IGVFFI2489GBMS.yaml.gz"
    "2025/10/13/1a500026-d608-4402-9985-1fa53812b4b7/IGVFFI1351YFFX.yaml.gz"
    "2025/10/13/221aed7f-ed7f-43ae-90e8-d6f0a8851dbc/IGVFFI9654TXXQ.yaml.gz"
    "2025/10/13/2bfbe16b-4fc7-47dc-8f29-09512706ddd5/IGVFFI3234CEPG.yaml.gz"
    "2025/10/13/3e4e9045-d11f-400a-9948-7e3a17f52ff4/IGVFFI8690XLVU.yaml.gz"
    "2025/10/13/3eb41290-ad82-48d8-ba28-79b11fa64d02/IGVFFI7385WYWW.yaml.gz"
    "2025/10/13/5944466b-57e8-4fb1-970b-e1ca76fd33ca/IGVFFI4245DFAH.yaml.gz"
    "2025/10/13/5c4e4389-da44-4de5-9b43-789dd88c8891/IGVFFI1285EAND.yaml.gz"
    "2025/10/13/691a7229-8874-4e16-a7ce-e0a00078d422/IGVFFI1457ZWBO.yaml.gz"
    "2025/10/13/9094a400-ea89-4395-94cc-98bb94627dcb/IGVFFI8130ZIAF.yaml.gz"
    "2025/10/13/94c430a8-9f74-4de1-b1a9-0a42cc88bc22/IGVFFI4097TYIU.yaml.gz"
    "2025/10/13/9adad163-d659-430e-81fb-c23ebb6eb06b/IGVFFI3503FEVJ.yaml.gz"
    "2025/10/13/9ae4f95d-72f7-451d-942b-3291628abc1c/IGVFFI0166PDSQ.yaml.gz"
    "2025/10/13/afb674fc-cd2c-425d-bc2c-141ff8ed73e7/IGVFFI1734XHHP.yaml.gz"
    "2025/10/13/b752fe0f-8e60-4751-8ee8-65f6ccdd9d8e/IGVFFI8589HBMB.yaml.gz"
    "2025/10/13/cf00746e-85da-4a76-bd6d-593b37eee295/IGVFFI1211IVIF.yaml.gz"
    "2025/10/13/d08979cf-0ce1-4306-acb8-84df293dccde/IGVFFI5563KUWN.yaml.gz"
    "2025/10/13/d31b1cc7-8ae9-4a40-920a-309dcd3071e5/IGVFFI9227FGVG.yaml.gz"
    "2025/10/13/d3999b5f-ae48-4629-95a9-ec612d053800/IGVFFI9000NROB.yaml.gz"
    "2025/10/13/d918e2f2-2728-4a1a-8369-b7e3c87e5065/IGVFFI4952AEFI.yaml.gz"
    "2025/10/13/e1ddc4b0-f4ab-4272-9ce5-8ffe79743edb/IGVFFI7159NAZZ.yaml.gz"
    "2025/10/13/e2bf8ea7-4fe0-4f7c-97fc-f2915ff22155/IGVFFI2224ZDCC.yaml.gz"
    "2025/10/13/e9cd81ef-c7e6-4e99-9465-14f744855ab5/IGVFFI6200CVTW.yaml.gz"
    "2025/10/13/ed30a605-2ec5-4551-b5b6-a6c5c376aca2/IGVFFI0267XCPU.yaml.gz"
    "2025/10/13/f06bdaf1-ea7b-4229-a4d6-5fe7bdd3fb3d/IGVFFI0790CZMT.yaml.gz"
    "2025/10/13/fa0445a2-efdb-4b81-b38b-e3ec45820602/IGVFFI0619ZBBU.yaml.gz"
    "2025/10/22/0f2f156e-ec9a-44c2-8750-47b28a57a1c7/IGVFFI0190XNRB.yaml.gz"
    "2025/10/22/141cb7c1-0816-48ec-bb8f-44172bbb1285/IGVFFI7422MEYL.yaml.gz"
    "2025/10/22/18b2e472-b697-4789-ada4-886c950f71da/IGVFFI6793QXXI.yaml.gz"
    "2025/10/22/25c1afac-51f6-44d5-b1b2-70636c633f57/IGVFFI0756NJJZ.yaml.gz"
    "2025/10/22/2739110b-b6b0-4eea-aeb9-d39d929f6f50/IGVFFI7399HNJF.yaml.gz"
    "2025/10/22/27d794c0-a4c4-4fd8-9132-3eb48a2b7626/IGVFFI3599FGRM.yaml.gz"
    "2025/10/22/390d8378-fcce-4889-b6bf-7565e4ed815c/IGVFFI8816BUIS.yaml.gz"
    "2025/10/22/3ab6e6cc-6187-410e-be1c-b9a40f922d24/IGVFFI7589EDXY.yaml.gz"
    "2025/10/22/3f6469f8-9648-41e5-a4d6-6f314507a8de/IGVFFI9078PBTG.yaml.gz"
    "2025/10/22/4163f4ff-e458-4413-b3f2-452b12c0a475/IGVFFI7175YGID.yaml.gz"
    "2025/10/22/48410694-db32-43bb-8404-ad41c759f229/IGVFFI6854UNOO.yaml.gz"
    "2025/10/22/49f616ac-1ff4-4dce-aeb7-c374f4e0ffbf/IGVFFI6213LVZQ.yaml.gz"
    "2025/10/22/4dff7a7c-2222-4c6c-a789-ee115e58b414/IGVFFI1802EPMA.yaml.gz"
    "2025/10/22/74917dd4-27ea-4263-acf9-3da3bfa0bcd6/IGVFFI3206IZRV.yaml.gz"
    "2025/10/22/762b91c7-74fb-4887-b8a3-3ab825732208/IGVFFI2822RAJY.yaml.gz"
    "2025/10/22/82752a12-d746-4ebc-bd9f-187af615eec4/IGVFFI6782VXTU.yaml.gz"
    "2025/10/22/86e6c96a-fe10-41ed-9989-b377f33c3d46/IGVFFI1935TYAP.yaml.gz"
    "2025/10/22/92a26d0d-2cab-4b50-962a-41a7c62c62c9/IGVFFI3849TLEO.yaml.gz"
    "2025/10/22/92f4aa03-4877-473a-a5d5-5b8907afa27a/IGVFFI3834PTZH.yaml.gz"
    "2025/10/22/9ab393e7-f55e-4a71-aa6e-d22562204c1d/IGVFFI3495VUNK.yaml.gz"
    "2025/10/22/9f9c779b-9ef1-482b-b75a-ffd53fe4217c/IGVFFI8147EZAD.yaml.gz"
    "2025/10/22/a50b7805-958b-4a51-a2ab-0760c1ad2fdd/IGVFFI9190GLBU.yaml.gz"
    "2025/10/22/ac4496ae-0302-4719-a847-bef0845d03ce/IGVFFI9154LOZE.yaml.gz"
    "2025/10/22/b967863f-8d95-4966-bcd9-75e13f2d8dd9/IGVFFI5709ZIJM.yaml.gz"
    "2025/10/22/b9f4cf0e-2a54-4282-b6cc-296b86a1e81c/IGVFFI7592OHYU.yaml.gz"
    "2025/10/22/ccdbdb96-1abd-4140-b10a-0f2af28ef44e/IGVFFI0013TXLI.yaml.gz"
    "2025/10/22/f28920d8-dec5-40d0-a880-c86bf9362349/IGVFFI7570LGFF.yaml.gz"
    "2025/10/22/fc8ef829-4673-4208-96d0-905cf858bca2/IGVFFI3050ZFNE.yaml.gz"
)

echo ""
echo "=== Decompressing ${#SEQSPEC_RELPATHS[@]} seqspec YAML files (stripping index reads) ==="
for rel in "${SEQSPEC_RELPATHS[@]}"; do
    fname=$(basename "${rel}" .gz)         # e.g. IGVFFI6782VXTU.yaml
    echo "  ${fname}"
    gsutil cat "${BUCKET}/${DATASET}/${DATE}/${rel}" | gunzip | python3 -c "$STRIP_PY" | gsutil cp - "${PATCH_DIR}/${fname}"
done

# =============================================================================
# BARCODE HASHTAG MAP (IGVFFI3028BJFI)
# =============================================================================
# Matches Hon benchmark: prepend "hash_id\tsequence" header, keep it gzipped.
# Output file name: IGVFFI3028BJFI_with_header.tsv.gz

HASHTAG_MAP_SRC="${BUCKET}/${DATASET}/${DATE}/2024/12/20/73fc75aa-9486-48c4-9b5a-0d16a32306b2/IGVFFI3028BJFI.tsv.gz"

echo "=== Patching hashtag_map (IGVFFI3028BJFI) with header ==="
( printf "hash_id\tsequence\n"; gsutil cat "${HASHTAG_MAP_SRC}" | gunzip ) \
    | gzip \
    | gsutil cp - "${PATCH_DIR}/IGVFFI3028BJFI_with_header.tsv.gz"

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "=== Done! Files uploaded to ${PATCH_DIR}/ ==="
gsutil ls "${PATCH_DIR}/" | wc -l
echo ""
echo "Next: generate patched sample metadata CSV with paths → patch/ directory"
