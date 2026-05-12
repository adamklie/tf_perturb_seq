#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Gersbach Hep — generate sample_metadata.csv from explicit measurement-set list
#
# Hep has no overarching analysis set, so we enumerate 47 sub-pool GEX
# measurement sets explicitly (with their paired CRISPR auxiliary sets via
# _samplesheets/hep_ms_aux_pairs.tsv) and call the from-ms-list generator.
# =============================================================================

BASE_DIR=/carter/users/aklie/projects/tf_perturb_seq
DATASET_NAME=Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq
DATASET_DIR=${BASE_DIR}/datasets/${DATASET_NAME}
SCRIPT=${BASE_DIR}/src/tf_perturb_seq/portal/generate_per_sample_from_ms_list.py

# Hep-specific portal references
GUIDE_DESIGN=IGVFFI8270UPKB        # harmonized_guide_file_poolabcdf (47-MS link is library IGVFDS3299AXST)
BARCODE_ONLIST_FALLBACK=IGVFFI9487JPEN  # 10x v2 onlist (737K-august-2016) — matches Sara's is_10x3v3=false canonical config; portal field is null on all 47 Hep MS

# Seqspec fallbacks (Hep portal has no seqspecs; using Hon CM's 10x 3' v3 YAMLs;
# confirm with Sara's actual seqspecs before submitting the run)
RNA_SEQSPEC=${BASE_DIR}/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/seqspec/rna_seqspec.yml
SGRNA_SEQSPEC=${BASE_DIR}/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/seqspec/guide_seqspec.yml

MS_LIST=${DATASET_DIR}/inputs/meta/hep_measurement_sets.txt
AUX_MAP=${DATASET_DIR}/inputs/meta/hep_ms_aux_pairs.tsv
OUTPUT=${DATASET_DIR}/samplesheets/sample_metadata.csv

echo "=========================================="
echo "Generate Hep Per-Sample Metadata (47 MS)"
echo "=========================================="
echo "MS list:     ${MS_LIST}"
echo "Aux map:     ${AUX_MAP}"
echo "Guide:       ${GUIDE_DESIGN}"
echo "Onlist:      ${BARCODE_ONLIST_FALLBACK}"
echo "RNA seqspec: ${RNA_SEQSPEC}"
echo "gRNA seqspec:${SGRNA_SEQSPEC}"
echo "Output:      ${OUTPUT}"
echo ""

python3 ${SCRIPT} \
    --measurement_sets ${MS_LIST} \
    --aux_map ${AUX_MAP} \
    --guide_design ${GUIDE_DESIGN} \
    --barcode_onlist ${BARCODE_ONLIST_FALLBACK} \
    --rna_seqspec ${RNA_SEQSPEC} \
    --sgrna_seqspec ${SGRNA_SEQSPEC} \
    --output ${OUTPUT}

echo ""
echo "Done. Output: ${OUTPUT}"
