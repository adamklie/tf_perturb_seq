#!/bin/bash

LOCAL_DIR=/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/fastq_files
ACCESSIONS=/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/bin/1_get_data/metadata/IGVF_accessions.txt

echo "Starting download of accessions..."

# Loop through each accession ID in the file
while IFS= read -r ACCESSION || [ -n "$ACCESSION" ]; do
    OUTPUT_PATH="${LOCAL_DIR}/${ACCESSION}.fastq.gz"

    echo "Downloading ${ACCESSION}..."

    curl -L -u "${IGVF_API_KEY}:${IGVF_SECRET_KEY}" \
        "https://api.data.igvf.org/sequence-files/${ACCESSION}/@@download/${ACCESSION}.fastq.gz" \
        -o "${OUTPUT_PATH}"

    echo "Saved to ${OUTPUT_PATH}"
done < "${ACCESSIONS}"
