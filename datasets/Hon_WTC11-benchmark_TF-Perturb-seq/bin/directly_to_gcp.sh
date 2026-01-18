# get s3_uri (property on file object on the portal) first i.e. SRC + INCLUDE in example below
#!/usr/bin/env bash
set -euo pipefail
gcloud config set project "igvf-pertub-seq-pipeline"
cat > role-arn.json <<'JSON'
{
  "roleArn": "arn:aws:iam::407227577691:role/S3toPertubSeqGoogleCloudTransfer"
}
JSON
# some files are on s3://igvf-public/ or s3://igvf-files/
SRC='s3://igvf-private/'
DEST='gs://igvf-pertub-seq-pipeline-data/WTC11_CM_TF_PerturbSeq/'
JOB='transferJobs/IGVFFI0115REVJ'
INCLUDE='2024/08/22/d352ce63-5b21-4db8-9192-c1855919a261/IGVFFI0115REVJ.fastq.gz'
gcloud transfer jobs create "$SRC" "$DEST" \
--name="$JOB" \
--description="Transfer IGVFFI0115REVJ" \
--source-creds-file=role-arn.json \
--include-prefixes="$INCLUDE" \
--overwrite-when=never