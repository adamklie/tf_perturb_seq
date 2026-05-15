#!/usr/bin/env bash
#SBATCH -J hon_cm_sync
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH --partition=carter-compute
#SBATCH -t 08:00:00
#SBATCH -o /cellar/users/aklie/projects/tf_perturb_seq/scratch/sync_logs/hon_cm_sync.%j.out
#SBATCH -e /cellar/users/aklie/projects/tf_perturb_seq/scratch/sync_logs/hon_cm_sync.%j.err

mkdir -p /cellar/users/aklie/projects/tf_perturb_seq/scratch/sync_logs

# Source bashrc to pick up SYNAPSE_AUTH_TOKEN (errors in module load OK)
set +e
source ~/.bashrc
set -e
set -uo pipefail

source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate

python <<'PYEOF'
import os, synapseclient
import synapseutils as su

syn = synapseclient.Synapse()
syn.login(authToken=os.environ["SYNAPSE_AUTH_TOKEN"])

LOCAL = "/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/2026_04_19_no_spacer/crispr_pipeline"
SYN_ID = "syn74520421"

print(f"Syncing {SYN_ID} -> {LOCAL}", flush=True)
files = su.syncFromSynapse(syn, SYN_ID, path=LOCAL)
print(f"DONE_SYNC. {len(files)} files synced.", flush=True)
PYEOF
