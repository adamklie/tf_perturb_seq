#!/usr/bin/env python3
"""
Generate a Carter-cluster SLURM script for sk-cNMF (CPU) Stage 1 inference.

CPU counterpart of make_torch_inference_slurm.py. Use sk-cNMF for the FINAL K
once selected (high tolerance, many iters). Use torch-cNMF for the K-selection
sweep. See PerturbNMF/.claude/skills/perturbNMF-runner/SKILL.md for guidance.

NOTE: requires the Engreitz fork of `cnmf` (extra args like --algo). Pip-install:
    pip install git+https://github.com/EngreitzLab/sk_cNMF.git
The public PyPI `cnmf` package does NOT accept these args and will fail.

Usage:
    python scripts/perturbnmf/make_skcnmf_inference_slurm.py \
        --counts_fn datasets/<dataset>/PerturbNMF/Data/<file>.h5ad \
        --out_dir   datasets/<dataset>/PerturbNMF/Result \
        --script_dir datasets/<dataset>/PerturbNMF/Script \
        --run_name MMDDYY_<dataset>_skcnmf_K<value>
"""
import argparse

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--dataset', required=True)
parser.add_argument('--counts_fn', required=True)
parser.add_argument('--out_dir', required=True)
parser.add_argument('--script_dir', required=True)
parser.add_argument('--run_name', required=True)
parser.add_argument('--K', default='50')
parser.add_argument('--numiter', type=int, default=10)
parser.add_argument('--numhvgenes', type=int, default=5000)
parser.add_argument('--sel_thresh', default='2.0')
parser.add_argument('--mem', default='192G')
parser.add_argument('--time', default='12:00:00')
parser.add_argument('--cpus', type=int, default=8)
parser.add_argument('--partition', default='carter-compute')
parser.add_argument('--email', default='aklie@ucsd.edu')
args = parser.parse_args()

PROJECT_ROOT = '/cellar/users/aklie/projects/tf_perturb_seq'
PERTURBNMF_SRC = f'{PROJECT_ROOT}/external/PerturbNMF/src'
PERTURBNMF_SHIMS = f'{PROJECT_ROOT}/external/PerturbNMF/src/_local_shims'
VENV = f'{PROJECT_ROOT}/.venv/bin/activate'
INFERENCE_SCRIPT = f'{PERTURBNMF_SRC}/Inference/sk-cNMF/Slurm_Version/sk-cNMF_batch_inference_pipeline.py'
LOG_DIR = f'{args.out_dir}/{args.run_name}/Inference/logs'

slurm = f'''#!/bin/bash
#SBATCH --job-name={args.run_name}
#SBATCH --partition={args.partition}
#SBATCH --time={args.time}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.cpus}
#SBATCH --mem={args.mem}
#SBATCH --output={LOG_DIR}/%j.out
#SBATCH --error={LOG_DIR}/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user={args.email}

set -euo pipefail
mkdir -p "{LOG_DIR}"

echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"
echo "Partition: $SLURM_JOB_PARTITION"

# Activate project venv (provides cnmf, scanpy, anndata, mudata, muon, torch_cnmf shim)
source "{VENV}"
export PYTHONPATH="{PERTURBNMF_SHIMS}:{PERTURBNMF_SRC}:${{PYTHONPATH:-}}"
echo "Python: $(which python)"
echo "PYTHONPATH: $PYTHONPATH"

python -u {INFERENCE_SCRIPT} \\
    --counts_fn "{args.counts_fn}" \\
    --output_directory "{args.out_dir}" \\
    --run_name "{args.run_name}" \\
    --species human \\
    --K {args.K} \\
    --numiter {args.numiter} \\
    --numhvgenes {args.numhvgenes} \\
    --sel_thresh {args.sel_thresh} \\
    --seed 14 \\
    --algo mu \\
    --init random \\
    --tol 1e-4 \\
    --max_NMF_iter 500 \\
    --categorical_key batch \\
    --guide_names_key guide_names \\
    --guide_targets_key guide_targets \\
    --guide_assignment_key guide_assignment \\
    --gene_names_key symbol \\
    --run_factorize \\
    --run_refit \\
    --run_complie_annotation \\
    --run_diagnostic_plots

echo "Done @ $(date)"
'''

import os
out = os.path.join(args.script_dir, f'{args.run_name}_inference.sh')
with open(out, 'w') as f:
    f.write(slurm)
os.chmod(out, 0o755)
print(out)
