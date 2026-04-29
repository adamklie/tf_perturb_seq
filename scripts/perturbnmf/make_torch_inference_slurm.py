#!/usr/bin/env python3
"""
Generate a Carter-cluster SLURM script for torch-cNMF (GPU) Stage 1 inference.

Wraps the PerturbNMF `torch_cnmf_inference_pipeline.py` with the right SBATCH
header, conda env activation, and PYTHONPATH, then writes a `.sh` to
`<script_dir>/<run_name>_inference.sh`. Does NOT submit — caller decides.

Carter-specific defaults:
  - partition = carter-gpu
  - GPU = a30 (24 GB) by default; rtx5000 (16 GB) also available
  - venv = /cellar/users/aklie/projects/tf_perturb_seq/.venv (has torch_cnmf, nmf-torch from github)
  - PYTHONPATH adds external/PerturbNMF/src

Usage:
    python scripts/perturbnmf/make_torch_inference_slurm.py \
        --counts_fn datasets/<dataset>/PerturbNMF/Data/<file>.h5ad \
        --out_dir   datasets/<dataset>/PerturbNMF/Result \
        --script_dir datasets/<dataset>/PerturbNMF/Script \
        --run_name MMDDYY_<dataset>_torchcnmf_K<vals>
"""
import argparse
import os

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--counts_fn', required=True)
p.add_argument('--out_dir', required=True)
p.add_argument('--script_dir', required=True)
p.add_argument('--run_name', required=True)
p.add_argument('--K', default='5 7 10')
p.add_argument('--numiter', type=int, default=10)
p.add_argument('--numhvgenes', type=int, default=2000)
p.add_argument('--sel_thresh', default='2.0')
p.add_argument('--mem', default='96G')
p.add_argument('--time', default='06:00:00')
p.add_argument('--cpus', type=int, default=4)
p.add_argument('--gpu_type', default='a30', help='rtx5000 (16GB) or a30 (24GB)')
p.add_argument('--email', default='aklie@ucsd.edu')
args = p.parse_args()

PROJECT_ROOT='/cellar/users/aklie/projects/tf_perturb_seq'
PNMF_SRC=f'{PROJECT_ROOT}/external/PerturbNMF/src'
VENV=f'{PROJECT_ROOT}/.venv/bin/activate'
SCRIPT=f'{PNMF_SRC}/Inference/torch-cNMF/Slurm_Version/torch_cnmf_inference_pipeline.py'
LOG_DIR=f'{args.out_dir}/{args.run_name}/Inference/logs'

slurm = f'''#!/bin/bash
#SBATCH --job-name={args.run_name}
#SBATCH --partition=carter-gpu
#SBATCH --time={args.time}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.cpus}
#SBATCH --mem={args.mem}
#SBATCH --gres=gpu:{args.gpu_type}:1
#SBATCH --output={LOG_DIR}/%j.out
#SBATCH --error={LOG_DIR}/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user={args.email}

set -euo pipefail
mkdir -p "{LOG_DIR}"
echo "Job: $SLURM_JOB_ID @ $(hostname) $(date)"
nvidia-smi || true

source "{VENV}"
export PYTHONPATH="{PNMF_SRC}:${{PYTHONPATH:-}}"

python -u {SCRIPT} \\
    --counts_fn "{args.counts_fn}" \\
    --output_directory "{args.out_dir}" \\
    --run_name "{args.run_name}" \\
    --species human \\
    --K {args.K} \\
    --numiter {args.numiter} \\
    --numhvgenes {args.numhvgenes} \\
    --sel_thresh {args.sel_thresh} \\
    --seed 14 \\
    --algo halsvar \\
    --mode batch \\
    --tol 1e-4 \\
    --use_gpu \\
    --categorical_key batch \\
    --gene_names_key symbol \\
    --run_factorize --run_refit --run_compile_annotation --run_diagnostic_plots

echo "Done @ $(date)"
'''
out = os.path.join(args.script_dir, f'{args.run_name}_inference.sh')
open(out, 'w').write(slurm)
os.chmod(out, 0o755)
print(out)
