# Handing off to the external PerturbNMF skill

For everything except input conversion and the TFP3 prep steps, the **external** PerturbNMF skill is authoritative. This file documents the handoff.

## Where the external skill lives

**Upstream (canonical):** https://github.com/EngreitzLab/PerturbNMF/tree/main/.claude/skills/perturbNMF-runner — always available, no clone required.

**Local clone path (when present):**

```
<REPO_ROOT>/external/PerturbNMF/.claude/skills/perturbNMF-runner/
├── SKILL.md
├── references/
│   ├── 01-inference.md
│   ├── 02-evaluation.md
│   ├── 03-calibration.md
│   ├── 04-visualization.md
│   ├── 05-annotation-summary.md
│   ├── data-format-spec.md
│   └── parameter-catalog.md
└── scripts/
    ├── generate_slurm.py
    ├── generate_h5mu_structure.py
    ├── validate_data.py
    └── prepare_guide_data.py
```

**Important:** `external/PerturbNMF/` is **gitignored** in this repo — it doesn't ship with a clone. You need to clone it manually before invoking the skill or running anything that needs the upstream code:

```bash
# Clone from adamklie's fork (includes the required fix/utest-oom-leak branch)
git clone https://github.com/adamklie/PerturbNMF external/PerturbNMF
cd external/PerturbNMF
git checkout fix/utest-oom-leak
```

(Note: `docs/analysis/cnmf/PerturbNMF.md` calls this a "submodule" but it isn't currently registered as one in `.gitmodules`. Treat it as a manually-maintained sibling clone.)

If you only need to read the external skill's docs (not run the code), fetch them directly from upstream:

```bash
# E.g. via gh CLI without cloning:
gh api repos/EngreitzLab/PerturbNMF/contents/.claude/skills/perturbNMF-runner/SKILL.md \
  -H "Accept: application/vnd.github.raw"
```

Read the SKILL.md for stage mechanics: each stage's `--stage` value, conda env, reference file, parameters, and SLURM resource estimation.

## When to invoke which skill

```
TFP3 cnmf-runner             →  External perturbNMF-runner
─────────────────────────────────────────────────────────
h5mu → h5ad conversion       →  (after conversion done)
                                Stage 1 inference (--stage inference-torch)
TFP3 prep (h5mu_for_eval,    →  Stage 2a evaluation
inject_umap, guide_annotation)   (--stage evaluation)
                             →  Stage 2b U-test calibration
                                (--stage u-test-calibration; requires our fix branch)
                             →  Stage 3a K-selection (--stage k-selection)
                             →  Stage 3c perturbed-gene (--stage perturbed-gene)
                             →  Stage 3e Excel summary (--stage excel-summary)
Synapse upload               ←  (after all stages complete)
```

Rule of thumb: anything involving `generate_slurm.py`, `--stage`, conda env activation, or upstream parameters → external skill. Anything involving `<DS>/<RUN>/cnmf/<cnmf_run_name>/Data/` setup, the fix branch, or `upload_to_synapse.py` → this skill.

## How to invoke the external skill

From inside a Claude Code session, the external skill is auto-discovered (`.claude/skills/` under repo root *and* under any submodule). Just say:

```
/perturbNMF-runner
```

Or use keyword triggers (any of: `perturbnmf`, `cnmf inference`, `k-selection`, `u-test calibration`, `submit slurm`, `--stage evaluation`, etc.). Claude loads `external/PerturbNMF/.claude/skills/perturbNMF-runner/SKILL.md` and walks the stage.

You can also pre-emptively read the external SKILL.md (`Read external/PerturbNMF/.claude/skills/perturbNMF-runner/SKILL.md`) to get the constants block (Sherlock paths, conda envs, GWAS data, reference GTF) and the stage table. **Note:** the upstream skill is calibrated for Stanford Sherlock; for TFP3 on UCSD nrnb, paths differ — see "Path overrides" below.

## Path overrides (Sherlock → nrnb)

| Sherlock constant | TFP3 (UCSD nrnb) equivalent |
|---|---|
| `PIPELINE_ROOT=/oak/stanford/groups/engreitz/Users/ymo/Tools/PerturbNMF` | `<REPO_ROOT>/external/PerturbNMF` |
| `SKILL_DIR=/oak/.../PerturbNMF/.claude/skills/perturbNMF-runner` | `<REPO_ROOT>/external/PerturbNMF/.claude/skills/perturbNMF-runner` |
| `CONDA_BASE=/oak/.../miniforge3` | conda not used on nrnb; the SLURM scripts source `.venv` instead |
| `DEFAULT_EMAIL=ymo@stanford.edu` | `aklie@ucsd.edu` |
| `GWAS_DATA=/oak/.../Stage2_Evaluation/Resources/OpenTargets_L2G_Filtered.csv.gz` | `<REPO_ROOT>/external/PerturbNMF/src/Stage2_Evaluation/Resources/OpenTargets_L2G_Filtered.csv.gz` |
| `REFERENCE_GTF=/oak/.../IGVFFI9573KOZR.gtf.gz` | `<REPO_ROOT>/ref/IGVFFI9573KOZR.gtf.gz` (gitignored; download from IGVF portal) |

When invoking `generate_slurm.py`, pass these as `--script_output_path`, `--gwas_data`, `--reference_gtf` etc. The external skill's `01-inference.md` etc. enumerate which stages need which.

Also: the upstream skill assumes SLURM partition + `-C "GPU_MEM:..."` constraints from Sherlock. On nrnb, override:
- Partition: `carter-gpu`
- GPU constraint: `--gres=gpu:a30:1` (don't pass `-C`)
- Memory: per `03-compute-budget.md`

## TFP3 conventions to pass through

When the external skill asks for inputs, give it:

| External skill arg | TFP3 value |
|---|---|
| Run name | `MMDDYY_<short_desc>_torchcnmf_<rationale>` (e.g. `042926_huangfu_de_torchcnmf_KskillA`) |
| `<out_dir>` (Result/) | `<REPO_ROOT>/datasets/<DS>/<RUN>/cnmf/<cnmf_run_name>/Result` |
| `--script_output_path` | `<REPO_ROOT>/datasets/<DS>/<RUN>/cnmf/<cnmf_run_name>/Script/<run_name>_<stage>.sh` |
| K sweep | `30 50 60 80 100 200 250 300` (TFP3 production default) |
| HVG | 5000 |
| sel_thresh | 2.0 |
| iterations | 20 |
| Stage 1 method flags | `--method halsvar --batch_correction` |

For the rationale on these defaults, see `docs/analysis/cnmf/PerturbNMF.md`.

## What about the per-dataset Script/ scripts?

Each cNMF run directory has many SLURM scripts under `Script/`:

```
<cnmf_run_name>/Script/
├── Convert_file_adata.py                       # (this skill) — Stage 1 input builder
├── prepare_h5mu_for_eval.sh                    # (this skill) — h5mu prep 3a
├── inject_umap_into_h5mu.sh                    # (this skill) — h5mu prep 3b
├── <date>_<short>_torchcnmf_<rationale>_inference.sh   # (external skill output) Stage 1
├── <date>_<short>_torchcnmf_<rationale>_finish.sh      # (external skill output) Stage 1 cleanup
├── cNMF_evaluation_pipeline.sh                 # (external skill output) Stage 2a
├── cNMF_evaluation_trait_only.sh               # (external skill output) Stage 2a trait follow-up
├── U-test_perturbation_calibration.sh          # (external skill output) Stage 2b — requires fix branch
├── cNMF_k_selection.sh                         # (external skill output) Stage 3a
├── cNMF_perturbed_gene_analysis_<K>_<sel>.sh   # (external skill output) Stage 3c
├── cNMF_program_analysis_<K>_<sel>.sh          # (external skill output, but Stage 3b deferred)
├── cNMF_compile_excel_summary.sh               # (external skill output) Stage 3e
└── upload_to_synapse.py                        # (this skill) — Stage 5 upload
```

The external skill's `generate_slurm.py` produces stages-1-through-3e scripts. Each gets a hashed name based on its parameters and goes into `Script/`. Re-running `generate_slurm.py` for the same stage overwrites in place.

## Quick sanity sequence

After Stage 1 finishes:

```bash
ls <cnmf_run_name>/Result/Inference/adata/cNMF_*.h5mu        # one per K
ls <cnmf_run_name>/Result/Inference/adata/cNMF_*_structure.txt
```

After Stage 2 finishes:

```bash
ls <cnmf_run_name>/Result/Evaluation/<K>_2_0/
head <cnmf_run_name>/Result/Evaluation/<K>_2_0/<K>_perturbation_association_results_all.txt
```

After Stage 3a finishes:

```bash
ls <cnmf_run_name>/Result/Plot/k_selection/K-selection_panel_2.0.png   # the one plot the team reviews together
```

After Stage 3e finishes:

```bash
ls <cnmf_run_name>/Result/Interpretation/Summary_table/<K>_2_0/cNMF_<K>_2_0.xlsx
```

Then proceed to Synapse upload (parent SKILL.md Step 5).
