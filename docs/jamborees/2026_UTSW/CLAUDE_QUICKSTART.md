# Claude Code Quickstart — 2026 UTSW Jamboree

This page is for jamboree participants new to this repo and/or to Claude Code. By the end you'll know how to ask Claude to walk you through the canonical workflows.

## 0. Setup (5 min)

```bash
# 1. Install Claude Code (CLI). See https://docs.claude.com/en/docs/claude-code/quickstart
# 2. Clone the repo
git clone --recurse-submodules https://github.com/adamklie/tf_perturb_seq.git
cd tf_perturb_seq

# 3. Install Python deps (optional — only needed if you want to run things, not just browse)
uv sync

# 4. Start Claude Code from the repo root
claude
```

When you start `claude` in this directory:
- It auto-loads `CLAUDE.md` (top-level orientation) and `.claude/skills/` (the workflow skills).
- The repo's `.claude/settings.json` pre-allows safe read-only commands (`ls`, `gsutil ls`, `git status`, `squeue`, `synapse list`, etc.) so you won't be prompted for routine inspections.

## 1. Pipeline overview

This project has a 5-stage pipeline. Each stage has a skill that walks Claude through running it:

```
[1] IGVF Portal → GCS    /igvf-portal-staging
[2] CRISPR Nextflow      /crispr-pipeline-runner
[3] Local QC              /qc-runner
[4] Energy distance      /energy-distance-runner
[5] cNMF / PerturbNMF    /cnmf-runner
```

Plus three cross-cutting skills:

```
DEG calibration          /deg-calibration       (calibrate PerTurbo per-element p-values)
Dataset scaffolder       /dataset-scaffolder    (onboard a new dataset)
Synapse packaging        /synapse-jamboree-packaging   (mirror to syn64423137/2026_UTSW/)
```

A complete list with descriptions: [`.claude/skills/README.md`](../../../.claude/skills/README.md).

## 2. How to invoke a skill

Two ways:

**Slash command** — explicit:

```
/qc-runner
```

**Keyword trigger** — Claude offers the right skill when you mention what you're trying to do:

```
You: "I want to run QC on the Hon CM seqspec_v3 inference mudata"
Claude: [offers to load qc-runner skill]
```

Every skill's `SKILL.md` lists its trigger keywords near the top.

## 3. Five canonical prompts

Copy-paste these to get rolling. Replace `<DATASET>` and `<RUN>` with your dataset/run.

### A. "Walk me through CRISPR pipeline for a new run"

```
Using crispr-pipeline-runner, scaffold a new run for <DATASET> with
RUN_LABEL=test_sweep and DATA_DATE=<YYYY_MM_DD>. Use Hon CM seqspec_v3 as the
base-config. Show me the generated files before I edit anything.
```

### B. "Run QC across all my mirrored runs"

```
Using qc-runner, build a samples.tsv covering every dataset/run with an
inference_mudata.h5mu present, then show me the sbatch command to submit.
```

### C. "Calibrate DEGs for one run"

```
Using deg-calibration with method t-fit, calibrate
datasets/<DATASET>/<RUN>/crispr_pipeline/pipeline_outputs/perturbo_trans_per_element_output.tsv.gz
and write outputs to datasets/<DATASET>/<RUN>/calibration/.
```

### D. "Package a dataset for Synapse"

```
Using synapse-jamboree-packaging, dry-run the CRISPR pipeline mirror for
<DATASET> to syn64423137/2026_UTSW/. Use the HPC staging workdir.
```

### E. "Audit current state"

```
Show me the jamboree packaging status: which datasets have which bundles on
Synapse, and what's pending. Pull from synapse_paths.tsv.
```

## 4. Common gotchas

| Symptom | Cause | Fix |
|---|---|---|
| "Permission denied" prompts for every `gsutil ls` | First time on this machine | They'll auto-allow after the first approval — or rely on the pre-built project allowlist in `.claude/settings.json` |
| `synapseclient.core.exceptions.SynapseAuthenticationError: 401` | `SYNAPSE_AUTH_TOKEN` not in current shell | Invoke via `zsh -ic '...'` (token is in `~/.zshrc`) |
| Stage 1 scripts fail with "no such directory" | `BASE_DIR` hard-codes Adam's path | Edit `BASE_DIR` at the top of each `setup/scripts/N_*.sh` for your machine |
| `module: command not found` | Not on UCSD HPC | These scripts only work on HPC; run from `nrnb-login.ucsd.edu` |
| Stage 2 (`-c $CONFIG`) missing | Older driver scripts omit `-c` | Use the canonical pattern from Hon cardio or `crispr-pipeline-runner`'s scaffolder |
| Stage 4 step 3 doesn't run | Default skips it (need cutoffs first) | Edit `config3.json`, pass `--run-step3` |
| Stage 5 U-test OOMs | Submodule isn't on the fix branch | `cd external/PerturbNMF && git checkout fix/utest-oom-leak` |

## 5. Where to dig deeper

| Question | File |
|---|---|
| What is this project? | [`CLAUDE.md`](../../../CLAUDE.md) |
| Pipeline overview | [`docs/analysis/ANALYSIS.md`](../analysis/ANALYSIS.md) |
| Per-dataset layout | [`docs/data/DATA.md`](../data/DATA.md) |
| Jamboree plan | [`docs/jamborees/2026_UTSW/README.md`](README.md) |
| Skill catalog + pattern | [`.claude/skills/README.md`](../../../.claude/skills/README.md) |
| Today's tasks | [`docs/TODAY.md`](../TODAY.md) (gitignored — local) |

## 6. Asking Claude for help

When you're stuck, try one of:

- `"What does <error message> mean and how do I fix it?"`
- `"Show me how the <stage> works for this dataset"`
- `"Audit the state of <X> and tell me what's missing"`
- `"Pull the schema for <bundle type> and compare against what's on Synapse"`

Claude has access to all the references files in each skill, plus the full `docs/` tree. It will read what it needs.

## 7. What to avoid (Claude's perspective)

By default Claude will pause and ask before:
- Anything destructive (`rm -rf`, `git reset --hard`, `gsutil rm`)
- Anything that affects shared state (`git push`, `sbatch`, Synapse uploads)
- Anything that hits a billable service (GCS transfers, Tower)

If you want it to be more autonomous within a specific scope, say so explicitly: `"Go ahead and sbatch — I've reviewed the script."`

If you want it to never do something, add the rule to your local `CLAUDE.md` or the user-global `~/.claude/CLAUDE.md`.
