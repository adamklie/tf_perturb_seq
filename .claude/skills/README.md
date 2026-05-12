# TFP3 Claude Skills

Project-scoped Claude Code skills for running the TFP3 pipeline. Each skill is an interactive guide that walks Claude (and you, via Claude) through a specific stage of the workflow with the project-specific paths, scripts, env vars, and known pitfalls baked in.

These live under `.claude/skills/` in this repo, so they ship to every collaborator who clones it. They are **not** under `~/.claude/skills/` (which is user-global).

## Available skills

| Skill | Status | Covers |
|---|---|---|
| [`igvf-portal-staging`](igvf-portal-staging/) | done | Stage 1: IGVF portal query, S3→GCS upload, patch samplesheet |
| [`crispr-pipeline-runner`](crispr-pipeline-runner/) | done | Stage 2: CRISPR Nextflow on GCP Batch (launch / resume / monitor / mirror outputs) |
| [`deg-calibration`](deg-calibration/) | done | Empirical-null calibration of PerTurbo per-element DEG results (parallel to QC / energy distance / cNMF) |
| `qc-runner` | planned | Stage 3: SLURM array QC (mapping_gene / mapping_guide / intended_target) |
| `energy-distance-runner` | planned | Stage 4: energy_dist_pipeline submodule |
| `cnmf-runner` | planned | Stage 5: project-specific wrapper around the external `perturbNMF-runner` skill |
| `dataset-scaffolder` | planned | Scaffold a new `datasets/<Lab>_..._TF-Perturb-seq/` from template |
| `synapse-jamboree-packaging` | planned | Package a run for `syn64423137/2026_UTSW/` upload |

## How to use a skill

You don't run a skill the way you run a script. You invoke it from Claude Code (the agent) and it then drives the workflow conversationally — asking questions, reading reference files, generating commands.

There are two ways to invoke:

**1. Slash command** (when the skill has `user_invocable: true` in its frontmatter):

```
/igvf-portal-staging
```

In a Claude Code session, this loads the skill's `SKILL.md` into context and tells Claude "drive this workflow now."

**2. Keyword-triggered**:

If you say something like `"can you stage Hon_WTC11-benchmark to GCS"`, Claude scans available skills' `description` fields and offers to load the matching one. The trigger words for each skill are listed in its `SKILL.md` `description` field (look for "Triggers on keywords like...").

You can also just say `"use the igvf-portal-staging skill"` to force it.

### What the skill actually does

A typical session:

1. You invoke `/igvf-portal-staging` or mention keywords.
2. Claude reads `SKILL.md` → asks which substep + which dataset.
3. Claude reads the matching `references/0N-*.md` for full step details.
4. Claude generates the actual commands to run (Bash invocations of the existing `setup/scripts/N_*.sh` drivers, not new code).
5. Claude shows you each command before executing. You approve, it runs, it reports output.
6. On error, Claude consults the reference file's "Common errors" section before trying anything destructive.

The skill **does not replace** the underlying driver scripts in `datasets/<ds>/setup/scripts/` — it orchestrates them, knows the right invocation order, and surfaces the per-step nuances.

## The skill pattern

Every skill in this repo follows this structure:

```
<skill-name>/
├── SKILL.md              # required — top-level orchestration + routing
├── references/           # optional — long-form per-step details
│   ├── 01-<step>.md
│   ├── 02-<step>.md
│   ├── ...
│   └── <topic>.md
└── scripts/              # optional — helper Python/shell tools
    └── *.py
```

### `SKILL.md` frontmatter

```yaml
---
name: <skill-name>             # must match directory name
description: <one-paragraph>   # used for keyword-triggered invocation — include explicit "Triggers on keywords like ..." sentence
user_invocable: true           # set to allow /<skill-name> slash command
---
```

### `SKILL.md` body conventions

Roughly in order:

1. One-paragraph statement of what the skill does and where it sits in the pipeline.
2. **Constants** block — paths, project IDs, env vars, URLs that don't change per run.
3. **Funnel / step table** — one row per substep with input → output.
4. **Step 0: identify** — first thing Claude should ask: which substep? which dataset?
5. **One short section per substep** — pointer to the reference file, the single command to run, the common per-dataset edits. Don't duplicate the reference here.
6. **Handoff** section — what feeds the next stage.
7. **Important notes** — non-obvious gotchas that apply across substeps.

Keep `SKILL.md` short — it's the routing layer. Long-form goes in `references/`.

### `references/0N-*.md` conventions

One file per substep. Numbered. Each one has:

- Prereqs (env vars, tools, portal state).
- The exact Run command (single block).
- CLI option reference for the underlying Python tool.
- Output schema (what the file/dir looks like after).
- **Nuances** — date-tag conventions, idempotency, no-op cases, streaming behavior.
- **Common errors** — error → fix mapping.
- "After Step N" — verification commands.

### `scripts/` conventions

- Idempotent by default. Add `--force` if overwriting is sometimes desired.
- Self-contained Python (no project dependencies if avoidable) or shell.
- `--help` works and is informative.
- Print "Next steps:" at the end pointing at the next user action.

## Adding a new skill

1. Copy the structure of `igvf-portal-staging/` to `<new-skill-name>/`.
2. Rewrite `SKILL.md` (name, description with trigger keywords, constants, step table).
3. Write one reference file per substep. Lift nuances from the actual driver scripts and `docs/analysis/<stage>/`.
4. Add a helper script under `scripts/` only if there's a repeated mechanical task worth automating (scaffolding, validation, generation). Otherwise skip — most skills should just orchestrate existing repo scripts.
5. Add an entry to the "Available skills" table at the top of this file.
6. Commit. Collaborators get the skill on their next pull.

### Picking trigger keywords

The `description` field is what Claude matches against your prompt to decide whether to offer the skill. Good trigger lists:

- Include the stage number (`Stage 1`, `Stage 2`, ...) and the noun phrases people actually use (`samplesheet`, `GCS upload`, `qc array`, `energy distance`, `cNMF`, `K-selection`).
- Include user-facing artifacts (`sample_metadata.csv`, `inference_mudata.h5mu`, etc.).
- Include verbs (`onboard a new dataset`, `submit qc`, `decompress`).
- Avoid keywords so generic they'd fire on unrelated requests (`run`, `pipeline` alone are too broad — pair with the stage noun).

## Cross-references

- Pipeline overview: [docs/analysis/ANALYSIS.md](../../docs/analysis/ANALYSIS.md)
- Per-stage docs: [docs/analysis/](../../docs/analysis/)
- Dataset conventions: [docs/data/DATA.md](../../docs/data/DATA.md)
- External example we modeled on: [github.com/EngreitzLab/PerturbNMF/.claude/skills](https://github.com/EngreitzLab/PerturbNMF/tree/main/.claude/skills)
