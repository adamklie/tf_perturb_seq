# Container + submodule quirks

The energy-distance pipeline runs inside an apptainer container with a specific numpy version, and the submodule pin has a few known idiosyncrasies. These quirks are responsible for most of the silent failures.

## Container

```
CONTAINER_IMAGE:  docker://docker.io/takechikara/energy_distance_env:latest
LOCAL_PATH:       /cellar/users/aklie/opt/containers/edist_pipeline.sif
PULL_ON_DEMAND:   yes — runner pulls automatically if missing
```

### Pre-installed packages (inside container)

- `muon` 0.1.7 in `/app/.venv` — **use this one**
- `numpy` 1.26.4 — **fixed; do not upgrade**
- PyTorch (CUDA-enabled) for the step 1 / step 2 GPU ops
- scanpy, anndata, mudata, scikit-learn, scipy, etc.

### The numpy / muon trap

> The runner has an inline comment: "Use container's pre-installed muon 0.1.7 (in /app/.venv). DO NOT pip install muon to /tmp/muon_deps — that pulls a newer numpy (>=2.0) whose pickled artifacts cannot be deserialized by the container's numpy 1.26.4 in steps 1/2/2.1."

**What this looks like in practice:** if you (or a future PR) try to `pip install muon` from inside the container to a side-load directory, pip resolves a newer muon that requires numpy ≥ 2.0. The PCA pickle is then written with numpy 2.0 dtypes. Step 1 / step 2 / step 2.1 fail to deserialize with errors like:

```
ModuleNotFoundError: No module named 'numpy._core'
```

or

```
ValueError: numpy.dtype size changed, may indicate binary incompatibility
```

The fix is always: use the pre-installed muon. The runner does this correctly by default.

### Bind mounts

The runner sets up three bind mounts:

```
--bind $(dirname "$MUDATA_PATH"):$(dirname "$MUDATA_PATH")
--bind "${OUTPUT_DIR_ABS}:${OUTPUT_DIR_ABS}"
--bind "${REPO_ROOT}:${REPO_ROOT}"
```

So the container sees the MuData, the output dir, and the repo (for the preprocess script + pipeline bin) under their actual host paths.

## Submodule (`external/energy_dist_pipeline`)

```
PIN:        Specific commit of github.com/Chikara-Takeuchi/energy_dist_pipeline
BIN_DIR:    external/energy_dist_pipeline/bin/
INIT:       git submodule update --init external/energy_dist_pipeline
```

### Step 1 filename typo

At the currently pinned commit, step 1 is `1_filtereing_gRNA.py` (sic — typo carried from upstream). The upstream `main` branch has it renamed to `1_filtering_gRNA.py`. The runner handles both:

```bash
STEP1_SCRIPT="${PIPELINE_BIN}/1_filtering_gRNA.py"
[[ ! -f "$STEP1_SCRIPT" ]] && STEP1_SCRIPT="${PIPELINE_BIN}/1_filtereing_gRNA.py"
```

When the submodule is bumped past the rename commit, both branches work without change.

### Other pipeline scripts in `bin/`

- `2_e_distance_nontargeting.py` — step 2
- `2_1_Plot_figure.py` — step 2.1
- `3_e_distance_among_regions.py` — step 3

These names haven't changed upstream.

### PYTHONPATH

The runner sets `PYTHONPATH=<PIPELINE_BIN>` when invoking each step. This is so the step scripts can import their helper modules from sibling files in `bin/`.

## Wrapper (`scripts/run_energy_distance_pipeline.sh`)

### Config heredoc clobber

The runner unconditionally regenerates `config1_2.json` and `config3.json` from heredocs on every invocation. **Your edits get clobbered.** Workarounds:

1. Patch the heredoc directly. Per `[[feedback_dataset_local_scripts]]`, scripts under `scripts/` are frozen — copy to `datasets/<ds>/bin/` and edit there.
2. Bypass the runner for the affected step (see `references/03-outputs-step3.md` §"Running step 3" for the apptainer-exec pattern).

### Idempotency

Skipped when outputs already exist:
- MuData download (if `<OUTPUT_DIR>/inference_mudata.h5mu` exists).
- Container pull (if `.sif` exists).
- Preprocess (if `pca_dataframe.pickle` AND `gRNA_dict.pickle` exist).

**Not skipped:**
- Config heredoc (regenerated every run).
- Steps 1 / 2 / 2.1 (no existence check — they always run).
- Step 3 (gated by `--run-step3`).

To re-run step 2 with different params: edit the heredoc patch (in your `datasets/<ds>/bin/` copy of the runner), then resubmit. Steps 0 + 1 will re-execute too; they're cheap compared to step 2's permutation test.

## Path conventions in flux

The runner doesn't care where `--output-dir` points — both layouts work:

- **Legacy (still in active wrappers):** `datasets/<DS>/results/energy_distance/<RUN>/`
- **Canonical (per `docs/data/DATA.md`):** `datasets/<DS>/<RUN>/energy_distance/`

New `5_run_energy_distance.sh` wrappers should use the canonical layout. Existing wrappers are not retroactively migrated — keeps git history clean.

## Cross-references

- Upstream pipeline source: https://github.com/Chikara-Takeuchi/energy_dist_pipeline
- Upstream wrapper source: https://github.com/Chikara-Takeuchi/energy_dist_TFperturb
- Container source: https://hub.docker.com/r/takechikara/energy_distance_env
- TFP3-side notes: `docs/analysis/energy_dist/ENERGY_DISTANCE.md`
