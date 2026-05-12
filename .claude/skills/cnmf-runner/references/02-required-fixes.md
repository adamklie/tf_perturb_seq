# Required fixes (local PerturbNMF branch)

The `external/PerturbNMF` submodule **must** be on the local `fix/utest-oom-leak` branch for Stage 2b U-test calibration to work on production-scale datasets.

## State as of 2026-05-12

Branch: `fix/utest-oom-leak` (in [adamklie/PerturbNMF](https://github.com/adamklie/PerturbNMF), 2 commits ahead of `EngreitzLab/PerturbNMF#main`).

| Commit | What it fixes | Upstream PR |
|---|---|---|
| `9665129` | `AttributeError: 'Namespace' object has no attribute 'reference_targets'` in U-test fake-test | [#4](https://github.com/EngreitzLab/PerturbNMF/pull/4) |
| `520cb15` | Memory leak across U-test fake-test iterations and Ks (adds `del + gc.collect()`) | [#5](https://github.com/EngreitzLab/PerturbNMF/pull/5) |

Both PRs pending upstream merge.

## Verify the submodule is on the right branch

```bash
cd <REPO_ROOT>/external/PerturbNMF
git status -b
git log --oneline -3
# expect:
#   520cb15 ...
#   9665129 ...
#   <upstream main HEAD>  ...
```

If you see a different HEAD or branch:

```bash
cd <REPO_ROOT>/external/PerturbNMF
git fetch adamklie
git checkout fix/utest-oom-leak
git log --oneline -3   # confirm 520cb15 and 9665129 on top
```

(Add `adamklie` as a remote if needed: `git remote add adamklie https://github.com/adamklie/PerturbNMF`.)

## What happens without these fixes

### Without `9665129`

Stage 2b crashes early with:

```
AttributeError: 'Namespace' object has no attribute 'reference_targets'
```

The fake-test path references `args.reference_targets` but the argparse setup doesn't define it. The fix adds it with a sensible default.

### Without `520cb15`

Stage 2b appears to run, but memory grows monotonically across iterations:

- 50 fake iterations × 8 Ks = 400 inner loops
- Each loop holds onto intermediate DataFrames not released by the garbage collector
- Production datasets (~270 k cells × full TF library) hit ~700 GB by iteration ~150 and SIGKILL on the SLURM node

The fix inserts `del` + `gc.collect()` after each iteration and after each K. With it, Stage 2b runs comfortably in 256 GB.

## When PRs land upstream

Bump the submodule:

```bash
cd <REPO_ROOT>/external/PerturbNMF
git fetch origin
git checkout origin/main
git log --oneline -5
# verify 520cb15-equivalent and 9665129-equivalent are now on main

cd <REPO_ROOT>
git add external/PerturbNMF
git commit -m "chore(submodule): bump PerturbNMF (upstream merged #4 #5)"
```

Then update this reference file to drop the "required" status.

## Other known issues (no fix yet)

Filed upstream, no immediate action required, but worth knowing for failures:

| Issue | Symptom | Workaround |
|---|---|---|
| [#6](https://github.com/EngreitzLab/PerturbNMF/issues/6) | Stage 3c can't find `<run>/adata` | Create symlink: `cd <run> && ln -s Inference/adata adata` |
| [#7](https://github.com/EngreitzLab/PerturbNMF/issues/7) | Stage 3b OOM at production scale (>700 GB) | Skip Stage 3b entirely; coverage gap filled by Stage 3e Excel + Stage 3c per-TF |
| [#8](https://github.com/EngreitzLab/PerturbNMF/issues/8) | `merge_pdfs_in_folder` hangs on thousands of PDFs | Use `pdfunite` instead for the per-TF merged PDF |
| [#9](https://github.com/EngreitzLab/PerturbNMF/issues/9) | Stage 3c HVG breaks on raw counts | Pre-inject UMAP via Step 3b in `01-h5mu-prep.md`, or use `--compute_umap` upstream of Stage 1 |

## Reverting to upstream main (don't, but how)

If you need to verify upstream's behavior (e.g., to write a reproducer for a PR):

```bash
cd <REPO_ROOT>/external/PerturbNMF
git checkout origin/main
# … reproduce the failure …
git checkout fix/utest-oom-leak
```

Keep notes — the failure mode for U-test without our fixes is the canonical reproducer for PRs #4 and #5.
