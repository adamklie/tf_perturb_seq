# WG6 (optional) — Predictive modeling

**Topic 3 / Figure 5.** Goal: design and prototype a model trained on the uniformly processed TFP3 outputs that can predict context-specific TF-perturbation impact.

## Questions

From [`../../WORKING_GROUPS.md`](../../WORKING_GROUPS.md):

- Brainstorm model architectures that capitalize on TFP3 data as a training source (e.g., predicting context-specific TF perturbation impact).
- Specifically define inputs, outputs, task, and validation metrics, and potentially begin implementation.
- Alternatively, apply existing tools from member labs to this data source.

## Artifacts in this folder

| ID | File | Status | What it answers | Source data |
|---|---|---|---|---|
| WG6-A | `model_task_spec.md` | ✅ ready (writeup) | Task statement, inputs, outputs, validation metrics — agreed at the jamboree | group discussion |
| WG6-B | `baseline_model_results.tsv` (+ figures) | 🟡 partial | Per-lineage CV R², feature-importance ranking, predictions vs observations | per-dataset `pval_edist_full.csv` + `tf_metadata.tsv` + (eventually) cNMF cell-state features |

## Notes

- WG6-A is a writeup deliverable agreed at the jamboree — not a derivative computation. The placeholder lives here so the folder has shape; participants fill it in during the session.
- WG6-B is a prototype baseline (likely elastic-net or random forest on TF + gene features). Not the final model; mostly a feasibility demo.

## Run the examples

No runnable examples yet — WG6 is in brainstorm stage. Start by filling out [`task_spec_template.md`](task_spec_template.md) with the group, then drop the prototype baseline as `wg6_baseline.py` next to it.
