# WG6 — Predictive model task spec (template)

Fill this in during the jamboree with WG6 participants. Keep it tight; the goal
is to define one well-scoped task + baseline + validation plan that the group
can prototype in 1–2 days.

## Task statement

> _Given <inputs>, predict <output>. Why this prediction is useful: <impact>._

## Inputs

- **From this jamboree's data**: which artifacts (CRISPR pipeline / cNMF / ED / TF metadata)? Which datasets?
- **External inputs (if any)**: multiome, motifs, reference annotations.
- **Featurization**: how raw artifacts become model-ready tensors.

## Outputs

- Shape, units, interpretation.
- Per-perturbation, per-cell, per-program, etc.

## Training / test split

- **Train**: which datasets?
- **Held-out**: which dataset(s) and/or which perturbations?
- **Why this split** controls for leakage / probes generalization.

## Metrics

- Headline metric (single number for slide-deck reporting).
- Secondary metrics (calibration, per-class breakdowns).
- Statistical baseline for comparison (random, perturbation-mean, etc.).

## Baseline implementation

- Architecture, library, compute footprint.
- Pointer to scratch notebook / script when started.

## Validation plan

- How do we know the prediction means something biologically?
- Sanity checks (e.g., predict known TF→target relationships from existing GRN).
- External cross-reference dataset.

## Open questions / decisions for the group

- ...

---

*Created from template `working_groups/wg6_predictive_modeling/task_spec_template.md`. Add a `task_spec.md` in this folder with the filled-in spec.*
