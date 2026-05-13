# Gersbach WTC11 Hepatocyte — cNMF

**Status**: ⏳ Not run yet. Awaiting **Sara** (Gersbach team) to deliver outputs in our schema-defined format. Sara likely already has a complete run; we just need it shaped to our curation rule.

See [Issue 3: Gersbach Hep deliverables](https://github.com/adamklie/tf_perturb_seq/issues/gersbach-hep-deliverables.md) §3b.

## What we want delivered

Per [`schemas/cnmf.json`](../../../schemas/cnmf.json) `bundle_inclusion_rule`:

1. **Selected-k full data** for downstream analysis — integrated MuData, all loading variants (score / tpm / consensus / starcat), cell usages, full `Eval/<sel>_<dt>/` TXT bundle, selected-k `Plot/Program_*/` + `Plot/Perturb_gene_*/` + `Annotation/<sel>_<dt>.xlsx` + `Interpretation/Summary_table/<sel>_<dt>/`.
2. **Sweep-as-provenance** — `k_selection.png` + raw stats, all-k clustering pngs, all-k `gene_spectra_score`, all-k `Eval/<k>_<dt>/` TXT bundles, the `Plot/k_selection_<run_id>/` figure folder, `README.txt` with the selection rationale.

If Sara's params differ from our reference (Hon's `030726_20iter_5KHVG_torch_halsvar_batch_e7`), that's fine — just capture them in the run's `README.txt`.

## Synapse target

`2026_UTSW/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/cnmf/<run_name>/`

## Schema + walkthrough

- Machine-readable schema: [`schemas/cnmf.json`](../../../schemas/cnmf.json)
- Analysis-level walkthrough: [`docs/analysis/cNMF_OUTPUTS.md`](../../../../../analysis/cNMF_OUTPUTS.md), [`docs/analysis/cNMF.md`](../../../../../analysis/cNMF.md)
