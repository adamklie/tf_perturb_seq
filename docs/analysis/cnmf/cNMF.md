# cNMF gene program discovery — index

Gene expression program (GEP) discovery for TFP3 is done with **consensus
non-negative matrix factorization (cNMF)** on the per-dataset
`inference_mudata.h5mu`. The wrapper that runs cNMF + downstream evaluation is
[PerturbNMF](https://github.com/EngreitzLab/PerturbNMF).

This file is an index. The two real docs are:

- **[PerturbNMF.md](PerturbNMF.md)** — how to run it: pipeline stages, project-specific
  h5mu prep convention, required upstream fixes, compute budget, what to look at
  when reviewing a run.
- **[cNMF_OUTPUTS.md](cNMF_OUTPUTS.md)** — what each output file is, where it
  lives in the per-dataset run directory, and which downstream working group
  consumes it.

The IGVF consortium publishes a shared cNMF/PerturbNMF handbook in the Engreitz
Lab's
[gene_network_evaluation](https://github.com/EngreitzLab/gene_network_evaluation)
and [PerturbNMF](https://github.com/EngreitzLab/PerturbNMF) repos — that's the
authoritative source for the tool itself. The two docs above record only what's
specific to *this* project (TFP3): which dataset directories live where, which
K values we run, what fixes we maintain on local branches, and which Synapse
folders the outputs go to.

> **Legacy note.** An earlier version of this file held the consortium-wide
> torch-cNMF instructions copy-pasted from a Google Doc, alongside instructions
> for the older `cNMF_benchmarking` tool. That content has moved upstream and
> the tool has been superseded by PerturbNMF; the per-run scripts in
> `external/cNMF_benchmarking/` are no longer the canonical entrypoint.
