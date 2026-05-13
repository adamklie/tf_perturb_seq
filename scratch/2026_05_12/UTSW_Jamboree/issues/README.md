# Open issues for jamboree prep

Each file here is a self-contained report on one current problem — designed to be the doc Adam reads to inspect the situation, and the artifact handed to a collaborator (Sara, Weizhou, Chikara, etc.) when starting the conversation.

| # | Issue | Owner of fix | Status |
|---|---|---|---|
| 1 | [Energy-distance p-value calibration](edistance-calibration/) | us (apply Chikara's recommended fix, re-run on HPC) ± Chikara (upstream guidance) | ⚠ Huangfu DE + ESC anti-conservative; Hon CM + HTv2 calibrated; cause unresolved |
| 2 | [Hon CM CRISPR pipeline bundle gap](hon-cm-crispr-bundle.md) | Weizhou (Hon team) | ⚠ partial bundle on Synapse, missing `pipeline_info/` |
| 3 | [Gersbach Hep deliverables (CRISPR + cNMF + energy distance)](gersbach-hep-deliverables.md) | Sara (Gersbach team) | ⚠ non-canonical bundle on Synapse, no production runs in our schema |
| 4 | [Engreitz endothelial data not on portal](engreitz-no-data.md) | Engreitz team | ☐ blocked: no data |
| 5 | [HTv2 cNMF testbed verification → production launch](htv2-cnmf-testbed.md) | us | ✅ resolved — testbed cleared and Huangfu DE/ESC cNMF runs landed on Synapse (`syn74893844`, `syn74893846`); see [`synapse_paths.tsv`](../synapse_paths.tsv) |
| 6 | [Per-dataset READMEs](per-dataset-readmes.md) | us | ✅ landed 2026-05-09 (refresh in place as status changes) |

Each report includes:

- **Problem** — what's broken / pending
- **Observations** — what we measure (numbers, plots, log excerpts), with paths to source files
- **What we've ruled out** — with the evidence that ruled each candidate out
- **Code & objects** — pointers to schemas, scripts, Synapse IDs, HPC paths
- **Open questions** — what's still unresolved (not asserting a working hypothesis — leave that for the collaborator we're handing the doc to)
- **The ask / next action** — exactly what we want from whom
- **Acceptance criteria** — how we'll know it's resolved

Status legend: ✅ done • 🔄 in progress • ⚠ partial / has caveat • ⏳ pending • ☐ blocked.
