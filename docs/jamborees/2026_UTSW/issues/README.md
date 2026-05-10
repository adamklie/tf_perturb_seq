# Open issues for jamboree prep

Each file here is a self-contained report on one current problem — designed to be the doc Adam reads to inspect the situation, and the artifact handed to a collaborator (Sara, Weizhou, Chikara, etc.) when starting the conversation.

| # | Issue | Owner of fix | Status |
|---|---|---|---|
| 1 | [Energy-distance p-value calibration](edistance-calibration.md) | us (re-run with HVG-subset PCA) ± Chikara (upstream guidance) | ⚠ Huangfu DE + ESC anti-conservative; runs already on Synapse |
| 2 | [Hon CM CRISPR pipeline bundle gap](hon-cm-crispr-bundle.md) | Weizhou (Hon team) | ⚠ partial bundle on Synapse, missing `pipeline_info/` |
| 3 | [Gersbach Hep deliverables (CRISPR + cNMF + energy distance)](gersbach-hep-deliverables.md) | Sara (Gersbach team) | ⚠ non-canonical bundle on Synapse, no production runs in our schema |
| 4 | [Engreitz endothelial data not on portal](engreitz-no-data.md) | Engreitz team | ☐ blocked: no data |
| 5 | [HTv2 cNMF testbed verification → production launch](htv2-cnmf-testbed.md) | us (job is running) | 🔄 in progress |
| 6 | [Per-dataset READMEs](per-dataset-readmes.md) | us | ✅ landed 2026-05-09 (refresh in place as status changes) |

Each report includes:

- **Problem** — what's broken / pending
- **Evidence** — what we observe (plots, tables, log excerpts), with paths to source files
- **Code & objects** — pointers to schemas, scripts, Synapse IDs, HPC paths
- **Hypotheses** — what we think is going on (where applicable)
- **The ask / next action** — exactly what we want from whom
- **Acceptance criteria** — how we'll know it's resolved

Status legend: ✅ done • 🔄 in progress • ⚠ partial / has caveat • ⏳ pending • ☐ blocked.
