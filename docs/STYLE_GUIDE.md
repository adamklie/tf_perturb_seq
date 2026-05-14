---
title: "Style guide"
subtitle: "Editorial conventions for jamboree documentation"
author:
  - name: Adam Klie
    affiliation: UC San Diego
date-modified: last-modified
---

Editorial conventions for jamboree documentation under [`docs/jamborees/2026_UTSW/`](jamborees/2026_UTSW/). The goal is a plain-language, consistent voice that reads cleanly on GitHub and on the rendered ReadTheDocs site.

## Voice and tone

- As a rule, simpler is better. If a sentence, paragraph, table, or section can be cut without losing what a reader needs, cut it.
- Direct and concise. Shorter sentences beat longer ones.
- Plain language over jargon; when jargon is needed, define it on first use and link to the glossary.
- No hyperbolic language. No "unprecedented," "revolutionary," "cutting-edge," or similar.
- Don't use the word "bundle." Say what you actually mean — "the CRISPR pipeline outputs," "the dataset directory," "the Synapse folder," "the agreed-on layout," depending on context. "Bundle" is vague and overloaded.
- Don't use the word "canonical." Same reason as "bundle": vague and overloaded. Say what you actually mean — "the agreed-on," "the shared," "the project-standard," or drop the qualifier entirely.

## Formatting

- Em dashes are allowed in H1 subtitles (`# Title — subtitle`) and for parenthetical asides where commas would be ambiguous. Don't sprinkle them in lists.
- No bullet points inside narrative sections. Prose is the default. Bullets are reserved for genuinely enumerable items (datasets, working groups, file outputs, command steps).
- Sentence case for headings. Never Title Case.
- Oxford comma.
- No emojis. Anywhere. Not in prose, not in status tables, not in headings. Use text labels for status: the canonical vocabulary is `ready` / `caveat` / `blocked`. Include a one-line legend caption under any status table if the column isn't already self-explanatory. Don't substitute alternative vocabularies (`done/wip/todo`, `complete/partial/missing`, etc.).
- Always include links when applicable. Internal links use relative paths (`[label](../guide/X.md)`), not GitHub URLs to the same repo — relative paths work both on GitHub and in the ReadTheDocs build.
- Synapse IDs are written as a backticked, linked accession: [`syn74834227`](https://www.synapse.org/Synapse:syn74834227). Link every mention, not just the first.
- Only link to a GitHub issue if a real issue exists (i.e., a numeric URL like `https://github.com/adamklie/tf_perturb_seq/issues/42`). Don't link to slug-style placeholders such as `.../issues/edistance-calibration.md` — those aren't real URLs. If the work is tracked but not yet ticketed, write `*[FILL IN issue link]*`.
- For known problems with a dataset or output, use the section heading `## Issues` (or `**Issues:**` as a subsection). Don't use "Caveats", "Concerns", "Warnings". Each issue entry must lead with the linked GitHub issue (`[Issue #42](https://github.com/.../issues/42): one-line summary`), or `*[FILL IN issue link]*` if not yet ticketed.
- Use `*[FILL IN]*` (italicized brackets) as the canonical placeholder for incomplete content. It's grep-able, scans visually as "missing," and gets cleared during cleanup passes. Don't substitute `TODO`, `TBD`, `???`, or empty parens.

## What not to do

- **Do not present analytical results or findings in these docs.** Jamboree docs are navigational and methodological — they tell readers where the data lives, what the columns mean, and how to load it. Findings, claims, interpretations, and figures-of-results belong in manuscripts, notebooks, or slides; link to those instead of restating them.
- **Do not make unvetted scientific claims.** This extends the no-findings rule to *general biology* assertions, not just analyses from this project. Don't write "SOX17 marks endoderm," "POU5F1/NANOG/SOX2 score highly in ESCs," "ISL1 is the cardiomyocyte master regulator," or similar. Even when they're well-known biology, these are claims that should be vetted, cited, or omitted — not casually dropped into shared participant-facing docs. When you need a concrete example to explain a method, use a generic placeholder like *"a lineage-defining program"* or *"the cell-cycle program"* without naming specific marker genes; or cite a primary reference.
- **Do not link to absolute GitHub URLs for files in this repo.** Use a relative path (`[label](../guide/X.md)`) — `https://github.com/adamklie/tf_perturb_seq/blob/main/...` breaks if the repo is renamed/moved, doesn't work in the ReadTheDocs build, and is harder to grep.
- **Do not reference personal HPC paths, SSH targets, or login hosts.** No `ssh aklie@nrnb-login.ucsd.edu`, no `/cellar/users/aklie/...` paths, no `sbatch /cellar/.../script.sh`. Jamboree docs are read by collaborators who don't have UCSD HPC access. For local-execution context, refer to repo-relative paths (`datasets/<ds>/...`, `scripts/...`) and let internal HPC details live in private notes or per-dataset README sections elsewhere.
- Do not paste in PR text or email content that has not been reviewed for confidentiality.
- Do not commit generated PDFs or HTML outputs; those are build artifacts.
