# Presentation Notes — Benchmark Update

Notes for revising the 03/20/2026 benchmark deck.
PDF of current version: `docs/2026_03_20_benchmark_update.pdf`
Google Slides: https://docs.google.com/presentation/d/1fyr9OrqZmboFOfFIwkGFrFS6yRxnXyzQrPNeZfhttt4/edit

---

## Slide 1: Title
- No changes needed

## Slide 2: Overview table (5 labs, technologies, cell counts)
- Slide is fine as-is
- Could consider adding guide library info (shared ~50-gene library, ~416 guides) to reinforce that the perturbation space is constant (**TODO**: add as a hidden slide, have a previously used version I can use)

## Slide 3: Methods
- Currently brief — tells what, not why
- Add a one-liner on what Sceptre vs CLEANSER actually differ on (guide assignment strategy)
- Flag the pipeline version mismatch (CLEANSER and Sceptre ran on slightly different versions) — important caveat
- Move the downstream analysis roadmap (QC → repression efficiency → calibration → TF enrichment) to its own slide

## NEW SLIDE: Analysis roadmap
- After methods, add a slide that maps out the analysis steps so the audience knows the arc of the talk
- Could be a simple flow diagram or bullet list

## NEW SLIDE(S): Dataset characteristics & filtering funnel (before metrics)

We extracted upstream metrics from the GCP pipeline outputs. These should be 1-3 new slides that set the stage before showing the QC bar charts.

### Data extracted (in `results/cross_tech_comparison/`):
- `per_ms_read_counts.tsv` — per measurement set, per modality
- `dataset_summary_metrics.tsv` — compiled summary
- `filtering_funnel.tsv` — cell counts at each filtering stage

### Slide: Dataset characteristics table

| Dataset | MS | scRNA Reads | Guide Reads | Guide % Aligned | Final Cells | Median UMI/cell | Reads/Cell |
|---------|:--:|:-----------:|:-----------:|:---:|:-----------:|:---:|:---------:|
| Hon | 3 | 1,947M | 568M | 62.5% | 73,458 | 4,889 | ~26K |
| Huangfu | 4 | 1,497M | 381M | 16.1%* | 141,251 | 3,856 | ~11K |
| Engreitz | 16 | 6,426M | 768M | 90.9% | 467,119 | 1,098 | ~14K |
| Gersbach GEM-Xv3 | 4 | 6,404M | 477M | 95.8% | 67,436 | 27,184 | ~95K |
| Gersbach HTv2 | 4 | 6,383M | 692M | 95.9% | 43,376 | 15,640 | ~147K |

*Huangfu guide alignment at 16.1% is a known aspect of 10x 3' v3 chemistry, already troubleshot.

**Key takeaway for audience:** The 5 datasets have vastly different sequencing depths and cell counts. Gersbach has ~6.4B reads across ~50-70K cells (very deep), Engreitz has ~6.4B reads across ~467K cells (shallow), Hon/Huangfu are moderate.

### Slide: Pipeline config differences

| Parameter | Hon | Huangfu | Engreitz | Gersbach GEM-Xv3 | Gersbach HTv2 |
|-----------|-----|---------|----------|-------------------|---------------|
| `QC_min_genes_per_cell` | 500 | 500 | 500 | **1000** | **1000** |
| `QC_pct_mito` | **15%** | 20% | 20% | 20% | 20% |
| `capture_method` | CROP-seq | CROP-seq | CROP-seq | **direct-capture** | CROP-seq |
| `ENABLE_DATA_HASHING` | **true** | false | false | false | false |
| `is_10x3v3` | false | **true** | false | false | false |
| `reverse_complement_guides` | true | **false** | true | true | true |

Source: GCP `pipeline_info/params_*.json` (actual configs used in runs)

### Slide: Filtering funnel (side-by-side bar chart or Sankey)

| Stage | Hon | Huangfu | Engreitz | GEM-Xv3 | HTv2 |
|-------|----:|--------:|---------:|--------:|-----:|
| Concatenated (unfiltered) | 580,607 | 1,413,616 | 1,275,903 | 1,275,090 | 742,959 |
| After knee filter | 146,509 | 141,885 | 582,627 | 67,688 | 79,368 |
| After mito filter | 146,484 | 141,251 | 467,119 | 67,436 | 43,376 |
| After guide intersection | 73,459 | 141,251 | 467,119 | 67,436 | 43,376 |
| After HTO filter | 73,458 | — | — | — | — |

Key stories:
- **Hon**: Loses 50% at guide intersection (HTO singlet filtering), almost nothing at mito
- **Engreitz**: Loses 20% at mito filter (583K → 467K)
- **Gersbach HTv2**: Loses 45% at mito (79K → 43K) — biggest mito loss
- **Gersbach GEM-Xv3**: Aggressive knee filter (1.3M → 68K) driven by `min_genes=1000`
- **Huangfu**: Clean — minimal loss at mito and guide stages

**TODO**: More upstream metrics to explore (saturation curves from BUS files, evaluation plots, seqspeccheck data). See `docs/analysis/CRISPR_PIPELINE_OUTPUTS.md`.

## Slide 4: Metrics overview (bar charts)
- With the new upstream slides, the audience will now understand WHY the bar charts look the way they do
- Gene expression metrics are identical between CLEANSER and Sceptre (same upstream) — could show Tx UMIs once and note this
- Engreitz low sgRNA assignment is likely a samplesheet issue with how measurement sets were split — we think this is fixed now, need to confirm
- **Talking point:** Engreitz has 16 measurement sets (8 pairs of _1/_2) due to how the CC-Perturb-seq data was split in the samplesheet. This means each "lane" in the per-lane plots is actually a merged pair, and the 8 Engreitz bars aren't directly comparable to the 3-4 true lanes in other datasets. This split may affect cell calling and filtering behavior.
- Consider whether this slide is even needed in its current form, or if the new dataset characteristics table replaces it

## Slides 5, 7, 9: CLEANSER vs Sceptre comparisons
- **Consider removing or demoting to appendix.** The Sceptre runs lack standard error / confidence intervals, making the comparison incomplete. Since we can't calibrate Sceptre trans results either (missing `log2_fc_std`), the CLEANSER-only story is cleaner.
- If we keep any CLEANSER vs Sceptre content, consolidate into one summary slide rather than three.
- Slide 5 (scatter: guide assignment & target detection) overlaps with the new dataset characteristics table.

## Slide 6: Repression efficiency (auROC/auPRC + fold change boxplots)
- Rework to show CLEANSER only, cross-dataset comparison
- Drop the Sceptre side-by-side
- Key question this slide should answer: "how well does each technology detect intended target knockdown?"

## Slide 8: Cross-dataset log2FC correlation matrices
- CLEANSER only
- This is a good slide — shows whether the same perturbations give similar effect sizes across technologies
- Could add annotation highlighting which pairs have highest/lowest correlation and why (e.g., sequencing depth, capture method)

## Slide 10: Guide representation
- Current version is too busy — distributions, overlap, variable guides all crammed in
- Rework: Flag that guide assignment distributions vary across technologies, and this depends on upstream factors (read depth, filtering, capture method, etc.)
- Use a worked example: pick one positive control guide (e.g., CD81#strong) and show cells/guide, fold change, DEG count across all 5 datasets
- Second panel: simple bar chart of median cells/guide per dataset as library-level context
- Move "most variable guides" and detailed overlap to appendix

## Slide 11: Inference calibration (NTC null distributions)
- CLEANSER only, already is
- Needs to walk through the methodology simply:
  1. **Problem**: Pipeline gives posterior p-values per perturbation × gene pair, but not calibrated for hit calling
  2. **Key insight**: Non-targeting controls define what "noise" looks like — they shouldn't have real effects
  3. **Method**: Compare each targeting element's test statistic against NTC distribution → empirical p-value → BH FDR
  4. **Show**: NTC distributions across datasets as evidence the null is well-behaved
- Could be 1 slide with a clear diagram or 2 slides (concept then results)
- This sets up everything that follows (DEG counts, overlap, TF enrichment)

## Slides 12-15: DEG calls and cross-dataset discrepancies — KEY OPEN QUESTION

These slides are where we see major discrepancies that need serious investigation before we can present confidently. The same perturbation in the same cell type gives wildly different DEG counts across technologies.

**What we see:**
- Slide 12: Trans DEG counts range from ~1.4K (Engreitz) to ~4.2K (Huangfu) at FDR < 0.1
- Slide 13: Positive controls show very different # DEGs per element across datasets (same guide, same target)
- Slide 14: UpSet plots show most DEGs are dataset-specific, not shared (BEX3: 151 unique to one dataset vs 11 shared across all)
- Slide 15: TF enrichment heatmaps show different patterns for the same TF across datasets

**Candidate explanations to investigate:**
1. **Power differences** — reads/cell ranges from 11K (Huangfu) to 147K (Gersbach HTv2). More depth = more power = more DEGs. Can we normalize for this?
2. **Cell counts per element** — varies by guide assignment and filtering. Directly affects statistical power.
3. **Pipeline config** — different QC thresholds remove different numbers of cells, changing the denominator
4. **Calibration method sensitivity** — does the t-fit null behave the same across datasets with very different cell/read counts?
5. **Guide assignment** — same guide, different capture methods → different cell assignments → different effect estimates

**What we need to do:**
- [ ] For positive controls: plot cells/element vs DEGs/element across datasets — is it purely a power issue?
- [ ] Downsample to matched cells/element and rerun calibration — do discrepancies disappear?
- [ ] Check if the discrepancies correlate with reads/cell or cells/element rather than technology per se
- [ ] Look at effect sizes (log2FC) rather than significance — are the directions and magnitudes consistent even if significance differs?
- [ ] Consider whether FDR < 0.1 is the right threshold, or if we need dataset-specific thresholds

**Analysis priorities when we come back to this:**
- Focus on effect size (log2FC) correlation rather than significance — more comparable across power levels
- Worked examples: pick positive controls, show log2FC for direct target + top trans DEGs across all 5 datasets
- Scatter plots: Dataset A vs Dataset B log2FC for all tested genes for a given perturbation
- TF enrichment discrepancies need specific investigation — the heatmaps show different patterns for the same TF across datasets, which is harder to explain by power alone
- If log2FCs correlate well but DEG counts don't → power story. If log2FCs also diverge → something deeper (guide assignment, technical artifact)

**For the presentation:** Frame honestly as "here's what we see, here's what we think is driving it, here's what we're investigating."

## Slide 16: TODOs

---

## Old FG meeting announcements

Our next FG meeting is scheduled for Thursday, March 12th @ 10-11am PT / 1-2pm ET. Here is the planned agenda:

Presentation from Anish Karpurapu, Gersbach Lab (talk title: DeltaNMF: A Two-Stage Neural NMF for Differential Gene Program Discovery) — 30m
TFP3 benchmark datasets update @Adam Klie — 15m
Discuss contributions to the flagship ("Data and Analysis Products" list and timelines) @Jesse Engreitz — 15m

Best,
Sara, Adam
