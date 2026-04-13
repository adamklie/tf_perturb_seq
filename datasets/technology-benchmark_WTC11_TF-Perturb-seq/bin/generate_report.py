#!/usr/bin/env python3
"""Generate an HTML report comparing CLEANSER and Sceptre pipeline runs."""

import base64
import os
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path("/cellar/users/aklie/projects/tf_perturb_seq")
BENCHMARK_DIR = PROJECT_ROOT / "datasets" / "technology-benchmark_WTC11_TF-Perturb-seq"
RESULTS_DIR = BENCHMARK_DIR / "results" / "cross_tech_comparison"
OUTPUT_HTML = BENCHMARK_DIR / "results" / "cross_tech_comparison" / "benchmark_report.html"

RUNS = ["cleanser_unified", "sceptre_v11"]
RUN_LABELS = {"cleanser_unified": "CLEANSER (unified)", "sceptre_v11": "Sceptre v11"}

DATASETS_SHORT = {
    "Hon_WTC11-benchmark_TF-Perturb-seq": "Hon",
    "Huangfu_WTC11-benchmark_TF-Perturb-seq": "Huangfu",
    "Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3": "Gersbach GEM-Xv3",
    "Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2": "Gersbach HTv2",
    "Engreitz_WTC11-benchmark_TF-Perturb-seq": "Engreitz",
}


def embed_pdf(pdf_path: Path) -> str:
    """Convert PDF to base64-embedded <img> via pdf2image, or return placeholder."""
    try:
        from pdf2image import convert_from_path
        images = convert_from_path(str(pdf_path), dpi=150, first_page=1, last_page=1)
        import io
        buf = io.BytesIO()
        images[0].save(buf, format="PNG")
        b64 = base64.b64encode(buf.getvalue()).decode()
        return f'<img src="data:image/png;base64,{b64}" style="max-width:100%; border:1px solid #ddd; border-radius:4px;">'
    except Exception:
        # Fallback: embed PDF directly
        if pdf_path.exists():
            b64 = base64.b64encode(pdf_path.read_bytes()).decode()
            return f'<object data="data:application/pdf;base64,{b64}" type="application/pdf" width="100%" height="500px"><p>PDF: {pdf_path.name}</p></object>'
        return f'<p style="color:gray;">Plot not available: {pdf_path.name}</p>'


def load_tsv(run: str, name: str) -> pd.DataFrame | None:
    path = RESULTS_DIR / run / name
    if path.exists():
        return pd.read_csv(path, sep="\t")
    return None


def df_to_html(df: pd.DataFrame, max_rows: int = 20) -> str:
    if df is None:
        return '<p style="color:gray;">Data not available</p>'
    return df.head(max_rows).to_html(index=False, classes="data-table", border=0, float_format="%.3f")


def section_plots(title: str, pdf_names: list, description: str = "") -> str:
    """Create a section with side-by-side plots for both runs."""
    html = f'<h2>{title}</h2>\n'
    if description:
        html += f'<p>{description}</p>\n'

    for pdf_name in pdf_names:
        html += f'<h3>{pdf_name.replace(".pdf", "").replace("_", " ").title()}</h3>\n'
        html += '<div class="comparison-grid">\n'
        for run in RUNS:
            pdf_path = RESULTS_DIR / run / pdf_name
            # Also check inference subdir
            if not pdf_path.exists():
                pdf_path = RESULTS_DIR / run / "inference" / pdf_name
            html += f'<div class="run-panel"><h4>{RUN_LABELS[run]}</h4>{embed_pdf(pdf_path)}</div>\n'
        html += '</div>\n'
    return html


def build_summary_tables() -> str:
    """Build the key metrics comparison tables."""
    html = '<h2>Summary Metrics</h2>\n'

    for run in RUNS:
        summary = load_tsv(run, "metrics_summary.tsv")
        if summary is not None:
            html += f'<h3>{RUN_LABELS[run]} — Metrics Summary</h3>\n'
            html += df_to_html(summary) + '\n'

    # Compare auROC/auPRC across runs
    html += '<h3>Intended Target Performance (auROC / auPRC)</h3>\n'
    html += '<table class="data-table"><thead><tr><th>Dataset</th>'
    for run in RUNS:
        html += f'<th>{RUN_LABELS[run]} auROC</th><th>{RUN_LABELS[run]} auPRC</th>'
    html += '</tr></thead><tbody>\n'

    run_metrics = {}
    for run in RUNS:
        df = load_tsv(run, "combined_intended_target_metrics.tsv")
        if df is not None:
            run_metrics[run] = df.set_index("dataset")

    for ds_full, ds_short in DATASETS_SHORT.items():
        html += f'<tr><td>{ds_short}</td>'
        for run in RUNS:
            if run in run_metrics and ds_full in run_metrics[run].index:
                row = run_metrics[run].loc[ds_full]
                auroc = row.get("auroc", float("nan"))
                auprc = row.get("auprc", float("nan"))
                html += f'<td>{auroc:.3f}</td><td>{auprc:.3f}</td>'
            else:
                html += '<td>—</td><td>—</td>'
        html += '</tr>\n'
    html += '</tbody></table>\n'

    return html


def build_takeaways() -> str:
    """Generate key takeaways by analyzing the data."""
    html = '<h2>Key Takeaways</h2>\n<div class="takeaways">\n'

    takeaways = []

    # Compare gene metrics
    for run in RUNS:
        summary = load_tsv(run, "metrics_summary.tsv")
        if summary is not None:
            if "umi_median" in summary.columns:
                best = summary.loc[summary["umi_median"].idxmax()]
                ds_name = DATASETS_SHORT.get(best.get("dataset", ""), "Unknown")
                takeaways.append(
                    f"<b>{RUN_LABELS[run]}</b>: {ds_name} has the highest median Tx UMIs/cell "
                    f"({best['umi_median']:,.0f})"
                )

    # Compare auROC across runs
    run_aurocs = {}
    for run in RUNS:
        df = load_tsv(run, "combined_intended_target_metrics.tsv")
        if df is not None and "auroc" in df.columns:
            run_aurocs[run] = df["auroc"].mean()

    if len(run_aurocs) == 2:
        better = max(run_aurocs, key=run_aurocs.get)
        worse = min(run_aurocs, key=run_aurocs.get)
        diff = run_aurocs[better] - run_aurocs[worse]
        takeaways.append(
            f"<b>Intended target detection</b>: {RUN_LABELS[better]} achieves higher mean auROC "
            f"({run_aurocs[better]:.3f}) vs {RUN_LABELS[worse]} ({run_aurocs[worse]:.3f}), "
            f"a difference of {diff:.3f}"
        )

    # Compare guide assignment
    for run in RUNS:
        guide_df = load_tsv(run, "combined_guide_metrics.tsv")
        if guide_df is not None and "frac_cells_with_guide" in guide_df.columns:
            agg = guide_df[guide_df.get("batch", "all") == "all"] if "batch" in guide_df.columns else guide_df
            if len(agg) > 0:
                mean_frac = agg["frac_cells_with_guide"].mean()
                takeaways.append(
                    f"<b>{RUN_LABELS[run]}</b>: Average {mean_frac:.1%} of cells have at least one guide assigned"
                )

    # Dataset-specific observations
    for run in RUNS:
        df = load_tsv(run, "combined_intended_target_metrics.tsv")
        if df is not None and "auroc" in df.columns:
            worst = df.loc[df["auroc"].idxmin()]
            best = df.loc[df["auroc"].idxmax()]
            worst_name = DATASETS_SHORT.get(worst.get("dataset", ""), "Unknown")
            best_name = DATASETS_SHORT.get(best.get("dataset", ""), "Unknown")
            takeaways.append(
                f"<b>{RUN_LABELS[run]}</b>: Best intended target auROC in {best_name} "
                f"({best['auroc']:.3f}), worst in {worst_name} ({worst['auroc']:.3f})"
            )

    # General observations
    takeaways.append(
        "<b>Analysis decisions</b>: Section 1 (Cell Counts per Guide) now uses "
        "<code>n_cells_assigned</code> from the guide_assignment layer rather than "
        "<code>n_cells_detected</code> (UMI > 0), which dramatically reduces counts "
        "(e.g., CD81#strong: 62k detected → 2.5k assigned in Hon/CLEANSER)"
    )
    takeaways.append(
        "<b>Analysis decisions</b>: Notebook 4 (DEG Comparison) requires trans QC outputs "
        "which are not yet generated by the QC pipeline. These analyses are deferred."
    )

    for t in takeaways:
        html += f'<div class="takeaway-item">{t}</div>\n'

    html += '</div>\n'
    return html


def main():
    html = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>TF Perturb-seq Cross-Technology Benchmark Report</title>
<style>
    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; margin: 0; padding: 20px; background: #f8f9fa; color: #333; max-width: 1400px; margin: 0 auto; }
    h1 { color: #1a1a2e; border-bottom: 3px solid #16213e; padding-bottom: 10px; }
    h2 { color: #16213e; border-bottom: 2px solid #e2e8f0; padding-bottom: 8px; margin-top: 40px; }
    h3 { color: #2d3748; }
    h4 { color: #4a5568; margin-bottom: 5px; text-align: center; }
    .comparison-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin: 15px 0; }
    .run-panel { background: white; padding: 15px; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
    .data-table { border-collapse: collapse; width: 100%; background: white; border-radius: 8px; overflow: hidden; box-shadow: 0 1px 3px rgba(0,0,0,0.1); margin: 10px 0; }
    .data-table th { background: #16213e; color: white; padding: 10px 12px; text-align: left; font-size: 13px; }
    .data-table td { padding: 8px 12px; border-bottom: 1px solid #e2e8f0; font-size: 13px; }
    .data-table tr:hover { background: #f7fafc; }
    .takeaways { background: white; padding: 20px; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
    .takeaway-item { padding: 10px 15px; border-left: 4px solid #3182ce; margin: 10px 0; background: #ebf8ff; border-radius: 0 4px 4px 0; }
    .toc { background: white; padding: 20px; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); margin: 20px 0; }
    .toc ul { list-style: none; padding-left: 0; }
    .toc li { padding: 4px 0; }
    .toc a { color: #3182ce; text-decoration: none; }
    .toc a:hover { text-decoration: underline; }
    .methodology { background: #fffbeb; padding: 15px; border-radius: 8px; border-left: 4px solid #d69e2e; margin: 15px 0; }
    code { background: #edf2f7; padding: 2px 6px; border-radius: 3px; font-size: 12px; }
    .meta { color: #718096; font-size: 14px; }
</style>
</head>
<body>
<h1>TF Perturb-seq Cross-Technology Benchmark Report</h1>
<p class="meta">
    Comparing 5 datasets across 2 inference pipelines: <b>CLEANSER (unified)</b> and <b>Sceptre v11</b><br>
    Datasets: Hon, Huangfu, Gersbach GEM-Xv3, Gersbach HTv2, Engreitz (all WTC11 cell line)<br>
    Generated from: <code>technology-benchmark_WTC11_TF-Perturb-seq</code>
</p>

<div class="toc">
<h3>Table of Contents</h3>
<ul>
    <li><a href="#takeaways">Key Takeaways</a></li>
    <li><a href="#methodology">Methodology & Analysis Decisions</a></li>
    <li><a href="#summary">Summary Metrics</a></li>
    <li><a href="#gene-qc">Gene Expression QC</a></li>
    <li><a href="#guide-qc">Guide QC</a></li>
    <li><a href="#combined-metrics">Combined Key Metrics</a></li>
    <li><a href="#repression">Repression Efficiency</a></li>
    <li><a href="#guide-capture">Guide Capture & Assignment</a></li>
    <li><a href="#guide-uniformity">Guide Uniformity</a></li>
    <li><a href="#intended-target">Intended Target Performance</a></li>
    <li><a href="#lane-metrics">Per-Lane Metrics</a></li>
    <li><a href="#inference">Calibration & Inference</a> (CLEANSER only)</li>
</ul>
</div>
"""

    # Key takeaways
    html += '<div id="takeaways">\n'
    html += build_takeaways()
    html += '</div>\n'

    # Methodology
    html += """
<div id="methodology">
<h2>Methodology & Analysis Decisions</h2>
<div class="methodology">
<p><b>Pipeline comparison:</b> Each dataset was processed through both CLEANSER (unified) and Sceptre v11
inference pipelines. QC metrics were computed locally using the project's QC pipeline
(<code>src/tf_perturb_seq/qc/</code>), which extracts gene mapping, guide mapping, and intended
target metrics from the MuData outputs.</p>

<p><b>Guide assignment vs detection:</b> A key distinction in this analysis is between guide
<em>detection</em> (any cell with UMI > 0 for a guide) and guide <em>assignment</em> (cells where
the guide was called as present by the assignment algorithm). Detection counts are inflated by
ambient RNA and barcode hopping — e.g., top guides can appear detected in 80%+ of cells despite
being assigned to only ~3%. The "Cell Counts per Guide" section uses assignment counts for accurate
representation.</p>

<p><b>Calibration:</b> Empirical null calibration (t-fit method) was applied to trans-effect results
to obtain calibrated p-values. This corrects for dataset-specific null distributions and enables
fair cross-dataset comparison of significant hits.</p>

<p><b>Trans QC:</b> Notebook 4 (DEG Comparison) is currently skipped because trans-effect QC metrics
have not been generated by the QC pipeline. The trans-effect analysis is covered in the 3_inference/
calibration and comparison notebooks instead.</p>
</div>
</div>
"""

    # Summary tables
    html += '<div id="summary">\n'
    html += build_summary_tables()
    html += '</div>\n'

    # Gene QC
    html += '<div id="gene-qc">\n'
    html += section_plots(
        "Gene Expression QC",
        ["gene_metrics_comparison.pdf"],
        "Median UMI, median genes, and mitochondrial fraction across datasets."
    )
    html += '</div>\n'

    # Guide QC
    html += '<div id="guide-qc">\n'
    html += section_plots(
        "Guide QC",
        ["guide_metrics_comparison.pdf"],
        "Guide UMI, assignment rates, and cells per guide across datasets."
    )
    html += '</div>\n'

    # Combined key metrics
    html += '<div id="combined-metrics">\n'
    html += section_plots(
        "Combined Key Metrics",
        ["key_combined_metrics_comparison.pdf"],
        "Single-row comparison of gene + guide metrics: Median Tx UMIs/cell, Median sgRNA UMIs/cell, "
        "Cells w/ sgRNA UMIs > 0, Mean assigned sgRNA/cell, Cells w/ exactly 1 assigned sgRNA."
    )
    html += '</div>\n'

    # Repression efficiency
    html += '<div id="repression">\n'
    html += section_plots(
        "Repression Efficiency",
        [
            "repression_efficiency_boxplot.pdf",
            "repression_efficiency_all_targeting_boxplot.pdf",
            "repression_efficiency_by_label_boxplot.pdf",
        ],
        "Fold change distributions for positive controls only, all targeting guides, "
        "and split by guide type. Red dashed line at FC=1.0 (no effect)."
    )
    html += '</div>\n'

    # Guide capture
    html += '<div id="guide-capture">\n'
    html += section_plots(
        "Guide Capture & Assignment",
        [
            "guide_cells_assigned_ranked_barplot.pdf",
            "guide_cells_assigned_violin.pdf",
        ],
        "Per-guide cell assignment counts (not UMI detection). "
        "Ranked waterfall and violin distributions."
    )
    # Guide capture is only in cleanser for now
    for pdf in ["guide_most_variable_heatmap.pdf", "guide_top_bottom_assigned.pdf",
                 "guide_overlap_summary.pdf", "guide_detection_patterns.pdf"]:
        pdf_path = RESULTS_DIR / "cleanser_unified" / pdf
        if pdf_path.exists():
            html += f'<h3>{pdf.replace(".pdf", "").replace("_", " ").title()} (CLEANSER)</h3>\n'
            html += f'<div class="run-panel">{embed_pdf(pdf_path)}</div>\n'
    html += '</div>\n'

    # Guide uniformity
    html += '<div id="guide-uniformity">\n'
    html += section_plots(
        "Guide Uniformity",
        ["guide_uniformity_lorenz.pdf"],
        "Lorenz curves showing how uniformly guides are represented. "
        "Closer to diagonal = more uniform library representation."
    )
    html += '</div>\n'

    # Intended target
    html += '<div id="intended-target">\n'
    html += section_plots(
        "Intended Target Performance",
        [
            "auroc_auprc_comparison.pdf",
            "log2fc_violin_by_label.pdf",
            "log2fc_correlation_matrix.pdf",
        ],
        "auROC/auPRC for intended target detection, log2FC distributions by label, "
        "and pairwise correlation of guide effects across datasets."
    )
    html += '</div>\n'

    # Per-lane metrics
    html += '<div id="lane-metrics">\n'
    html += section_plots(
        "Per-Lane Metrics",
        ["metrics_by_lane_comparison.pdf"],
        "Metrics broken down by individual sequencing lanes within each dataset."
    )
    html += '</div>\n'

    # Inference (calibration)
    html += '<div id="inference">\n'
    html += '<h2>Calibration & Inference</h2>\n'
    html += '<p>Empirical null calibration results. Currently available for CLEANSER unified.</p>\n'

    inference_pdfs = [
        "ntc_null_distributions.pdf",
        "significant_hits_by_category.pdf",
        "direct_target_fc_boxplot.pdf",
        "direct_target_correlation_matrix.pdf",
        "direct_target_volcano.pdf",
    ]
    for run in RUNS:
        inf_dir = RESULTS_DIR / run / "inference"
        if inf_dir.exists():
            html += f'<h3>{RUN_LABELS[run]}</h3>\n'
            for pdf_name in inference_pdfs:
                pdf_path = inf_dir / pdf_name
                if pdf_path.exists():
                    html += f'<h4>{pdf_name.replace(".pdf", "").replace("_", " ").title()}</h4>\n'
                    html += f'<div class="run-panel">{embed_pdf(pdf_path)}</div>\n'
    html += '</div>\n'

    html += """
<hr>
<p class="meta" style="text-align:center; margin-top:40px;">
    Generated by <code>generate_report.py</code> | TF Perturb-seq Cross-Technology Benchmark
</p>
</body>
</html>"""

    OUTPUT_HTML.write_text(html)
    print(f"Report written to: {OUTPUT_HTML}")
    print(f"File size: {OUTPUT_HTML.stat().st_size / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    main()
