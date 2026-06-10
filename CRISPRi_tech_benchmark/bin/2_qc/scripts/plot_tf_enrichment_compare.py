"""TF enrichment comparison across benchmark runs.

For each dataset × filter config, the pipeline (with ENABLE_BENCHMARK=true) runs
a Fisher-exact enrichment test asking whether the trans DEGs of each TF
perturbation overlap with that TF's known ENCODE ChIP-seq targets.

The pipeline tests 5 promoter window widths around each gene's TSS:
    [500, 1000, 2500, 5000, 10000] bp
(hard-coded in `bin/tf_benchmark.py`; default for single-window mode is 5 kb).

Each `tf/benchmark_output/benchmark_tables/enrichment_all.tsv` has 5 TFs ×
5 windows = 25 rows (one ChIP dataset per TF).

Two views in this script:
  HEADLINE: "best across windows" — for each (dataset × run × TF), take the
    most generous window (min p, max OR). Shows whether the TF can be enriched
    AT ALL with permissive window choice.
  PER-WINDOW: shows how stable each TF's enrichment is to window choice.

Outputs to: results/cross_tech_comparison/tf_enrichment/
"""
from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings('ignore')

PROJECT_ROOT = Path('/cellar/users/aklie/projects/tf_perturb_seq')
BASE = PROJECT_ROOT / 'datasets' / 'technology-benchmark_WTC11_TF-Perturb-seq'
RESULTS = BASE / 'results' / 'cross_tech_comparison'
OUT = RESULTS / 'tf_enrichment'
OUT.mkdir(parents=True, exist_ok=True)

sys.path.append(str(PROJECT_ROOT / 'config'))
from loader import load_colors
cfg = load_colors('technology-benchmark_WTC11_TF-Perturb-seq')
DATASET_COLORS = cfg['dataset_colors']
DATASET_ORDER = cfg['dataset_order']

RUNS = [
    ('cleanser_extremes_200_mito_15pc',  '200 UMI'),
    ('cleanser_500_mito_15pc',           '500 UMI'),
    ('cleanser_800_mito_15pc',           '800 UMI'),
    ('cleanser_extremes_2000_mito_15pc', '2000 UMI'),
    ('cleanser_knee2_mito_15pc',         'Knee2'),
    ('scrublet_off_cleanser_800_mito_15pc', 'scrub OFF\n(sceptre 800)'),
    ('scrublet_on_sceptre_800_mito_15pc',   'scrub ON\n(sceptre 800)'),
]
RUN_LABELS = [lbl for _, lbl in RUNS]
PROMOTER_WINDOWS = [500, 1000, 2500, 5000, 10000]
WINDOW_LABELS = ['500 bp', '1 kb', '2.5 kb', '5 kb', '10 kb']

SHORT = {
    'Hon_WTC11-benchmark_TF-Perturb-seq': 'Hon',
    'Huangfu_WTC11-benchmark_TF-Perturb-seq': 'Huangfu',
    'Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3': 'Gersbach\nGEM-Xv3',
    'Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2': 'Gersbach\nHTv2',
    'Engreitz_WTC11-benchmark_TF-Perturb-seq': 'Engreitz',
}


def load_one(dataset, run):
    p = BASE.parent / dataset / 'runs' / run / 'tf' / 'benchmark_output' / 'benchmark_tables' / 'enrichment_all.tsv'
    return pd.read_csv(p, sep='\t').assign(dataset=dataset, run=run) if p.exists() else None


# === Load all ===
frames = [load_one(d, r) for d in DATASET_ORDER for r, _ in RUNS]
frames = [f for f in frames if f is not None]
all_df = pd.concat(frames, ignore_index=True)
all_df['neg_log10_p'] = -np.log10(all_df['pvalue'].clip(lower=1e-300))
all_df['run_label'] = all_df['run'].map(dict(RUNS))
print(f'Loaded {len(frames)} (dataset × run) files; {len(all_df)} total rows')
print(f'Promoter windows tested: {sorted(all_df["promoter_window_width"].unique())} bp')

ds_in = [d for d in DATASET_ORDER if d in all_df['dataset'].values]
TFs = sorted(all_df['TF_display'].unique())


# === HEADLINE VIEW: best across windows per TF ===
agg_best = all_df.groupby(['dataset', 'run', 'run_label', 'TF_display']).agg(
    best_pvalue=('pvalue', 'min'),
    best_odds_ratio=('odds_ratio', 'max'),
    best_window=('promoter_window_width', 'first'),  # placeholder; will recompute below
).reset_index()
# best window per (dataset, run, TF) = window with min p
idx_min = all_df.loc[all_df.groupby(['dataset', 'run', 'TF_display'])['pvalue'].idxmin()]
agg_best = agg_best.merge(
    idx_min[['dataset', 'run', 'TF_display', 'promoter_window_width']].rename(
        columns={'promoter_window_width': 'best_window_bp'}),
    on=['dataset', 'run', 'TF_display'], how='left'
).drop(columns='best_window')
agg_best['neg_log10_p'] = -np.log10(agg_best['best_pvalue'].clip(lower=1e-300))
agg_best['is_significant'] = agg_best['best_pvalue'] < 0.05
agg_best['enriched'] = agg_best['is_significant'] & (agg_best['best_odds_ratio'] > 1)


def plot_count_heatmap(df, value_col, title, fname, cmap='YlGnBu', vmax=5, label_fmt='{:.0f}'):
    counts = df.groupby(['dataset', 'run_label']).agg(n=(value_col, 'sum')).reset_index()
    piv = counts.pivot(index='dataset', columns='run_label', values='n').reindex(ds_in).reindex(columns=RUN_LABELS)
    fig, axes = plt.subplots(1, 2, figsize=(15, 4.5))
    ax = axes[0]
    im = ax.imshow(piv.values, aspect='auto', cmap=cmap, vmin=0, vmax=vmax)
    ax.set_xticks(range(len(RUN_LABELS))); ax.set_xticklabels(RUN_LABELS, rotation=20, ha='right', fontsize=8)
    ax.set_yticks(range(len(ds_in))); ax.set_yticklabels([SHORT.get(d, d) for d in ds_in], fontsize=9)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.iloc[i, j]
            if pd.notna(v):
                color = 'white' if v >= vmax * 0.6 else 'black'
                ax.text(j, i, label_fmt.format(v), ha='center', va='center', fontsize=10, color=color)
            else:
                ax.text(j, i, '—', ha='center', va='center', fontsize=10, color='gray')
    ax.set_title(title, fontsize=11, fontweight='bold')
    plt.colorbar(im, ax=ax, fraction=0.04)

    ax = axes[1]
    for ds in ds_in:
        sub = counts[counts['dataset'] == ds].set_index('run_label').reindex(RUN_LABELS)
        ax.plot(RUN_LABELS, sub['n'], 'o-', color=DATASET_COLORS[ds], label=SHORT.get(ds, ds), linewidth=2, markersize=7)
    ax.set_ylabel(value_col)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.set_ylim(-0.3, vmax + 0.5)
    ax.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.02, 1.0))
    ax.tick_params(axis='x', rotation=20)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(OUT / f'{fname}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(OUT / f'{fname}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved: {fname}.pdf / .png')


# 01: # TFs enriched (best across windows)
plot_count_heatmap(agg_best, 'enriched', '# TFs enriched (best across windows; p<0.05 & OR>1)',
                   '01_n_enriched_TFs_best_across_windows', vmax=5)


def plot_per_TF_heatmap(df, value_col, title, fname, cmap='RdPu', vmin=0, vmax=None,
                        annot_thresh=None, value_fmt='{:.1f}', box_thresh=None):
    if vmax is None:
        vmax = df[value_col].quantile(0.95)
    fig, axes = plt.subplots(1, len(ds_in), figsize=(3.2 * len(ds_in), 4.5), sharey=True)
    for ax, ds in zip(axes, ds_in):
        sub = df[df['dataset'] == ds].pivot_table(
            index='TF_display', columns='run_label', values=value_col, aggfunc='max'
        ).reindex(index=TFs, columns=RUN_LABELS)
        ax.imshow(sub.values, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xticks(range(len(RUN_LABELS))); ax.set_xticklabels(RUN_LABELS, rotation=30, ha='right', fontsize=7)
        ax.set_yticks(range(len(TFs))); ax.set_yticklabels(TFs, fontsize=9)
        ax.set_title(SHORT.get(ds, ds).replace('\n', ' '), fontsize=10, color=DATASET_COLORS[ds], fontweight='bold')
        for i in range(sub.shape[0]):
            for j in range(sub.shape[1]):
                v = sub.iloc[i, j]
                if pd.notna(v):
                    color = 'white' if v > vmax * 0.6 else 'black'
                    ax.text(j, i, value_fmt.format(v), ha='center', va='center', fontsize=7, color=color)
                else:
                    ax.text(j, i, '—', ha='center', va='center', fontsize=7, color='gray')
                if box_thresh is not None and pd.notna(v) and v >= box_thresh:
                    ax.add_patch(plt.Rectangle((j-0.45, i-0.45), 0.9, 0.9, fill=False, edgecolor='black', linewidth=1.5))
    fig.suptitle(title, fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(OUT / f'{fname}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(OUT / f'{fname}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved: {fname}.pdf / .png')


# 02 + 03: per-TF best across windows
plot_per_TF_heatmap(agg_best, 'neg_log10_p',
                    'Per-TF best -log10(p) — best across 5 promoter windows (boxed = p<0.05)',
                    '02_per_TF_neg_log10_p_best_across_windows',
                    box_thresh=-np.log10(0.05))
plot_per_TF_heatmap(agg_best, 'best_odds_ratio',
                    'Per-TF best odds ratio — best across 5 promoter windows (>1 = enriched)',
                    '03_per_TF_odds_ratio_best_across_windows',
                    cmap='Blues', vmin=0.5)


# 04: distribution boxplot (best across windows)
fig, ax = plt.subplots(figsize=(13, 5))
positions = []; data = []; colors_list = []; labels = []
xx = 0
for run_lbl in RUN_LABELS:
    for ds in ds_in:
        sub = agg_best[(agg_best['dataset'] == ds) & (agg_best['run_label'] == run_lbl)]
        if len(sub):
            positions.append(xx); data.append(sub['neg_log10_p'].values)
            colors_list.append(DATASET_COLORS[ds]); labels.append((run_lbl, ds))
        xx += 1
    xx += 1
bp = ax.boxplot(data, positions=positions, widths=0.7, patch_artist=True, showfliers=False)
for patch, c in zip(bp['boxes'], colors_list):
    patch.set_facecolor(c); patch.set_alpha(0.7)
ax.axhline(-np.log10(0.05), color='red', linestyle='--', linewidth=1, alpha=0.7)
midpoints = []
for run_lbl in RUN_LABELS:
    starts = [j for j, lbl in enumerate(labels) if lbl[0] == run_lbl]
    if starts:
        midpoints.append((min(positions[s] for s in starts) + max(positions[s] for s in starts)) / 2)
ax.set_xticks(midpoints); ax.set_xticklabels(RUN_LABELS, rotation=20, ha='right', fontsize=9)
ax.set_ylabel('-log10(best p) per TF (5 TFs / box, best across 5 windows)')
ax.set_title('Distribution of TF enrichment significance (best across windows)', fontsize=11, fontweight='bold')
legend_handles = [Patch(facecolor=DATASET_COLORS[d], alpha=0.7, label=SHORT.get(d, d).replace('\n', ' ')) for d in ds_in]
legend_handles.append(plt.Line2D([0], [0], color='red', linestyle='--', label='p=0.05'))
ax.legend(handles=legend_handles, fontsize=8, loc='upper left', bbox_to_anchor=(1.01, 1.0))
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(OUT / '04_neg_log10_p_distribution_best_across_windows.pdf', dpi=300, bbox_inches='tight')
plt.savefig(OUT / '04_neg_log10_p_distribution_best_across_windows.png', dpi=300, bbox_inches='tight')
plt.close()
print('Saved: 04_neg_log10_p_distribution_best_across_windows.pdf / .png')


# === PER-WINDOW VIEW ===
# 05: # TFs enriched per (dataset × run) for EACH window separately
all_df['enriched_at_this_window'] = (all_df['pvalue'] < 0.05) & (all_df['odds_ratio'] > 1)
counts_w = all_df.groupby(['promoter_window_width', 'dataset', 'run_label']).agg(
    n=('enriched_at_this_window', 'sum')
).reset_index()

fig, axes = plt.subplots(1, len(PROMOTER_WINDOWS), figsize=(3.0 * len(PROMOTER_WINDOWS), 5), sharey=True)
for ax, w, lbl in zip(axes, PROMOTER_WINDOWS, WINDOW_LABELS):
    sub = counts_w[counts_w['promoter_window_width'] == w]
    piv = sub.pivot(index='dataset', columns='run_label', values='n').reindex(ds_in).reindex(columns=RUN_LABELS)
    im = ax.imshow(piv.values, aspect='auto', cmap='YlGnBu', vmin=0, vmax=5)
    ax.set_xticks(range(len(RUN_LABELS))); ax.set_xticklabels(RUN_LABELS, rotation=30, ha='right', fontsize=7)
    ax.set_yticks(range(len(ds_in))); ax.set_yticklabels([SHORT.get(d, d) for d in ds_in], fontsize=8)
    ax.set_title(f'window = {lbl}', fontsize=10, fontweight='bold')
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.iloc[i, j]
            if pd.notna(v):
                color = 'white' if v >= 3 else 'black'
                ax.text(j, i, f'{int(v)}', ha='center', va='center', fontsize=8, color=color)
            else:
                ax.text(j, i, '—', ha='center', va='center', fontsize=8, color='gray')
fig.suptitle('# TFs enriched (p<0.05 & OR>1) per (dataset × config) at each TSS window',
             fontsize=12, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(OUT / '05_n_enriched_TFs_per_window.pdf', dpi=300, bbox_inches='tight')
plt.savefig(OUT / '05_n_enriched_TFs_per_window.png', dpi=300, bbox_inches='tight')
plt.close()
print('Saved: 05_n_enriched_TFs_per_window.pdf / .png')


# 06: per-TF -log10(p) heatmap, faceted by window (one row of dataset panels per window)
fig, axes = plt.subplots(len(PROMOTER_WINDOWS), len(ds_in),
                         figsize=(2.8 * len(ds_in), 2.4 * len(PROMOTER_WINDOWS)),
                         sharey=True, sharex=True)
vmax6 = all_df['neg_log10_p'].quantile(0.95)
for row, (w, lbl) in enumerate(zip(PROMOTER_WINDOWS, WINDOW_LABELS)):
    sub_w = all_df[all_df['promoter_window_width'] == w]
    for col, ds in enumerate(ds_in):
        ax = axes[row, col]
        sub = sub_w[sub_w['dataset'] == ds].pivot_table(
            index='TF_display', columns='run_label', values='neg_log10_p', aggfunc='max'
        ).reindex(index=TFs, columns=RUN_LABELS)
        ax.imshow(sub.values, aspect='auto', cmap='RdPu', vmin=0, vmax=vmax6)
        ax.set_xticks(range(len(RUN_LABELS)))
        ax.set_xticklabels(RUN_LABELS if row == len(PROMOTER_WINDOWS) - 1 else [], rotation=45, ha='right', fontsize=6)
        ax.set_yticks(range(len(TFs)))
        ax.set_yticklabels(TFs if col == 0 else [], fontsize=8)
        if row == 0:
            ax.set_title(SHORT.get(ds, ds).replace('\n', ' '), fontsize=9, color=DATASET_COLORS[ds], fontweight='bold')
        if col == 0:
            ax.set_ylabel(lbl, fontsize=10, fontweight='bold')
        for i in range(sub.shape[0]):
            for j in range(sub.shape[1]):
                v = sub.iloc[i, j]
                if pd.notna(v):
                    color = 'white' if v > vmax6 * 0.6 else 'black'
                    ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=5.5, color=color)
                if pd.notna(v) and v >= -np.log10(0.05):
                    ax.add_patch(plt.Rectangle((j-0.45, i-0.45), 0.9, 0.9, fill=False, edgecolor='black', linewidth=0.8))
fig.suptitle('Per-TF -log10(p) faceted by TSS window (rows) × dataset (cols); boxed = p<0.05',
             fontsize=11, fontweight='bold', y=0.99)
plt.tight_layout()
plt.savefig(OUT / '06_per_TF_neg_log10_p_per_window.pdf', dpi=300, bbox_inches='tight')
plt.savefig(OUT / '06_per_TF_neg_log10_p_per_window.png', dpi=300, bbox_inches='tight')
plt.close()
print('Saved: 06_per_TF_neg_log10_p_per_window.pdf / .png')


# 07: per-TF ODDS RATIO heatmap, faceted by window (rows = window widths, cols = datasets)
fig, axes = plt.subplots(len(PROMOTER_WINDOWS), len(ds_in),
                         figsize=(2.8 * len(ds_in), 2.4 * len(PROMOTER_WINDOWS)),
                         sharey=True, sharex=True)
vmax7 = all_df['odds_ratio'].quantile(0.95)
for row, (w, lbl) in enumerate(zip(PROMOTER_WINDOWS, WINDOW_LABELS)):
    sub_w = all_df[all_df['promoter_window_width'] == w]
    for col, ds in enumerate(ds_in):
        ax = axes[row, col]
        sub = sub_w[sub_w['dataset'] == ds].pivot_table(
            index='TF_display', columns='run_label', values='odds_ratio', aggfunc='max'
        ).reindex(index=TFs, columns=RUN_LABELS)
        ax.imshow(sub.values, aspect='auto', cmap='Blues', vmin=0.5, vmax=vmax7)
        ax.set_xticks(range(len(RUN_LABELS)))
        ax.set_xticklabels(RUN_LABELS if row == len(PROMOTER_WINDOWS) - 1 else [], rotation=45, ha='right', fontsize=6)
        ax.set_yticks(range(len(TFs)))
        ax.set_yticklabels(TFs if col == 0 else [], fontsize=8)
        if row == 0:
            ax.set_title(SHORT.get(ds, ds).replace('\n', ' '), fontsize=9, color=DATASET_COLORS[ds], fontweight='bold')
        if col == 0:
            ax.set_ylabel(lbl, fontsize=10, fontweight='bold')
        for i in range(sub.shape[0]):
            for j in range(sub.shape[1]):
                v = sub.iloc[i, j]
                if pd.notna(v):
                    color = 'white' if v > vmax7 * 0.6 else 'black'
                    ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=5.5, color=color)
                if pd.notna(v) and v > 1:
                    ax.add_patch(plt.Rectangle((j-0.45, i-0.45), 0.9, 0.9, fill=False, edgecolor='black', linewidth=0.5))
fig.suptitle('Per-TF odds ratio faceted by TSS window (rows) × dataset (cols); boxed = OR>1',
             fontsize=11, fontweight='bold', y=0.99)
plt.tight_layout()
plt.savefig(OUT / '07_per_TF_odds_ratio_per_window.pdf', dpi=300, bbox_inches='tight')
plt.savefig(OUT / '07_per_TF_odds_ratio_per_window.png', dpi=300, bbox_inches='tight')
plt.close()
print('Saved: 07_per_TF_odds_ratio_per_window.pdf / .png')


# 08: COMBINED DOT PLOT — color = log2(OR), size = -log10(p)
# Faceted: rows = windows, cols = datasets. x = run, y = TF.
def or_to_color(or_vals, vmin=-1.5, vmax=1.5):
    """Map log2(OR) to a diverging colormap (blue=depleted, red=enriched)."""
    return np.log2(or_vals.clip(lower=0.01))

fig, axes = plt.subplots(len(PROMOTER_WINDOWS), len(ds_in),
                         figsize=(3.0 * len(ds_in), 2.6 * len(PROMOTER_WINDOWS)),
                         sharey=False, sharex=True)
# Common scaling
all_log2_or = np.log2(all_df['odds_ratio'].clip(lower=0.01))
all_neg_log_p = all_df['neg_log10_p']
or_vmin, or_vmax = -1.5, 1.5  # symmetric diverging
size_min, size_max = all_neg_log_p.quantile(0.05), all_neg_log_p.quantile(0.97)
import matplotlib.colors as mcolors
norm = mcolors.Normalize(vmin=or_vmin, vmax=or_vmax)
cmap = plt.cm.RdBu_r

def neg_log_p_to_size(values):
    """Map -log10(p) to dot size (scatter `s` in points^2)."""
    s_min, s_max = 8, 200
    v = np.clip(values, 0, size_max)
    return s_min + (s_max - s_min) * (v / size_max if size_max > 0 else 0)

x_pos = {lbl: i for i, lbl in enumerate(RUN_LABELS)}
y_pos = {tf: i for i, tf in enumerate(TFs)}

for row, (w, lbl) in enumerate(zip(PROMOTER_WINDOWS, WINDOW_LABELS)):
    sub_w = all_df[all_df['promoter_window_width'] == w]
    for col, ds in enumerate(ds_in):
        ax = axes[row, col]
        sub = sub_w[sub_w['dataset'] == ds]
        for _, r in sub.iterrows():
            x = x_pos.get(r['run_label'])
            y = y_pos.get(r['TF_display'])
            if x is None or y is None:
                continue
            log2_or = np.log2(max(r['odds_ratio'], 0.01))
            size = neg_log_p_to_size(r['neg_log10_p'])
            edge = 'black' if r['pvalue'] < 0.05 else 'none'
            ax.scatter(x, y, s=size, c=[log2_or], cmap=cmap, norm=norm,
                       edgecolors=edge, linewidths=1.0, alpha=0.95)
        ax.set_xlim(-0.5, len(RUN_LABELS) - 0.5)
        # Invert before setting ticks so the labels stick.
        ax.set_ylim(len(TFs) - 0.5, -0.5)  # inverted: top-down TFs
        ax.set_xticks(range(len(RUN_LABELS)))
        ax.set_xticklabels(RUN_LABELS if row == len(PROMOTER_WINDOWS) - 1 else [], rotation=45, ha='right', fontsize=6)
        ax.set_yticks(range(len(TFs)))
        ax.set_yticklabels(TFs if col == 0 else [''] * len(TFs), fontsize=8)
        ax.grid(True, linestyle=':', alpha=0.3)
        if row == 0:
            ax.set_title(SHORT.get(ds, ds).replace('\n', ' '), fontsize=9, color=DATASET_COLORS[ds], fontweight='bold')
        if col == 0:
            ax.set_ylabel(lbl, fontsize=10, fontweight='bold')

# Colorbar (log2 OR)
cbar_ax = fig.add_axes([1.005, 0.5, 0.012, 0.3])
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
cb = fig.colorbar(sm, cax=cbar_ax)
cb.set_label('log2(odds ratio)', fontsize=8)
cb.ax.tick_params(labelsize=7)

# Size legend
size_legend_handles = []
for v in [1, 2, 4, 6]:
    s = neg_log_p_to_size(np.array([v]))[0]
    size_legend_handles.append(plt.scatter([], [], s=s, c='gray', alpha=0.7, label=f'-log10(p) = {v}'))
fig.legend(handles=size_legend_handles, loc='center right', bbox_to_anchor=(1.10, 0.20),
           fontsize=8, frameon=True, title='Dot size', title_fontsize=8)

# Significance edge legend
fig.text(1.005, 0.05, "Black ring = p<0.05", fontsize=8)

fig.suptitle('TF enrichment dotplot: color = log2(OR), size = -log10(p), ring = p<0.05',
             fontsize=11, fontweight='bold', y=1.00)
plt.tight_layout()
plt.savefig(OUT / '08_dotplot_per_window.pdf', dpi=300, bbox_inches='tight')
plt.savefig(OUT / '08_dotplot_per_window.png', dpi=300, bbox_inches='tight')
plt.close()
print('Saved: 08_dotplot_per_window.pdf / .png')


# Save tables
agg_best.to_csv(OUT / 'tf_enrichment_per_TF_best_across_windows.tsv', sep='\t', index=False)
all_df.to_csv(OUT / 'tf_enrichment_all_windows.tsv', sep='\t', index=False)
counts_w.to_csv(OUT / 'tf_enrichment_n_enriched_per_window.tsv', sep='\t', index=False)
print('\nSaved tables: tf_enrichment_per_TF_best_across_windows.tsv, tf_enrichment_all_windows.tsv, tf_enrichment_n_enriched_per_window.tsv')
print(f'\nDone. Outputs in: {OUT}')
