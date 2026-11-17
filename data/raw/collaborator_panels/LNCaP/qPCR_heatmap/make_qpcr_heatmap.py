#!/usr/bin/env python3
"""
LNCaP qPCR time-course heatmap: ATR, BRCA1, BRCA2, ERCC4
Timepoints: 0, 4, 8, 12 hr after CDK12 inhibition
Fold change vs 0hr; stats = t-test vs 0hr, Bonferroni-corrected (3 comparisons)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import seaborn as sns
from scipy.stats import ttest_ind
import pickle
import os, sys

# Usage: python make_qpcr_heatmap.py [qPCR_DNA_repair_analysis_dir] [outdir]
# Defaults: the folder above this script, and <that folder>/outputs/figures
BASE = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else os.path.join(BASE, 'outputs', 'figures')
os.makedirs(OUTDIR, exist_ok=True)
INTER = os.path.join(BASE, 'scripts', '_intermediate')

# ---------- font ----------
try:
    prop = fm.FontProperties(family='Liberation Sans')
    matplotlib.rcParams['font.family'] = 'Liberation Sans'
except:
    matplotlib.rcParams['font.family'] = 'DejaVu Sans'

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

# ---------- load data ----------
with open(os.path.join(INTER, 'fc_summary.pkl'), 'rb') as f:
    fc_summary = pickle.load(f)

with open(os.path.join(INTER, 'pivot_results.pkl'), 'rb') as f:
    piv = pickle.load(f)

genes    = ['hATR', 'hBRCA1', 'hBRCA2', 'hERCC4']
timepoints     = [0, 4, 8, 12]   # all timepoints (needed for stats vs 0hr)
plot_timepoints = [0, 4, 8, 12]   # all timepoints
labels_genes = ['ATR', 'BRCA1', 'BRCA2', 'ERCC4']
labels_tp    = ['0 hr', '4 hr', '8 hr', '12 hr']
labels_tp_plot = ['0 hr', '4 hr', '8 hr', '12 hr']

# ---------- build FC matrix & SEM matrix ----------
fc_mat  = np.zeros((len(genes), len(timepoints)))
sem_mat = np.zeros((len(genes), len(timepoints)))

for i, gene in enumerate(genes):
    for j, tp in enumerate(timepoints):
        row = fc_summary[(fc_summary['Gene'] == gene) & (fc_summary['Timepoint'] == str(tp))]
        fc_mat[i, j]  = row['Mean_FC'].values[0]
        sem_mat[i, j] = row['SEM_FC'].values[0]

# ---------- statistics: t-test each tp vs 0hr on dCq, Bonferroni x3 ----------
def stars(p):
    if   p < 0.0001: return '****'
    elif p < 0.001:  return '***'
    elif p < 0.01:   return '**'
    elif p < 0.05:   return '*'
    else:            return 'ns'

n_comparisons = 3  # 4hr, 8hr, 12hr vs 0hr

pval_mat = np.ones((len(genes), len(timepoints)))
star_mat = [['' for _ in timepoints] for _ in genes]

for i, gene in enumerate(genes):
    dCq_col = f'dCq_{gene}'
    ref = piv[piv['timepoint_num'] == 0][dCq_col].dropna().values
    for j, tp in enumerate(timepoints):
        if tp == 0:
            star_mat[i][j] = ''
            continue
        grp = piv[piv['timepoint_num'] == tp][dCq_col].dropna().values
        _, p = ttest_ind(ref, grp)
        p_corr = min(p * n_comparisons, 1.0)
        pval_mat[i, j] = p_corr
        star_mat[i][j] = stars(p_corr)

# ---------- plot ----------
# ---------- subset to plot timepoints only (drop 0hr reference row) ----------
plot_idx   = [timepoints.index(tp) for tp in plot_timepoints]
fc_plot    = fc_mat[:, plot_idx]   # shape (genes, 3)
sem_plot   = sem_mat[:, plot_idx]
star_plot  = [[star_mat[g][t] for t in plot_idx] for g in range(len(genes))]

# genes as rows, timepoints as columns
fig, ax = plt.subplots(figsize=(2.8, 3.2))

hm = sns.heatmap(
    fc_plot,
    ax=ax,
    cmap='RdBu_r',
    center=1.0,
    vmin=0.0,
    vmax=2.0,
    annot=False,
    linewidths=0.8,
    linecolor='black',
    cbar=True,
    cbar_kws={'shrink': 0.7, 'pad': 0.02},
    xticklabels=labels_tp_plot,
    yticklabels=labels_genes,
)

# ---------- annotate cells: mean ± SEM and significant stars only ----------
for i in range(len(genes)):           # rows
    for j in range(len(plot_timepoints)):  # cols
        mean = fc_plot[i, j]
        sem  = sem_plot[i, j]
        sig  = star_plot[i][j]

        txt_color = 'white' if (mean < 0.55 or mean > 1.6) else 'black'

        cell_text = f'{mean:.2f}\n±{sem:.2f}'
        if sig and sig != 'ns':
            cell_text += f'\n{sig}'

        ax.text(j + 0.5, i + 0.5, cell_text,
                ha='center', va='center',
                fontsize=6.5, color=txt_color,
                fontproperties=fm.FontProperties(family='Liberation Sans'))

# ---------- axes styling ----------
ax.set_xticklabels(labels_tp_plot, fontsize=8, rotation=90)
ax.set_yticklabels(labels_genes, fontsize=8, rotation=90, ha='center', va='bottom')
ax.tick_params(axis='both', length=0)
ax.set_xlabel('')
ax.set_ylabel('')

# colorbar label
cbar = ax.collections[0].colorbar
cbar.set_label('Fold Change', fontsize=7)
cbar.ax.tick_params(labelsize=7)

plt.tight_layout(pad=0.4)

outbase = os.path.join(OUTDIR, 'LNCaP_CDK12_qPCR_DNArepair_heatmap')
fig.savefig(outbase + '.pdf', bbox_inches='tight', dpi=300)
fig.savefig(outbase + '.svg', bbox_inches='tight', dpi=300)
fig.savefig(outbase + '.png', bbox_inches='tight', dpi=300)
print(f"Saved {outbase}.[pdf/svg/png]")
print()
print("FC matrix:")
print(pd.DataFrame(fc_mat, index=labels_genes, columns=labels_tp).round(3))
print()
print("Stars:")
for i, g in enumerate(labels_genes):
    print(f"  {g}: {star_mat[i]}")
