"""
UpSet plot — CDK12 variant DE gene overlap.
Layout follows upset.app convention; aesthetic follows Sam Gould / lab style.

Sets:  G822E (confirmed hypomorph), G731K (confirmed hypomorph), I733V (unconfirmed)
Input: existing DESeq2 CSVs in pc3_rnaseq_results/deseq2_results/
"""

import os
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

# ── Style ──────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':      'sans-serif',
    'font.sans-serif':  ['Liberation Sans', 'Arial', 'Helvetica Neue',
                         'Helvetica', 'DejaVu Sans'],
    'pdf.fonttype':     42,
    'ps.fonttype':      42,
    'axes.linewidth':   0.8,
    'xtick.major.width': 0.7,
    'ytick.major.width': 0.7,
})

# ── Paths ──────────────────────────────────────────────────────────────────────
_HERE  = os.path.dirname(os.path.abspath(__file__))
_ROOT  = os.path.dirname(_HERE)
RESDIR = os.path.join(_ROOT, 'pc3_rnaseq_results', 'deseq2_results')
OUTDIR = os.path.join(_ROOT, 'pc3_rnaseq_results', 'figures', '03_hypomorphic')
os.makedirs(OUTDIR, exist_ok=True)

FDR = 0.05

# ── Load DE gene sets ──────────────────────────────────────────────────────────
def _sig_genes(fname):
    df = pd.read_csv(os.path.join(RESDIR, fname))
    df = df.dropna(subset=['padj', 'gene_name']).drop_duplicates('gene_name')
    return set(df.loc[df['padj'] < FDR, 'gene_name'])

SETS = ['G822E', 'G731K', 'I733V']
SET_COLORS = {
    'G822E': '#c0392b',
    'G731K': '#8e0152',
    'I733V': '#2980b9',
}

sig = {
    'G822E': _sig_genes('NT_vs_G822E_Control.csv'),
    'G731K': _sig_genes('NT_vs_G731K_Control.csv'),
    'I733V': _sig_genes('NT_vs_I733V_Control.csv'),
}

# ── Compute exclusive intersections ───────────────────────────────────────────
intersections = []
for r in range(1, len(SETS) + 1):
    for combo in combinations(SETS, r):
        included = frozenset(combo)
        excluded = frozenset(SETS) - included
        members  = set.intersection(*[sig[s] for s in included])
        for s in excluded:
            members -= sig[s]
        if members:
            intersections.append({'sets': included, 'count': len(members)})

intersections.sort(key=lambda x: -x['count'])

for d in intersections:
    print(f"  {'+'.join(sorted(d['sets'])): <30}  n = {d['count']:,}")

# ── Bar colors ─────────────────────────────────────────────────────────────────
# Hypomorph-only → red; I733V-only → blue; mixed (involves I733V) → gray
def _bar_color(included):
    if 'I733V' not in included:
        return '#c0392b'
    if included == frozenset({'I733V'}):
        return '#2980b9'
    return '#95a5a6'

# ── Layout ────────────────────────────────────────────────────────────────────
N_INTER = len(intersections)
N_SETS  = len(SETS)

FIG_W = max(9, N_INTER * 1.15 + 2.8)
FIG_H = 5.8

fig = plt.figure(figsize=(FIG_W, FIG_H))

# Outer GridSpec: 2 rows (bars | matrix), 2 cols (set-sizes | main)
outer = gridspec.GridSpec(
    2, 2,
    height_ratios=[3.0, 1.8],
    width_ratios=[1.0, N_INTER + 0.5],
    hspace=0.0,
    wspace=0.06,
    left=0.05, right=0.97,
    top=0.87, bottom=0.17,
)

ax_corner = fig.add_subplot(outer[0, 0])   # top-left: hidden
ax_bars   = fig.add_subplot(outer[0, 1])   # top-right: intersection size bars
ax_sets   = fig.add_subplot(outer[1, 0])   # bottom-left: set size bars
ax_mat    = fig.add_subplot(outer[1, 1])   # bottom-right: dot matrix

ax_corner.set_visible(False)

# ── Top-right: intersection size bars ─────────────────────────────────────────
x     = np.arange(N_INTER)
cnts  = [d['count'] for d in intersections]
bcols = [_bar_color(d['sets']) for d in intersections]

bars = ax_bars.bar(x, cnts, color=bcols, edgecolor='black',
                   linewidth=0.7, width=0.55, zorder=3)

# Count labels above bars
for xi, n in zip(x, cnts):
    ax_bars.text(xi, n + max(cnts) * 0.025, f'{n:,}',
                 ha='center', va='bottom', fontsize=8.5, color='#111111')

ax_bars.set_xlim(-0.65, N_INTER - 0.35)
ax_bars.set_ylim(0, max(cnts) * 1.22)
ax_bars.set_xticks([])
ax_bars.set_ylabel('Intersection size', fontsize=10, labelpad=6)
ax_bars.spines[['top', 'right', 'bottom']].set_visible(False)
ax_bars.tick_params(axis='y', labelsize=9, length=3, width=0.7)
ax_bars.yaxis.grid(True, lw=0.35, color='#DDDDDD', zorder=0)
ax_bars.set_axisbelow(True)

# Legend (top right of bar panel)
_handles = [
    mpatches.Patch(facecolor='#c0392b', edgecolor='black', lw=0.7,
                   label='G822E / G731K (confirmed hypomorphs)'),
    mpatches.Patch(facecolor='#2980b9', edgecolor='black', lw=0.7,
                   label='I733V only (unconfirmed)'),
    mpatches.Patch(facecolor='#95a5a6', edgecolor='black', lw=0.7,
                   label='Mixed (includes I733V)'),
]
ax_bars.legend(handles=_handles, frameon=False, fontsize=7.5,
               loc='upper right', handlelength=1.0, handletextpad=0.5,
               borderpad=0, labelspacing=0.35)

# ── Bottom-right: dot matrix ───────────────────────────────────────────────────
# Row order: G822E at top, I733V at bottom  (y=2, 1, 0)
set_y = {s: i for i, s in enumerate(reversed(SETS))}
# → G822E:2, G731K:1, I733V:0

DOT_ON  = 130
DOT_OFF = 55
CONN_LW = 3.0

for xi, inter in enumerate(intersections):
    included = inter['sets']
    # Background (excluded) dots first
    for s, yi in set_y.items():
        if s not in included:
            ax_mat.scatter(xi, yi, s=DOT_OFF, color='#DEDEDE',
                           zorder=3, linewidths=0)
    # Connecting line between included dots
    ys_on = sorted(set_y[s] for s in included)
    if len(ys_on) > 1:
        ax_mat.plot([xi, xi], [ys_on[0], ys_on[-1]],
                    color='#444444', lw=CONN_LW, zorder=4, solid_capstyle='round')
    # Filled dots on top
    for s in included:
        ax_mat.scatter(xi, set_y[s], s=DOT_ON,
                       color=SET_COLORS[s], zorder=5,
                       linewidths=0.6, edgecolors='white')

ax_mat.set_xlim(-0.65, N_INTER - 0.35)
ax_mat.set_ylim(-0.65, N_SETS - 0.35)
ax_mat.set_xticks([])

# Set name labels as y-tick labels, colored by variant
ax_mat.set_yticks(list(set_y.values()))
ax_mat.set_yticklabels(list(reversed(SETS)), fontsize=11, fontweight='bold')
for lbl, s in zip(ax_mat.get_yticklabels(), reversed(SETS)):
    lbl.set_color(SET_COLORS[s])
ax_mat.tick_params(axis='y', length=0, pad=4)
ax_mat.spines[:].set_visible(False)

# Subtle row separator lines
for yi in set_y.values():
    ax_mat.axhline(yi, color='#F0F0F0', lw=0.6, zorder=1)

# ── Bottom-left: set size bars (point LEFT = reversed x) ──────────────────────
set_sizes = {s: len(sig[s]) for s in SETS}
max_set   = max(set_sizes.values())

for s in SETS:
    yi = set_y[s]
    n  = set_sizes[s]
    ax_sets.barh(yi, n, color=SET_COLORS[s], edgecolor='black',
                 linewidth=0.7, height=0.45, zorder=3)
    # Label to the left of bar in the inverted axis
    ax_sets.text(n + max_set * 0.06, yi, f'{n:,}',
                 va='center', ha='left', fontsize=8, color='#333333')

ax_sets.set_xlim(0, max_set * 1.35)
ax_sets.invert_xaxis()    # bars now point left toward the matrix
ax_sets.set_ylim(-0.65, N_SETS - 0.35)
ax_sets.set_yticks([])
ax_sets.set_xlabel('Set size', fontsize=9, labelpad=4)
ax_sets.spines[['top', 'left', 'right']].set_visible(False)
ax_sets.tick_params(axis='x', labelsize=8, length=3, width=0.7)
ax_sets.xaxis.grid(True, lw=0.35, color='#DDDDDD', zorder=0)
ax_sets.set_axisbelow(True)

# Subtle row separators (same y as matrix)
for yi in set_y.values():
    ax_sets.axhline(yi, color='#F0F0F0', lw=0.6, zorder=1)

# ── Title ─────────────────────────────────────────────────────────────────────
fig.suptitle(
    'CDK12 Variant DE Gene Overlap — Control Condition (FDR < 0.05)',
    fontsize=12, fontweight='bold', x=0.55, y=0.95
)

# ── Save ──────────────────────────────────────────────────────────────────────
_out = os.path.join(OUTDIR, 'de_gene_upset')
plt.savefig(f'{_out}.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{_out}.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"\nSaved: {_out}.pdf / .png")
