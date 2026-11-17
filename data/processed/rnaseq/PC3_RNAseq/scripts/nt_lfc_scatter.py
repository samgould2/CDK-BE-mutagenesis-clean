"""
nt_lfc_scatter.py
-----------------
Genome-wide LFC scatter plots: variant drug response vs NT drug response.
3 rows (variants) x 2 cols (drugs), ordered by resistance phenotype.

Points:
  - All shared genes: grey
  - NT DEGs (padj<0.05, |LFC|>1): colored by direction (up=coral, dn=blue)
  - NT DEGs recovered in variant (same direction, padj<0.05 |LFC|>1): filled
  - NT DEGs lost in variant: open circles

Outputs
-------
    figures/04_gsea_cdk12/NT_LFC_scatter.{pdf,png,svg}

Run from the project root
    cd ~/Desktop/"Comprehensive PC3 RNA-seq Analysis"
    python3 scripts/nt_lfc_scatter.py
"""

import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

plt.rcParams['font.family']     = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'Arial', 'Helvetica Neue',
                                    'Helvetica', 'DejaVu Sans', 'sans-serif']
plt.rcParams['pdf.fonttype']    = 42
plt.rcParams['ps.fonttype']     = 42

_ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(_ROOT, 'pc3_rnaseq_results', 'deseq2_results')
FIGDIR = os.path.join(_ROOT, 'pc3_rnaseq_results', 'figures', '04_gsea_cdk12')
os.makedirs(FIGDIR, exist_ok=True)

def load_csv(fname):
    path = os.path.join(RESDIR, f'{fname}.csv')
    df = pd.read_csv(path, index_col=0)
    df = df.dropna(subset=['log2FoldChange', 'padj', 'gene_name'])
    df = df.drop_duplicates(subset='gene_name').set_index('gene_name')
    return df

print("Loading DESeq2 results...")
nt_hq461 = load_csv('NT_HQ461_vs_NT_Control')
nt_zly   = load_csv('NT_ZLY025_vs_NT_Control')

def de_set(df, direction):
    mask = (df['padj'] < 0.05) & (
        (df['log2FoldChange'] > 1) if direction == 'up' else (df['log2FoldChange'] < -1)
    )
    return set(df[mask].index)

nt_up = {'HQ461': de_set(nt_hq461, 'up'),   'ZLY-025': de_set(nt_zly, 'up')}
nt_dn = {'HQ461': de_set(nt_hq461, 'down'), 'ZLY-025': de_set(nt_zly, 'down')}

# Ordered: most resistant top row, most sensitive bottom row
COMPARISONS = [
    ('G731K_HQ461_vs_G731K_Control',  'G731K', 'HQ461',   nt_hq461),
    ('I733V_HQ461_vs_I733V_Control',  'I733V', 'HQ461',   nt_hq461),
    ('G822E_HQ461_vs_G822E_Control',  'G822E', 'HQ461',   nt_hq461),
    ('G731K_ZLY025_vs_G731K_Control', 'G731K', 'ZLY-025', nt_zly),
    ('G822E_ZLY025_vs_G822E_Control', 'G822E', 'ZLY-025', nt_zly),
    ('I733V_ZLY025_vs_I733V_Control', 'I733V', 'ZLY-025', nt_zly),
]

C_GREY = '#CCCCCC'
C_UP   = '#D6604D'
C_DN   = '#4393C3'

fig, axes = plt.subplots(3, 2, figsize=(7.5, 10.0))
fig.subplots_adjust(hspace=0.42, wspace=0.38)

for idx, (fname, variant, drug, ref_df) in enumerate(COMPARISONS):
    row_i = idx % 3
    col_i = idx // 3
    ax    = axes[row_i, col_i]

    var_df  = load_csv(fname)
    shared  = ref_df.index.intersection(var_df.index)
    x = ref_df.loc[shared, 'log2FoldChange'].dropna()
    y = var_df.loc[shared, 'log2FoldChange'].dropna()
    shared2 = x.index.intersection(y.index)
    x, y    = x.loc[shared2], y.loc[shared2]

    r, pval = pearsonr(x, y)

    up_genes     = nt_up[drug] & set(shared2)
    dn_genes     = nt_dn[drug] & set(shared2)
    var_up_set   = de_set(var_df, 'up')
    var_dn_set   = de_set(var_df, 'down')
    recovered_up = list(up_genes & var_up_set)
    recovered_dn = list(dn_genes & var_dn_set)
    lost_up      = list(up_genes - var_up_set)
    lost_dn      = list(dn_genes - var_dn_set)
    background   = list(set(shared2) - up_genes - dn_genes)

    total_nt  = len(up_genes) + len(dn_genes)
    recovered = len(recovered_up) + len(recovered_dn)
    recov_pct = 100 * recovered / total_nt if total_nt > 0 else 0

    # Background
    ax.scatter(x.loc[background], y.loc[background],
               s=3, color=C_GREY, alpha=0.3, linewidths=0, rasterized=True, zorder=1)
    # Lost NT DEGs — same color as background, no highlight
    for g, c in [(lost_up, C_UP), (lost_dn, C_DN)]:
        ax.scatter(x.loc[g], y.loc[g], s=16,
                   facecolors='none', edgecolors=c, linewidths=0.6, alpha=0.55, zorder=2)
    # Recovered NT DEGs — filled with thick black outline
    for g, c in [(recovered_up, C_UP), (recovered_dn, C_DN)]:
        ax.scatter(x.loc[g], y.loc[g], s=22,
                   color=c, alpha=0.9, edgecolors='black', linewidths=0.7, zorder=3)

    lim = max(abs(x).max(), abs(y).max()) * 1.06
    ax.axhline(0, color='#aaaaaa', lw=0.4, zorder=0)
    ax.axvline(0, color='#aaaaaa', lw=0.4, zorder=0)
    ax.plot([-lim, lim], [-lim, lim], color='#555555', lw=0.8, linestyle='--', zorder=0)

    pstr = f'p = {pval:.1e}' if pval >= 1e-300 else 'p < 1e-300'
    ax.text(0.04, 0.97,
            f'r = {r:.3f},  {pstr}\n{recov_pct:.0f}% of NT DEGs recovered',
            transform=ax.transAxes, va='top', ha='left', fontsize=7.5,
            color='#222222',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='none', alpha=0.85))

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(f'NT {drug}  log₂FC', fontsize=8.5)
    ax.set_ylabel(f'{variant} {drug}  log₂FC', fontsize=8.5)
    ax.set_title(f'{variant} — {drug}', fontsize=9.5, fontweight='bold', pad=5)
    ax.tick_params(labelsize=7.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    print(f"  {variant} {drug:7s}  r={r:.3f}  recovered={recov_pct:.0f}%")

legend_elements = [
    mlines.Line2D([0],[0], marker='o', color='w', markerfacecolor=C_UP,
                  markeredgecolor='black', markeredgewidth=0.7,
                  markersize=7, label='NT UP DEG — recovered in variant'),
    mlines.Line2D([0],[0], marker='o', color='w', markerfacecolor=C_DN,
                  markeredgecolor='black', markeredgewidth=0.7,
                  markersize=7, label='NT DN DEG — recovered in variant'),
    mlines.Line2D([0],[0], marker='o', color='w', markerfacecolor=C_GREY,
                  markersize=6, label='All other genes'),
    mlines.Line2D([0],[0], color='#555555', lw=0.8, linestyle='--',
                  label='y = x  (perfect preservation)'),
]
fig.legend(handles=legend_elements, loc='lower center', ncol=3, fontsize=7.5,
           bbox_to_anchor=(0.5, -0.03), frameon=False)

fig.suptitle('Transcriptional Drug Response Preservation — PC3 CDK12 variants\n'
             'Genome-wide LFC correlation; NT DEGs (padj<0.05, |LFC|>1) highlighted',
             fontsize=10, y=1.01)

outpath = os.path.join(FIGDIR, 'NT_LFC_scatter')
plt.savefig(f'{outpath}.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{outpath}.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{outpath}.svg', dpi=300, bbox_inches='tight')
plt.close()
print(f"\nSaved: {outpath}.{{pdf,png,svg}}")
