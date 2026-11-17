"""
nt_deg_counts.py
────────────────
Bar chart of total DEG counts per genotype × drug.
Two panels (HQ461, ZLY-025); bars ordered most resistant → sensitive.
DEG threshold: padj < 0.05, |LFC| > 1

Output: figures/04_gsea_cdk12/NT_DEG_counts.{pdf,png,svg}
"""

import os, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family']     = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'Arial', 'DejaVu Sans']
plt.rcParams['pdf.fonttype']    = 42
plt.rcParams['ps.fonttype']     = 42

_ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(_ROOT, 'pc3_rnaseq_results', 'deseq2_results', 'drug_response')
FIGDIR = os.path.join(_ROOT, 'pc3_rnaseq_results', 'figures', '04_gsea_cdk12')
os.makedirs(FIGDIR, exist_ok=True)

COMPARISONS = {
    'HQ461':   [('NT',    'NT_HQ461_vs_NT_Control'),
                ('G731K', 'G731K_HQ461_vs_G731K_Control'),
                ('I733V', 'I733V_HQ461_vs_I733V_Control'),
                ('G822E', 'G822E_HQ461_vs_G822E_Control')],
    'ZLY-025': [('NT',    'NT_ZLY025_vs_NT_Control'),
                ('G731K', 'G731K_ZLY025_vs_G731K_Control'),
                ('I733V', 'I733V_ZLY025_vs_I733V_Control'),
                ('G822E', 'G822E_ZLY025_vs_G822E_Control')],
}

def count_degs(fname):
    df = pd.read_csv(os.path.join(RESDIR, f'{fname}.csv'), index_col=0)
    df = df.dropna(subset=['log2FoldChange', 'padj'])
    return ((df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 1)).sum()

C_NT  = '#888888'
COLORS = {'NT': C_NT, 'G731K': '#2166AC', 'G822E': '#4DAC26', 'I733V': '#D01C8B'}

fig, axes = plt.subplots(1, 2, figsize=(7, 4), sharey=False)
fig.subplots_adjust(wspace=0.35)

for ax, (drug, comps) in zip(axes, COMPARISONS.items()):
    labels = [g for g, _ in comps]
    counts = [count_degs(f) for _, f in comps]
    colors = [COLORS[g] for g in labels]
    x = np.arange(len(labels))
    bars = ax.bar(x, counts, color=colors, width=0.6, edgecolor='white', linewidth=0.5)
    for bar, n in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 40,
                f'{n:,}', ha='center', va='bottom', fontsize=8, color='#333333')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('DEGs (padj < 0.05, |LFC| > 1)', fontsize=8.5)
    ax.set_title(drug, fontsize=10, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.set_ylim(0, max(counts) * 1.18)

fig.suptitle('DEG counts per genotype — drug response\nPC3 CDK12 variants',
             fontsize=10, y=1.02)

out = os.path.join(FIGDIR, 'NT_DEG_counts')
plt.savefig(f'{out}.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{out}.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{out}.svg', dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {out}.{{pdf,png,svg}}")
