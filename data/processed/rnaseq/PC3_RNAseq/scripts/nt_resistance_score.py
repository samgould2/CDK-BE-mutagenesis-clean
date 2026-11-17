"""
nt_resistance_score.py
──────────────────────
Computes a genome-wide transcriptional response preservation score for
each CDK12 variant under each drug, defined as the Pearson r between
the variant's LFC vector and the NT LFC vector for the same drug.

    r = 1.0  →  responds exactly like NT (fully sensitive)
    r = 0    →  no response preserved (fully resistant)

Also computes a concordant DE gene recovery rate as a complementary
metric (% of NT drug response DEGs that are also DE in the same
direction in the variant; padj < 0.05, |LFC| > 1 in both).

Outputs
-------
    figures/04_gsea_cdk12/NT_resistance_score.{pdf,png,svg}

Run from the project root
-------------------------
    cd ~/Desktop/"Comprehensive PC3 RNA-seq Analysis"
    python3 scripts/nt_resistance_score.py
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ── Style ─────────────────────────────────────────────────────────────────────
plt.rcParams['font.family']      = 'sans-serif'
plt.rcParams['font.sans-serif']  = ['Liberation Sans', 'Arial', 'Helvetica Neue',
                                     'Helvetica', 'DejaVu Sans', 'sans-serif']
plt.rcParams['pdf.fonttype']     = 42
plt.rcParams['ps.fonttype']      = 42

# ── Paths ─────────────────────────────────────────────────────────────────────
_ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(_ROOT, 'pc3_rnaseq_results', 'deseq2_results')
FIGDIR = os.path.join(_ROOT, 'pc3_rnaseq_results', 'figures', '04_gsea_cdk12')
os.makedirs(FIGDIR, exist_ok=True)

# ── Load DESeq2 CSVs ──────────────────────────────────────────────────────────
def load_csv(fname):
    path = os.path.join(RESDIR, f'{fname}.csv')
    df = pd.read_csv(path, index_col=0)
    df = df.dropna(subset=['log2FoldChange', 'gene_name'])
    df = df.drop_duplicates(subset='gene_name').set_index('gene_name')
    return df

print("Loading DESeq2 results...")
nt_hq461 = load_csv('NT_HQ461_vs_NT_Control')
nt_zly   = load_csv('NT_ZLY025_vs_NT_Control')
print(f"  NT HQ461 response: {len(nt_hq461)} genes")
print(f"  NT ZLY-025 response: {len(nt_zly)} genes")

# ── Variant comparisons to score ──────────────────────────────────────────────
# (csv_fname, display_label, variant, drug, NT_reference_df)
COMPARISONS = [
    ('I733V_HQ461_vs_I733V_Control', 'I733V', 'HQ461',   nt_hq461),
    ('G822E_HQ461_vs_G822E_Control', 'G822E', 'HQ461',   nt_hq461),
    ('G731K_HQ461_vs_G731K_Control', 'G731K', 'HQ461',   nt_hq461),
    ('I733V_ZLY025_vs_I733V_Control','I733V', 'ZLY-025', nt_zly),
    ('G822E_ZLY025_vs_G822E_Control','G822E', 'ZLY-025', nt_zly),
    ('G731K_ZLY025_vs_G731K_Control','G731K', 'ZLY-025', nt_zly),
]

# ── NT DE gene sets (for recovery rate) ───────────────────────────────────────
def de_genes(df, direction='down'):
    mask = (df['padj'] < 0.05) & (
        (df['log2FoldChange'] < -1) if direction == 'down' else (df['log2FoldChange'] > 1)
    )
    return set(df[mask].index)

nt_hq461_dn = de_genes(nt_hq461, 'down');  nt_hq461_up = de_genes(nt_hq461, 'up')
nt_zly_dn   = de_genes(nt_zly,   'down');  nt_zly_up   = de_genes(nt_zly,   'up')
print(f"\n  NT HQ461 DEGs:  {len(nt_hq461_dn)} down, {len(nt_hq461_up)} up")
print(f"  NT ZLY-025 DEGs:{len(nt_zly_dn)} down, {len(nt_zly_up)} up")

# ── Score each variant ────────────────────────────────────────────────────────
print("\nComputing resistance scores...")
records = []
for fname, variant, drug, ref_df in COMPARISONS:
    var_df  = load_csv(fname)
    ref_lfc = ref_df['log2FoldChange']
    var_lfc = var_df['log2FoldChange']

    # Genome-wide Pearson r on shared genes
    shared = ref_lfc.index.intersection(var_lfc.index)
    x, y   = ref_lfc.loc[shared], var_lfc.loc[shared]
    shared2 = x.dropna().index.intersection(y.dropna().index)
    r, pval = pearsonr(x.loc[shared2], y.loc[shared2])

    # DE gene recovery rate
    if drug == 'HQ461':
        nt_dn, nt_up = nt_hq461_dn, nt_hq461_up
    else:
        nt_dn, nt_up = nt_zly_dn, nt_zly_up

    var_dn = de_genes(var_df, 'down')
    var_up = de_genes(var_df, 'up')
    total_nt_de = len(nt_dn) + len(nt_up)
    recovered   = len(nt_dn & var_dn) + len(nt_up & var_up)
    recovery_pct = 100 * recovered / total_nt_de if total_nt_de > 0 else 0

    records.append({
        'variant':      variant,
        'drug':         drug,
        'r':            r,
        'pval':         pval,
        'n_genes':      len(shared2),
        'recovery_pct': recovery_pct,
        'recovered_n':  recovered,
        'nt_de_n':      total_nt_de,
    })
    print(f"  {variant} {drug:7s}  r={r:.3f}  recovery={recovery_pct:.1f}%  (n={len(shared2):,})")

scores = pd.DataFrame(records)

# ── Plot ──────────────────────────────────────────────────────────────────────
VARIANT_COLORS = {
    'I733V': '#4393C3',   # blue
    'G822E': '#F4A582',   # orange
    'G731K': '#D6604D',   # red
}
VARIANTS  = ['G731K', 'G822E', 'I733V']
DRUGS     = ['HQ461', 'ZLY-025']

fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.2), sharey=False)
fig.subplots_adjust(wspace=0.38)

for ax, drug in zip(axes, DRUGS):
    sub = scores[scores['drug'] == drug].set_index('variant').loc[VARIANTS]

    x_pos  = np.arange(len(VARIANTS))
    colors = [VARIANT_COLORS[v] for v in VARIANTS]
    bars   = ax.bar(x_pos, sub['r'], color=colors, width=0.55,
                    edgecolor='black', linewidth=0.7, zorder=3)

    # Recovery rate as text above each bar
    for xi, (v, row) in zip(x_pos, sub.iterrows()):
        ax.text(xi, row['r'] + 0.015, f"{row['recovery_pct']:.0f}%",
                ha='center', va='bottom', fontsize=7.5, color='#333333')

    # NT reference line
    ax.axhline(1.0, color='#555555', lw=1.0, linestyle='--', zorder=2, label='NT (r = 1)')

    ax.set_xlim(-0.5, len(VARIANTS) - 0.5)
    ax.set_ylim(0, 1.18)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(VARIANTS, fontsize=9)
    ax.set_title(drug, fontsize=10, fontweight='bold', pad=6)
    ax.set_ylabel('Pearson r\n(vs NT drug response)', fontsize=8.5)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax.tick_params(axis='both', labelsize=8)

    # Remove top/right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Subtle grid
    ax.yaxis.grid(True, color='#e0e0e0', linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)

# Shared legend annotation
fig.text(0.5, 0.01,
         'Pearson r: genome-wide LFC correlation (variant drug response vs NT drug response)\n'
         'Percentages: NT drug response DEGs recovered in variant (padj < 0.05, |LFC| > 1, same direction)',
         ha='center', va='bottom', fontsize=6.5, color='#555555', style='italic')

fig.suptitle('Transcriptional Drug Response Preservation — PC3 CDK12 variants',
             fontsize=10, y=1.02)

outpath = os.path.join(FIGDIR, 'NT_resistance_score')
plt.savefig(f'{outpath}.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{outpath}.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{outpath}.svg', dpi=300, bbox_inches='tight')
plt.close()
print(f"\nSaved: {outpath}.{{pdf,png,svg}}")
print("\nScores summary:")
print(scores[['variant','drug','r','recovery_pct','n_genes']].to_string(index=False))
