"""
nt_gsea_enrichment_plots.py
───────────────────────────
GSEA enrichment plots for the core NT CDK inhibition UP signature
(intersection of HQ461 and ZLY-025 NT DEGs, padj<0.05, LFC>1)
scored against each CDK12 variant drug response.

Layout: 3 rows (variants, most resistant → most sensitive) × 2 cols (drugs)
Each panel: enrichment curve + hit rug + NES / FDR annotation

Outputs
-------
    figures/04_gsea_cdk12/NT_GSEA_enrichment.{pdf,png,svg}

Run from project root
    cd ~/Desktop/"Comprehensive PC3 RNA-seq Analysis"
    python3 scripts/nt_gsea_enrichment_plots.py
"""

import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import gseapy as gs

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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
    df = df.dropna(subset=['stat', 'gene_name'])
    df = df.drop_duplicates(subset='gene_name').set_index('gene_name')
    return df

print("Loading DESeq2 results...")
nt_hq461 = load_csv('NT_HQ461_vs_NT_Control')
nt_zly   = load_csv('NT_ZLY025_vs_NT_Control')

def de_set(df, direction):
    mask = (df['padj'] < 0.05) & (
        (df['log2FoldChange'] > 1) if direction == 'up' else (df['log2FoldChange'] < -1)
    )
    return sorted(set(df[mask].index))

core_up = de_set(nt_hq461, 'up')
core_up = sorted(set(core_up) & set(de_set(nt_zly, 'up')))
print(f"  Core UP signature: {len(core_up)} genes")

GENE_SETS = {'NT_CORE_UP': core_up}

def run_prerank(df, seed=42):
    ranking = df[['stat']].sort_values('stat', ascending=False)
    pre = gs.prerank(
        rnk=ranking, gene_sets=GENE_SETS,
        permutation_num=1000, seed=seed,
        min_size=5, max_size=5000,
        verbose=False, outdir=None,
    )
    return pre

# Most resistant → most sensitive per drug
COMPARISONS = [
    ('G731K_HQ461_vs_G731K_Control',  'G731K', 'HQ461'),
    ('I733V_HQ461_vs_I733V_Control',  'I733V', 'HQ461'),
    ('G822E_HQ461_vs_G822E_Control',  'G822E', 'HQ461'),
    ('G731K_ZLY025_vs_G731K_Control', 'G731K', 'ZLY-025'),
    ('G822E_ZLY025_vs_G822E_Control', 'G822E', 'ZLY-025'),
    ('I733V_ZLY025_vs_I733V_Control', 'I733V', 'ZLY-025'),
]

VARIANT_COLORS = {
    'I733V': '#4393C3',
    'G822E': '#F4A582',
    'G731K': '#D6604D',
}

# ── Precompute NT reference enrichment curves ─────────────────────────────────
print("\nComputing NT reference curves...")
nt_ref_curves = {}
for drug, ref_df in [('HQ461', nt_hq461), ('ZLY-025', nt_zly)]:
    ranking  = ref_df[['stat']].sort_values('stat', ascending=False)
    genes    = list(ranking.index)
    n        = len(genes)
    sig_set  = set(core_up)
    hits     = np.array([1 if g in sig_set else 0 for g in genes])
    n_hits   = hits.sum()
    n_miss   = n - n_hits
    weights  = np.abs(ranking['stat'].values) * hits
    weight_sum = weights.sum()
    miss_val   = 1.0 / n_miss if n_miss > 0 else 0
    running  = np.zeros(n)
    cumw = cumm = 0.0
    for i in range(n):
        if hits[i]:
            cumw += weights[i] / weight_sum
            running[i] = cumw - cumm
        else:
            cumm += miss_val
            running[i] = cumw - cumm
    nt_ref_curves[drug] = {'curve': running, 'NES': None}
    print(f"  NT {drug}: peak ES = {running.max():.3f}")

# Get NT NES values from self-prerank
print("\nComputing NT self NES (for ΔNES)...")
for drug, ref_df in [('HQ461', nt_hq461), ('ZLY-025', nt_zly)]:
    pre = run_prerank(ref_df)
    nt_ref_curves[drug]['NES'] = float(pre.res2d.set_index('Term').loc['NT_CORE_UP']['NES'])
    print(f"  NT {drug}: NES = {nt_ref_curves[drug]['NES']:.3f}")

print("\nRunning GSEA (1000 permutations × 6 comparisons)...")
results = []
for fname, variant, drug in COMPARISONS:
    df  = load_csv(fname)
    pre = run_prerank(df)
    res = pre.res2d.set_index('Term').loc['NT_CORE_UP']
    nes = float(res['NES'])
    fdr = float(res['FDR q-val'])
    nom = float(res['NOM p-val'])

    # Enrichment curve data
    ranking = df[['stat']].sort_values('stat', ascending=False)
    genes   = list(ranking.index)
    n       = len(genes)
    sig_set = set(core_up)
    hits    = np.array([1 if g in sig_set else 0 for g in genes])
    n_hits  = hits.sum()
    n_miss  = n - n_hits

    # Running enrichment score (weighted, p=1)
    weights = np.abs(ranking['stat'].values) * hits
    weight_sum = weights.sum()
    miss_val   = 1.0 / n_miss if n_miss > 0 else 0

    running = np.zeros(n)
    cumw    = 0.0
    cumm    = 0.0
    for i in range(n):
        if hits[i]:
            cumw    += weights[i] / weight_sum
            running[i] = cumw - cumm
        else:
            cumm    += miss_val
            running[i] = cumw - cumm

    results.append({
        'variant': variant, 'drug': drug,
        'NES': nes, 'FDR': fdr, 'NOM': nom,
        'running': running, 'hits': hits, 'n': n,
    })
    print(f"  {variant} {drug:7s}  NES={nes:+.3f}  FDR={fdr:.3f}  NOM={nom:.3f}")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(8.0, 10.5))
outer = gridspec.GridSpec(3, 2, figure=fig, hspace=0.50, wspace=0.38)

for idx, rec in enumerate(results):
    row_i = idx % 3
    col_i = idx // 3

    inner = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer[row_i, col_i],
        height_ratios=[5, 1], hspace=0.06
    )
    ax_es  = fig.add_subplot(inner[0])
    ax_rug = fig.add_subplot(inner[1], sharex=ax_es)

    color   = VARIANT_COLORS[rec['variant']]
    running = rec['running']
    hits    = rec['hits']
    n       = rec['n']
    x       = np.arange(n)

    # NT reference curve (grey background)
    nt_curve = nt_ref_curves[rec['drug']]['curve']
    nt_x     = np.arange(len(nt_curve))
    ax_es.plot(nt_x, nt_curve, color='#AAAAAA', lw=1.2, zorder=2, label='NT reference')
    ax_es.fill_between(nt_x, nt_curve, 0,
                        where=(nt_curve > 0), color='#AAAAAA', alpha=0.12, zorder=1)

    # Variant enrichment curve
    ax_es.plot(x, running, color=color, lw=1.5, zorder=3, label=rec['variant'])
    ax_es.axhline(0, color='#888888', lw=0.5, zorder=1)
    ax_es.fill_between(x, running, 0,
                        where=(running > 0), color=color, alpha=0.15, zorder=2)

    # NES / FDR annotation
    fdr_str = f'FDR = {rec["FDR"]:.3f}' if rec["FDR"] >= 0.001 else 'FDR < 0.001'
    nom_str = f'NOM p = {rec["NOM"]:.3f}' if rec["NOM"] >= 0.001 else 'NOM p < 0.001'
    delta_nes = rec['NES'] - nt_ref_curves[rec['drug']]['NES']
    ax_es.text(0.97, 0.97,
               f'NES = {rec["NES"]:+.2f}  (NT: {nt_ref_curves[rec["drug"]]["NES"]:+.2f})\n'
               f'ΔNES = {delta_nes:+.2f}\n{fdr_str}\n{nom_str}',
               transform=ax_es.transAxes, va='top', ha='right', fontsize=7.5,
               color='#222222',
               bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='none', alpha=0.85))

    ax_es.set_xlim(0, n)
    ax_es.set_ylabel('Enrichment score', fontsize=7.5)
    ax_es.tick_params(labelsize=7, labelbottom=False)
    ax_es.spines['top'].set_visible(False)
    ax_es.spines['right'].set_visible(False)
    ax_es.set_title(f'{rec["variant"]} — {rec["drug"]}',
                    fontsize=9.5, fontweight='bold', pad=4)

    # Hit rug
    hit_pos = np.where(hits)[0]
    ax_rug.vlines(hit_pos, 0, 1, color=color, lw=0.4, alpha=0.6)
    ax_rug.set_xlim(0, n)
    ax_rug.set_ylim(0, 1)
    ax_rug.set_yticks([])
    ax_rug.set_xlabel('Gene rank (high → low Wald stat)', fontsize=7.5)
    ax_rug.tick_params(labelsize=7)
    ax_rug.spines['top'].set_visible(False)
    ax_rug.spines['right'].set_visible(False)
    ax_rug.spines['left'].set_visible(False)

fig.suptitle(
    'Core NT CDK Inhibition UP Signature Enrichment — PC3 CDK12 variants\n'
    'Gene set: intersection of NT HQ461 ∩ ZLY-025 upregulated DEGs (padj<0.05, LFC>1)',
    fontsize=10, y=1.01
)

outpath = os.path.join(FIGDIR, 'NT_GSEA_enrichment')
plt.savefig(f'{outpath}.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{outpath}.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{outpath}.svg', dpi=300, bbox_inches='tight')
plt.close()
print(f"\nSaved: {outpath}.{{pdf,png,svg}}")


# ── Fold-change of NES bar plot ───────────────────────────────────────────────
import matplotlib.ticker as mticker

VARIANT_ORDER_HQ  = sorted([r for r in results if r['drug']=='HQ461'],
                             key=lambda r: r['NES'])
VARIANT_ORDER_ZLY = sorted([r for r in results if r['drug']=='ZLY-025'],
                             key=lambda r: r['NES'])

fig2, axes2 = plt.subplots(1, 2, figsize=(7.5, 3.2))
fig2.subplots_adjust(wspace=0.38)

def fdr_stars(fdr):
    return '***' if fdr < 0.001 else '**' if fdr < 0.01 else '*' if fdr < 0.05 else ''

for ax, drug, ordered in zip(axes2, ['HQ461','ZLY-025'],
                              [VARIANT_ORDER_HQ, VARIANT_ORDER_ZLY]):
    nt_nes = nt_ref_curves[drug]['NES']
    x_pos  = np.arange(len(ordered))
    for xi, rec in enumerate(ordered):
        fc    = rec['NES'] / nt_nes      # 0–1, NT = 1.0
        color = VARIANT_COLORS[rec['variant']]
        ax.bar(xi, fc, width=0.55, color=color,
               edgecolor='black', linewidth=0.7, zorder=3)
        s = fdr_stars(rec['FDR'])
        if s:
            ax.text(xi, fc + 0.02, s, ha='center', va='bottom',
                    fontsize=8, color='#333333')

    ax.axhline(1.0, color='#555555', lw=1.0, linestyle='--', zorder=2)
    ax.set_xlim(-0.5, len(ordered) - 0.5)
    ax.set_ylim(0, 1.25)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([r['variant'] for r in ordered], fontsize=9)
    ax.set_title(drug, fontsize=10, fontweight='bold', pad=6)
    ax.set_ylabel('NES / NT NES\n(NT = 1.0)', fontsize=8.5)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax.tick_params(labelsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, color='#e0e0e0', linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)

fig2.text(0.5, 0.01,
          'NES fold-change = NES_variant / NES_NT for core UP signature.  '
          'NT = 1.0 (dashed).  Stars: GSEA FDR (***<0.001, **<0.01, *<0.05).',
          ha='center', va='bottom', fontsize=6.5, color='#555555', style='italic')

fig2.suptitle('Core NT CDK UP Signature — NES Fold-Change vs NT',
              fontsize=10, y=1.02)

outpath2 = os.path.join(FIGDIR, 'NT_GSEA_NES_foldchange')
fig2.savefig(f'{outpath2}.pdf', dpi=300, bbox_inches='tight')
fig2.savefig(f'{outpath2}.png', dpi=300, bbox_inches='tight')
fig2.savefig(f'{outpath2}.svg', dpi=300, bbox_inches='tight')
plt.close(fig2)
print(f"Saved: {outpath2}.{{pdf,png,svg}}")
