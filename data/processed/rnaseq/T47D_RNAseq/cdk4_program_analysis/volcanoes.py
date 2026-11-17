"""
Volcano plots in the repository style (see t47d_cdk4_variants_20260927_deseq2_gsea.ipynb),
including the W198* vs D196N contrast, which was not in the original notebook.

Refits the shared ~condition DESeq2 model over all 12 DMSO samples so that the
variant-vs-variant contrast comes from the same dispersion estimates as the
variant-vs-control ones.
"""
import sys
sys.path.insert(0, '/private/tmp/claude-501/pydeseq2_lib')

from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
OUT  = HERE / 'results'; OUT.mkdir(exist_ok=True)
FIG  = HERE / 'figures'; FIG.mkdir(exist_ok=True)

# ── repo conventions ────────────────────────────────────────────────────────
LFC_THRESH, PADJ_THRESH = 1.0, 0.05
UP_COLOR, DOWN_COLOR, NS_COLOR, THRESH_COLOR = '#d6604d', '#4393c3', '#d0d0d0', '#999999'
X_LIM, N_LABELS = 6, 12
ALWAYS_MARK = ['CDK4']
MIN_COUNT, MIN_SAMPLES = 10, 3

counts = pd.read_csv(BASE / 'salmon.merged.gene_counts_scaled.tsv', sep='\t', index_col=0)
gene_name = counts.pop('gene_name')
meta = pd.DataFrame(index=counts.columns)
meta['condition'] = counts.columns.str.extract(r'^T47D_(.+)_veh_\d+$')[0].values
# pydeseq2 splits design-matrix column names on '_', so levels must not contain one
LEVEL = {'CBE_NT': 'NT', 'CDK4_D196N': 'D196N', 'CDK4_E56K': 'E56K', 'CDK4_W198': 'W198'}
meta['condition'] = meta['condition'].map(LEVEL)

CONTRASTS = [
    ('D196N', 'NT',    'CDK4_D196N_vs_CBE_NT_veh',    'CDK4 D196N vs CBE NT (DMSO)'),
    ('E56K',  'NT',    'CDK4_E56K_vs_CBE_NT_veh',     'CDK4 E56K vs CBE NT (DMSO)'),
    ('W198',  'NT',    'CDK4_W198_vs_CBE_NT_veh',     'CDK4 W198* vs CBE NT (DMSO)'),
    ('W198',  'D196N', 'CDK4_W198_vs_CDK4_D196N_veh', 'CDK4 W198* vs CDK4 D196N (DMSO)'),
]

csvs = {n: OUT / f'{n}.csv' for _, _, n, _ in CONTRASTS}
if not all(p.exists() for p in csvs.values()):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    X = counts.T
    keep = (X >= MIN_COUNT).sum(axis=0) >= MIN_SAMPLES
    X = X.loc[:, keep].round().astype(int)
    print(f'{X.shape[1]} genes kept of {counts.shape[0]}')
    dds = DeseqDataSet(counts=X, metadata=meta, design_factors='condition', refit_cooks=True, quiet=True)
    dds.deseq2()
    for a, b, name, _ in CONTRASTS:
        st = DeseqStats(dds, contrast=['condition', a, b], quiet=True)
        st.summary()
        r = st.results_df.copy()
        r['gene_name'] = gene_name.reindex(r.index).values
        r.index.name = 'gene_id'
        r.to_csv(csvs[name])
        print(f'{name}: {int(((r.padj < PADJ_THRESH) & (r.log2FoldChange.abs() > LFC_THRESH)).sum())} DEGs')

results = {n: pd.read_csv(p).set_index('gene_id') for n, p in csvs.items()}

try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False

HALLMARK = json.loads((BASE / 'T47D_CDK4_20260927_results' / 'hallmark_gene_sets_cache.json').read_text())
_E2F, _G2M = set(HALLMARK.get('E2F Targets', [])), set(HALLMARK.get('G2-M Checkpoint', []))
CC_STYLE = {'E2F Targets': ('#e6850d', _E2F - _G2M), 'G2-M Checkpoint': ('#c51b8a', _G2M - _E2F),
            'E2F + G2-M': ('#000000', _E2F & _G2M)}
CC_ALL = _E2F | _G2M


def plot_volcano_grid(comparison_list, ncols, outpath, label_mode='top', panel=4.3):
    """Repo-standard volcano grid (ported from the 20260927 notebook)."""
    nrows = int(np.ceil(len(comparison_list) / ncols))
    ys = [-np.log10(results[n]['padj'].dropna().clip(lower=1e-300)) for _, _, n, _ in comparison_list]
    y_lim = np.ceil(max(y.max() for y in ys)) * 1.05 + 1
    fig, axes = plt.subplots(nrows, ncols, figsize=(panel * ncols, panel * nrows), squeeze=False)
    axes = axes.flatten()
    for k, (ax, (_, _, name, title)) in enumerate(zip(axes, comparison_list)):
        rc = results[name].dropna(subset=['padj', 'log2FoldChange']).copy()
        rc['y'] = -np.log10(rc['padj'].clip(lower=1e-300))
        rc['x'] = rc['log2FoldChange'].clip(-X_LIM * 0.98, X_LIM * 0.98)
        up   = (rc['padj'] < PADJ_THRESH) & (rc['log2FoldChange'] >  LFC_THRESH)
        down = (rc['padj'] < PADJ_THRESH) & (rc['log2FoldChange'] < -LFC_THRESH)
        ns = ~(up | down)
        ax.scatter(rc.loc[ns, 'x'], rc.loc[ns, 'y'], s=4, color=NS_COLOR, alpha=0.5, lw=0, rasterized=True)
        up_pts = ax.scatter(rc.loc[up, 'x'], rc.loc[up, 'y'], s=10, color=UP_COLOR, alpha=0.85, lw=0, rasterized=True)
        dn_pts = ax.scatter(rc.loc[down, 'x'], rc.loc[down, 'y'], s=10, color=DOWN_COLOR, alpha=0.85, lw=0, rasterized=True)

        if label_mode == 'cell_cycle':
            for lbl, (col, gset) in CC_STYLE.items():
                sel = rc[rc['gene_name'].isin(gset)]
                ax.scatter(sel['x'], sel['y'], s=20, color=col, edgecolors='black', lw=0.3, zorder=5,
                           label=(lbl if k == 0 else None))
            pool = rc[rc['gene_name'].isin(CC_ALL) & (up | down)]
        else:
            pool = rc[up | down]
        to_label = pool.nsmallest(N_LABELS, 'padj')
        mark = rc[rc['gene_name'].isin(ALWAYS_MARK)]
        ax.scatter(mark['x'], mark['y'], s=40, facecolor='none', edgecolor='black', lw=1.0, zorder=6)
        to_label = pd.concat([to_label, mark])
        to_label = to_label[~to_label.index.duplicated()]
        texts = [ax.text(r['x'], r['y'], r['gene_name'], fontsize=7, style='italic')
                 for _, r in to_label.iterrows()]
        if HAS_ADJUSTTEXT and texts:
            try:
                adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='black', lw=0.4),
                            expand=(1.2, 1.4), time_lim=5, iter_lim=300)
            except (AttributeError, TypeError):      # older adjustText signature
                adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='black', lw=0.4),
                            lim=300)

        ax.set_xlim(-X_LIM, X_LIM); ax.set_ylim(0, y_lim)
        for xv in (LFC_THRESH, -LFC_THRESH):
            ax.axvline(xv, color=THRESH_COLOR, lw=1.0, ls='--', alpha=0.8)
        ax.axhline(-np.log10(PADJ_THRESH), color=THRESH_COLOR, lw=1.0, ls='--', alpha=0.8)
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_xlabel(r'log$_2$ fold change', fontsize=10)
        ax.set_ylabel(r'-log$_{10}$(padj)', fontsize=10)
        ax.set_title(title, fontsize=10)
        leg = ax.legend(handles=[up_pts, dn_pts],
                        labels=[f'{int(up.sum())} DE UP', f'{int(down.sum())} DE DOWN'],
                        loc='upper right', fontsize=8, frameon=True, framealpha=1, edgecolor='#cccccc')
        ax.add_artist(leg)
    for ax in axes[len(comparison_list):]:
        ax.axis('off')
    if label_mode == 'cell_cycle':
        axes[0].legend(fontsize=7, frameon=False, loc='upper left', title='Cell cycle', title_fontsize=7)
    plt.tight_layout()
    plt.savefig(f'{outpath}.svg', dpi=300, bbox_inches='tight')
    plt.savefig(f'{outpath}.png', dpi=300, bbox_inches='tight')
    plt.close()


plot_volcano_grid(CONTRASTS, ncols=2, outpath=FIG / 'volcano_all_contrasts_DMSO', label_mode='top')
plot_volcano_grid(CONTRASTS, ncols=2, outpath=FIG / 'volcano_all_contrasts_DMSO_cellcycle', label_mode='cell_cycle')
plot_volcano_grid(CONTRASTS[3:], ncols=1, outpath=FIG / 'volcano_W198_vs_D196N', label_mode='cell_cycle', panel=4.8)

print('\nDEG counts (padj<0.05, |log2FC|>1):')
for _, _, n, t in CONTRASTS:
    r = results[n].dropna(subset=['padj', 'log2FoldChange'])
    u = int(((r.padj < PADJ_THRESH) & (r.log2FoldChange > LFC_THRESH)).sum())
    d = int(((r.padj < PADJ_THRESH) & (r.log2FoldChange < -LFC_THRESH)).sum())
    print(f'  {t:36s} up {u:5d}  down {d:5d}  total {u + d:5d}')
cdk4 = results['CDK4_W198_vs_CDK4_D196N_veh']
cdk4 = cdk4[cdk4.gene_name == 'CDK4']
print('\nCDK4 in W198* vs D196N:', cdk4[['log2FoldChange', 'padj']].to_string(index=False))
