"""Shared figure style, palette and loaders for the LNCaP CDK12 re-analysis.

Deliberately mirrors the PC3 `_style.py` (same palette, same rcParams, same helper
semantics) so figures from the two cell lines are directly comparable. Kept as a separate
self-contained copy rather than imported across datasets: the designs differ (2x2 here vs
4x3 there) and each analysis directory should stand alone.
"""
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
RES = os.path.join(BASE, 'results')
FIG = os.path.join(BASE, 'figures')
os.makedirs(FIG, exist_ok=True)

GENOS = ['NT', 'I733V']
VARS = GENOS[1:]
DRUGS = ['HQ461']
TREATS = ['Control', 'HQ461']

GCOL = {'NT': '#3D3D3D', 'I733V': '#1B6CA8'}
TLAB = {'Control': 'DMSO', 'HQ461': 'HQ461'}
TMARK = {'Control': 'o', 'HQ461': '^'}
DIVERGING = 'RdBu_r'
INK, MUTED, GRID = '#0B0B0B', '#898781', '#E1E0D9'
BAR_EDGE = dict(edgecolor='black', linewidth=1.0)
PADJ, LFC = 0.05, 1.0

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7, 'axes.titlesize': 7, 'axes.labelsize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 6,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'legend.frameon': False, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
})


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + '.pdf'), transparent=True)
    fig.savefig(os.path.join(FIG, name + '.png'), dpi=300)
    plt.close(fig)
    print('saved', name)


def load(name):
    return pd.read_csv(os.path.join(RES, 'contrasts', name + '.csv'), index_col='gene_id')


def drug_con(g, d='HQ461'):
    return load(f'{g}_{d}_vs_{g}_Control')


def meta():
    return pd.read_csv(os.path.join(RES, 'sample_metadata.csv')).set_index('sample')


def features():
    return pd.read_csv(os.path.join(RES, 'annotation', 'gene_features.tsv'), sep='\t', index_col='gene_id')


def expressed_pc(min_base_mean=50):
    """Protein-coding genes expressed across the experiment (DESeq2 baseMean, all samples)."""
    f = features()
    r = drug_con('NT')
    idx = r.index[r.baseMean > min_base_mean]
    f = f.reindex(idx).dropna(subset=['span_bp'])
    return f[f.biotype == 'protein_coding']


def replicate_profiles(v, bins, test, ref):
    """Per-replicate median log2 ratio within each bin (rows = replicates)."""
    r = v[ref].mean(axis=1)
    return pd.DataFrame([(v[s] - r).groupby(bins).median().rename(s) for s in test])
