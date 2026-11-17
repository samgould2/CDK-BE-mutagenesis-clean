"""Shared figure style, palette and loaders for the PC3 CDK12 re-analysis."""
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
RES = os.path.join(BASE, 'results')
FIG = os.path.join(BASE, 'figures')
os.makedirs(FIG, exist_ok=True)

GENOS = ['NT', 'I733V', 'G822E', 'G731K']
VARS = GENOS[1:]
DRUGS = ['HQ461', 'ZLY025']
TREATS = ['Control', 'HQ461', 'ZLY025']

# Palette: neutral ink for the NT reference + blue / gold / plum for the variants.
# Validated colourblind-safe for the three variant hues; gold sits below 3:1 contrast
# on white, which the black bar outlines (BAR_EDGE) provide relief for.
GCOL = {'NT': '#3D3D3D', 'I733V': '#1B6CA8', 'G822E': '#E1A730', 'G731K': '#8E4585'}
TLAB = {'Control': 'DMSO', 'HQ461': 'HQ461', 'ZLY025': 'ZLY-025'}
TMARK = {'Control': 'o', 'HQ461': '^', 'ZLY025': 's'}
TFILL = {'Control': '#F0EFEC', 'HQ461': '#898781', 'ZLY025': '#3D3D3D'}
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
    # PDF transparent for compositing in Illustrator; PNG kept on white for quick viewing
    fig.savefig(os.path.join(FIG, name + '.pdf'), transparent=True)
    fig.savefig(os.path.join(FIG, name + '.png'), dpi=300)
    plt.close(fig)
    print('saved', name)


def load(name, sub='contrasts'):
    return pd.read_csv(os.path.join(RES, sub, name + '.csv'), index_col='gene_id')


def drug_con(g, d):
    return load(f'{g}_{d}_vs_{g}_Control')


def meta():
    return pd.read_csv(os.path.join(RES, 'sample_metadata.csv')).set_index('sample')


def features():
    """Per-gene structure from the Ensembl 113 canonical transcript."""
    return pd.read_csv(os.path.join(RES, 'annotation', 'gene_features.tsv'), sep='\t', index_col='gene_id')


def expressed_pc(min_base_mean=50):
    """Protein-coding genes expressed across the experiment (baseMean in the NT contrast)."""
    f = features()
    r = drug_con('NT', 'HQ461')
    idx = r.index[r.baseMean > min_base_mean]
    f = f.reindex(idx).dropna(subset=['span_bp'])
    return f[f.biotype == 'protein_coding']


# ── Elongation helpers (shared by 04_elongation.py and 05_elongation_length.py) ──
def decile_profile(values, lfc, n=10):
    """Bin assignment, bin centres, and median/IQR/size of lfc within n equal-count bins."""
    import pandas as _pd
    b = _pd.qcut(values.rank(method='first'), n, labels=False)
    g = lfc.groupby(b)
    return b, values.groupby(b).median(), g.median(), g.quantile(0.25), g.quantile(0.75), g.size()


def gap_ci(top, allgenes, rng, n_boot=2000):
    """median(log2FC | top decile) - median(log2FC | all genes), bootstrap 95% CI.

    Preferred over an OLS slope because some profiles are non-monotonic, where a linear
    slope is a poor summary. Referenced to the overall median rather than the bottom
    decile because single-intron genes behave atypically under ZLY-025 and would
    otherwise dominate the contrast.
    """
    import numpy as _np
    top = top[_np.isfinite(top)]
    allgenes = allgenes[_np.isfinite(allgenes)]
    g = _np.median(top) - _np.median(allgenes)
    boots = [_np.median(top[rng.integers(0, len(top), len(top))]) -
             _np.median(allgenes[rng.integers(0, len(allgenes), len(allgenes))]) for _ in range(n_boot)]
    return g, *_np.percentile(boots, [2.5, 97.5])


def replicate_profiles(v, bins, test, ref):
    """Per-replicate median log2 ratio within each bin (rows = replicates).

    v    VST matrix restricted to the genes of interest
    bins bin assignment per gene (e.g. gene-length decile)
    test replicate sample names to profile
    ref  replicate sample names whose mean is the denominator
    """
    import pandas as _pd
    r = v[ref].mean(axis=1)
    return _pd.DataFrame([(v[s] - r).groupby(bins).median().rename(s) for s in test])
