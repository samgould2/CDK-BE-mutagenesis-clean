"""
Is E56K's weaker 9p21 loss a property of the locus, or of the pool?

Plots, for each variant line, the residual expression of two independent marker
classes as a fraction of CBE NT: genes inside the 9p21 deleted interval, and
ER/luminal-identity genes lost in the de-differentiated clone. If E56K simply
contains more unswept parental cells, both classes rise together by the same factor.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

mpl.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})

HERE = Path(__file__).resolve().parent
BASE, OUT, FIG = HERE.parent, HERE / 'results', HERE / 'figures'

VARIANTS = ['CDK4_D196N', 'CDK4_E56K', 'CDK4_W198']
LABEL = {'CDK4_D196N': 'D196N', 'CDK4_E56K': 'E56K', 'CDK4_W198': 'W198*'}
DEL = ['FOCAD', 'HACD4', 'KLHL9', 'MTAP', 'CDKN2A']
LUM = ['ESR1', 'PGR', 'GREB1', 'CDH1', 'FOXA1', 'GATA3', 'AGR2', 'EPCAM']

tpm = pd.read_csv(BASE / 'salmon.merged.gene_tpm.tsv', sep='\t', index_col=0)
gn = tpm.pop('gene_name')
meta = pd.Series(tpm.columns.str.extract(r'^T47D_(.+)_veh_\d+$')[0].values, index=tpm.columns)
gm = tpm.T.groupby(meta).mean().T
order = tpm.mean(axis=1).sort_values(ascending=False).index
s2i = pd.Series(order, index=gn.reindex(order).values); s2i = s2i[~s2i.index.duplicated()]

rows = []
for cls, genes in [('9p21 interval', DEL), ('ER / luminal', LUM)]:
    for s in genes:
        g = s2i[s]
        for v in VARIANTS:
            rows.append(dict(marker_class=cls, gene=s, line=v,
                             residual_pct=100 * gm.loc[g, v] / gm.loc[g, 'CBE_NT']))
df = pd.DataFrame(rows)
df.to_csv(OUT / 'pool_purity_residual_fractions.csv', index=False)

fig, ax = plt.subplots(figsize=(5.2, 3.6))
style = {'9p21 interval': ('#d6604d', 'o'), 'ER / luminal': ('#4393c3', 's')}
rng = np.random.default_rng(0)
for cls, (c, m) in style.items():
    sub = df[df.marker_class == cls]
    for gene, d in sub.groupby('gene'):
        d = d.set_index('line').loc[VARIANTS]
        xs = np.arange(3) + (0.1 if cls == '9p21 interval' else -0.1)
        ax.plot(xs, d.residual_pct.clip(lower=0.05), color=c, lw=.6, alpha=.35, zorder=1)
        ax.scatter(xs, d.residual_pct.clip(lower=0.05), color=c, marker=m, s=22,
                   edgecolor='black', lw=.4, zorder=3)
ax.set_yscale('log')
ax.set_xticks(range(3)); ax.set_xticklabels([LABEL[v] for v in VARIANTS], fontsize=9)
ax.set_ylabel('residual expression\n(% of CBE NT)', fontsize=9)
ax.legend(handles=[Line2D([], [], marker=m, ls='', color=c, mec='black', mew=.4, ms=6, label=k)
                   for k, (c, m) in style.items()], fontsize=8, frameon=False, loc='lower right')
ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=8)
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'pool_purity.{e}', dpi=300, bbox_inches='tight')
plt.close()
print(df.pivot_table(index=['marker_class', 'gene'], columns='line', values='residual_pct').round(1).to_string())
