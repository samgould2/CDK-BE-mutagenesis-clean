"""
Does the 9p21 (p16) loss plus CDK6 up-regulation look like compensation for CDK4 loss?

Tests the prediction that compensation should scale with how much CDK4 function is lost:
the CDK4-null line should carry the most CDK6 / cyclin D and the least INK4 inhibitor.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
OUT, FIG = HERE / 'results', HERE / 'figures'

GENOTYPES = ['CBE_NT', 'CDK4_D196N', 'CDK4_E56K', 'CDK4_W198']
LABEL = {'CBE_NT': 'CBE NT', 'CDK4_D196N': 'D196N', 'CDK4_E56K': 'E56K', 'CDK4_W198': 'W198*'}
COLOR = {'CBE_NT': '#9e9e9e', 'CDK4_D196N': '#4393c3', 'CDK4_E56K': '#66bd63', 'CDK4_W198': '#d6604d'}

tpm = pd.read_csv(BASE / 'salmon.merged.gene_tpm.tsv', sep='\t', index_col=0)
gn = tpm.pop('gene_name')
meta = pd.Series(tpm.columns.str.extract(r'^T47D_(.+)_veh_\d+$')[0].values, index=tpm.columns)
gm = tpm.T.groupby(meta).mean().T[GENOTYPES]
sd = tpm.T.groupby(meta).std().T[GENOTYPES]
order = tpm.mean(axis=1).sort_values(ascending=False).index
s2i = pd.Series(order, index=gn.reindex(order).values); s2i = s2i[~s2i.index.duplicated()]

fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.2), gridspec_kw={'width_ratios': [1.75, 1]})

# (a) individual pathway genes
genes = ['CDK4', 'CDK6', 'CCND1', 'CDKN2A', 'CDKN2D', 'CDKN1B']
ax = axes[0]; w = 0.2
for k, g in enumerate(GENOTYPES):
    v = [gm.loc[s2i[s], g] for s in genes]
    e = [sd.loc[s2i[s], g] for s in genes]
    ax.bar(np.arange(len(genes)) + (k - 1.5) * w, v, width=w, color=COLOR[g], edgecolor='black',
           lw=.4, label=LABEL[g], yerr=e, error_kw=dict(lw=.6, capsize=1.5))
ax.set_yscale('log'); ax.set_ylim(0.3, 1200)
ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, fontsize=8, style='italic')
ax.set_ylabel('TPM (log scale)', fontsize=8.5); ax.tick_params(labelsize=7.5)
ax.legend(fontsize=7, frameon=False, ncol=4, loc='upper center', bbox_to_anchor=(.5, 1.16))
ax.spines[['top', 'right']].set_visible(False)

# (b) budgets: is the net cyclin D-CDK capacity actually restored?
ax = axes[1]
tot = lambda lst, g: sum(gm.loc[s2i[s], g] for s in lst)
budgets = [('CDK4 + CDK6', ['CDK4', 'CDK6']),
           ('INK4 total\n(p16+p15+p18+p19)', ['CDKN2A', 'CDKN2B', 'CDKN2C', 'CDKN2D'])]
w = 0.2
for k, g in enumerate(GENOTYPES):
    v = [tot(lst, g) / tot(lst, 'CBE_NT') for _, lst in budgets]
    ax.bar(np.arange(len(budgets)) + (k - 1.5) * w, v, width=w, color=COLOR[g], edgecolor='black', lw=.4)
ax.axhline(1, color='black', lw=.8, ls='--')
ax.set_xticks(range(len(budgets))); ax.set_xticklabels([n for n, _ in budgets], fontsize=8)
ax.set_ylabel('total mRNA relative to CBE NT', fontsize=8.5)
ax.tick_params(labelsize=7.5); ax.spines[['top', 'right']].set_visible(False)
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'cdk6_compensation_test.{e}', dpi=300, bbox_inches='tight')
plt.close()

rows = []
for s in ['CDK4', 'CDK6', 'CCND1', 'CCND3', 'CDKN2A', 'CDKN2C', 'CDKN2D', 'CDKN1B', 'RB1', 'CCNE1', 'CDK2']:
    rows.append(dict(gene=s, **{LABEL[g]: gm.loc[s2i[s], g] for g in GENOTYPES}))
for n, lst in budgets + [('cyclin D total', ['CCND1', 'CCND2', 'CCND3'])]:
    rows.append(dict(gene=n.replace('\n', ' '), **{LABEL[g]: tot(lst, g) for g in GENOTYPES}))
pd.DataFrame(rows).set_index('gene').round(1).to_csv(OUT / 'cdk6_compensation_budget.csv')
print(pd.DataFrame(rows).set_index('gene').round(1).to_string())
