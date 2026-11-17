"""
T47D-specific context for the ER/luminal collapse: the variant lines lost the entire
oestrogen-response programme without paying any proliferative cost, which is why an
ER-null subclone can sweep a T47D pool in serum-containing medium.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, str(Path(__file__).resolve().parent))
import genesets as GS

mpl.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})

HERE = Path(__file__).resolve().parent
BASE, OUT, FIG = HERE.parent, HERE / 'results', HERE / 'figures'
C = ['CBE_NT', 'CDK4_D196N', 'CDK4_E56K', 'CDK4_W198']
LAB = {'CBE_NT': 'CBE NT', 'CDK4_D196N': 'D196N', 'CDK4_E56K': 'E56K', 'CDK4_W198': 'W198*'}
COL = {'CBE_NT': '#9e9e9e', 'CDK4_D196N': '#4393c3', 'CDK4_E56K': '#66bd63', 'CDK4_W198': '#d6604d'}

tpm = pd.read_csv(BASE / 'salmon.merged.gene_tpm.tsv', sep='\t', index_col=0)
gn = tpm.pop('gene_name')
meta = pd.Series(tpm.columns.str.extract(r'^T47D_(.+)_veh_\d+$')[0].values, index=tpm.columns)
gm = tpm.T.groupby(meta).mean().T
sdv = tpm.T.groupby(meta).std().T
log_tpm = np.log2(tpm + 1)
order = tpm.mean(axis=1).sort_values(ascending=False).index
s2i = pd.Series(order, index=gn.reindex(order).values); s2i = s2i[~s2i.index.duplicated()]

ER = ['ESR1', 'PGR', 'GREB1', 'TFF1', 'TFF3', 'STC2', 'CA12', 'AREG', 'SIAH2', 'NRIP1']
PROLIF = {'S-phase': GS.S_SPECIFIC, 'G2/M': GS.G2M_SPECIFIC, 'core E2F': GS.CORE_E2F}

fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.3), gridspec_kw={'width_ratios': [1.55, 1]})

# (a) ER programme genes
ax = axes[0]; w = 0.2
for k, c in enumerate(C):
    v = [gm.loc[s2i[g], c] for g in ER]
    e = [sdv.loc[s2i[g], c] for g in ER]
    ax.bar(np.arange(len(ER)) + (k - 1.5) * w, np.maximum(v, 0.02), width=w, color=COL[c],
           edgecolor='black', lw=.4, label=LAB[c], yerr=e, error_kw=dict(lw=.5, capsize=1.2))
ax.set_yscale('log'); ax.set_ylim(0.015, 900)
ax.set_xticks(range(len(ER))); ax.set_xticklabels(ER, rotation=45, ha='right',
                                                  rotation_mode='anchor', fontsize=8, style='italic')
ax.set_ylabel('TPM (log scale)', fontsize=8.5); ax.tick_params(labelsize=7.5)
handles, labels = ax.get_legend_handles_labels()
ax.spines[['top', 'right']].set_visible(False)
ax.set_title('ER / oestrogen-response programme', fontsize=9, pad=6)

# (b) the trade-off that wasn't
def geo(sy):
    g = [s2i[s] for s in sy if s in s2i.index and log_tpm.loc[s2i[s]].mean() > 1]
    return log_tpm.loc[g].mean(axis=0)

ax = axes[1]
tracks = [('ER programme', ER)] + list(PROLIF.items())
w = 0.2
rows = []
for k, c in enumerate(C):
    vals = []
    for nm, sy in tracks:
        s = geo(sy)
        vals.append(2 ** (s[meta[meta == c].index].mean() - s[meta[meta == 'CBE_NT'].index].mean()))
        rows.append(dict(genotype=LAB[c], track=nm, fold_vs_NT=vals[-1]))
    ax.bar(np.arange(len(tracks)) + (k - 1.5) * w, vals, width=w, color=COL[c], edgecolor='black', lw=.4)
ax.axhline(1, color='black', lw=.8, ls='--')
ax.set_xticks(range(len(tracks))); ax.set_xticklabels([n for n, _ in tracks], rotation=20,
                                                      ha='right', rotation_mode='anchor', fontsize=8)
ax.set_ylabel('abundance relative to CBE NT', fontsize=8.5)
ax.set_ylim(0, 1.35); ax.tick_params(labelsize=7.5)
ax.spines[['top', 'right']].set_visible(False)
ax.set_title('ER collapses; proliferation does not', fontsize=9, pad=6)
pd.DataFrame(rows).to_csv(OUT / 'er_vs_proliferation.csv', index=False)
fig.legend(handles, labels, fontsize=8, frameon=False, ncol=4,
           loc='upper center', bbox_to_anchor=(.5, 1.07))
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'ER_context_T47D.{e}', dpi=300, bbox_inches='tight')
plt.close()
print(pd.DataFrame(rows).pivot(index='track', columns='genotype', values='fold_vs_NT').round(2).to_string())
