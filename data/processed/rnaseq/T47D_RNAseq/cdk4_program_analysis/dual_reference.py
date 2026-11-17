"""
Justifying CBE NT as the comparator for the cell-cycle quantification.

E56K carries the full background drift (9p21.3 loss, ER/luminal collapse) but near-WT
CDK4 mRNA, so E56K vs NT isolates what the background alone does to each axis. Where
that is ~1.0, CBE NT is a legitimate reference; where it is not, it is not.

Panel A: background effect (E56K vs CBE NT), with a bootstrap CI over genes.
Panel B: D196N and W198* under BOTH references - agreement is the robustness argument.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats as sps
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, str(Path(__file__).resolve().parent))
import genesets as GS

mpl.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})

HERE = Path(__file__).resolve().parent
BASE, OUT, FIG = HERE.parent, HERE / 'results', HERE / 'figures'

tpm = pd.read_csv(BASE / 'salmon.merged.gene_tpm.tsv', sep='\t', index_col=0)
gn = tpm.pop('gene_name')
meta = pd.Series(tpm.columns.str.extract(r'^T47D_(.+)_veh_\d+$')[0].values, index=tpm.columns)
lt = np.log2(tpm + 1)
order = tpm.mean(axis=1).sort_values(ascending=False).index
s2i = pd.Series(order, index=gn.reindex(order).values); s2i = s2i[~s2i.index.duplicated()]
grp = lambda g: list(meta.index[meta == g])
ALL = list(tpm.columns)
VARONLY = [c for c in ALL if meta[c] != 'CBE_NT']

hall = json.loads((BASE / 'T47D_CDK4_20260927_results' / 'hallmark_gene_sets_cache.json').read_text())
SETS = {'core RB-E2F': GS.CORE_E2F, 'Hallmark\nE2F targets': hall['E2F Targets'],
        'S-phase\nmarkers': GS.S_SPECIFIC, 'G2/M\nmarkers': GS.G2M_SPECIFIC,
        'Hallmark\nG2M checkpoint': hall['G2-M Checkpoint']}
stars = lambda p: '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < .05 else 'ns'
rng = np.random.default_rng(0)

def score(sy, cols):
    g = [s2i[s] for s in sy if s in s2i.index and lt.loc[s2i[s], cols].mean() > 1]
    return lt.loc[g, cols].mean(axis=0)

def contrast(sy, a, b, cols):
    s = score(sy, cols)
    fc = 2 ** (s[grp(a)].mean() - s[grp(b)].mean())
    p = sps.ttest_ind(s[grp(a)], s[grp(b)], equal_var=False)[1]
    # bootstrap CI over genes, on per-gene differences of group means
    g = [s2i[x] for x in sy if x in s2i.index and lt.loc[s2i[x], cols].mean() > 1]
    d = (lt.loc[g, grp(a)].mean(axis=1) - lt.loc[g, grp(b)].mean(axis=1)).to_numpy()
    bs = np.array([np.mean(rng.choice(d, len(d))) for _ in range(4000)])
    return fc, p, 2 ** np.percentile(bs, 2.5), 2 ** np.percentile(bs, 97.5)

names = list(SETS)
fig, axes = plt.subplots(1, 2, figsize=(11.6, 3.5), gridspec_kw={'width_ratios': [1, 1.5]})

# ── A: what the background alone does ───────────────────────────────────────
ax = axes[0]
rows = []
for i, n in enumerate(names):
    fc, p, lo, hi = contrast(SETS[n], 'CDK4_E56K', 'CBE_NT', ALL)
    rows.append(dict(panel='background', gene_set=n.replace('\n', ' '), fc=fc, p=p, lo=lo, hi=hi))
    c = '#9e9e9e' if p >= .05 else '#d6604d'
    ax.barh(i, fc - 1, left=1, height=.6, color=c, edgecolor='black', lw=.4)
    ax.plot([lo, hi], [i, i], color='black', lw=1)
    ax.text(1.085, i, stars(p), va='center', fontsize=8 if p < .05 else 6.5,
            color='#a33' if p < .05 else '#777')
ax.axvline(1, color='black', lw=.9, ls='--')
ax.axvspan(0.97, 1.03, color='#cfe8cf', alpha=.5, zorder=0)
ax.set_yticks(range(len(names))); ax.set_yticklabels([n.replace('\n', ' ') for n in names], fontsize=8)
ax.invert_yaxis(); ax.set_xlim(0.88, 1.11)
ax.set_xlabel('E56K / CBE NT', fontsize=8.5)
ax.set_title('A. Background drift alone\n(E56K carries it; CDK4 near-WT)', fontsize=9)
ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=7.5)

# ── B: both references agree where it matters ───────────────────────────────
ax = axes[1]
w = 0.19
cols = {'CDK4_D196N': '#4393c3', 'CDK4_W198': '#d6604d'}
lab = {'CDK4_D196N': 'D196N', 'CDK4_W198': 'W198*'}
for k, v in enumerate(['CDK4_D196N', 'CDK4_W198']):
    for j, (ref, refcols, hatch) in enumerate([('CBE NT', ALL, ''), ('E56K', VARONLY, '////')]):
        vals = []
        for n in names:
            b = 'CBE_NT' if ref == 'CBE NT' else 'CDK4_E56K'
            fc, p, lo, hi = contrast(SETS[n], v, b, refcols)
            vals.append(fc)
            rows.append(dict(panel='effect', variant=lab[v], reference=ref,
                             gene_set=n.replace('\n', ' '), fc=fc, p=p, lo=lo, hi=hi))
        off = (k * 2 + j - 1.5) * w
        ax.bar(np.arange(len(names)) + off, vals, width=w, color=cols[v], edgecolor='black', lw=.4,
               hatch=hatch, label=f'{lab[v]} vs {ref}')
ax.axhline(1, color='black', lw=.9, ls='--')
ax.set_xticks(range(len(names))); ax.set_xticklabels(names, fontsize=7.5)
ax.set_ylim(0.8, 1.14)
ax.set_ylabel('abundance relative to reference', fontsize=8.5)
ax.set_title('B. Effect of CDK4 loss under both references', fontsize=9)
ax.legend(fontsize=7, frameon=False, ncol=2, loc='upper left')
ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=7.5)
for i, n in enumerate(names):
    if 'G2/M' in n:
        ax.axvspan(i - .45, i + .45, color='#f8d7d3', alpha=.5, zorder=0)
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'dual_reference_cell_cycle.{e}', dpi=300, bbox_inches='tight')
plt.close()
pd.DataFrame(rows).to_csv(OUT / 'dual_reference_stats.csv', index=False)

d = pd.DataFrame(rows)
print('Background effect (E56K vs CBE NT):')
print(d[d.panel == 'background'][['gene_set', 'fc', 'lo', 'hi', 'p']].round(3).to_string(index=False))
print('\nEffect under both references:')
piv = d[d.panel == 'effect'].pivot_table(index='gene_set', columns=['variant', 'reference'], values='fc')
print(piv.round(3).to_string())
