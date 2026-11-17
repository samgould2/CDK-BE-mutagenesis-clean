"""
Compact 2 x 2 summary figure, emitted at two aspect ratios.

  a  CDK4 transcript level (NMD of the W198* allele)
  b  Hallmark E2F targets score
  c  volcano, W198* vs CBE NT
  d  volcano, W198* vs D196N  (matched background)

Both volcanoes use contrasts refit with DESeq2 independent filtering disabled, so the
two panels share an identical gene universe and their DEG counts are comparable.
Order throughout: CBE NT, E56K, D196N, W198*.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats as sps
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': 0.7,
    'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
    'xtick.major.size': 2.4, 'ytick.major.size': 2.4,
})

HERE = Path(__file__).resolve().parent
BASE, OUT, FIG = HERE.parent, HERE / 'results', HERE / 'figures'

BASEMEAN_MIN = 50
LFC_T, PADJ_T = 1.0, 0.05
X_LIM = 6

CTRL = 'CBE_NT'
VARIANTS = ['CDK4_E56K', 'CDK4_D196N', 'CDK4_W198']
GENOTYPES = [CTRL] + VARIANTS
LABEL = {'CBE_NT': 'CBE NT', 'CDK4_E56K': 'E56K', 'CDK4_D196N': 'D196N', 'CDK4_W198': 'W198*'}
COLOR = {'CBE_NT': '#9E9E9E', 'CDK4_E56K': '#92C5DE', 'CDK4_D196N': '#2166AC', 'CDK4_W198': '#B2182B'}
UP, DOWN, NS = '#B2182B', '#2166AC', '#D9D9D9'

# (name, figsize, bracket offset, bracket step, headroom, volcano y-ticks, panel-letter x)
LAYOUTS = [('pub_fig_square_4panel',      (6.6, 6.4), .13, .17, .78, 6, -.20),
           ('pub_fig_square_4panel_mid',  (7.6, 5.6), .11, .15, .68, 5, -.17),
           ('pub_fig_square_4panel_wide', (8.8, 4.9), .10, .13, .58, 4, -.13)]

tpm = pd.read_csv(BASE / 'salmon.merged.gene_tpm.tsv', sep='\t', index_col=0)
gn = tpm.pop('gene_name')
meta = pd.Series(tpm.columns.str.extract(r'^T47D_(.+)_veh_\d+$')[0].values, index=tpm.columns)
samples = [s for g in GENOTYPES for s in meta.index[meta == g]]
tpm, meta = tpm[samples], meta[samples]
log_tpm = np.log2(tpm + 1)
order = tpm.mean(axis=1).sort_values(ascending=False).index
s2i = pd.Series(order, index=gn.reindex(order).values); s2i = s2i[~s2i.index.duplicated()]
grp = lambda g: list(meta.index[meta == g])
res_nt = {v: pd.read_csv(OUT / f'{v}_vs_{CTRL}_veh.csv').set_index('gene_id') for v in VARIANTS}
volc_nt = pd.read_csv(OUT / 'CDK4_W198_vs_CBE_NT_veh_nofilt.csv').set_index('gene_id')
volc_wd = pd.read_csv(OUT / 'CDK4_W198_vs_CDK4_D196N_veh_nofilt.csv').set_index('gene_id')
hall = json.loads((BASE / 'T47D_CDK4_20260927_results' / 'hallmark_gene_sets_cache.json').read_text())

stars = lambda p: '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < .05 else 'n.s.'

def holm(p):
    p = np.asarray(p, float); o = np.empty(len(p)); r = 0.
    for k, i in enumerate(np.argsort(p)):
        r = max(r, min(1., (len(p) - k) * p[i])); o[i] = r
    return o

E2F_IDS = [s2i[s] for s in hall['E2F Targets'] if s in s2i.index and log_tpm.loc[s2i[s]].mean() > 1]
_m = log_tpm.loc[E2F_IDS]
E2F_SCORE = _m.sub(_m.mean(axis=1), axis=0).div(_m.std(axis=1).replace(0, np.nan), axis=0).mean(axis=0)
_ps = []
for v in VARIANTS:
    t, p2 = sps.ttest_ind(E2F_SCORE[grp(v)], E2F_SCORE[grp(CTRL)], equal_var=False)
    _ps.append(p2 / 2 if t < 0 else 1 - p2 / 2)
E2F_ADJ = holm(_ps)


def volcano(ax, r, title, label_genes, yn):
    d = r.dropna(subset=['padj', 'log2FoldChange'])
    d = d[d.baseMean >= BASEMEAN_MIN].copy()
    d['y'] = -np.log10(d['padj'].clip(lower=1e-300))
    d['x'] = d['log2FoldChange'].clip(-X_LIM * .98, X_LIM * .98)
    up = (d.padj < PADJ_T) & (d.log2FoldChange > LFC_T)
    dn = (d.padj < PADJ_T) & (d.log2FoldChange < -LFC_T)
    ax.scatter(d.loc[~(up | dn), 'x'], d.loc[~(up | dn), 'y'], s=2.5, color=NS, lw=0, rasterized=True)
    ax.scatter(d.loc[dn, 'x'], d.loc[dn, 'y'], s=5, color=DOWN, lw=0, rasterized=True)
    ax.scatter(d.loc[up, 'x'], d.loc[up, 'y'], s=5, color=UP, lw=0, rasterized=True)
    for xv in (LFC_T, -LFC_T):
        ax.axvline(xv, color='#999999', lw=.7, ls='--')
    ax.axhline(-np.log10(PADJ_T), color='#999999', lw=.7, ls='--')
    for sym in label_genes:
        hit = d[d.gene_name == sym]
        if not len(hit):
            continue
        rr = hit.iloc[0]
        ax.scatter(rr['x'], rr['y'], s=26, facecolor='none', edgecolor='black', lw=.9, zorder=5)
        ax.annotate(sym, (rr['x'], rr['y']), fontsize=6.5, style='italic', zorder=6,
                    xytext=(4, 3), textcoords='offset points')
    ax.set_xlim(-X_LIM, X_LIM)
    ax.set_xlabel('log$_2$ fold change', fontsize=8.5)
    ax.set_ylabel('$-$log$_{10}$ $P_{adj}$', fontsize=8.5)
    ax.set_title(title, fontsize=8.5, pad=3)
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(yn))
    ax.tick_params(labelsize=7.5); ax.spines[['top', 'right']].set_visible(False)
    ax.text(.03, .055, f'{int(dn.sum())} down', transform=ax.transAxes, fontsize=7,
            color=DOWN, va='bottom', ha='left')
    ax.text(.97, .055, f'{int(up.sum())} up', transform=ax.transAxes, fontsize=7,
            color=UP, va='bottom', ha='right')
    return int(up.sum()), int(dn.sum()), len(d)


def build(figsize, br0, brs, head, yn, lx):
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    for ax, L in zip(axes.flat, 'abcd'):
        ax.text(lx, 1.07, L, transform=ax.transAxes, fontsize=11, fontweight='bold', va='top')

    # (a) CDK4 transcript
    ax = axes[0, 0]
    g = s2i['CDK4']
    rng = np.random.default_rng(1)
    for xi, geno in enumerate(GENOTYPES):
        v = tpm.loc[g, grp(geno)].to_numpy()
        ax.bar(xi, v.mean(), width=.64, color=COLOR[geno], edgecolor='black', lw=.6, zorder=2)
        ax.errorbar(xi, v.mean(), yerr=v.std(), color='black', lw=.7, capsize=2.2, zorder=3)
        ax.scatter(xi + rng.uniform(-.1, .1, 3), v, s=9, color='black', zorder=4, lw=0)
    for xi, geno in enumerate(VARIANTS, start=1):
        p = res_nt[geno].loc[g, 'padj']
        ax.text(xi, tpm.loc[g, grp(geno)].mean() + tpm.loc[g, grp(geno)].std() + 16,
                stars(p), ha='center', fontsize=8 if p < .05 else 6.2)
    ax.set_ylim(0, 480)
    ax.set_ylabel('CDK4 (TPM)', fontsize=9)
    ax.set_xticks(range(4))
    ax.set_xticklabels([LABEL[x] for x in GENOTYPES], rotation=45, ha='right',
                       rotation_mode='anchor', fontsize=8)
    ax.tick_params(labelsize=7.5); ax.spines[['top', 'right']].set_visible(False)

    # (b) Hallmark E2F targets
    ax = axes[0, 1]
    rng = np.random.default_rng(0)
    vals = {geno: E2F_SCORE[grp(geno)].to_numpy() for geno in GENOTYPES}
    for xi, geno in enumerate(GENOTYPES):
        ax.scatter(xi + rng.uniform(-.1, .1, 3), vals[geno], s=28, facecolor=COLOR[geno],
                   edgecolor='black', lw=.5, zorder=3, clip_on=False)
        ax.plot([xi - .27, xi + .27], [vals[geno].mean()] * 2, color='black', lw=1.3, zorder=4)
    top = max(v.max() for v in vals.values())
    span = top - min(v.min() for v in vals.values())
    for k, v in enumerate(VARIANTS):
        xv = GENOTYPES.index(v)
        yb = top + span * (br0 + brs * k)
        ax.plot([0, 0, xv, xv], [yb - span * .04, yb, yb, yb - span * .04], color='black', lw=.7)
        s = stars(E2F_ADJ[k])
        ax.text(xv / 2, yb + span * .01, s, ha='center', va='bottom',
                fontsize=6.4 if s == 'n.s.' else 8)
    ax.set_ylim(top=top + span * head)
    ax.set_ylabel(f'Hallmark E2F targets\n(mean $z$-score, $n$ = {len(E2F_IDS)})', fontsize=8.5)
    ax.set_xticks(range(4))
    ax.set_xticklabels([LABEL[x] for x in GENOTYPES], rotation=45, ha='right',
                       rotation_mode='anchor', fontsize=8)
    ax.tick_params(labelsize=7.5); ax.spines[['top', 'right']].set_visible(False)

    # (c, d) volcanoes
    u1, d1, n1 = volcano(axes[1, 0], volc_nt, 'W198* vs CBE NT', ['CDK4'], yn)
    u2, d2, n2 = volcano(axes[1, 1], volc_wd, 'W198* vs D196N', ['CDK4', 'PIP', 'MGP'], yn)

    plt.tight_layout()
    return fig, (u1, d1, n1, u2, d2, n2)


for name, figsize, br0, brs, head, yn, lx in LAYOUTS:
    fig, counts = build(figsize, br0, brs, head, yn, lx)
    for e in ('png', 'svg', 'pdf'):
        fig.savefig(FIG / f'{name}.{e}', dpi=400, bbox_inches='tight')
    plt.close(fig)
    print(f'{name:30s} {figsize[0]} x {figsize[1]} in')
u1, d1, n1, u2, d2, n2 = counts
print(f'\nbaseMean >= {BASEMEAN_MIN}, identical universe ({n1} genes):')
print(f'  W198* vs CBE NT   {u1} up / {d1} down')
print(f'  W198* vs D196N    {u2} up / {d2} down')
