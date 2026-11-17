"""
NT-referenced E2F / cell-cycle figures with ONE-SIDED statistics.

Pre-specified direction: CDK4 phosphorylates RB and releases E2F, so loss of CDK4
function is expected to LOWER E2F / S-phase output. All tests are lower-tail
(variant < CBE NT), Holm-adjusted across the three variants within each score.
Gene-level stats are recomputed one-sided from the DESeq2 Wald statistic and
re-FDR'd (Benjamini-Hochberg) across all tested genes in that contrast.

Outputs `*_onesided.{png,svg}` in figures/. Two-sided versions from figures_vs_NT.py
are left in place for comparison.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats as sps
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import genesets as GS

mpl.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
OLD  = BASE / 'T47D_CDK4_20260927_results' / 'deseq2_results'
OUT  = HERE / 'results'; OUT.mkdir(exist_ok=True)
FIG  = HERE / 'figures'; FIG.mkdir(exist_ok=True)

CTRL      = 'CBE_NT'
VARIANTS  = ['CDK4_D196N', 'CDK4_E56K', 'CDK4_W198']
GENOTYPES = [CTRL] + VARIANTS
LABEL = {'CBE_NT': 'CBE NT', 'CDK4_D196N': 'D196N', 'CDK4_E56K': 'E56K', 'CDK4_W198': 'W198*'}
COLOR = {'CBE_NT': '#9e9e9e', 'CDK4_D196N': '#4393c3', 'CDK4_E56K': '#66bd63', 'CDK4_W198': '#d6604d'}
FOOT = ('one-sided Welch $t$-test (lower tail, pre-specified from RB-E2F direction), '
        'Holm-adjusted across variants')

tpm = pd.read_csv(BASE / 'salmon.merged.gene_tpm.tsv', sep='\t', index_col=0)
gene_name = tpm.pop('gene_name')
meta = pd.Series(tpm.columns.str.extract(r'^T47D_(.+)_veh_\d+$')[0].values, index=tpm.columns)
samples = [s for g in GENOTYPES for s in meta.index[meta == g]]
tpm, meta = tpm[samples], meta[samples]
log_tpm = np.log2(tpm + 1)
order = tpm.mean(axis=1).sort_values(ascending=False).index
sym2id = pd.Series(order, index=gene_name.reindex(order).values)
sym2id = sym2id[~sym2id.index.duplicated()]
res = {v: pd.read_csv(OLD / f'{v}_vs_{CTRL}_veh.csv').set_index('gene_id') for v in VARIANTS}
grp = lambda g: list(meta.index[meta == g])

stars = lambda p: '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < .05 else 'ns'

def holm(p):
    p = np.asarray(p, float); out = np.empty(len(p)); run = 0.0
    for rank, i in enumerate(np.argsort(p)):
        run = max(run, min(1.0, (len(p) - rank) * p[i])); out[i] = run
    return out

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); run = 1.0
    for k, i in enumerate(o[::-1]):
        run = min(run, p[i] * n / (n - k)); q[i] = run
    return q

hall = json.loads((BASE / 'T47D_CDK4_20260927_results' / 'hallmark_gene_sets_cache.json').read_text())
SETS = {'core_E2F': GS.CORE_E2F, 'Hallmark_E2F_Targets': hall['E2F Targets'],
        'S_specific': GS.S_SPECIFIC, 'G2M_specific': GS.G2M_SPECIFIC,
        'S_Tirosh': GS.TIROSH_S, 'G2M_Tirosh': GS.TIROSH_G2M}

def ids(symbols):
    g = [sym2id[s] for s in symbols if s in sym2id.index]
    return [i for i in g if log_tpm.loc[i].mean() > 1.0]

def zscore(symbols):
    m = log_tpm.loc[ids(symbols)]
    return m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0).mean(axis=0)

def geo(symbols):
    return log_tpm.loc[ids(symbols)].mean(axis=0)

scores = pd.DataFrame({'genotype': meta})
for n, g in SETS.items():
    scores[n + '_z'] = zscore(g)
    scores[n + '_geo'] = geo(g)

# ── one-sided (and two-sided, for the comparison panel) score stats ──────────
rows = []
for n in SETS:
    cz = scores.loc[grp(CTRL), n + '_z']; cg = scores.loc[grp(CTRL), n + '_geo']
    p1s, p2s = [], []
    for v in VARIANTS:
        xz, xg = scores.loc[grp(v), n + '_z'], scores.loc[grp(v), n + '_geo']
        t, p2 = sps.ttest_ind(xz, cz, equal_var=False)
        p1 = p2 / 2 if t < 0 else 1 - p2 / 2
        p1s.append(p1); p2s.append(p2)
        rows.append(dict(score=n, variant=v, fold_change=2 ** (xg.mean() - cg.mean()),
                         delta_z=xz.mean() - cz.mean(), p_one=p1, p_two=p2))
    for r, a1, a2 in zip(rows[-3:], holm(p1s), holm(p2s)):
        r['p_one_holm'], r['p_two_holm'] = a1, a2
stats = pd.DataFrame(rows)
stats.to_csv(OUT / 'program_score_stats_vs_NT_onesided.csv', index=False)
sig1 = {(r.score, r.variant): r.p_one_holm for r in stats.itertuples()}

# ── 1. score strip plots, one-sided stars ───────────────────────────────────
def strip(ax, score, ylab):
    rng = np.random.default_rng(0)
    vals = {g: scores.loc[grp(g), score + '_z'].to_numpy() for g in GENOTYPES}
    for xi, g in enumerate(GENOTYPES):
        ax.scatter(xi + rng.uniform(-.11, .11, 3), vals[g], s=30, facecolor=COLOR[g],
                   edgecolor='black', lw=.5, zorder=3)
        ax.plot([xi - .26, xi + .26], [vals[g].mean()] * 2, color='black', lw=1.3, zorder=4)
    top = max(v.max() for v in vals.values()); span = top - min(v.min() for v in vals.values())
    for k, v in enumerate(VARIANTS):
        xv = GENOTYPES.index(v); yb = top + span * (.10 + .145 * k)
        ax.plot([0, 0, xv, xv], [yb - span * .035, yb, yb, yb - span * .035], color='black', lw=.8)
        s = stars(sig1[(score, v)])
        ax.text(xv / 2, yb, s, ha='center', va='bottom', fontsize=7 if s == 'ns' else 9)
    ax.set_ylim(top=top + span * (.2 + .145 * len(VARIANTS)))
    ax.set_xticks(range(len(GENOTYPES)))
    ax.set_xticklabels([LABEL[g] for g in GENOTYPES], rotation=45, ha='right',
                       rotation_mode='anchor', fontsize=8)
    ax.set_ylabel(ylab, fontsize=8); ax.tick_params(labelsize=7.5)
    ax.spines[['top', 'right']].set_visible(False)

fig, axes = plt.subplots(1, 4, figsize=(10.4, 3.3))
for ax, (s, l) in zip(axes, [('core_E2F', 'Core RB-E2F targets\n(mean z, log$_2$TPM)'),
                             ('Hallmark_E2F_Targets', 'Hallmark E2F Targets\n(mean z)'),
                             ('S_specific', 'S-phase markers\n(mean z)'),
                             ('G2M_specific', 'G2/M markers\n(mean z)')]):
    strip(ax, s, l)
fig.text(0.5, -0.04, FOOT, ha='center', fontsize=6.5, style='italic', color='#555555')
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'program_scores_vs_NT_onesided.{e}', dpi=300, bbox_inches='tight')
plt.close()

# ── 2. phase-marker fold change, one-sided stars ────────────────────────────
rng = np.random.default_rng(0)
fig, ax = plt.subplots(figsize=(4.8, 3.3))
width, gap = 0.26, 0.04
setlist = [('S_specific', 'S-phase markers'), ('G2M_specific', 'G2/M markers')]
for si, (sname, slab) in enumerate(setlist):
    gids = ids(SETS[sname])
    for vi, v in enumerate(VARIANTS):
        g = scores[sname + '_geo']
        fc = 2 ** (g[grp(v)].mean() - g[grp(CTRL)].mean())
        l = res[v].reindex(gids)['log2FoldChange'].dropna().to_numpy()
        boot = np.array([np.median(rng.choice(l, len(l))) for _ in range(4000)])
        lo, hi = 2 ** np.percentile(boot, [2.5, 97.5]); med = 2 ** np.median(l)
        x = si * (3 * (width + gap) + 0.3) + vi * (width + gap)
        ax.bar(x, fc, width=width, color=COLOR[v], edgecolor='black', lw=.5,
               label=LABEL[v] if si == 0 else None)
        ax.errorbar(x, med, yerr=[[med - lo], [hi - med]], fmt='D', color='black',
                    ms=3.2, lw=1, capsize=2.5, zorder=4)
        s = stars(sig1[(sname, v)])
        ax.text(x, max(fc, hi) + .012, s, ha='center', va='bottom', fontsize=7 if s == 'ns' else 9)
ax.axhline(1, color='black', lw=.8, ls='--')
ax.set_xticks([(3 * (width + gap) + 0.3) * si + (width + gap) for si in range(2)])
ax.set_xticklabels([l for _, l in setlist], fontsize=9)
ax.set_ylabel('phase-marker abundance relative to CBE NT\n(bars: geometric mean of gene set)', fontsize=8)
ax.set_ylim(0.8, 1.18)
h, lb = ax.get_legend_handles_labels()
h.append(Line2D([], [], marker='D', color='black', ms=3.2, lw=1, ls='-'))
lb.append('median per-gene FC (95% CI)')
ax.legend(h, lb, fontsize=7, frameon=False, loc='upper left', ncol=2, columnspacing=1.0)
ax.tick_params(labelsize=8); ax.spines[['top', 'right']].set_visible(False)
fig.text(0.5, -0.04, FOOT, ha='center', fontsize=6.5, style='italic', color='#555555')
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'phase_marker_fold_change_vs_NT_onesided.{e}', dpi=300, bbox_inches='tight')
plt.close()

# ── 3. core E2F heatmap, one-sided gene-level FDR ───────────────────────────
qmap = {}
for v in VARIANTS:
    d = res[v].dropna(subset=['padj', 'stat']).copy()
    d['p_one'] = sps.norm.cdf(d['stat'])          # lower tail from the Wald statistic
    d['q_one'] = bh(d['p_one'].to_numpy())
    qmap[v] = d['q_one']
rows = []
for sym in GS.CORE_E2F:
    if sym not in sym2id.index:
        continue
    gid = sym2id[sym]; r = {'gene': sym}
    for v in VARIANTS:
        if gid in res[v].index:
            r[v] = res[v].loc[gid, 'log2FoldChange']
            r[v + '_q_onesided'] = qmap[v].get(gid, np.nan)
            r[v + '_padj_twosided'] = res[v].loc[gid, 'padj']
    rows.append(r)
core = pd.DataFrame(rows).set_index('gene')
core.to_csv(OUT / 'core_E2F_target_log2FC_vs_NT_onesided.csv')
M = core[VARIANTS].astype(float)
M = M.loc[M.mean(axis=1).sort_values().index]

fig, ax = plt.subplots(figsize=(3.1, 6.0))
vmax = np.ceil(np.nanmax(np.abs(M.values)) * 10) / 10
im = ax.imshow(M.values, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
for i, gene in enumerate(M.index):
    for j, v in enumerate(VARIANTS):
        q = core.loc[gene, v + '_q_onesided']
        if pd.notna(q) and q < .05:
            ax.text(j, i, stars(q), ha='center', va='center', fontsize=7,
                    color='white' if abs(M.values[i, j]) > vmax * .6 else 'black')
ax.set_xticks(range(len(VARIANTS)))
ax.set_xticklabels([LABEL[v] for v in VARIANTS], rotation=45, ha='right', rotation_mode='anchor', fontsize=8)
ax.set_yticks(range(len(M))); ax.set_yticklabels(M.index, fontsize=7.5, style='italic')
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, shrink=.45, pad=.04, aspect=14)
cb.set_label('log$_2$FC vs CBE NT', fontsize=8); cb.ax.tick_params(labelsize=7)
fig.text(0.5, -0.02, 'one-sided Wald test (lower tail),\nBenjamini-Hochberg across all tested genes',
         ha='center', fontsize=6.5, style='italic', color='#555555')
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'core_E2F_targets_vs_NT_onesided.{e}', dpi=300, bbox_inches='tight')
plt.close()

# ── 4. what the tail choice actually changes ────────────────────────────────
show = ['core_E2F', 'Hallmark_E2F_Targets', 'S_specific', 'G2M_specific']
d = stats[stats.score.isin(show)].copy()
d['key'] = d.score.map({'core_E2F': 'core RB-E2F', 'Hallmark_E2F_Targets': 'Hallmark E2F',
                        'S_specific': 'S-phase', 'G2M_specific': 'G2/M'}) + '  ·  ' + d.variant.map(LABEL)
d = d.iloc[::-1].reset_index(drop=True)
fig, ax = plt.subplots(figsize=(5.6, 4.6))
lp = lambda p: -np.log10(np.clip(p, 1e-6, 1))
for i, r in d.iterrows():
    ax.plot([lp(r.p_two_holm), lp(r.p_one_holm)], [i, i], color='#bbbbbb', lw=1, zorder=1)
    ax.scatter(lp(r.p_two_holm), i, s=34, facecolor='white', edgecolor='black', lw=.9, zorder=3)
    ax.scatter(lp(r.p_one_holm), i, s=34, facecolor=COLOR[r.variant], edgecolor='black', lw=.5, zorder=3)
ax.axvline(-np.log10(.05), color='#d6604d', lw=1, ls='--')
ax.text(-np.log10(.05), len(d) - .3, ' p = 0.05', color='#d6604d', fontsize=7, va='top')
ax.set_yticks(range(len(d))); ax.set_yticklabels(d.key, fontsize=7.5)
ax.set_xlabel('$-\\log_{10}$(Holm-adjusted p)', fontsize=8.5)
ax.tick_params(labelsize=8); ax.spines[['top', 'right']].set_visible(False)
ax.legend(handles=[Line2D([], [], marker='o', ls='', mfc='white', mec='black', ms=6, label='two-sided'),
                   Line2D([], [], marker='o', ls='', mfc='#777777', mec='black', ms=6, label='one-sided (lower tail)')],
          fontsize=7.5, frameon=False, loc='lower right')
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'onesided_vs_twosided_comparison.{e}', dpi=300, bbox_inches='tight')
plt.close()

print('score-level, one-sided Holm (vs CBE NT):\n')
print(f"{'score':24s} " + ' '.join(f'{LABEL[v]:>16s}' for v in VARIANTS))
for n in SETS:
    s = stats[stats.score == n].set_index('variant')
    print(f'{n:24s} ' + ' '.join(f'{s.loc[v, "fold_change"]:7.2f} {stars(s.loc[v, "p_one_holm"]):<8s}'
                                 for v in VARIANTS))
print('\ncore RB-E2F targets passing one-sided BH q<0.05 (two-sided padj<0.05 in brackets):')
for v in VARIANTS:
    a = int((core[v + '_q_onesided'] < .05).sum()); b = int((core[v + '_padj_twosided'] < .05).sum())
    dn = int(((core[v + '_q_onesided'] < .05) & (core[v] < 0)).sum())
    print(f'  {LABEL[v]:7s} {a:2d} of {len(core)}  [{b}]   of which down: {dn}')
print('\nwrote figures/*_onesided.{png,svg}')
