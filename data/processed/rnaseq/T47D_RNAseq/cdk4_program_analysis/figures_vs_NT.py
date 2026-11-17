"""
E2F / cell-cycle figures referenced to the CBE NT control (rather than to the
missense lines as in analyze.py).

Outputs `*_vs_NT.{png,svg}` in figures/ and the matching stats in results/.
See FINDINGS.md section 0: CBE NT is not isogenic with the variant lines, so these
contrasts carry the cell-line background in addition to the CDK4 genotype effect.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats as sps
import matplotlib as mpl
import matplotlib.pyplot as plt

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

hall = json.loads((BASE / 'T47D_CDK4_20260927_results' / 'hallmark_gene_sets_cache.json').read_text())
SETS = {'core_E2F': GS.CORE_E2F, 'Hallmark_E2F_Targets': hall['E2F Targets'],
        'S_specific': GS.S_SPECIFIC, 'G2M_specific': GS.G2M_SPECIFIC,
        'S_Tirosh': GS.TIROSH_S, 'G2M_Tirosh': GS.TIROSH_G2M,
        'Hallmark_G2M_Checkpoint': hall['G2-M Checkpoint']}

def ids(symbols, min_log_tpm=1.0):
    g = [sym2id[s] for s in symbols if s in sym2id.index]
    return [i for i in g if log_tpm.loc[i].mean() > min_log_tpm]

def zscore(symbols):
    m = log_tpm.loc[ids(symbols)]
    return m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0).mean(axis=0)

def geo(symbols):
    return log_tpm.loc[ids(symbols)].mean(axis=0)

scores = pd.DataFrame({'genotype': meta})
for n, g in SETS.items():
    scores[n + '_z'] = zscore(g)
    scores[n + '_geo'] = geo(g)
scores.to_csv(OUT / 'cell_cycle_scores_vs_NT.csv')

# ── stats: each variant vs CBE NT, Holm over the three variants, per score ────
rows = []
for n in SETS:
    c_z, c_g = scores.loc[grp(CTRL), n + '_z'], scores.loc[grp(CTRL), n + '_geo']
    ps = []
    for v in VARIANTS:
        x_z, x_g = scores.loc[grp(v), n + '_z'], scores.loc[grp(v), n + '_geo']
        p = sps.ttest_ind(x_z, c_z, equal_var=False)[1]
        ps.append(p)
        rows.append(dict(score=n, variant=v, ctrl_mean_z=c_z.mean(), variant_mean_z=x_z.mean(),
                         delta_z=x_z.mean() - c_z.mean(),
                         fold_change=2 ** (x_g.mean() - c_g.mean()), p=p))
    for r, pa in zip(rows[-3:], holm(ps)):
        r['p_holm'] = pa
stats = pd.DataFrame(rows)
stats.to_csv(OUT / 'program_score_stats_vs_NT.csv', index=False)
sig = {(r.score, r.variant): r.p_holm for r in stats.itertuples()}

# ── 1. strip plots of the scores, all four genotypes ─────────────────────────
def strip(ax, score, ylab):
    rng = np.random.default_rng(0)
    vals = {g: scores.loc[grp(g), score + '_z'].to_numpy() for g in GENOTYPES}
    for xi, g in enumerate(GENOTYPES):
        ax.scatter(xi + rng.uniform(-.11, .11, 3), vals[g], s=30, facecolor=COLOR[g],
                   edgecolor='black', lw=.5, zorder=3)
        ax.plot([xi - .26, xi + .26], [vals[g].mean()] * 2, color='black', lw=1.3, zorder=4)
    top = max(v.max() for v in vals.values())
    span = top - min(v.min() for v in vals.values())
    for k, v in enumerate(VARIANTS):
        xv = GENOTYPES.index(v); yb = top + span * (.10 + .145 * k)
        ax.plot([0, 0, xv, xv], [yb - span * .035, yb, yb, yb - span * .035], color='black', lw=.8)
        s = stars(sig[(score, v)])
        ax.text(xv / 2, yb, s, ha='center', va='bottom', fontsize=7 if s == 'ns' else 9)
    ax.set_ylim(top=top + span * (.2 + .145 * len(VARIANTS)))
    ax.set_xticks(range(len(GENOTYPES)))
    ax.set_xticklabels([LABEL[g] for g in GENOTYPES], rotation=45, ha='right',
                       rotation_mode='anchor', fontsize=8)
    ax.set_ylabel(ylab, fontsize=8); ax.tick_params(labelsize=7.5)
    ax.spines[['top', 'right']].set_visible(False)

panels = [('core_E2F', 'Core RB-E2F targets\n(mean z, log$_2$TPM)'),
          ('Hallmark_E2F_Targets', 'Hallmark E2F Targets\n(mean z)'),
          ('S_specific', 'S-phase markers\n(mean z)'),
          ('G2M_specific', 'G2/M markers\n(mean z)')]
fig, axes = plt.subplots(1, 4, figsize=(10.4, 3.1))
for ax, (s, l) in zip(axes, panels):
    strip(ax, s, l)
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'program_scores_vs_NT.{e}', dpi=300, bbox_inches='tight')
plt.close()

# ── 2. fold change in phase-marker abundance vs NT ───────────────────────────
rng = np.random.default_rng(0)
fig, ax = plt.subplots(figsize=(4.6, 3.1))
width, gap = 0.26, 0.04
setlist = [('S_specific', 'S-phase markers'), ('G2M_specific', 'G2/M markers')]
for si, (sname, slab) in enumerate(setlist):
    gids = [i for i in ids(SETS[sname])]
    for vi, v in enumerate(VARIANTS):
        g = scores[sname + '_geo']
        fc = 2 ** (g[grp(v)].mean() - g[grp(CTRL)].mean())
        l = res[v].reindex(gids)['log2FoldChange'].dropna().to_numpy()
        boot = [np.median(rng.choice(l, len(l))) for _ in range(4000)]
        lo, hi = 2 ** np.percentile(boot, [2.5, 97.5])
        med = 2 ** np.median(l)
        x = si * (3 * (width + gap) + 0.3) + vi * (width + gap)
        ax.bar(x, fc, width=width, color=COLOR[v], edgecolor='black', lw=.5,
               label=LABEL[v] if si == 0 else None)
        ax.errorbar(x, med, yerr=[[med - lo], [hi - med]], fmt='D', color='black',
                    ms=3.2, lw=1, capsize=2.5, zorder=4)
        s = stars(sig[(sname, v)])
        ax.text(x, max(fc, hi) + .012, s, ha='center', va='bottom', fontsize=7 if s == 'ns' else 9)
ax.axhline(1, color='black', lw=.8, ls='--')
ax.set_xticks([(3 * (width + gap) + 0.3) * si + (width + gap) for si in range(2)])
ax.set_xticklabels([l for _, l in setlist], fontsize=9)
ax.set_ylabel('phase-marker abundance relative to CBE NT\n(bars: geometric mean of gene set)', fontsize=8)
ax.set_ylim(0.8, 1.18)
from matplotlib.lines import Line2D
h, lb = ax.get_legend_handles_labels()
h.append(Line2D([], [], marker='D', color='black', ms=3.2, lw=1, ls='-'))
lb.append('median per-gene FC (95% CI)')
ax.legend(h, lb, fontsize=7, frameon=False, loc='upper left', ncol=2, columnspacing=1.0)
ax.tick_params(labelsize=8); ax.spines[['top', 'right']].set_visible(False)
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'phase_marker_fold_change_vs_NT.{e}', dpi=300, bbox_inches='tight')
plt.close()

# ── 3. estimated phase proportions anchored on CBE NT ────────────────────────
ANCHOR = dict(G1=0.55, S=0.27, G2M=0.18)
prop = {}
for g in GENOTYPES:
    fS = 2 ** (scores.loc[grp(g), 'S_specific_geo'].mean() - scores.loc[grp(CTRL), 'S_specific_geo'].mean())
    fM = 2 ** (scores.loc[grp(g), 'G2M_specific_geo'].mean() - scores.loc[grp(CTRL), 'G2M_specific_geo'].mean())
    s, m = ANCHOR['S'] * fS, ANCHOR['G2M'] * fM
    g1 = ANCHOR['G1'] + (ANCHOR['S'] - s) + (ANCHOR['G2M'] - m)
    t = s + m + g1
    prop[g] = dict(G1=g1 / t, S=s / t, G2M=m / t)
prop = pd.DataFrame(prop).T
prop.to_csv(OUT / 'estimated_phase_proportions_vs_NT.csv')

fig, ax = plt.subplots(figsize=(3.7, 3.1))
bottom = np.zeros(len(GENOTYPES))
for phase, c in [('G1', '#cfd8dc'), ('S', '#4393c3'), ('G2M', '#d6604d')]:
    v = prop.loc[GENOTYPES, phase].values
    ax.bar([LABEL[g] for g in GENOTYPES], v, bottom=bottom, color=c, edgecolor='black', lw=.5, label=phase)
    for xi, (a, b) in enumerate(zip(v, bottom)):
        ax.text(xi, b + a / 2, f'{a:.0%}', ha='center', va='center', fontsize=7)
    bottom += v
ax.set_ylabel('estimated fraction of cells', fontsize=8.5)
ax.set_title('CBE NT anchored at 55 / 27 / 18', fontsize=8)
ax.legend(fontsize=7.5, frameon=False, bbox_to_anchor=(1.0, .5), loc='center left')
ax.tick_params(labelsize=8); ax.spines[['top', 'right']].set_visible(False)
plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'estimated_phase_proportions_vs_NT.{e}', dpi=300, bbox_inches='tight')
plt.close()

# ── 4. core RB-E2F target heatmap, log2FC vs CBE NT ──────────────────────────
rows = []
for sym in GS.CORE_E2F:
    if sym not in sym2id.index:
        continue
    gid = sym2id[sym]
    r = {'gene': sym}
    for v in VARIANTS:
        if gid in res[v].index:
            r[v] = res[v].loc[gid, 'log2FoldChange']
            r[v + '_padj'] = res[v].loc[gid, 'padj']
    rows.append(r)
core = pd.DataFrame(rows).set_index('gene')
core.to_csv(OUT / 'core_E2F_target_log2FC_vs_NT.csv')
M = core[VARIANTS].astype(float)
M = M.loc[M.mean(axis=1).sort_values().index]

fig, ax = plt.subplots(figsize=(3.0, 6.0))
vmax = np.ceil(np.nanmax(np.abs(M.values)) * 10) / 10
im = ax.imshow(M.values, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
for i, gene in enumerate(M.index):
    for j, v in enumerate(VARIANTS):
        p = core.loc[gene, v + '_padj']
        if pd.notna(p) and p < .05:
            ax.text(j, i, stars(p), ha='center', va='center', fontsize=7,
                    color='white' if abs(M.values[i, j]) > vmax * .6 else 'black')
ax.set_xticks(range(len(VARIANTS)))
ax.set_xticklabels([LABEL[v] for v in VARIANTS], rotation=45, ha='right', rotation_mode='anchor', fontsize=8)
ax.set_yticks(range(len(M))); ax.set_yticklabels(M.index, fontsize=7.5, style='italic')
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, shrink=.45, pad=.04, aspect=14)
cb.set_label('log$_2$FC vs CBE NT', fontsize=8); cb.ax.tick_params(labelsize=7)
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'core_E2F_targets_vs_NT.{e}', dpi=300, bbox_inches='tight')
plt.close()

# ── console summary ──────────────────────────────────────────────────────────
print('Scores vs CBE NT (mean z over 12 samples; FC = geometric-mean ratio of the set):\n')
print(f"{'score':24s} " + ' '.join(f'{LABEL[g]:>14s}' for g in GENOTYPES))
for n in ['core_E2F', 'Hallmark_E2F_Targets', 'S_specific', 'S_Tirosh',
          'G2M_specific', 'G2M_Tirosh', 'Hallmark_G2M_Checkpoint']:
    cells = [f"{scores.loc[grp(CTRL), n + '_z'].mean():+14.2f}"]
    for v in VARIANTS:
        r = stats.query('score == @n and variant == @v').iloc[0]
        cells.append(f"{r.variant_mean_z:+6.2f} {r.fold_change:.2f}{stars(r.p_holm):>3s}")
    print(f'{n:24s} ' + ' '.join(cells))
print('\nEstimated phase proportions (CBE NT anchored at 55/27/18):')
print((prop.loc[GENOTYPES] * 100).round(1).to_string())
print('\nwrote figures/*_vs_NT.{png,svg} and results/*_vs_NT.csv')
