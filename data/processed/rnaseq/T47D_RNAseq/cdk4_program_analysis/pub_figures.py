"""
Publication figures, referenced to the CBE NT control.

GENE SETS: published only.
  MSigDB Hallmark 2020               (Liberzon et al. 2015)   - via Enrichr, cached
  WHITFIELD_CELL_CYCLE_*             (Whitfield et al. 2002)  - MSigDB C2 CGP
  REACTOME_*                         - MSigDB C2 CP:REACTOME
  FISCHER_DREAM_TARGETS              (Fischer et al. 2016)    - MSigDB C2 CGP
  ISHIDA_E2F_TARGETS                 (Ishida et al. 2001)     - MSigDB C2 CGP
GMTs are cached in results/genesets_published/. The hand-curated panels in
genesets.py are exploratory only and are NOT used here.

Palette: nonsense allele (W198*) red; missense alleles blue, E56K the lighter shade.
Order: CBE NT, E56K, D196N, W198*  (increasing severity of CDK4 disruption).
"""
from pathlib import Path
import glob
import json
import numpy as np
import pandas as pd
from scipy import stats as sps
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

mpl.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': 0.7,
    'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
    'xtick.major.size': 2.6, 'ytick.major.size': 2.6,
})

HERE = Path(__file__).resolve().parent
BASE, OUT, FIG = HERE.parent, HERE / 'results', HERE / 'figures'

CTRL = 'CBE_NT'
VARIANTS = ['CDK4_E56K', 'CDK4_D196N', 'CDK4_W198']        # requested order
GENOTYPES = [CTRL] + VARIANTS
LABEL = {'CBE_NT': 'CBE NT', 'CDK4_E56K': 'E56K', 'CDK4_D196N': 'D196N', 'CDK4_W198': 'W198*'}
COLOR = {'CBE_NT': '#9E9E9E', 'CDK4_E56K': '#92C5DE', 'CDK4_D196N': '#2166AC', 'CDK4_W198': '#B2182B'}

tpm = pd.read_csv(BASE / 'salmon.merged.gene_tpm.tsv', sep='\t', index_col=0)
gn = tpm.pop('gene_name')
meta = pd.Series(tpm.columns.str.extract(r'^T47D_(.+)_veh_\d+$')[0].values, index=tpm.columns)
samples = [s for g in GENOTYPES for s in meta.index[meta == g]]
tpm, meta = tpm[samples], meta[samples]
log_tpm = np.log2(tpm + 1)
order = tpm.mean(axis=1).sort_values(ascending=False).index
s2i = pd.Series(order, index=gn.reindex(order).values); s2i = s2i[~s2i.index.duplicated()]
grp = lambda g: list(meta.index[meta == g])
res = {v: pd.read_csv(OUT / f'{v}_vs_{CTRL}_veh.csv').set_index('gene_id') for v in VARIANTS}

hall = json.loads((BASE / 'T47D_CDK4_20260927_results' / 'hallmark_gene_sets_cache.json').read_text())
SETS = {'Hallmark E2F targets': hall['E2F Targets'],
        'Hallmark G2-M checkpoint': hall['G2-M Checkpoint']}
for f in sorted(glob.glob(str(OUT / 'genesets_published' / '*.gmt'))):
    p = open(f).read().strip().split('\t')
    SETS[p[0]] = [g for g in p[2:] if g]

def holm(p):
    p = np.asarray(p, float); o = np.empty(len(p)); r = 0.
    for k, i in enumerate(np.argsort(p)):
        r = max(r, min(1., (len(p) - k) * p[i])); o[i] = r
    return o

stars = lambda p: '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < .05 else 'n.s.'

def members(name):
    return [s2i[s] for s in SETS[name] if s in s2i.index and log_tpm.loc[s2i[s]].mean() > 1]

def zscore(name):
    m = log_tpm.loc[members(name)]
    return m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0).mean(axis=0)

def onesided(name):
    """Welch, lower tail (pre-specified from CDK4/RB/E2F direction), Holm over the 3 variants."""
    sc = zscore(name)
    ps = []
    for v in VARIANTS:
        t, p2 = sps.ttest_ind(sc[grp(v)], sc[grp(CTRL)], equal_var=False)
        ps.append(p2 / 2 if t < 0 else 1 - p2 / 2)
    return sc, holm(ps)

FOOT = ('one-sided Welch $t$-test (lower tail, pre-specified from the CDK4/RB/E2F direction), '
        'Holm-adjusted across the three variants; CDK4 panel: DESeq2 two-sided $P_{adj}$')

# ═══════════════ Figure 1 ════════════════════════════════════════════════════
PANELS = [('Hallmark E2F targets', 'Hallmark\nE2F targets'),
          ('WHITFIELD_CELL_CYCLE_G1_S', 'Whitfield\nG1/S'),
          ('WHITFIELD_CELL_CYCLE_S', 'Whitfield\nS phase'),
          ('WHITFIELD_CELL_CYCLE_G2_M', 'Whitfield\nG2/M')]

fig = plt.figure(figsize=(11.4, 2.9))
gs = fig.add_gridspec(1, 5, wspace=0.60, left=.05, right=.99, bottom=.30, top=.76)

ax = fig.add_subplot(gs[0])
g = s2i['CDK4']
rng = np.random.default_rng(1)
for xi, geno in enumerate(GENOTYPES):
    v = tpm.loc[g, grp(geno)].to_numpy()
    ax.bar(xi, v.mean(), width=.62, color=COLOR[geno], edgecolor='black', lw=.6, zorder=2)
    ax.errorbar(xi, v.mean(), yerr=v.std(), color='black', lw=.7, capsize=2.2, zorder=3)
    ax.scatter(xi + rng.uniform(-.1, .1, 3), v, s=9, color='black', zorder=4, lw=0)
for xi, geno in enumerate(VARIANTS, start=1):
    p = res[geno].loc[g, 'padj']
    ax.text(xi, tpm.loc[g, grp(geno)].mean() + tpm.loc[g, grp(geno)].std() + 18,
            stars(p), ha='center', fontsize=7.5 if p < .05 else 6)
ax.set_ylim(0, 500); ax.set_ylabel('CDK4 (TPM)', fontsize=8)
ax.set_title('CDK4 transcript', fontsize=8.5, pad=4)
ax.set_xticks(range(4)); ax.set_xticklabels([LABEL[x] for x in GENOTYPES], rotation=45,
                                            ha='right', rotation_mode='anchor', fontsize=7.5)
ax.tick_params(labelsize=7); ax.spines[['top', 'right']].set_visible(False)

stats_rows = []
for pi, (name, title) in enumerate(PANELS):
    ax = fig.add_subplot(gs[pi + 1])
    sc, adj = onesided(name)
    n = len(members(name))
    rng = np.random.default_rng(0)
    vals = {geno: sc[grp(geno)].to_numpy() for geno in GENOTYPES}
    for xi, geno in enumerate(GENOTYPES):
        ax.scatter(xi + rng.uniform(-.1, .1, 3), vals[geno], s=26, facecolor=COLOR[geno],
                   edgecolor='black', lw=.5, zorder=3, clip_on=False)
        ax.plot([xi - .27, xi + .27], [vals[geno].mean()] * 2, color='black', lw=1.3, zorder=4)
    top = max(v.max() for v in vals.values()); span = top - min(v.min() for v in vals.values())
    for k, v in enumerate(VARIANTS):
        xv = GENOTYPES.index(v); yb = top + span * (.13 + .17 * k)
        ax.plot([0, 0, xv, xv], [yb - span * .04, yb, yb, yb - span * .04], color='black', lw=.7)
        s = stars(adj[k])
        ax.text(xv / 2, yb + span * .01, s, ha='center', va='bottom',
                fontsize=6.2 if s == 'n.s.' else 8)
        stats_rows.append(dict(gene_set=name, n_genes=n, variant=LABEL[v], p_one_holm=adj[k]))
    ax.set_ylim(top=top + span * .78)
    ax.set_ylabel('mean $z$-score', fontsize=8)
    ax.set_title(f'{title}\n($n$ = {n} genes)', fontsize=8.5, pad=4)
    ax.set_xticks(range(4)); ax.set_xticklabels([LABEL[x] for x in GENOTYPES], rotation=45,
                                                ha='right', rotation_mode='anchor', fontsize=7.5)
    ax.tick_params(labelsize=7); ax.spines[['top', 'right']].set_visible(False)

fig.text(0.5, -0.02, FOOT, ha='center', fontsize=6.2, style='italic', color='#555555')
for e in ('png', 'svg', 'pdf'):
    fig.savefig(FIG / f'pub_fig1_CDK4_and_programs.{e}', dpi=400, bbox_inches='tight')
plt.close(fig)

# ═══════════════ Figure 2: cyclins / CDKs / CKIs ═════════════════════════════
CLASSES = [('D-type cyclins', ['CCND1', 'CCND2', 'CCND3']),
           ('E-type cyclins', ['CCNE1', 'CCNE2']),
           ('A/B cyclins', ['CCNA2', 'CCNB1', 'CCNB2']),
           ('CDKs', ['CDK4', 'CDK6', 'CDK2', 'CDK1']),
           ('G1/S effectors', ['E2F1', 'CDC25A', 'RB1']),
           ('CDK inhibitors', ['CDKN1A', 'CDKN1B', 'CDKN2D'])]
rows, seps, blocks = [], [], []
for cname, gs_ in CLASSES:
    blocks.append((len(rows) + len(gs_) / 2 - .5, cname)); rows += gs_; seps.append(len(rows) - .5)

fig, ax = plt.subplots(figsize=(4.4, 6.0))
off = {'CDK4_E56K': +.24, 'CDK4_D196N': 0.0, 'CDK4_W198': -.24}
tab = []
for yi, sym in enumerate(rows):
    g = s2i[sym]
    for v in VARIANTS:
        r = res[v].loc[g]
        sig = pd.notna(r['padj']) and r['padj'] < .05
        ax.errorbar(r['log2FoldChange'], yi + off[v], xerr=r['lfcSE'], fmt='none',
                    ecolor='black', elinewidth=.7, capsize=1.6, zorder=2)
        ax.scatter(r['log2FoldChange'], yi + off[v], s=30, zorder=3,
                   facecolor=COLOR[v] if sig else 'white',
                   edgecolor=COLOR[v] if not sig else 'black', lw=.8 if sig else 1.1)
        tab.append(dict(gene=sym, variant=LABEL[v], log2FC=r['log2FoldChange'],
                        lfcSE=r['lfcSE'], padj=r['padj'], significant=bool(sig)))
ax.axvline(0, color='black', lw=.8)
for s in seps[:-1]:
    ax.axhline(s, color='#cccccc', lw=.7)
ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=8, style='italic')
ax.invert_yaxis(); ax.set_xlim(-2.45, 1.35)
ax.set_xlabel('log$_2$ fold change vs CBE NT', fontsize=8.5)
ax.tick_params(labelsize=7.5); ax.spines[['top', 'right']].set_visible(False)
for y, cname in blocks:
    ax.text(1.33, y, cname, fontsize=7, ha='right', va='center', color='#666666', style='italic')
h = [Line2D([], [], marker='o', ls='', mfc=COLOR[v], mec='black', ms=6, label=LABEL[v]) for v in VARIANTS]
h.append(Line2D([], [], marker='o', ls='', mfc='white', mec='#555555', ms=6, label='$P_{adj}\\geq$ 0.05'))
ax.legend(handles=h, fontsize=7, frameon=False, loc='lower left', handletextpad=.3)
plt.tight_layout()
for e in ('png', 'svg', 'pdf'):
    fig.savefig(FIG / f'pub_fig2_cyclins.{e}', dpi=400, bbox_inches='tight')
plt.close(fig)
pd.DataFrame(tab).to_csv(OUT / 'cyclin_panel_vs_NT.csv', index=False)

# ═══════════════ Supplementary: every published set ══════════════════════════
ALL = ['Hallmark E2F targets', 'ISHIDA_E2F_TARGETS', 'FISCHER_DREAM_TARGETS',
       'REACTOME_G1_S_TRANSITION', 'WHITFIELD_CELL_CYCLE_G1_S', 'REACTOME_DNA_REPLICATION',
       'WHITFIELD_CELL_CYCLE_S', 'WHITFIELD_CELL_CYCLE_G2', 'WHITFIELD_CELL_CYCLE_G2_M',
       'Hallmark G2-M checkpoint', 'WHITFIELD_CELL_CYCLE_M_G1']
fig, ax = plt.subplots(figsize=(5.6, 4.6))
sup = []
for yi, name in enumerate(ALL):
    sc, adj = onesided(name)
    n = len(members(name))
    for v, pa in zip(VARIANTS, adj):
        fc = 2 ** (log_tpm.loc[members(name), grp(v)].mean(axis=1).mean()
                   - log_tpm.loc[members(name), grp(CTRL)].mean(axis=1).mean())
        sig = pa < .05
        ax.scatter(fc, yi + {'CDK4_E56K': .24, 'CDK4_D196N': 0., 'CDK4_W198': -.24}[v], s=34,
                   facecolor=COLOR[v] if sig else 'white',
                   edgecolor=COLOR[v] if not sig else 'black', lw=.8 if sig else 1.1, zorder=3)
        sup.append(dict(gene_set=name, n_genes=n, variant=LABEL[v], fold_change=fc, p_one_holm=pa))
ax.axvline(1, color='black', lw=.8)
ax.set_yticks(range(len(ALL)))
ax.set_yticklabels([f'{a.replace("_"," ").replace("WHITFIELD CELL CYCLE","Whitfield").replace("REACTOME","Reactome")}'
                    f'  ({len(members(a))})' for a in ALL], fontsize=7.5)
ax.invert_yaxis()
ax.set_xlabel('abundance relative to CBE NT (geometric mean of set)', fontsize=8.5)
ax.tick_params(labelsize=7.5); ax.spines[['top', 'right']].set_visible(False)
ax.legend(handles=h, fontsize=7, frameon=False, loc='lower left', handletextpad=.3)
fig.text(0.5, -0.03, FOOT.split(';')[0], ha='center', fontsize=6.2, style='italic', color='#555555')
plt.tight_layout()
for e in ('png', 'svg', 'pdf'):
    fig.savefig(FIG / f'pub_figS_all_published_sets.{e}', dpi=400, bbox_inches='tight')
plt.close(fig)
pd.DataFrame(sup).to_csv(OUT / 'published_geneset_stats.csv', index=False)
pd.DataFrame(stats_rows).to_csv(OUT / 'pub_fig1_onesided_stats.csv', index=False)

print('Fig 1 (published sets, one-sided Holm):')
print(pd.DataFrame(stats_rows).pivot(index='gene_set', columns='variant',
                                     values='p_one_holm').round(4).to_string())
print('\nWritten: pub_fig1_CDK4_and_programs, pub_fig2_cyclins, pub_figS_all_published_sets')
