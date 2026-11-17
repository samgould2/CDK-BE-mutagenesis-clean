"""The CDK12/13 elongation defect: gene length and intron number, plus targeted GSEA.

Expectation under CDK12/13 inhibition: loss of processive RNA Pol II elongation causes
premature cleavage and polyadenylation at intronic poly(A) sites. Transcripts that require
the most elongation — long genes, intron-rich genes — are lost preferentially, while short
intron-poor transcripts rise. This script tests that directly rather than inferring it from
pathway enrichment.

Outputs:
  Fig7_elongation_defect          length- and intron-decile response curves + slope metrics
  Fig6_targeted_GSEA              mechanistically-motivated gene sets instead of a 50-set sweep
  results/elongation_slopes.csv   slope metrics with 95% CI
  results/gsea/targeted_sets.gmt  the custom sets, for reuse

Run:  cd reanalysis/scripts && python 04_elongation.py
"""
import os
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from _style import (RES, FIG, GENOS, VARS, DRUGS, GCOL, TLAB, DIVERGING, MUTED, GRID,
                    BAR_EDGE, save, load, drug_con, features, expressed_pc)

rng = np.random.default_rng(0)
pc = expressed_pc()
print(f'{len(pc):,} expressed protein-coding genes with annotation')

AXES = [
    ('span_bp', 'Gene length (bp)', 'length', True),
    ('n_introns', 'Introns per gene', 'intron', False),
]


def decile_profile(values, lfc, n=10):
    """Median and IQR of lfc within n equal-count bins of values."""
    b = pd.qcut(values.rank(method='first'), n, labels=False)
    centre = values.groupby(b).median()
    g = lfc.groupby(b)
    return b, centre, g.median(), g.quantile(0.25), g.quantile(0.75), g.size()


def slope_ci(x, y, n_boot=1000):
    """OLS slope of y on x with a bootstrap 95% CI."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    b = np.polyfit(x, y, 1)[0]
    boots = [np.polyfit(x[k], y[k], 1)[0] for k in (rng.integers(0, len(x), len(x)) for _ in range(n_boot))]
    return b, *np.percentile(boots, [2.5, 97.5])


def gap_ci(top, allgenes, n_boot=2000):
    """median(log2FC | top decile) - median(log2FC | all genes), bootstrap 95% CI.

    Preferred over an OLS slope because the NT ZLY-025 profile is non-monotonic across
    intron deciles, where a linear slope is a poor summary. Referenced to the overall
    median rather than the bottom decile because single-intron genes behave atypically
    under ZLY-025, which would otherwise dominate the contrast.
    """
    top = top[np.isfinite(top)]
    allgenes = allgenes[np.isfinite(allgenes)]
    g = np.median(top) - np.median(allgenes)
    boots = [np.median(top[rng.integers(0, len(top), len(top))]) -
             np.median(allgenes[rng.integers(0, len(allgenes), len(allgenes))]) for _ in range(n_boot)]
    return g, *np.percentile(boots, [2.5, 97.5])


# ═════════════════════════════════════════════════════════════════════════════
# Fig 7 — elongation defect
# ═════════════════════════════════════════════════════════════════════════════
rows, prof_rows = [], []
fig, axes = plt.subplots(2, 3, figsize=(6.8, 3.9),
                         gridspec_kw={'width_ratios': [1, 1, 0.9], 'wspace': 0.42, 'hspace': 0.85})

for r_i, (col, xlabel, key, logx) in enumerate(AXES):
    xvals = pc[col]
    xfit = np.log10(xvals.values) if logx else np.log10(xvals.values + 1)
    for c_i, d in enumerate(DRUGS):
        ax = axes[r_i, c_i]
        for g in GENOS:
            lfc = drug_con(g, d).log2FC_shrunk.reindex(pc.index)
            b, centre, med, q1, q3, nbin = decile_profile(xvals, lfc)
            ax.plot(centre, med, color=GCOL[g], lw=1.3, marker='o', ms=3.0,
                    mec='black', mew=0.5, label=g, zorder=3, clip_on=False)
            for bb in centre.index:
                prof_rows.append({'axis': key, 'drug': d, 'genotype': g, 'decile': int(bb),
                                  'x_median': centre[bb], 'median_log2FC': med[bb],
                                  'q25': q1[bb], 'q75': q3[bb], 'n_genes': int(nbin[bb])})
            sl, slo, shi = slope_ci(xfit, lfc.values)
            gap, glo, ghi = gap_ci(lfc[b == 9].values, lfc.values)
            rows.append({'axis': key, 'drug': d, 'genotype': g, 'slope': sl,
                         'slope_ci_lo': slo, 'slope_ci_hi': shi,
                         'gap_top_decile_vs_median': gap, 'gap_ci_lo': glo, 'gap_ci_hi': ghi})
        ax.axhline(0, color=MUTED, lw=0.6, ls='--')
        ax.set_xscale('log')
        if not logx:
            ax.set_xticks([1, 3, 10, 30], ['1', '3', '10', '30'])
        ax.set_ylim(-0.85, 1.05)
        ax.set_xlabel(xlabel)
        ax.set_title(f'{TLAB[d]} vs DMSO', fontsize=6.5)
        if c_i == 0:
            ax.set_ylabel('median log$_2$FC per decile')

slopes = pd.DataFrame(rows)
pd.DataFrame(prof_rows).to_csv(os.path.join(FIG, 'FigS_elongation_decile_profiles.csv'), index=False)
slopes.to_csv(os.path.join(RES, 'elongation_slopes.csv'), index=False)

for r_i, (col, xlabel, key, logx) in enumerate(AXES):
    ax = axes[r_i, 2]
    s = slopes[slopes.axis == key]
    for j, d in enumerate(DRUGS):
        for i, g in enumerate(GENOS):
            x = j * 5 + i
            row = s[(s.drug == d) & (s.genotype == g)].iloc[0]
            v = row.gap_top_decile_vs_median
            ax.bar(x, v, width=0.78, color=GCOL[g], **BAR_EDGE)
            ax.errorbar(x, v, yerr=[[v - row.gap_ci_lo], [row.gap_ci_hi - v]],
                        color='black', lw=0.7, capsize=1.5, capthick=0.7)
    ax.set_xticks([j * 5 + i for j in range(2) for i in range(4)], GENOS * 2, rotation=90, fontsize=5.5)
    ax.axhline(0, color='black', lw=0.6)
    ax.set_ylim(-1.15, 0.30)
    for j, d in enumerate(DRUGS):
        ax.text(j * 5 + 1.5, 0.19, TLAB[d], ha='center', va='bottom', fontsize=6.5)
    noun = 'Long-gene' if logx else 'Intron-rich'
    ax.set_ylabel(f'{noun} response gap\n(top decile − median log$_2$FC)', fontsize=6.5)

axes[0, 0].legend(loc='upper right', handlelength=1.2, fontsize=5.5, ncol=2, columnspacing=0.8)
for ax, lab in [(axes[0, 0], 'Gene length'), (axes[1, 0], 'Intron number')]:
    ax.text(-0.30, 1.34, lab, transform=ax.transAxes, fontsize=8, fontweight='bold', va='bottom')
save(fig, 'FigS_elongation_length_and_introns')


# ═════════════════════════════════════════════════════════════════════════════
# Targeted gene sets
# ═════════════════════════════════════════════════════════════════════════════
name_of = load('NT_HQ461_vs_NT_Control').gene_name


def syms(idx):
    return sorted(set(name_of.reindex(idx).dropna()))


q_len = pc.span_bp.rank(pct=True)
q_int = pc.n_introns.rank(pct=True)
targeted = {
    'LONG_GENES (top 10% span)': syms(pc.index[q_len > 0.9]),
    'SHORT_GENES (bottom 10% span)': syms(pc.index[q_len < 0.1]),
    'INTRON_RICH (top 10%)': syms(pc.index[q_int > 0.9]),
    'INTRON_POOR (bottom 10%)': syms(pc.index[q_int < 0.1]),
}


def read_gmt(path):
    out = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            out[p[0]] = [g for g in p[2:] if g]
    return out


SIG = os.path.join(RES, '..', '..', 'pc3_rnaseq_results', 'gsea_signatures')
curated = read_gmt(os.path.join(SIG, 'CDK12_13_custom_sets.gmt'))
ddr = read_gmt(os.path.join(SIG, 'DDR_sets.gmt'))
targeted.update({
    'HR / DDR (CDK12-dependent)': curated['CDK12_LOSS_DN'],
    'cGAS-STING / interferon': curated['CDK12_LOSS_UP'],
    'CDK12+CDK13 shared targets': curated['CDK12_CDK13_SHARED_DN'],
    'DSB repair via HR': ddr['GOBP_DSB_REPAIR_VIA_HR'],
    'Fanconi anaemia pathway': ddr['KEGG_FANCONI_ANEMIA_PATHWAY'],
})
# Immediate-early genes: short, intron-poor, induced by stress — the expected "up" class
targeted['Immediate-early genes'] = ['FOS', 'FOSB', 'JUN', 'JUNB', 'EGR1', 'EGR2', 'EGR3',
                                     'ATF3', 'IER2', 'IER3', 'DUSP1', 'NR4A1', 'KLF2', 'KLF4', 'ZFP36']
hallmark = read_gmt(os.path.join(RES, 'gsea', 'MSigDB_Hallmark_2020.gmt'))
for h in ['E2F Targets', 'G2-M Checkpoint', 'Myc Targets V1', 'DNA Repair',
          'mTORC1 Signaling', 'TNF-alpha Signaling via NF-kB', 'Hypoxia']:
    targeted[f'HALLMARK {h}'] = hallmark[h]

gmt_path = os.path.join(RES, 'gsea', 'targeted_sets.gmt')
with open(gmt_path, 'w') as f:
    for k, v in targeted.items():
        f.write('\t'.join([k, 'targeted'] + list(v)) + '\n')
print('targeted sets:', {k: len(v) for k, v in targeted.items()})

contrasts = [f'{g}_{d}_vs_{g}_Control' for d in DRUGS for g in GENOS] + \
            [f'{g}_Control_vs_NT_Control' for g in VARS]
out = []
for nm in contrasts:
    df = load(nm).dropna(subset=['stat', 'gene_name'])
    df = df.loc[df['stat'].abs().groupby(df['gene_name']).idxmax()]
    rnk = df.set_index('gene_name')['stat'].sort_values(ascending=False)
    pre = gp.prerank(rnk=rnk, gene_sets=targeted, permutation_num=1000, seed=42,
                     min_size=10, max_size=3000, threads=6, outdir=None, verbose=False)
    r = pre.res2d.copy()
    r['contrast'] = nm
    out.append(r)
    print(nm, 'done')
tg = pd.concat(out, ignore_index=True).rename(columns={'Term': 'term', 'NOM p-val': 'nom_p', 'FDR q-val': 'fdr'})
tg['NES'] = pd.to_numeric(tg.NES, errors='coerce')
tg['fdr'] = pd.to_numeric(tg.fdr, errors='coerce')
tg.to_csv(os.path.join(RES, 'gsea', 'gsea_targeted.csv'), index=False)


# ═════════════════════════════════════════════════════════════════════════════
# Fig 6 — targeted GSEA dot plot
# ═════════════════════════════════════════════════════════════════════════════
term_order = list(targeted.keys())
drug_cons = [f'{g}_{d}_vs_{g}_Control' for d in DRUGS for g in GENOS]
base_cons = [f'{g}_Control_vs_NT_Control' for g in VARS]
cols = drug_cons + base_cons
nes = tg.pivot(index='term', columns='contrast', values='NES').reindex(index=term_order, columns=cols)
fdr = tg.pivot(index='term', columns='contrast', values='fdr').reindex(index=term_order, columns=cols).clip(lower=1e-3)

fig, ax = plt.subplots(figsize=(3.5, 0.19 * len(term_order) + 1.1))
yy, xx = np.meshgrid(range(len(term_order)), range(len(cols)), indexing='ij')
sizes = (-np.log10(fdr.values)).clip(0, 3) / 3 * 48 + 2
sig = fdr.values < 0.05
s_ = ax.scatter(xx.ravel(), yy.ravel(), s=sizes.ravel(), c=nes.values.ravel(), cmap=DIVERGING,
                norm=TwoSlopeNorm(0, -3, 3), edgecolor=np.where(sig.ravel(), 'black', 'none'), linewidth=0.5)
ax.set_xticks(range(len(cols)), [g for d in DRUGS for g in GENOS] + VARS, rotation=90)
for t, g in zip(ax.get_xticklabels(), [g for d in DRUGS for g in GENOS] + VARS):
    t.set_color(GCOL[g])
ax.set_yticks(range(len(term_order)), term_order)
ax.set_xlim(-0.6, len(cols) - 0.4)
ax.set_ylim(len(term_order) - 0.5, -0.5)
ax.grid(color=GRID, lw=0.4)
ax.set_axisbelow(True)
for x in (3.5, 7.5):
    ax.axvline(x, color=MUTED, lw=0.6)
for x, lab in [(1.5, TLAB['HQ461']), (5.5, TLAB['ZLY025']), (9, 'DMSO,\nvs NT')]:
    ax.text(x, -1.0, lab, ha='center', va='bottom', fontsize=6.5)
ax.axhline(3.5, color=MUTED, lw=0.6, ls=':')
ax.axhline(8.5, color=MUTED, lw=0.6, ls=':')
cb = fig.colorbar(s_, ax=ax, shrink=0.22, anchor=(0, 1), pad=0.03)
cb.set_label('NES', fontsize=6)
cb.outline.set_linewidth(0.4)
cb.ax.tick_params(labelsize=5.5, width=0.4)
hs = [Line2D([], [], ls='', marker='o', ms=np.sqrt(v / 3 * 48 + 2), mfc=MUTED, mec='none', label=l)
      for v, l in [(-np.log10(0.25), '0.25'), (-np.log10(0.05), '0.05'), (2, '0.01'), (3, '≤0.001')]]
ax.legend(handles=hs, title='FDR', loc='upper left', bbox_to_anchor=(1.02, 0.60),
          labelspacing=0.85, title_fontsize=6)
ax.set_title('Targeted GSEA: elongation-dependent\nand CDK12-pathway gene sets', fontsize=6.5, pad=18)
save(fig, 'Fig6_targeted_GSEA')

print(slopes.round(3).to_string(index=False))
