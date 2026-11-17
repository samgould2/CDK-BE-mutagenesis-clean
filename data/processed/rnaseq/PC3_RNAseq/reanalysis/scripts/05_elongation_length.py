"""Gene-length-only elongation figures.

Fig7_elongation_length   drug response vs gene length (the length half of Fig7, standalone)
Fig9_baseline_elongation baseline DMSO: variant vs NT, same gene-length axis
                         -> tests whether the variants carry a chronic elongation defect
                            without drug, i.e. whether they are hypomorphic at steady state

Run:  cd reanalysis/scripts && python 05_elongation_length.py
"""
import os
import numpy as np
import pandas as pd
from _style import (RES, FIG, GENOS, VARS, DRUGS, GCOL, TLAB, INK, MUTED, BAR_EDGE,
                    save, load, drug_con, expressed_pc, decile_profile, gap_ci)
import matplotlib.pyplot as plt

rng = np.random.default_rng(0)
pc = expressed_pc()
span = pc.span_bp
print(f'{len(pc):,} expressed protein-coding genes')
print('decile median span (kb):',
      (span.groupby(pd.qcut(span.rank(method="first"), 10, labels=False)).median() / 1e3).round(1).tolist())

rows, prof_rows = [], []


def add(lfc, label, group):
    b, centre, med, q1, q3, n = decile_profile(span, lfc)
    gap, lo, hi = gap_ci(lfc[b == 9].values, lfc.values, rng)
    rows.append({'group': group, 'contrast': label, 'gap_top_decile_vs_median': gap,
                 'gap_ci_lo': lo, 'gap_ci_hi': hi})
    for i in centre.index:
        prof_rows.append({'group': group, 'contrast': label, 'decile': int(i),
                          'span_median_bp': centre[i], 'median_log2FC': med[i],
                          'q25': q1[i], 'q75': q3[i], 'n_genes': int(n[i])})
    return centre, med


def decile_axis(ax):
    ax.axhline(0, color=MUTED, lw=0.6, ls='--')
    ax.set_xscale('log')
    ax.set_xlabel('Gene length (bp)')


def gap_bars(ax, sub, labels, groups_per_block, block_labels, ylim, label_y):
    for j, blk in enumerate(block_labels):
        for i, g in enumerate(groups_per_block):
            x = j * (len(groups_per_block) + 1) + i
            r = sub[(sub.contrast == labels[(blk, g)])].iloc[0]
            v = r.gap_top_decile_vs_median
            ax.bar(x, v, width=0.78, color=GCOL[g], **BAR_EDGE)
            ax.errorbar(x, v, yerr=[[v - r.gap_ci_lo], [r.gap_ci_hi - v]],
                        color=INK, lw=0.7, capsize=1.5, capthick=0.7)
    step = len(groups_per_block) + 1
    ax.set_xticks([j * step + i for j in range(len(block_labels)) for i in range(len(groups_per_block))],
                  list(groups_per_block) * len(block_labels), rotation=90, fontsize=5.5)
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_ylim(*ylim)
    for j, blk in enumerate(block_labels):
        ax.text(j * step + (len(groups_per_block) - 1) / 2, label_y, blk,
                ha='center', va='bottom', fontsize=6.5)
    ax.set_ylabel('Long-gene response gap\n(top decile − median log$_2$FC)', fontsize=6.5)


# ═════════════════════════════════════════════════════════════════════════════
# Fig 7 (length only) — drug response
# ═════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(6.4, 1.95),
                         gridspec_kw={'width_ratios': [1, 1, 0.85], 'wspace': 0.38})
for c_i, d in enumerate(DRUGS):
    ax = axes[c_i]
    for g in GENOS:
        centre, med = add(drug_con(g, d).log2FC_shrunk.reindex(pc.index), f'{g}_{d}', 'drug')
        ax.plot(centre, med, color=GCOL[g], lw=1.3, marker='o', ms=3.0,
                mec=INK, mew=0.5, label=g, zorder=3, clip_on=False)
    decile_axis(ax)
    ax.set_ylim(-0.85, 1.05)
    ax.set_title(f'{TLAB[d]} vs DMSO', fontsize=6.5)
    if c_i == 0:
        ax.set_ylabel('median log$_2$FC per decile')
axes[0].legend(loc='upper right', handlelength=1.2, fontsize=5.5, ncol=2, columnspacing=0.8)

allrows = pd.DataFrame(rows)
gap_bars(axes[2], allrows[allrows.group == 'drug'],
         {(d, g): f'{g}_{d}' for d in DRUGS for g in GENOS},
         GENOS, DRUGS, (-0.55, 0.16), 0.08)
axes[2].set_xticklabels(GENOS * 2, rotation=90, fontsize=5.5)
for j, d in enumerate(DRUGS):
    axes[2].texts[j].set_text(TLAB[d])
save(fig, 'Fig7_elongation_length')


# ═════════════════════════════════════════════════════════════════════════════
# Fig 9 — baseline (DMSO) elongation defect: variant vs NT, no drug
# ═════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(4.3, 1.95),
                         gridspec_kw={'width_ratios': [1.25, 0.75], 'wspace': 0.42})
ax = axes[0]
# Axis is scaled to the baseline effect itself; the much larger acute-drug effect is
# compared on-scale in the bar panel instead of being overlaid here, which would
# flatten these curves to invisibility.
for g in VARS:
    centre, med = add(load(f'{g}_Control_vs_NT_Control').log2FC_shrunk.reindex(pc.index),
                      f'{g}_Control', 'baseline')
    ax.plot(centre, med, color=GCOL[g], lw=1.3, marker='o', ms=3.0,
            mec=INK, mew=0.5, label=g, zorder=3, clip_on=False)
decile_axis(ax)
ax.set_ylim(-0.13, 0.17)
ax.set_ylabel('median log$_2$FC per decile\n(variant vs NT, DMSO)')
ax.set_title('Baseline, no drug', fontsize=6.5)
ax.legend(loc='lower left', handlelength=1.2, fontsize=5.5)

allrows = pd.DataFrame(rows)
base = allrows[allrows.group == 'baseline']
ax = axes[1]
for i, g in enumerate(VARS):
    r = base[base.contrast == f'{g}_Control'].iloc[0]
    v = r.gap_top_decile_vs_median
    ax.bar(i, v, width=0.78, color=GCOL[g], **BAR_EDGE)
    ax.errorbar(i, v, yerr=[[v - r.gap_ci_lo], [r.gap_ci_hi - v]],
                color=INK, lw=0.7, capsize=1.5, capthick=0.7)
ax.set_xticks(range(len(VARS)), VARS, rotation=90)
ax.axhline(0, color=INK, lw=0.6)
# NT's acute HQ461 gap on the same scale: the baseline defect is ~4x smaller
nt_gap = allrows[allrows.contrast == 'NT_HQ461'].gap_top_decile_vs_median.iloc[0]
ax.axhline(nt_gap, color=MUTED, lw=0.8, ls=':')
ax.text(len(VARS) - 0.45, nt_gap + 0.008, 'NT + HQ461', ha='right', va='bottom',
        fontsize=5.5, color=MUTED)
ax.set_ylim(-0.47, 0.05)
ax.set_ylabel('Long-gene response gap\n(top decile − median log$_2$FC)', fontsize=6.5)
ax.set_title('Baseline gap', fontsize=6.5)
save(fig, 'Fig9_baseline_elongation')

allrows.to_csv(os.path.join(RES, 'elongation_length_gaps.csv'), index=False)
pd.DataFrame(prof_rows).to_csv(os.path.join(FIG, 'Fig7_Fig9_length_decile_profiles.csv'), index=False)
print(allrows.round(3).to_string(index=False))
