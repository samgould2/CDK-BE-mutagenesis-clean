"""Combined elongation figure: drug response and baseline on one gene-length axis.

Fig10_elongation_combined — three panels sharing the gene-length axis:
    HQ461 vs DMSO | ZLY-025 vs DMSO | baseline (variant vs NT, no drug)

Error bars are mean +/- SD across the three biological replicates, not the gene-to-gene
IQR. Replicate-level log2 ratios are built from the VST matrix (each replicate against the
mean of its reference group) because a DESeq2 shrunken LFC is a single contrast-level
estimate and has no per-replicate decomposition. VST ratios track the shrunken LFCs in
shape and sign but are compressed in magnitude (VST shrinks low-count genes toward the
mean), so values here run slightly smaller than in Fig7/Fig9 — the panel is internally
consistent and should not be read against those absolute numbers.

Run:  cd reanalysis/scripts && python 06_elongation_combined.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from _style import (RES, FIG, GENOS, VARS, DRUGS, GCOL, TLAB, INK, MUTED,
                    save, expressed_pc, meta as load_meta, replicate_profiles)

meta = load_meta()
vst = pd.read_csv(os.path.join(RES, 'vst.csv'), index_col='gene_id')
vst.pop('gene_name')
pc = expressed_pc()
span = pc.span_bp
v = vst.loc[pc.index, meta.index]

N_BINS = 10
bins = pd.qcut(span.rank(method='first'), N_BINS, labels=False)
centre = span.groupby(bins).median()
print(f'{len(pc):,} expressed protein-coding genes; decile median span (kb):',
      (centre / 1e3).round(1).tolist())


def samples(g, t):
    return list(meta.index[(meta.genotype == g) & (meta.treatment == t)])


rows = []


def plot_condition(ax, g, test, ref, label, group):
    R = replicate_profiles(v, bins, test, ref)
    mu, sd = R.mean(), R.std(ddof=1)
    ax.errorbar(centre, mu, yerr=sd, color=GCOL[g], lw=1.3, marker='o', ms=3.0,
                mec=INK, mew=0.5, capsize=1.6, capthick=0.6, elinewidth=0.7,
                label=g, zorder=3, clip_on=False)
    for i in centre.index:
        rows.append({'group': group, 'contrast': label, 'genotype': g, 'decile': int(i),
                     'span_median_bp': centre[i], 'mean_log2ratio': mu[i], 'sd': sd[i],
                     'n_replicates': len(R),
                     **{f'rep{k + 1}': R.iloc[k, i] for k in range(len(R))}})


fig, axes = plt.subplots(1, 3, figsize=(6.7, 2.05),
                         gridspec_kw={'width_ratios': [1, 1, 1], 'wspace': 0.3})

DRUG_YLIM = (-0.45, 0.85)
BASE_YLIM = (-0.12, 0.16)

for c_i, d in enumerate(DRUGS):
    ax = axes[c_i]
    for g in GENOS:
        plot_condition(ax, g, samples(g, d), samples(g, 'Control'), f'{g}_{d}', 'drug')
    ax.set_ylim(*DRUG_YLIM)
    ax.set_title(f'{TLAB[d]} vs DMSO\n ', fontsize=6.5)
    if c_i == 0:
        ax.set_ylabel('median log$_2$ ratio per decile\n(mean ± SD, n = 3)')
    else:
        ax.tick_params(labelleft=False)

ax = axes[2]
for g in VARS:
    plot_condition(ax, g, samples(g, 'Control'), samples('NT', 'Control'), f'{g}_Control', 'baseline')
ax.set_ylim(*BASE_YLIM)
ax.set_title('No drug: variant vs NT\n(note expanded y-scale)', fontsize=6.5)
ax.set_ylabel('median log$_2$ ratio per decile', fontsize=6.5)

for ax in axes:
    ax.axhline(0, color=MUTED, lw=0.6, ls='--', zorder=1)
    ax.set_xscale('log')
    ax.set_xlabel('Gene length (bp)')
    ax.set_axisbelow(True)
axes[0].legend(loc='upper right', handlelength=1.2, fontsize=5.5, ncol=2, columnspacing=0.8)
axes[2].legend(loc='lower left', handlelength=1.2, fontsize=5.5)

save(fig, 'Fig10_elongation_combined')

out = pd.DataFrame(rows)
out.to_csv(os.path.join(FIG, 'Fig10_replicate_decile_profiles.csv'), index=False)
print('\nSD across replicates of the top-decile median:')
print(out[out.decile == 9][['group', 'contrast', 'mean_log2ratio', 'sd']].round(3).to_string(index=False))
