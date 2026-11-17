"""LNCaP equivalent of the PC3 Fig10a elongation figure, on a shared y-axis.

Fig10a_elongation_shared_y_LNCaP
    panel 1  HQ461 vs DMSO, NT and I733V
    panel 2  no drug, I733V vs NT

The PC3 version has three panels (HQ461 / ZLY-025 / baseline); LNCaP is a 2x2 design with
HQ461 only, so the ZLY-025 panel has no counterpart here.

Method matches the PC3 analysis exactly: gene-length deciles over expressed protein-coding
genes, per-replicate log2 ratios from the VST matrix (each replicate against the mean of its
reference group), mean +/- SD across n = 3.

NOTE ON LIBRARY TYPE: LNCaP is 3' DGE (UMI-deduplicated 3'-end reads, length correction
disabled upstream), whereas PC3 is full-length RNA-seq. Effect *sizes* are therefore not
directly comparable between the two cell lines; shape and genotype ordering are.
See README_elongation.md.

Run:  cd reanalysis/scripts && python 02_elongation_shared_y.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from _style import (RES, FIG, GENOS, VARS, GCOL, TLAB, INK, MUTED,
                    save, load, drug_con, expressed_pc, meta as load_meta, replicate_profiles)

meta = load_meta()
vst = pd.read_csv(os.path.join(RES, 'vst.csv'), index_col='gene_id')
vst.pop('gene_name')
pc = expressed_pc()
span = pc.span_bp
v = vst.loc[pc.index, meta.index]

bins = pd.qcut(span.rank(method='first'), 10, labels=False)
centre = span.groupby(bins).median()
print(f'{len(pc):,} expressed protein-coding genes; decile median span (kb):',
      (centre / 1e3).round(1).tolist())

SHARED_YLIM = (-0.45, 0.85)
ZOOM_YLIM = (-0.06, 0.10)
rows = []


def samples(g, t):
    return list(meta.index[(meta.genotype == g) & (meta.treatment == t)])


def draw(ax, g, test, ref, label, group):
    R = replicate_profiles(v, bins, test, ref)
    mu, sd = R.mean(), R.std(ddof=1)
    ax.errorbar(centre, mu, yerr=sd, color=GCOL[g], lw=1.3, marker='o', ms=3.0,
                mec=INK, mew=0.5, capsize=1.6, capthick=0.6, elinewidth=0.7,
                label=g, zorder=3, clip_on=False)
    for i in centre.index:
        rows.append({'group': group, 'contrast': label, 'genotype': g, 'decile': int(i),
                     'span_median_bp': centre[i], 'mean_log2ratio': mu[i], 'sd': sd[i],
                     'n_replicates': len(R)})


def build(ylim, tag):
    rows.clear()
    fig, axes = plt.subplots(1, 2, figsize=(4.5, 2.05), sharey=True,
                             gridspec_kw={'wspace': 0.14})
    ax = axes[0]
    for g in GENOS:
        draw(ax, g, samples(g, 'HQ461'), samples(g, 'Control'), f'{g}_HQ461', 'drug')
    ax.set_title('HQ461 vs DMSO', fontsize=6.5)
    ax.set_ylabel('median log$_2$ ratio per decile\n(mean ± SD, n = 3)')
    ax.legend(loc='upper right', handlelength=1.2, fontsize=5.5)

    ax = axes[1]
    for g in VARS:
        draw(ax, g, samples(g, 'Control'), samples('NT', 'Control'), f'{g}_Control', 'baseline')
    ax.set_title('No drug: variant vs NT', fontsize=6.5)
    ax.legend(loc='lower left', handlelength=1.2, fontsize=5.5)

    for a in axes:
        a.axhline(0, color=MUTED, lw=0.6, ls='--', zorder=1)
        a.set_xscale('log')
        a.set_xlabel('Gene length (bp)')
        a.set_axisbelow(True)
    axes[0].set_ylim(*ylim)
    save(fig, tag)


# Matched to the PC3 Fig10a y-scale, so the two cell lines sit side by side.
build(SHARED_YLIM, 'Fig10a_elongation_shared_y_LNCaP')
# Zoomed ~10x: shows the flat profile is genuinely structureless, not a trend
# hidden by the PC3-matched scale.
build(ZOOM_YLIM, 'Fig10a_elongation_zoom_LNCaP')

out = pd.DataFrame(rows)
out.to_csv(os.path.join(FIG, 'Fig10a_LNCaP_decile_profiles.csv'), index=False)

# Top-decile response gap, matching the PC3 summary metric
rng = np.random.default_rng(0)


def gap_ci(top, allg, n_boot=2000):
    top, allg = top[np.isfinite(top)], allg[np.isfinite(allg)]
    g = np.median(top) - np.median(allg)
    b = [np.median(top[rng.integers(0, len(top), len(top))]) -
         np.median(allg[rng.integers(0, len(allg), len(allg))]) for _ in range(n_boot)]
    return g, *np.percentile(b, [2.5, 97.5])


gaps = []
for g in GENOS:
    l = drug_con(g).log2FC_shrunk.reindex(pc.index)
    gaps.append({'contrast': f'{g} HQ461 vs DMSO', **dict(zip(
        ['gap', 'ci_lo', 'ci_hi'], gap_ci(l[bins == 9].values, l.values)))})
l = load('I733V_Control_vs_NT_Control').log2FC_shrunk.reindex(pc.index)
gaps.append({'contrast': 'I733V vs NT (DMSO)', **dict(zip(
    ['gap', 'ci_lo', 'ci_hi'], gap_ci(l[bins == 9].values, l.values)))})
gdf = pd.DataFrame(gaps)
gdf.to_csv(os.path.join(RES, 'elongation_length_gaps.csv'), index=False)
print('\nTop-decile response gap (shrunken LFC, bootstrap 95% CI):')
print(gdf.round(3).to_string(index=False))
