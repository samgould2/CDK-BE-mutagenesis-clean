"""Two alternative views of the gene-length response.

Fig11_length_decile_heatmap  length decile x condition, colour = median log2FC.
                             Compact summary of all 11 contrasts on one set of rows.
                             Drug and baseline columns carry separate colour scales
                             because the baseline effect is ~4x smaller.

FigS_length_scatter          per-gene density (hexbin) of log2FC vs gene length, with the
                             decile-median line overlaid. Its job is to show the spread the
                             binned plots hide, and that the trend is not outlier-driven.
                             A raw point scatter of 8,547 genes is an overplotted blob --
                             density + overlaid trend is the readable form.

Run:  cd reanalysis/scripts && python 09_length_heatmap_scatter.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patheffects import withStroke
from _style import (RES, FIG, GENOS, VARS, DRUGS, GCOL, TLAB, DIVERGING, INK, MUTED,
                    save, load, drug_con, expressed_pc)

pc = expressed_pc()
span = pc.span_bp
bins = pd.qcut(span.rank(method='first'), 10, labels=False)
centre = span.groupby(bins).median()
row_lab = [f'{c / 1e3:.0f}' if c >= 1e4 else f'{c / 1e3:.1f}' for c in centre]

drug_cols = [(g, d, f'{g}_{d}_vs_{g}_Control') for d in DRUGS for g in GENOS]
base_cols = [(g, 'Control', f'{g}_Control_vs_NT_Control') for g in VARS]


def decile_medians(contrast):
    return load(contrast).log2FC_shrunk.reindex(pc.index).groupby(bins).median()


M_drug = pd.DataFrame({c: decile_medians(c) for _, _, c in drug_cols})
M_base = pd.DataFrame({c: decile_medians(c) for _, _, c in base_cols})
M_drug.index = M_base.index = row_lab
pd.concat([M_drug, M_base], axis=1).to_csv(os.path.join(FIG, 'Fig11_length_decile_matrix.csv'))


# ═════════════════════════════════════════════════════════════════════════════
# Fig 11 — decile x condition heatmap
# ═════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.5),
                         gridspec_kw={'width_ratios': [8, 3], 'wspace': 0.10})


def panel(ax, M, cols, vlim, title, show_rows):
    im = ax.imshow(M.values, aspect='auto', cmap=DIVERGING,
                   norm=TwoSlopeNorm(0, -vlim, vlim), interpolation='nearest')
    labs = [g for g, _, _ in cols]
    ax.set_xticks(range(len(cols)), labs, rotation=90)
    for t, g in zip(ax.get_xticklabels(), labs):
        t.set_color(GCOL[g])
    if show_rows:
        ax.set_yticks(range(len(M)), M.index)
        ax.set_ylabel('Gene length decile (median, kb)')
    else:
        ax.set_yticks([])
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(title, fontsize=6.5, pad=14)
    return im


im1 = panel(axes[0], M_drug, drug_cols, 1.0, 'Drug vs DMSO', True)
axes[0].axvline(3.5, color='white', lw=2)
for j, d in enumerate(DRUGS):
    axes[0].text(j * 4 + 1.5, -0.85, TLAB[d], ha='center', va='bottom', fontsize=6)
im2 = panel(axes[1], M_base, base_cols, 0.12, 'DMSO, vs NT', False)
axes[1].text(1.0, -0.85, 'no drug', ha='center', va='bottom', fontsize=6)

for im, ax, lab in [(im1, axes[0], 'median log$_2$FC'), (im2, axes[1], 'median log$_2$FC\n(note ×8 scale)')]:
    cb = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.27, shrink=0.92, aspect=17)
    cb.set_label(lab, fontsize=5.5)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=5, width=0.4)
save(fig, 'Fig11_length_decile_heatmap')


# ═════════════════════════════════════════════════════════════════════════════
# Fig S — per-gene density scatter with decile-median overlay
# ═════════════════════════════════════════════════════════════════════════════
YL = 4.0
# Genes called out by name. MCL1 is a canonical short-lived CDK9/CDK12-13 target
# (18.8 kb, length decile 4 of 10), so it sits mid-axis rather than in the long tail.
HIGHLIGHT = ['MCL1']
name_of = drug_con('NT', 'HQ461').gene_name
hi_ids = {n: i for i, n in name_of.items() if n in HIGHLIGHT and i in pc.index}
missing = set(HIGHLIGHT) - set(hi_ids)
if missing:
    print('WARNING: highlight genes absent from the filtered universe:', sorted(missing))

fig, axes = plt.subplots(2, 4, figsize=(6.4, 3.2), sharex=True, sharey=True,
                         gridspec_kw={'wspace': 0.1, 'hspace': 0.3})
x = np.log10(span.values)
for r_i, d in enumerate(DRUGS):
    for c_i, g in enumerate(GENOS):
        ax = axes[r_i, c_i]
        res = drug_con(g, d)
        y = res.log2FC_shrunk.reindex(pc.index).values
        ok = np.isfinite(y)
        ax.hexbin(x[ok], np.clip(y[ok], -YL, YL), gridsize=34, cmap='Greys', bins='log',
                  mincnt=1, linewidths=0, rasterized=True)
        med = pd.Series(y, index=pc.index).groupby(bins).median()
        ax.plot(np.log10(centre.values), med.values, color=GCOL[g], lw=1.4, marker='o',
                ms=2.8, mec=INK, mew=0.4, zorder=5)
        ax.axhline(0, color=MUTED, lw=0.6, ls='--', zorder=4)
        for nm, gi in hi_ids.items():
            hx, hy = np.log10(span[gi]), res.log2FC_shrunk[gi]
            if not np.isfinite(hy):
                continue
            hy = np.clip(hy, -YL, YL)
            ax.scatter(hx, hy, s=16, marker='D', facecolor='white',
                       edgecolor=INK, linewidth=0.7, zorder=6)
            # place the label on the far side of the decile-median line and give it a
            # white halo, so it never sits under the trend line or a dark hexbin
            dy = 4.0 if hy >= med[bins[gi]] else -8.0
            ax.annotate(nm, (hx, hy), textcoords='offset points', xytext=(5.0, dy),
                        fontsize=5, fontstyle='italic', color=INK, zorder=7,
                        path_effects=[withStroke(linewidth=1.8, foreground='white')])
        ax.set_ylim(-YL, YL)
        if r_i == 0:
            ax.set_title(g, color=GCOL[g], fontweight='bold', fontsize=6.5)
        if c_i == 0:
            ax.set_ylabel(f'{TLAB[d]} vs DMSO\nlog$_2$FC per gene')
        if r_i == 1:
            ax.set_xlabel('log$_{10}$ gene length (bp)')
fig.text(0.5, -0.02, f'Grey: per-gene density (n = {len(pc):,}, log count scale). '
                     'Coloured line: decile median, as plotted in Fig7. '
                     'White diamond: ' + ', '.join(HIGHLIGHT) + '.',
         ha='center', fontsize=5.5)
save(fig, 'FigS_length_scatter')

print('\nDrug panel decile medians:')
print(M_drug.round(2).to_string())
print('\nBaseline panel decile medians:')
print(M_base.round(3).to_string())
