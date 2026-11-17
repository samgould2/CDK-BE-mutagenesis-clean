#!/usr/bin/env python3
"""
Granular view of the pocket signal: every CDK, every KLIFS position, every variant.

`peak_significance.pdf` panel d collapses each pocket column to one number (the mean
depletion across CDKs). These figures open that up:

  klifs_position_detail.pdf   the 9 x 85 KLIFS heatmap, with per-residue significance
                              drawn on the cells, a recurrence track above it and the
                              column-level permutation test marked
  klifs_position_strips.pdf   the ranked columns with all nine CDK scores shown as
                              points, so the mean is never read without its spread
  dfg_motif_detail.pdf        the DFG motif and activation-loop start, residue by
                              residue: which variant each editor actually made, how
                              well it edited, and where each CDK landed

Run after `t0_plasmid_essentiality.py` and `essentiality_peaks.py`.
"""
import os, sys, warnings, collections
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)

from .essentiality_figures import (CMAP, NORM, VLIM, DEPLETED, ENRICHED, NO_DATA,
                                  INK, INK2, MUTED, POCKET_FILL, GENES, res, KLIFS, FIG, TAB)
from .essentiality_peaks import (alignment_maps, column_annotations, GENES as PGENES,
                                SUBDOMAIN_NAME, DFG_BLOCK, FDR_CUT)
from .._chimerax_config import CDK_COLORS
from .. import paths

var = pd.read_csv(f'{TAB}/variant_essentiality.csv')
# every guide x residue observation behind the scores, for the un-averaged versions
guide_rows = pd.read_csv(f'{TAB}/residue_guide_rows.csv').rename(
    columns={'%Reads': 'pct'})   # positional attribute access is fragile after merges
# the within-pocket null is the one the result is quoted from: it permutes each CDK
# among pocket positions, so the domain's own background cancels
cols_t = pd.read_csv(f'{TAB}/klifs_within_pocket_tests.csv')
sub_t = pd.read_csv(f'{TAB}/klifs_subdomain_tests.csv')

SIG = 0.05


def ink_on(rgba):
    """Readable text colour for a filled cell: relative luminance decides."""
    r, g, b = rgba[:3]
    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return INK if lum > 0.55 else 'white'
say = lambda *a: print(*a, flush=True)

# the window shown in the DFG figure: the motif plus one position either side
DFG_WINDOW = ['VII:78', 'VIII:79', 'xDFG:80', 'xDFG:81', 'xDFG:82', 'xDFG:83',
              'a.l:84', 'a.l:85']
MOTIF_ROLE = {'xDFG:80': 'DFG−1', 'xDFG:81': 'DFG Asp', 'xDFG:82': 'DFG Phe',
              'xDFG:83': 'DFG Gly', 'a.l:84': 'DFG+1', 'a.l:85': 'DFG+2'}


# ------------------------------------------------------------------ data
def build():
    """
    One row per (KLIFS position, CDK): the residue, what was made there, where it landed.

    The gene -> codon mapping comes straight from `KLIFS_compiled.csv`, not via the Clustal
    alignment. Going through the alignment loses any KLIFS index whose column is already
    claimed by a different index for another gene — `linker:50` dropped out that way and
    left a blank column in the figure.
    """
    rows = []
    for _, row in KLIFS.iterrows():
        k = row.KLIFS_idx
        for g in GENES:
            codon, aa = row[g], row.get(f'{g}_seq')
            if pd.isna(codon):
                rows.append(dict(klifs=k, Gene=g, aa='-', codon=np.nan, status='absent'))
                continue
            codon = int(codon)
            d = res[(res.Gene == g) & (res.Codon == codon)]
            v = var[(var.Gene == g) & (var.Codon == codon)]
            rows.append(dict(
                klifs=k, Gene=g, aa=(aa if isinstance(aa, str) else '?'), codon=codon,
                status='measured' if len(d) else 'not tested',
                ess=float(d.ess_early.iloc[0]) if len(d) else np.nan,
                FDR=float(d.FDR_early.iloc[0]) if len(d) else np.nan,
                n_guides=int(d.n_guides.iloc[0]) if len(d) else 0,
                max_edit=float(d.max_edit.iloc[0]) if len(d) else np.nan,
                variants=';'.join(sorted(set(v.SAV))) if len(v) else ''))
    t = pd.DataFrame(rows)
    t['subdomain'] = t.klifs.str.split(':').str[0]
    # one-sided: significantly *depleted*
    t['sig'] = (t.FDR < SIG) & (t.ess < 0)
    t.to_csv(f'{TAB}/klifs_position_by_cdk.csv', index=False)
    say(f'{len(t)} (KLIFS position x CDK) rows, {t.klifs.nunique()} positions '
        f'-> {paths.rel(TAB)}/klifs_position_by_cdk.csv')
    return t, None


def order_positions():
    return [k for k in KLIFS.KLIFS_idx]


# ------------------------------------------- 1. the pocket as an alignment
# subdomain blocks and colours, identical to CDK_resistance_patterns.ipynb so this panel
# and KLIFS_resistance_MSA.pdf can sit next to each other
KLIFS_ELEMENTS = [((1, 3), 'I', 'goldenrod'), ((4, 9), 'g.l', 'green'),
                  ((10, 13), 'II', 'goldenrod'), ((14, 19), 'III', 'goldenrod'),
                  ((20, 30), 'aC', 'red'), ((31, 37), 'b.l', 'green'),
                  ((38, 41), 'IV', 'goldenrod'), ((42, 44), 'V', 'goldenrod'),
                  ((45, 45), 'GK', 'orange'), ((46, 48), 'Hinge', 'purple'),
                  ((49, 52), 'Linker', 'blue'), ((53, 59), 'aD', 'red'),
                  ((60, 64), 'aE', 'red'), ((65, 67), 'VI', 'goldenrod'),
                  ((68, 75), 'c.l', 'orange'), ((76, 78), 'VII', 'goldenrod'),
                  ((79, 79), 'VIII', 'goldenrod'), ((80, 83), 'xDFG', 'blue'),
                  ((84, 85), 'a.l', 'blue')]
NOT_ASSAYED = '#c9c9c9'
ABSENT = '#ffffff'


def fig_heatmap(t):
    """
    The pocket drawn as an alignment, in the house style of
    `figures/resistance_patterns/*/KLIFS_resistance_MSA.pdf`: one row per CDK, the residue
    letter in every cell, the cell shaded by its essentiality score, a KLIFS subdomain strip
    above, recurrence markers, and a conservation track with the consensus sequence below.
    """
    order = order_positions()
    pos = {k: i + 1 for i, k in enumerate(order)}      # KLIFS index -> 1..85
    npos = len(order)

    grid = {}
    for r in t.itertuples():
        if r.klifs in pos:
            grid[(r.Gene, pos[r.klifs])] = r

    # conservation and consensus over the CDKs present at each column
    cons, consensus = np.zeros(npos + 1), [''] * (npos + 1)
    for j in range(1, npos + 1):
        aas = [grid[(g, j)].aa for g in GENES
               if (g, j) in grid and grid[(g, j)].aa not in ('-', '')]
        if aas:
            c = collections.Counter(aas).most_common(1)[0]
            cons[j] = c[1] / len(GENES)
            consensus[j] = c[0] if c[1] > len(GENES) / 2 else 'X'
    n_sig = np.array([sum(1 for g in GENES if (g, j) in grid and grid[(g, j)].sig)
                      for j in range(npos + 1)])

    fig, axes = plt.subplots(3, 1, figsize=(17.2, 6.2),
                             height_ratios=[0.5, 3.1, 1.15], sharex=True)
    fig.subplots_adjust(left=0.055, right=0.995, top=0.875, bottom=0.21, hspace=0.06)

    # ---------------------------------------------------- subdomain strip
    # which sites and columns are significant at the family level
    sd_sig = set(sub_t[sub_t.FDR < 0.05].subdomain)
    col_sig = set(cols_t[cols_t.FDR < 0.05].klifs)
    SD_KEY = {'aC': 'αC', 'aD': 'αD', 'aE': 'αE'}

    ax = axes[0]
    for (lo, hi), label, col in KLIFS_ELEMENTS:
        key = SD_KEY.get(label, label)
        hit = key in sd_sig
        ax.plot([lo - 0.42, hi + 0.42], [0.12, 0.12], color=col,
                lw=4.0 if hit else 2.2, solid_capstyle='butt', clip_on=False)
        ax.text((lo + hi) / 2, 0.44, label + (' *' if hit else ''),
                ha='center', va='bottom', fontsize=9.5 if hit else 8.5,
                color=col, fontweight='bold')
    for j in range(1, npos + 1):
        k = order[j - 1]
        if k in col_sig:
            ax.plot(j, -0.45, marker='*', ms=9, color=INK, clip_on=False, zorder=5)
        elif n_sig[j] >= 3:
            ax.plot(j, -0.45, marker='v', ms=6, color=INK, clip_on=False)
        elif n_sig[j] == 2:
            ax.plot(j, -0.45, marker='v', ms=5, color='#9a9a9a', clip_on=False)
    ax.set_ylim(-0.8, 1.05)
    ax.set_yticks([]); ax.spines[:].set_visible(False)

    # ---------------------------------------------------- the alignment
    ax = axes[1]
    for i, g in enumerate(GENES):
        y = len(GENES) - 1 - i
        for j in range(1, npos + 1):
            r = grid.get((g, j))
            if r is None or r.status == 'absent':
                fc, txt, tc = ABSENT, '-', MUTED
            elif r.status == 'not tested':
                fc, txt, tc = NOT_ASSAYED, r.aa, '#4d4d4d'
            else:
                fc = CMAP(NORM(r.ess)); txt, tc = r.aa, ink_on(fc)
            ax.add_patch(Rectangle((j - .5, y - .5), 1, 1, fc=fc, ec='white', lw=.45))
            if r is not None and r.status == 'measured' and r.sig:
                ax.add_patch(Rectangle((j - .46, y - .46), .92, .92, fc='none',
                                       ec=INK, lw=1.15, zorder=4))
            ax.text(j, y, txt, ha='center', va='center', fontsize=6.4, color=tc, zorder=5)
    ax.set_yticks(range(len(GENES)))
    ax.set_yticklabels(GENES[::-1], fontsize=9)
    for lbl, g in zip(ax.get_yticklabels(), GENES[::-1]):
        lbl.set_color(CDK_COLORS[g])
    ax.set_ylim(-0.5, len(GENES) - 0.5)
    ax.tick_params(axis='y', length=0)
    ax.spines[:].set_visible(False)

    # ---------------------------------------------------- consensus track
    ax = axes[2]
    xs = np.arange(1, npos + 1)
    ax.bar(xs, cons[1:], width=0.82, color=[DEPLETED if c >= 0.78 else '#9ecae1'
                                            for c in cons[1:]], lw=0)
    for j in xs:
        ax.text(j, -0.14, consensus[j], ha='center', va='top', fontsize=6.4, color=INK2)
    ax.set_ylim(0, 1.05)
    ax.set_yticks([0, 0.5, 1.0]); ax.set_yticklabels(['0', '', '1'], fontsize=7)
    ax.set_ylabel('conservation', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlim(0.2, npos + 0.8)
    ax.set_xticks(list(range(5, npos + 1, 5)))
    ax.set_xticklabels(list(range(5, npos + 1, 5)), fontsize=7.5)
    ax.tick_params(axis='x', pad=14)
    ax.set_xlabel('KLIFS pocket position', fontsize=10, labelpad=12)

    # ---------------------------------------------------- legend
    sm = mpl.cm.ScalarMappable(norm=NORM, cmap=CMAP)
    cax = fig.add_axes([0.055, 0.055, 0.13, 0.022])
    cb = fig.colorbar(sm, cax=cax, orientation='horizontal', extend='both',
                      ticks=[-VLIM, 0, VLIM])
    cb.ax.tick_params(labelsize=6.5, width=.6, size=2.5)
    cb.outline.set_linewidth(.5)
    fig.text(0.055, 0.105, 'essentiality score (early LFC, corrected)',
             fontsize=7.5, color=INK2)
    handles = [Rectangle((0, 0), 1, 1, fc='none', ec=INK, lw=1.15,
                         label='significantly depleted (FDR < 0.05)'),
               Rectangle((0, 0), 1, 1, fc=NOT_ASSAYED, ec='white',
                         label='no guide edits this residue'),
               Rectangle((0, 0), 1, 1, fc=ABSENT, ec='#dddddd',
                         label='position absent in this CDK'),
               Line2D([], [], marker='v', ls='none', ms=5, color='#9a9a9a',
                      label='depleted in 2 CDKs'),
               Line2D([], [], marker='*', ls='none', ms=9, color=INK,
                      label='column depleted family-wide (FDR < 0.05)'),
               Line2D([], [], ls='-', lw=4, color='#777777',
                      label='subdomain depleted family-wide (FDR < 0.05), marked *')]
    fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.215, 0.028),
               ncol=5, frameon=False, fontsize=7.5, handlelength=1.5,
               columnspacing=1.5, labelcolor=INK2)

    fig.suptitle('Early essentiality across the KLIFS pocket, one row per CDK',
                 x=0.055, ha='left', y=0.975, fontsize=11.5, color=INK)
    fig.text(0.055, 0.945,
             'Per-residue recurrence is thin — only one position is significantly depleted '
             'in more than one CDK — so the family-level calls (starred column, bold '
             'subdomain) are the ones that carry weight; see peak_significance.pdf.',
             fontsize=7.2, color=INK2, va='top')
    fig.savefig(f'{FIG}/klifs_position_detail.pdf')
    say('klifs_position_detail.pdf written (alignment style)')
    # returned rather than closed, so a notebook can render it inline;
    # the pdf above is still written exactly as before
    return fig


# -------------------------------------- 2. ranked columns, all CDKs shown
def fig_strips(t, n_show=26):
    rank = cols_t.sort_values('mean_depletion', ascending=False).head(n_show)
    fig, ax = plt.subplots(figsize=(7.6, 7.4))
    ys = np.arange(len(rank))[::-1]
    rng = np.random.RandomState(0)

    for y, r in zip(ys, rank.itertuples()):
        s = t[(t.klifs == r.klifs) & (t.status == 'measured')]
        is_dfg = r.klifs.startswith(DFG_BLOCK)
        ax.axhspan(y - .5, y + .5, color=DEPLETED if is_dfg else 'white',
                   alpha=.07 if is_dfg else 0, lw=0, zorder=0)
        ax.plot([-r.mean_depletion, -r.mean_depletion], [y - .34, y + .34],
                color=INK, lw=1.8, zorder=4)
        for v in s.itertuples():
            jitter = rng.uniform(-.17, .17)
            ax.scatter(v.ess, y + jitter, s=34 if v.sig else 24,
                       facecolor=CDK_COLORS[v.Gene] if v.sig else 'white',
                       edgecolor=CDK_COLORS[v.Gene], linewidths=1.0, zorder=3)
        n_abs = int((t[t.klifs == r.klifs].status != 'measured').sum())
        if n_abs:
            ax.text(2.55, y, f'{n_abs}', fontsize=5.5, color=MUTED, va='center')

    ax.axvline(0, color=INK2, lw=.7)
    ax.set_yticks(ys)
    ax.set_yticklabels([f'{r.klifs}  ({r.modal_aa})' +
                        (f'  {MOTIF_ROLE[r.klifs]}' if r.klifs in MOTIF_ROLE else '')
                        for r in rank.itertuples()], fontsize=6.6)
    for lbl, r in zip(ax.get_yticklabels(), rank.itertuples()):
        if r.klifs.startswith(DFG_BLOCK):
            lbl.set_color(DEPLETED); lbl.set_fontweight('bold')
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim()); ax2.set_yticks(ys)
    ax2.set_yticklabels([('p < 0.001' if r.p < 0.001 else f'p = {r.p:.3f}')
                         + ('  *' if r.FDR < 0.05 else '') for r in rank.itertuples()],
                        fontsize=5.8, color=MUTED)
    ax2.tick_params(length=0); ax2.spines[['top', 'right', 'left']].set_visible(False)
    ax.set_xlim(-4.6, 2.8)
    ax.set_xlabel('residue essentiality score (one point per CDK)', fontsize=8.5)
    ax.set_title(f'The {n_show} most depleted pocket positions, with every CDK shown\n'
                 'black bar = family mean; filled point = that CDK significantly depleted '
                 '(FDR < 0.05); grey number = CDKs with no measurement.  '
                 'p from the within-pocket null; * = column FDR < 0.05',
                 loc='left', fontsize=8.5, color=INK)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(handles=[Line2D([], [], marker='o', ls='none', mfc=c, mec=c, ms=5, label=g)
                       for g, c in [(g, CDK_COLORS[g]) for g in GENES]],
              frameon=False, fontsize=6.2, ncol=5, loc='lower left',
              bbox_to_anchor=(0, -0.145), columnspacing=.9, handletextpad=.1)
    fig.tight_layout()
    fig.savefig(f'{FIG}/klifs_position_strips.pdf')
    plt.close(fig)
    say('klifs_position_strips.pdf written')


# ------------------------------------------------ 3. the DFG motif itself
def _guide_points(t, w):
    """One row per guide x residue, carrying its KLIFS position."""
    key = t[(t.status == 'measured') & t.klifs.isin(w)][['klifs', 'Gene', 'codon']]
    key = key.rename(columns={'codon': 'Codon'})
    key['Codon'] = key.Codon.astype(int)
    gp = guide_rows.merge(key, on=['Gene', 'Codon'], how='inner')
    return gp


def fig_dfg(t, level='residue'):
    """
    level='residue' : one point per (CDK, position) — the score as it is used elsewhere
    level='guide'   : one point per guide x residue — nothing averaged away
    """
    w = [k for k in DFG_WINDOW if (t.klifs == k).any()]
    sub = t[t.klifs.isin(w)].copy()
    ct = cols_t.set_index('klifs')
    gp = _guide_points(t, w) if level == 'guide' else None
    suffix = '' if level == 'residue' else '_guides'
    unit = 'CDK' if level == 'residue' else 'guide'
    n_unit = len(sub[sub.status == "measured"]) if level == "residue" else len(gp)

    fig = plt.figure(figsize=(11.2, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], width_ratios=[1.0, 1.05],
                          hspace=0.42, wspace=0.22, left=0.07, right=0.975,
                          top=0.885, bottom=0.08)

    # ---- (a) per-position points, one per CDK
    ax = ax_a = fig.add_subplot(gs[0, 0])
    rng = np.random.RandomState(1)
    for i, k in enumerate(w):
        s = sub[(sub.klifs == k) & (sub.status == 'measured')]
        if k in MOTIF_ROLE:
            ax.axvspan(i - .5, i + .5, color=DEPLETED, alpha=.06, lw=0, zorder=0)
        if level == 'guide':
            d = gp[gp.klifs == k]
            for v in d.itertuples():
                ax.scatter(i + rng.uniform(-.3, .3), v.ess_early,
                           s=6 + 22 * v.pct / 100.0, facecolor='none',
                           edgecolor=CDK_COLORS[v.Gene], linewidths=.8, alpha=.85, zorder=3)
            if len(d):
                ax.plot([i - .33, i + .33], [d.ess_early.mean()] * 2,
                        color=INK, lw=2, zorder=4)
        else:
            for v in s.itertuples():
                ax.scatter(i + rng.uniform(-.2, .2), v.ess, s=42 if v.sig else 30,
                           facecolor=CDK_COLORS[v.Gene] if v.sig else 'white',
                           edgecolor=CDK_COLORS[v.Gene], linewidths=1.1, zorder=3)
            if len(s):
                ax.plot([i - .33, i + .33], [s.ess.mean()] * 2, color=INK, lw=2, zorder=4)
        p = ct.p.get(k, np.nan)
        q = ct.FDR.get(k, np.nan)
        if pd.notna(p):
            ax.text(i, 2.05, ('p<0.001' if p < 0.001 else f'p={p:.3f}'), ha='center',
                    fontsize=5.6, color=INK if pd.notna(q) and q < SIG else MUTED,
                    fontweight='bold' if pd.notna(q) and q < SIG else 'normal')
        n_abs = int((sub[sub.klifs == k].status != 'measured').sum())
        if n_abs:
            ax.text(i, -4.62, f'{n_abs} not\ntested', ha='center', fontsize=5.0,
                    color=MUTED, linespacing=.95)
    ax.axhline(0, color=INK2, lw=.7)
    ax.set_xticks(range(len(w)))
    ax.set_xticklabels([f'{k}\n{MOTIF_ROLE.get(k, "")}' for k in w], fontsize=6.4)
    ax.set_xlim(-.6, len(w) - .4); ax.set_ylim(-4.9, 2.4)
    ax.set_ylabel('essentiality score', fontsize=8)
    ax.set_title(f'a   every {unit} at every position of the motif  '
                 f'(n = {n_unit} points)', loc='left', fontsize=8.5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(handles=[Line2D([], [], marker='o', ls='none', mfc=CDK_COLORS[g],
                              mec=CDK_COLORS[g], ms=4.5, label=g) for g in GENES]
                     + [Line2D([], [], marker='o', ls='none', mfc='white', mec=INK2,
                               ms=4.5, label='FDR ≥ 0.05')],
              frameon=False, fontsize=5.9, ncol=10, loc='upper center',
              bbox_to_anchor=(0.5, -0.17), columnspacing=.7, handletextpad=.1)
    if level == 'guide':
        ax.text(0.0, -0.265, 'one circle per guide; circle size = % canonical reads '
                'making that variant; bar = mean of the guides',
                transform=ax.transAxes, fontsize=5.8, color=INK2)

    # ---- (b) the same points, regrouped by CDK, on the same score axis
    ax = fig.add_subplot(gs[0, 1], sharey=ax_a)
    meas = sub[sub.status == 'measured']
    per = (meas.groupby('Gene').agg(mean=('ess', 'mean'), n=('ess', 'size'),
                                    nsig=('sig', 'sum')).reindex(GENES))
    rng_b = np.random.RandomState(2)
    for i, g in enumerate(GENES):
        if level == 'guide':
            d = gp[gp.Gene == g]
            for v in d.itertuples():
                ax.scatter(i + rng_b.uniform(-.3, .3), v.ess_early,
                           s=6 + 22 * v.pct / 100.0, facecolor='none',
                           edgecolor=CDK_COLORS[g], linewidths=.8, alpha=.85, zorder=3)
            if len(d):
                ax.plot([i - .33, i + .33], [d.ess_early.mean()] * 2,
                        color=INK, lw=2, zorder=4)
            ax.text(i, 2.05, f'{len(d)}', ha='center', va='top', fontsize=5.4, color=INK2)
        else:
            d = meas[meas.Gene == g]
            for v in d.itertuples():
                ax.scatter(i + rng_b.uniform(-.21, .21), v.ess, s=42 if v.sig else 30,
                           facecolor=CDK_COLORS[g] if v.sig else 'white',
                           edgecolor=CDK_COLORS[g], linewidths=1.1, zorder=3)
            if len(d):
                ax.plot([i - .33, i + .33], [d.ess.mean()] * 2, color=INK, lw=2, zorder=4)
            r = per.loc[g]
            ax.text(i, 2.05, f"{int(r['n'])}" + (f"\n{int(r.nsig)} sig" if r.nsig else ''),
                    ha='center', va='top', fontsize=5.4, color=INK2, linespacing=.95)
    ax.axhline(0, color=INK2, lw=.7)
    ax.set_xticks(range(len(GENES)))
    ax.set_xticklabels(GENES, fontsize=6.6, rotation=90)
    for lbl, g in zip(ax.get_xticklabels(), GENES):
        lbl.set_color(CDK_COLORS[g])
    ax.set_xlim(-.6, len(GENES) - .4)
    ax.tick_params(labelleft=False)
    ax.set_title('b   the same points, regrouped by CDK  '
                 f'(number above = {"guides" if level == "guide" else "positions"} '
                 'measured)', loc='left', fontsize=8.5)
    ax.spines[['top', 'right']].set_visible(False)

    # ---- (c) the variant grid: what was actually made
    ax = fig.add_subplot(gs[1, :])
    for i, k in enumerate(w):
        for j, g in enumerate(GENES):
            r = sub[(sub.klifs == k) & (sub.Gene == g)]
            x, y = i, len(GENES) - 1 - j
            if not len(r):
                continue
            r = r.iloc[0]
            if r.status == 'absent':
                ax.add_patch(Rectangle((x - .5, y - .5), 1, 1, fc='white',
                                       ec='#E6E6E6', lw=.5))
                ax.text(x, y, '—', ha='center', va='center', fontsize=6, color=MUTED)
                continue
            if r.status == 'not tested':
                ax.add_patch(Rectangle((x - .5, y - .5), 1, 1, fc=NO_DATA, alpha=.5,
                                       ec='white', lw=.5, hatch='///'))
                ax.text(x, y, f'{r.aa}{int(r.codon)}', ha='center', va='center',
                        fontsize=5.4, color=INK2)
                continue
            ax.add_patch(Rectangle((x - .5, y - .5), 1, 1,
                                   fc=CMAP(NORM(r.ess)), ec='white', lw=.5))
            if r.sig:
                ax.add_patch(Rectangle((x - .46, y - .46), .92, .92, fc='none',
                                       ec=INK, lw=1.2))
            muts = '/'.join(sorted({v[-1] for v in r.variants.split(';') if v}))
            tc = ink_on(CMAP(NORM(r.ess)))
            ax.text(x, y + .17, f'{r.aa}{int(r.codon)}{muts}', ha='center', va='center',
                    fontsize=5.4, color=tc)
            ax.text(x, y - .22, f'{r.ess:+.1f}  {r.max_edit:.0f}%  n={int(r.n_guides)}',
                    ha='center', va='center', fontsize=4.6, color=tc)
    ax.set_xlim(-.5, len(w) - .5); ax.set_ylim(-.5, len(GENES) - .5)
    ax.set_xticks(range(len(w)))
    ax.set_xticklabels([f'{k}\n{MOTIF_ROLE.get(k, "")}' for k in w], fontsize=6.4)
    ax.set_yticks(range(len(GENES)))
    ax.set_yticklabels(GENES[::-1], fontsize=7)
    for lbl, g in zip(ax.get_yticklabels(), GENES[::-1]):
        lbl.set_color(CDK_COLORS[g])
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title('c   the variants the editors actually made  '
                 '(cell: the variants assayed; below: score, best editing, guides; '
                 'hatched = no guide reached 10% editing)', loc='left', fontsize=8.5)

    fig.suptitle('The DFG motif and the start of the activation loop, '
                 + ('residue by residue' if level == 'residue'
                    else 'every guide measurement shown'),
                 x=0.07, ha='left', y=0.975, fontsize=11, color=INK)
    fig.savefig(f'{FIG}/dfg_motif_detail{suffix}.pdf')
    plt.close(fig)
    say(f'dfg_motif_detail{suffix}.pdf written  ({n_unit} points)')


def report(t):
    sub = t[t.klifs.isin(DFG_WINDOW) & (t.status == 'measured')]
    per = (sub.groupby('Gene').agg(n=('ess', 'size'), mean=('ess', 'mean'),
                                   n_sig=('sig', 'sum')).reindex(GENES).round(3))
    say('\nDFG window, per CDK:')
    say(per.to_string())
    rest = res[~res.Codon.isin([])]
    say('\nmotif positions, family level:')
    say(cols_t[cols_t.klifs.isin(DFG_WINDOW)]
        [['klifs', 'modal_aa', 'n_cdks', 'mean_depletion', 'p', 'rank']]
        .round(3).to_string(index=False))
    a = sub[sub.Gene.isin(['CDK2', 'CDK4'])].ess.dropna()
    b = sub[~sub.Gene.isin(['CDK2', 'CDK4'])].ess.dropna()
    u = scipy.stats.mannwhitneyu(a, b)
    say(f'\nCDK2+CDK4 across the motif: median {a.median():+.2f} (n={len(a)})  vs '
        f'the other seven {b.median():+.2f} (n={len(b)})   p = {u.pvalue:.4f}')


if __name__ == '__main__':
    t, _ = build()
    fig_heatmap(t)
    fig_strips(t)
    fig_dfg(t, level='residue')
    fig_dfg(t, level='guide')
    report(t)
    say(f'\nfigures -> {FIG}')
