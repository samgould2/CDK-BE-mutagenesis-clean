#!/usr/bin/env python3
"""
Publication figures for the early (T0 vs plasmid) residue essentiality score.

Run `t0_plasmid_essentiality.py` first; this reads its tables.

Writes to figures/T0_plasmid_essentiality/ :
    essentiality_tracks.pdf          the 9-CDK linear track - pairs 1:1 with the
                                     ChimeraX panels (same palette, same limits)
    essentiality_klifs_heatmap.pdf   the same score on the 85 KLIFS pocket positions,
                                     structurally aligned across the family
    essentiality_colorbar.pdf        standalone colour key for the structure figure
    benchmark_early_vs_late.pdf      what the early contrast can and cannot support
    track_<GENE>.pdf                 one track per CDK, for figure assembly
"""
import os, json, warnings
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
# paths resolved from the repository root -- see _paths.py
from ._paths import (LIBRARY, KLIFS_DIR, DISORDER, ALPHAMISSENSE, DBNSFP, LFC_DIR, EDITING, FIG_DIR, RENDERS)  # noqa: E402
LEVEL = os.environ.get('ESS_LEVEL', 'L2').upper()
SUFFIX = '' if LEVEL == 'L2' else f'_{LEVEL}'
TAB = f'{LFC_DIR}/T0_plasmid_essentiality{SUFFIX}'
FIG = f'{FIG_DIR}/T0_plasmid_essentiality{SUFFIX}'
os.makedirs(FIG, exist_ok=True)

mpl.rcParams.update({
    'pdf.fonttype': 42, 'ps.fonttype': 42,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'axes.linewidth': 0.7, 'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'legend.fontsize': 7, 'axes.titlesize': 8,
    'savefig.bbox': 'tight', 'savefig.transparent': True,
})

# ----------------------------------------------------------------- palette
# diverging, two hues with a neutral (not white, not a hue) midpoint:
#   depleted / essential  <- cool  |  neutral  |  warm -> enriched
DEPLETED, MID, ENRICHED = '#2166AC', '#F2F2F2', '#B2182B'
NO_DATA = '#9E9E9E'          # residues with no guide above the editing threshold
CMAP = LinearSegmentedColormap.from_list('ess', [DEPLETED, '#67A9CF', '#D1E5F0',
                                                 MID, '#FDDBC7', '#EF8A62', ENRICHED])
CMAP.set_bad(NO_DATA)         # residues with no measurement
VLIM = 1.0                    # ~1.6x the permutation-null sd of a residue score,
                              # i.e. full saturation ≈ where FDR < 0.05 starts
NORM = Normalize(-VLIM, VLIM)

INK, INK2, MUTED = '#1A1A1A', '#4D4D4D', '#B3B3B3'
POCKET_FILL = '#E8E8E8'

GENES = ['CDK7', 'CDK8', 'CDK19', 'CDK9', 'CDK12', 'CDK13', 'CDK2', 'CDK4', 'CDK6']

PARAMS = json.load(open(f'{TAB}/parameters.json'))
NULL_SD = PARAMS['null_sd']['early']

res = pd.read_csv(f'{TAB}/residue_essentiality.csv')
hold = pd.read_csv(f'{TAB}/residue_essentiality_SY5609.csv')
guide = pd.read_csv(f'{TAB}/guide_essentiality.csv', low_memory=False)
KLIFS = pd.read_csv(f'{KLIFS_DIR}/KLIFS_compiled.csv')
def _read_qc(path):
    # at L1/L2 some QC tables are deliberately empty (the step they describe is skipped)
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


qc = {f[3:-4]: _read_qc(f'{TAB}/{f}') for f in os.listdir(TAB) if f.startswith('qc_')}

import sys
from .._chimerax_config import STRUCTURES, CDK_COLORS


def find_ranges(nums):
    s = sorted({int(n) for n in nums})
    if not s: return []
    out, start = [], s[0]
    for i in range(1, len(s)):
        if s[i] != s[i - 1] + 1:
            out.append((start, s[i - 1])); start = s[i]
    out.append((start, s[-1]))
    return out


def rolling(x, y, window=21):
    """Median in a sliding window over *residue position*, so gaps do not distort it."""
    xs, ys = [], []
    for c in x:
        m = (x >= c - window // 2) & (x <= c + window // 2)
        if m.sum() >= 5:
            xs.append(c); ys.append(np.median(y[m]))
    return np.array(xs), np.array(ys)


# ========================================================== 1. linear tracks
def draw_track(ax, gene, show_xlabel=True, label_top=None):
    d = res[res.Gene == gene].sort_values('Codon')
    x, y = d.Codon.values, d.ess_early.values
    lo, hi = 1, max(d.Codon.max(), STRUCTURES[gene]['kd'][1]) + 5

    for a, b in find_ranges(KLIFS[gene].dropna()):
        ax.axvspan(a - .5, b + .5, color=POCKET_FILL, lw=0, zorder=0)

    ax.axhline(0, color=INK2, lw=.6, zorder=1)
    for s in (-1, 1):
        ax.axhline(s * NULL_SD, color=MUTED, lw=.5, ls=(0, (2, 2)), zorder=1)

    ax.scatter(x, y, c=y, cmap=CMAP, norm=NORM, s=7, linewidths=.15,
               edgecolors='white', zorder=3, rasterized=True)
    rx, ry = rolling(x, y)
    ax.plot(rx, ry, color=INK, lw=1.3, zorder=4, solid_capstyle='round')

    # significance is one-sided: FDR < 0.05 *and* depleted. The two-sided list is not
    # usable at L2 (most of its calls are GC-driven enrichments), and depletion is the
    # direction this assay is powered for.
    sig = d[(d.FDR_early < 0.05) & (d.ess_early < 0)]
    if len(sig):
        ax.scatter(sig.Codon, np.full(len(sig), 3.35), marker='|', s=16,
                   color=INK, linewidths=.8, zorder=5, clip_on=False)

    kd = STRUCTURES[gene]['kd']
    ax.plot([kd[0], kd[1]], [-3.45, -3.45], color=INK2, lw=2.2,
            solid_capstyle='butt', zorder=5, clip_on=False)

    ax.set_xlim(lo, hi); ax.set_ylim(-3.6, 3.6)
    ax.set_yticks([-3, 0, 3])
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='y', pad=1)
    ax.set_ylabel(gene, rotation=0, ha='right', va='center', labelpad=6,
                  fontsize=9, color=INK)
    if label_top:
        ax.set_title(label_top, loc='left', color=INK2, pad=2)
    if show_xlabel:
        ax.set_xlabel('residue', labelpad=1)


def fig_tracks():
    fig, axes = plt.subplots(len(GENES), 1, figsize=(7.2, 9.4))
    for i, gene in enumerate(GENES):
        draw_track(axes[i], gene, show_xlabel=(i == len(GENES) - 1))
    fig.subplots_adjust(hspace=0.55, right=0.86, top=0.955, bottom=0.05, left=0.10)

    cax = fig.add_axes([0.885, 0.36, 0.016, 0.26])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax,
                      extend='both', ticks=[-VLIM, 0, VLIM])
    cb.set_label('essentiality score\n(early LFC, corrected)', fontsize=7.5, labelpad=4)
    cb.ax.tick_params(labelsize=7, width=.6, size=2.5)
    cb.outline.set_linewidth(.5)

    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=POCKET_FILL, ec='none', label='ATP pocket (KLIFS)'),
        Line2D([], [], color=INK, lw=1.3, label='21-residue rolling median'),
        Line2D([], [], color=MUTED, lw=.5, ls=(0, (2, 2)), label='± 1 null s.d.'),
        Line2D([], [], color=INK2, lw=2.2, label='catalytic domain'),
        Line2D([], [], color=INK, marker='|', ls='none', ms=5,
               label='depleted, FDR < 0.05'),
    ]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 1.0),
               ncol=5, frameon=False, handlelength=1.4, columnspacing=1.4,
               handletextpad=0.5, labelcolor=INK2)
    fig.suptitle('Residue essentiality across the CDK family — early dropout (T0 vs plasmid)',
                 y=1.025, fontsize=10, color=INK)
    fig.savefig(f'{FIG}/essentiality_tracks.pdf')
    plt.close(fig)

    for gene in GENES:
        f, a = plt.subplots(figsize=(6.4, 1.35))
        draw_track(a, gene)
        f.savefig(f'{FIG}/track_{gene}.pdf'); plt.close(f)
    print('tracks written')


# ================================================== 2. KLIFS-aligned heatmap
def fig_klifs():
    idx = KLIFS.KLIFS_idx.tolist()
    M = np.full((len(GENES), len(idx)), np.nan)
    for i, g in enumerate(GENES):
        s = KLIFS[['KLIFS_idx', g]].dropna()
        m = dict(zip(s[g].astype(int), s.KLIFS_idx))
        d = res[res.Gene == g]
        look = {m[c]: v for c, v in zip(d.Codon, d.ess_early) if c in m}
        for j, k in enumerate(idx):
            if k in look: M[i, j] = look[k]

    fig, ax = plt.subplots(figsize=(10.5, 2.5))
    ax.imshow(M, cmap=CMAP, norm=NORM, aspect='auto', interpolation='nearest')
    ax.set_yticks(range(len(GENES))); ax.set_yticklabels(GENES, fontsize=8)
    ax.set_xticks(range(len(idx)))
    ax.set_xticklabels(idx, rotation=90, fontsize=4.6)
    for j in range(len(idx) + 1):
        ax.axvline(j - .5, color='white', lw=.35)
    for i in range(len(GENES) + 1):
        ax.axhline(i - .5, color='white', lw=.35)
    # subdomain boundaries
    sub = [k.split(':')[0] for k in idx]
    for j in range(1, len(sub)):
        if sub[j] != sub[j - 1]:
            ax.axvline(j - .5, color=INK, lw=.8)
    ax.tick_params(length=0)
    ax.set_xlabel('KLIFS pocket position', labelpad=3)
    ax.set_title('Early essentiality on the structurally aligned ATP pocket  '
                 f'(grey = no guide editing that residue at ≥{PARAMS["MIN_EDIT"]:.0f}%)',
                 loc='left', pad=5)
    cax = fig.add_axes([0.92, 0.25, 0.008, 0.5])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax,
                      extend='both', ticks=[-VLIM, 0, VLIM])
    cb.set_label('essentiality score', fontsize=7.5)
    cb.ax.tick_params(labelsize=7, width=.6, size=2.5); cb.outline.set_linewidth(.5)
    fig.savefig(f'{FIG}/essentiality_klifs_heatmap.pdf'); plt.close(fig)
    print('klifs heatmap written')


# =============================================== 3. standalone colour key
def fig_colorbar():
    fig = plt.figure(figsize=(1.5, 2.4))
    cax = fig.add_axes([0.08, 0.08, 0.2, 0.84])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax,
                      extend='both', ticks=[-VLIM, -0.75, 0, 0.75, VLIM])
    cb.set_label('essentiality score (early LFC, corrected)', fontsize=7.5, labelpad=6)
    cb.ax.tick_params(labelsize=7, width=.6, size=2.5); cb.outline.set_linewidth(.5)
    fig.text(0.62, 0.95, 'depleted\n(essential)', fontsize=7, color=INK2, va='top')
    fig.text(0.62, 0.08, 'enriched', fontsize=7, color=INK2, va='bottom')
    fig.savefig(f'{FIG}/essentiality_colorbar.pdf'); plt.close(fig)
    print('colour key written')


# ====================================================== 4. benchmark panel
def fig_benchmark():
    fig = plt.figure(figsize=(9.6, 7.4))
    gs = fig.add_gridspec(3, 3, hspace=0.62, wspace=0.38,
                          left=0.08, right=0.97, top=0.9, bottom=0.07)
    EARLY_C, LATE_C = '#2166AC', '#C8742A'
    SHORT = {'Subpool1': 'S1', 'Subpool2': 'S2', 'Subpool3': 'S3', 'SY-5609': 'SY'}

    # (a) composition dependence
    ax = fig.add_subplot(gs[0, 0])
    g = qc['gc_dependence']
    g = g[g.guide_class == 'targeting']
    xs = np.arange(len(g))
    ax.scatter(xs, g.R_GC_raw, s=22, color=MUTED, label='raw', zorder=3)
    ax.scatter(xs, g.R_GC_final, s=22, color=EARLY_C, label='corrected', zorder=3)
    for i, r in enumerate(g.itertuples()):
        ax.plot([i, i], [r.R_GC_raw, r.R_GC_final], color=MUTED, lw=.7, zorder=2)
    ax.axhline(0, color=INK2, lw=.6)
    ax.set_xticks(xs)
    ax.set_xticklabels([f'{r.Editor}·{SHORT[r.Screen_ID]}' for r in g.itertuples()],
                       fontsize=6, rotation=90)
    ax.set_ylabel('r (protospacer GC, score)')
    ax.set_title('a   representation artifact removed', loc='left')
    ax.legend(frameon=False, loc='upper right')
    ax.spines[['top', 'right']].set_visible(False)

    # (b) editing-level trend — only an L3 step; at L1/L2 the trend is left in
    ax = fig.add_subplot(gs[0, 1])
    e = qc.get('editing_detrend', pd.DataFrame())
    if len(e):
        e = e[e.contrast == 'early']
        xs = np.arange(len(e))
        ax.scatter(xs, e.rho_before, s=22, color=MUTED, label='raw', zorder=3)
        ax.scatter(xs, e.rho_after, s=22, color=EARLY_C, label='corrected', zorder=3)
        for i, r in enumerate(e.itertuples()):
            ax.plot([i, i], [r.rho_before, r.rho_after], color=MUTED, lw=.7, zorder=2)
        ax.axhline(0, color=INK2, lw=.6)
        ax.set_xticks(xs)
        ax.set_xticklabels([f'{r.Editor}·{SHORT[r.Screen_ID]}' for r in e.itertuples()],
                           fontsize=6, rotation=90)
        ax.set_ylabel('ρ (editing level, score)')
        ax.set_title('b   editing-efficiency trend removed', loc='left')
        ax.legend(frameon=False, loc='upper right')
        ax.spines[['top', 'right']].set_visible(False)
    else:
        ax.axis('off')
        ax.set_title(f'b   no editing detrend at {LEVEL}', loc='left')
        ax.text(0.0, 0.55, 'The editing-level detrend is an L3 step and is not\n'
                           'applied here. Thresholding on editing level was tested\n'
                           'as an alternative and is worse — see threshold_sweep.pdf.',
                transform=ax.transAxes, fontsize=6.4, color=INK2, va='center')

    # (c) truncating dose response
    ax = fig.add_subplot(gs[0, 2])
    t = qc['truncating_dose']
    xs = np.arange(len(t))
    ax.plot(xs, t.early_raw, '-o', color=MUTED, lw=1.4, ms=5, label='early, raw')
    ax.plot(xs, t.early, '-o', color=EARLY_C, lw=1.6, ms=5, label='early, corrected')
    ax.plot(xs, t.late, '-o', color=LATE_C, lw=1.6, ms=5, label='late (DMSO/T0)')
    ax.axhline(0, color=INK2, lw=.6)
    ax.set_xticks(xs); ax.set_xticklabels(t.bin, fontsize=7)
    ax.set_xlabel('% canonical reads truncating')
    ax.set_ylabel('median score')
    ax.set_title('c   truncating variants must deplete', loc='left')
    ax.legend(frameon=False, fontsize=6)
    ax.spines[['top', 'right']].set_visible(False)

    # (d) positive control per CDK
    ax = fig.add_subplot(gs[1, :2])
    pc = qc['positive_control'].set_index('Gene').reindex(GENES).dropna(how='all').reset_index()
    xs = np.arange(len(pc)); w = 0.38
    ax.bar(xs - w / 2, pc.early, w, color=EARLY_C, label='early (T0/plasmid), corrected')
    ax.bar(xs + w / 2, pc.late, w, color=LATE_C, label='late (DMSO/T0)')
    for i, r in enumerate(pc.itertuples()):
        for off, v, p in [(-w / 2, r.early, r.p_early), (w / 2, r.late, r.p_late)]:
            if p < 0.01:
                ax.text(i + off, v + (0.05 if v >= 0 else -0.12), '∗',
                        ha='center', fontsize=8, color=INK)
    ax.axhline(0, color=INK2, lw=.6)
    ax.set_xticks(xs); ax.set_xticklabels(pc.Gene, fontsize=7)
    ax.set_ylabel('median score,\ntruncating guides (≥50%)')
    ax.set_title('d   the internal positive control — and where each contrast fails it', loc='left')
    ax.legend(frameon=False, ncol=2, loc='upper left')
    ax.spines[['top', 'right']].set_visible(False)
    i9 = list(pc.Gene).index('CDK9')
    ax.annotate('CDK9 — the one CDK with a measurable loss-of-function\n'
                'cost in this line — is the one the early contrast inverts',
                xy=(i9 + w / 2, float(pc.late[i9])), xycoords='data',
                xytext=(0.985, 0.04), textcoords='axes fraction', fontsize=6.2,
                color=INK2, ha='right', va='bottom',
                arrowprops=dict(arrowstyle='-', lw=.5, color=MUTED,
                                connectionstyle='arc3,rad=-0.2'))
    ax.set_ylim(min(pc.early.min(), pc.late.min()) - 0.35,
                max(pc.early.max(), pc.late.max()) + 0.28)

    # (e) pocket contrast
    ax = fig.add_subplot(gs[1, 2])
    pk = qc['pocket_contrast'].pivot_table(index='Gene', columns='contrast',
                                           values=['delta', 'p']).reindex(GENES).dropna(how='all')
    xs = np.arange(len(pk))
    ax.barh(xs - 0.19, pk[('delta', 'early')], 0.38, color=EARLY_C)
    ax.barh(xs + 0.19, pk[('delta', 'late')], 0.38, color=LATE_C)
    ax.axvline(0, color=INK2, lw=.6)
    ax.set_yticks(xs); ax.set_yticklabels(pk.index, fontsize=7); ax.invert_yaxis()
    ax.set_xlabel('pocket − non-pocket\n(median score)')
    ax.set_title('e   ATP pocket depletion', loc='left')
    ax.spines[['top', 'right']].set_visible(False)

    # (f) independent infection
    ax = fig.add_subplot(gs[2, 0])
    j = res.merge(hold, on=['Gene', 'Codon'], suffixes=('', '_h')).dropna(
        subset=['ess_early', 'ess_early_h'])
    ax.scatter(j.ess_early, j.ess_early_h, s=4, color=EARLY_C, alpha=.3,
               linewidths=0, rasterized=True)
    r = scipy.stats.pearsonr(j.ess_early, j.ess_early_h).statistic
    ax.plot([-4, 2.5], [-4, 2.5], color=MUTED, lw=.7, ls=(0, (2, 2)))
    ax.set_xlabel('score, screens 1–3'); ax.set_ylabel('score, SY-5609 infection')
    ax.set_title(f'f   independent infection   r = {r:.2f}', loc='left')
    ax.spines[['top', 'right']].set_visible(False)

    # (g) cross-editor at residue level
    ax = fig.add_subplot(gs[2, 1])
    rp = qc['reproducibility']
    rp = rp[rp.comparison.isin(['early: ABE vs CBE', 'late: ABE vs CBE',
                                'early: primary vs SY-5609'])]
    xs = np.arange(len(rp))
    ax.bar(xs, rp.R.astype(float), 0.55,
           color=[EARLY_C if c.startswith('early') else LATE_C for c in rp.comparison])
    ax.axhline(0, color=INK2, lw=.6)
    ax.set_xticks(xs)
    ax.set_xticklabels([c.replace(': ', '\n') for c in rp.comparison], fontsize=6)
    ax.set_ylabel('r, residue score')
    ax.set_title('g   does the score reproduce?', loc='left')
    ax.spines[['top', 'right']].set_visible(False)

    # (h) score distribution vs the matched permutation null
    ax = fig.add_subplot(gs[2, 2])
    nullv = qc['residue_null'].ess_early_permuted
    bins = np.linspace(-3.5, 2.5, 70)
    ax.hist(nullv, bins=bins, density=True, color=MUTED, alpha=.8, lw=0,
            label='permuted null')
    ax.hist(res.ess_early.dropna(), bins=bins, density=True, histtype='step',
            color=EARLY_C, lw=1.3, label='observed')
    ax.axvline(0, color=INK2, lw=.6)
    n_sig = int((res.FDR_early < 0.05).sum())
    ax.set_xlabel('essentiality score'); ax.set_ylabel('density')
    ax.set_title(f'h   excess over a matched null  ({n_sig} residues, FDR < 0.05)',
                 loc='left', fontsize=7.4)
    ax.legend(frameon=False, fontsize=6)
    ax.spines[['top', 'right']].set_visible(False)

    fig.suptitle('Benchmarking the early (T0 vs plasmid) essentiality readout',
                 y=0.965, fontsize=11, color=INK, x=0.08, ha='left')
    fig.savefig(f'{FIG}/benchmark_early_vs_late.pdf'); plt.close(fig)
    print('benchmark written')


# ============================== 5. distributions by mutational class
CLASSES = ['Non-targeting', 'Intron', 'Unedited', 'Silent', 'Missense', 'Truncating']
# controls run light -> dark grey (nothing should happen in any of them); the two
# classes that change the protein take the cool pole of the score palette
CLASS_COLOR = {'Non-targeting': '#C9C9C9', 'Intron': '#A8A8A8', 'Unedited': '#878787',
               'Silent': '#5F5F5F', 'Missense': '#67A9CF', 'Truncating': '#2166AC'}
CLASS_NOTE = {'Non-targeting': 'no target\nsite',
              'Intron': 'edited,\nintronic',
              'Unedited': '≥90% reads\nunedited',
              'Silent': 'edited, protein\nunchanged',
              'Missense': '≥25% reads\nmissense',
              'Truncating': '≥25% reads\nstop'}


def _violin(ax, data, order, colors, ylim=(-4, 3.2)):
    parts = ax.violinplot([data[k] for k in order], positions=range(len(order)),
                          widths=.82, showextrema=False, showmedians=False)
    for pc, k in zip(parts['bodies'], order):
        pc.set_facecolor(colors[k]); pc.set_alpha(.55); pc.set_linewidth(0)
    bp = ax.boxplot([data[k] for k in order], positions=range(len(order)),
                    widths=.17, showfliers=False, patch_artist=True, zorder=4,
                    medianprops=dict(color='white', lw=1.1),
                    whiskerprops=dict(color=INK2, lw=.7),
                    capprops=dict(color=INK2, lw=.7))
    for b, k in zip(bp['boxes'], order):
        b.set_facecolor(colors[k]); b.set_edgecolor('none')
    ax.set_ylim(*ylim)
    ax.spines[['top', 'right']].set_visible(False)


def fig_class_distributions():
    """Where every class of mutation — and every control — sits on the score."""
    gd = guide.dropna(subset=['GuideClass'])
    nt = gd[gd.GuideClass == 'Non-targeting']

    fig = plt.figure(figsize=(9.6, 7.0))
    gs = fig.add_gridspec(2, 2, hspace=0.52, wspace=0.26,
                          left=0.075, right=0.975, top=0.90, bottom=0.075)

    for col, (key, title) in enumerate([
            ('ess_early', 'a   early dropout (T0 vs plasmid), corrected'),
            ('ess_late', 'b   late (DMSO vs T0)')]):
        ax = fig.add_subplot(gs[0, col])
        data = {k: gd[gd.GuideClass == k][key].dropna().values for k in CLASSES}
        band = nt[key].dropna()
        ax.axhspan(band.quantile(.25), band.quantile(.75), color=MUTED, alpha=.22,
                   lw=0, zorder=0)
        ax.axhline(0, color=INK2, lw=.6, zorder=1)
        _violin(ax, data, CLASSES, CLASS_COLOR)
        for i, k in enumerate(CLASSES):
            v = data[k]
            ax.text(i, 3.05, f'{len(v):,}', ha='center', fontsize=5.8, color=INK2)
            if k not in ('Non-targeting',) and len(v) > 5:
                pv = scipy.stats.mannwhitneyu(v, band).pvalue
                mark = '***' if pv < 0.001 else ('*' if pv < 0.05 else 'n.s.')
                ax.text(i, -3.78, mark, ha='center',
                        fontsize=6.5 if mark == 'n.s.' else 8, color=INK)
        ax.set_xticks(range(len(CLASSES)))
        ax.set_xticklabels([f'{k}\n{CLASS_NOTE[k]}' for k in CLASSES],
                           fontsize=5.6, color=INK2, linespacing=1.15)
        ax.set_ylabel('guide score (LFC)')
        ax.set_title(title, loc='left')
        if col == 0:
            ax.text(0.015, 0.055, 'grey band = non-targeting IQR\n'
                    '*** p < 0.001 vs non-targeting (Mann–Whitney)',
                    transform=ax.transAxes, fontsize=5.6, color=INK2, va='bottom')

    # (c) cumulative distributions - the shift is easier to read than the violin
    ax = fig.add_subplot(gs[1, 0])
    for k in CLASSES:
        v = np.sort(gd[gd.GuideClass == k].ess_early.dropna().values)
        if len(v) < 10: continue
        ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=CLASS_COLOR[k], lw=1.5,
                label=f'{k} (n={len(v):,})')
    ax.axvline(0, color=INK2, lw=.6)
    ax.set_xlim(-3.5, 3); ax.set_ylim(0, 1)
    ax.set_xlabel('early score (LFC)'); ax.set_ylabel('cumulative fraction of guides')
    ax.set_title('c   cumulative distribution, early score', loc='left')
    ax.legend(frameon=False, fontsize=6, loc='lower right')
    ax.spines[['top', 'right']].set_visible(False)

    # (d) residue-level: ATP pocket vs the rest, per CDK
    ax = fig.add_subplot(gs[1, 1])
    PK = {'ATP pocket': '#2166AC', 'rest of protein': '#A8A8A8'}
    w = 0.38
    for i, gene in enumerate(GENES):
        d = res[res.Gene == gene]
        for j, (lab, c) in enumerate(PK.items()):
            v = d[d.in_pocket if lab == 'ATP pocket' else ~d.in_pocket].ess_early.dropna()
            if len(v) < 5: continue
            pos = i + (j - .5) * w
            ax.add_patch(plt.Rectangle((pos - w / 2 + .02, v.quantile(.25)), w - .04,
                                       v.quantile(.75) - v.quantile(.25),
                                       fc=c, ec='none', alpha=.8))
            ax.plot([pos - w / 2 + .02, pos + w / 2 - .02], [v.median()] * 2,
                    color='white', lw=1.1, zorder=3)
            ax.plot([pos, pos], [v.quantile(.05), v.quantile(.95)], color=c, lw=.7, zorder=1)
    ax.axhline(0, color=INK2, lw=.6)
    ax.set_xticks(range(len(GENES))); ax.set_xticklabels(GENES, fontsize=6.5, rotation=45)
    ax.set_ylim(-2.4, 1.5)
    ax.set_ylabel('residue score (LFC)')
    ax.set_title('d   missense residues, ATP pocket vs the rest', loc='left')
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=c, ec='none', label=l)
                       for l, c in PK.items()], frameon=False, fontsize=6.5,
              loc='lower right', ncol=2)
    ax.spines[['top', 'right']].set_visible(False)

    fig.suptitle('Score distributions by mutational class, with every control',
                 y=0.965, fontsize=11, color=INK, x=0.075, ha='left')
    fig.savefig(f'{FIG}/class_distributions.pdf'); plt.close(fig)
    print('class distributions written')


# ======================================= 6. composite: structures + tracks
# The renders ship in data/reference/, not in figures/: they cannot be
# regenerated without ChimeraX and the structure files, so they are input here
# rather than output. Regenerate with scripts/chimerax/generate_essentiality_cxc.py
# plus a ChimeraX run if the score changes.
CXDIR = os.path.join(RENDERS, 'essentiality')


def common_crop(paths, pad=20):
    """
    One crop box for every panel.

    All nine renders come from the same camera, so they are already on a common
    scale in pixels; cropping each to its own content would silently rescale them
    against one another. The union of the alpha bounding boxes keeps the scale.
    """
    import matplotlib.image as mpimg
    imgs = {p: mpimg.imread(p) for p in paths}
    r0 = c0 = 10 ** 9; r1 = c1 = 0
    for im in imgs.values():
        a = im[..., 3] > 0.02 if im.shape[2] > 3 else (im[..., :3] < 0.99).any(-1)
        if not a.any():
            continue
        rr, cc = np.where(a.any(1))[0], np.where(a.any(0))[0]
        r0, r1 = min(r0, rr[0]), max(r1, rr[-1])
        c0, c1 = min(c0, cc[0]), max(c1, cc[-1])
    h, w = next(iter(imgs.values())).shape[:2]
    box = (max(r0 - pad, 0), min(r1 + pad + 1, h), max(c0 - pad, 0), min(c1 + pad + 1, w))
    return {p: im[box[0]:box[1], box[2]:box[3]] for p, im in imgs.items()}


def fig_composite(variant=''):   # '' = the cyclin-bearing panels, as in the
                                 # rest of the repo's structural figures
    """The supp-16 replacement: each CDK's structure beside its linear track."""
    missing = [g for g in GENES if not os.path.exists(f'{CXDIR}/{g}{variant}.png')]
    if missing:
        print(f'composite skipped - no ChimeraX renders for {missing}; '
              f'run notebooks/11_ChimeraX/generate_essentiality_cxc.py and render_all.cxc')
        return
    panels = common_crop([f'{CXDIR}/{g}{variant}.png' for g in GENES])
    n = len(GENES)
    fig = plt.figure(figsize=(7.5, 1.14 * n + 0.9))
    gs = fig.add_gridspec(n, 2, width_ratios=[1.35, 4.2], wspace=0.01,
                          hspace=0.58, left=0.02, right=0.875,
                          top=1 - 0.75 / (1.14 * n + 0.9), bottom=0.05)
    for i, gene in enumerate(GENES):
        axs = fig.add_subplot(gs[i, 0])
        axs.imshow(panels[f'{CXDIR}/{gene}{variant}.png'])
        axs.axis('off')
        axs.set_title(gene, loc='center', fontsize=9, color=INK, pad=1)
        axt = fig.add_subplot(gs[i, 1])
        draw_track(axt, gene, show_xlabel=(i == n - 1))
        axt.set_ylabel('')
        axt.tick_params(axis='y', labelsize=6)

    cax = fig.add_axes([0.9, 0.36, 0.015, 0.25])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax,
                      extend='both', ticks=[-VLIM, 0, VLIM])
    cb.set_label('essentiality score\n(early LFC, corrected)', fontsize=7)
    cb.ax.tick_params(labelsize=6.5, width=.6, size=2.5); cb.outline.set_linewidth(.5)
    fig.text(0.9, 0.33, f'no data', fontsize=6.5, color=INK2, ha='left')
    fig.patches.append(plt.Rectangle((0.875, 0.325), 0.018, 0.013, fc=NO_DATA,
                                     ec='none', transform=fig.transFigure))

    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=POCKET_FILL, ec='none', label='ATP pocket (KLIFS)'),
        Line2D([], [], color=INK, lw=1.3, label='21-residue rolling median'),
        Line2D([], [], color=MUTED, lw=.5, ls=(0, (2, 2)), label='± 1 null s.d.'),
        Line2D([], [], color=INK2, lw=2.2, label='catalytic domain'),
        Line2D([], [], color=INK, marker='|', ls='none', ms=5,
               label='depleted, FDR < 0.05'),
    ]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.985),
               ncol=5, frameon=False, handlelength=1.4, columnspacing=1.3,
               handletextpad=0.5, labelcolor=INK2, fontsize=6.8)
    fig.suptitle('Early dropout (T0 vs plasmid) residue essentiality, mapped onto the '
                 'CDK structures', y=0.998, fontsize=10, color=INK)
    fig.savefig(f'{FIG}/supp16_essentiality_composite.pdf', dpi=400)
    plt.close(fig)
    print('composite written')


def fig_depleted_montage(variant='_depleted', cols=3, figsize=(8.6, 9.4),
                         fname='depleted_residues_montage.pdf'):
    """
    All nine CDKs in one figure, significantly depleted residues drawn as atoms.

    Every panel comes from the same camera, so one shared crop box keeps them on a common
    scale — cropping each to its own content would silently rescale them against each other.
    When the grid leaves a cell free (5 across puts 5 on top and 4 below), the colour key
    goes in it rather than stealing width from the panels.
    """
    missing = [g for g in GENES if not os.path.exists(f'{CXDIR}/{g}{variant}.png')]
    if missing:
        print(f'depleted montage skipped - no renders for {missing}')
        return
    panels = common_crop([f'{CXDIR}/{g}{variant}.png' for g in GENES], pad=14)
    counts = {g: int(((res.Gene == g) & (res.FDR_early < 0.05) & (res.ess_early < 0)).sum())
              for g in GENES}
    nrows = int(np.ceil(len(GENES) / cols))
    free = nrows * cols - len(GENES)          # empty cells in the last row

    fig = plt.figure(figsize=figsize)
    right = 0.86 if free == 0 else 0.995
    gs = fig.add_gridspec(nrows, cols, hspace=0.06, wspace=0.02,
                          left=0.012, right=right, top=0.90, bottom=0.015)
    for i, gene in enumerate(GENES):
        ax = fig.add_subplot(gs[i // cols, i % cols])
        ax.imshow(panels[f'{CXDIR}/{gene}{variant}.png'])
        ax.axis('off')
        ax.set_title(gene, loc='center', fontsize=11, color=CDK_COLORS[gene], pad=1)
        ax.text(0.5, -0.02, f'{counts[gene]} residues', transform=ax.transAxes,
                ha='center', va='top', fontsize=8, color=INK2)

    if free:                                   # colour key in the spare cell
        ax = fig.add_subplot(gs[nrows - 1, cols - free])
        ax.axis('off')
        cax = ax.inset_axes([0.12, 0.42, 0.28, 0.40])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax,
                          extend='both', ticks=[-VLIM, 0])
        cb.set_label('essentiality score', fontsize=9)
        cb.ax.tick_params(labelsize=8.5, width=.6, size=2.5); cb.outline.set_linewidth(.5)
    else:
        cax = fig.add_axes([0.885, 0.40, 0.016, 0.22])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax,
                          extend='both', ticks=[-VLIM, 0])
        cb.set_label('essentiality score', fontsize=7.5)
        cb.ax.tick_params(labelsize=7, width=.6, size=2.5); cb.outline.set_linewidth(.5)

    fig.suptitle('Significantly depleted residues across the CDK family',
                 x=0.012, ha='left', y=0.985, fontsize=13, color=INK)
    fig.text(0.012, 0.952,
             f'FDR < 0.05 and score < 0, drawn as atoms on a pale cartoon; cyclin in wheat. '
             f'{sum(counts.values())} residues in all. Every panel shares one camera and one '
             f'crop, so the structures are on a common scale.',
             fontsize=9, color=INK2, va='top')
    fig.savefig(f'{FIG}/{fname}', dpi=400)
    # returned rather than closed, so a notebook can render it inline;
    # the pdf above is still written exactly as before
    print(f'{fname} written ({cols} across, {sum(counts.values())} residues)')
    return fig

if __name__ == '__main__':
    fig_tracks(); fig_klifs(); fig_colorbar(); fig_benchmark()
    fig_class_distributions(); fig_composite()
    fig_depleted_montage()
    fig_depleted_montage(cols=5, figsize=(16.5, 7.6),
                         fname='depleted_residues_montage_wide.pdf')
    print('figures ->', FIG)
