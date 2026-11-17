"""Treatment-grouped endpoint panel for the in vivo xenograft cohorts.

Moved verbatim from `combined_panel.py`
-- it was already a module with its reasoning written down, so only this header
and the import of `cdkbe.plotting` changed. It draws figures 2m, 2n, 4n and 4o.

Original docstring follows.


Treatment-grouped endpoint panel for the in vivo xenograft cohorts, merging the
three plots that were previously separate (drug effect within each genotype, and
WT-vs-mutant within each of Control and treated) into one figure.

Layout: x is the treatment condition, so the WT and mutant tumors of the same
mouse sit next to each other and grey lines can join matched flanks.
Light fill = WT, solid fill = mutant; grey = Control, drug colour = treated.

Two families of comparison are annotated, and they need different tests:
  - inner brackets: WT vs mutant *within* a treatment group. Same mouse's two
    flanks where both grew, so this is paired when the design is complete.
  - outer brackets: Control vs treated *within* a genotype. Different mice, so
    always unpaired Welch.

Test choice for the inner bracket is made from the data, not hardcoded: the
paired t-test is used only when every tumor in the group belongs to a complete
pair (and there are >= MIN_PAIRS of them). The CDK6/palbociclib cohort is fully
paired; the CDK9 and CDK7 cohorts are not (several mice grew only one flank),
so those fall back to unpaired Welch, matching what their original figures
report. The bracket is labelled with whichever test was used, and both p-values
are always returned in the stats table.

Directional hypotheses are the pre-specified ones used throughout this project:
treatment lowers the measure (outer brackets, one-sided), and mutant > WT within
the treated group (inner bracket, one-sided). WT vs mutant within Control has no
predicted direction and is two-sided.
"""
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox

MIN_PAIRS = 3

CONTROL_SHADES = ('#DCDCDC', '#9A9A9A')
DRUG_SHADES = {
    'Palbociclib': ('#F3C7A4', '#CC5500'),   # burnt orange
    'Atuveciclib': ('#D5C6E8', '#9467BD'),   # purple, as in invivo_xenograft.ipynb
    'SY-5609':     ('#AEE8EF', '#0FA3B1'),   # cyan/teal, as in invivo_xenograft.ipynb
}


def stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


def pval_text(p):
    txt = f'p = {p:.2g}' if p >= 0.001 else f'p = {p:.1e}'
    return f'{txt} {stars(p)}'


def combined_panel(d, cohort, drug, y='FoldChange', ylabel='Fold change', title='',
                   wt_label='WT', ref_line=None, bottom=None, figsize=(5.4, 6.0),
                   box_w=0.42, fontscale=1.3, wrap_labels=False, force_square=False,
                   outfile=None):
    """`d` needs columns MouseID, Genotype_label, Treatment ('Control' or `drug`), and `y`."""
    genos = [wt_label, cohort]
    groups = ['Control', drug]
    shades = {'Control': CONTROL_SHADES, drug: DRUG_SHADES[drug]}
    fillcol = {(grp, g): shades[grp][i] for grp in groups for i, g in enumerate(genos)}
    pos = {(groups[0], genos[0]): 0.0, (groups[0], genos[1]): 0.52,
           (groups[1], genos[0]): 1.48, (groups[1], genos[1]): 2.0}

    fig, ax = plt.subplots(figsize=figsize)
    vals = {}
    for grp in groups:
        for g in genos:
            v = d[(d['Treatment'] == grp) & (d['Genotype_label'] == g)].set_index('MouseID')[y].dropna()
            vals[(grp, g)] = v
            ax.boxplot([v.values], positions=[pos[(grp, g)]], widths=box_w, patch_artist=True,
                       showfliers=False, medianprops=dict(color='black', lw=1.4),
                       whiskerprops=dict(color='black', lw=1.2), capprops=dict(color='black', lw=1.2),
                       boxprops=dict(facecolor=fillcol[(grp, g)], edgecolor='black', lw=1.2, alpha=0.6))

    # fixed jitter so the pair lines land on the points
    rng = np.random.default_rng(0)
    xs = {(grp, g): pos[(grp, g)] + rng.uniform(-0.21 * box_w, 0.21 * box_w, len(vals[(grp, g)]))
          for grp in groups for g in genos}

    for grp in groups:
        a, b = vals[(grp, genos[0])], vals[(grp, genos[1])]
        xa, xb = xs[(grp, genos[0])], xs[(grp, genos[1])]
        for m in a.index.intersection(b.index):   # join matched flanks
            ia, ib = a.index.get_loc(m), b.index.get_loc(m)
            ax.plot([xa[ia], xb[ib]], [a.iloc[ia], b.iloc[ib]], color='grey', lw=0.8, alpha=0.6, zorder=2)
        for g in genos:
            ax.scatter(xs[(grp, g)], vals[(grp, g)].values, s=64, color=fillcol[(grp, g)],
                       edgecolor='black', linewidth=1, zorder=3)

    if ref_line is not None:
        ax.axhline(ref_line, color='black', ls='--', lw=0.8, zorder=0)

    # bracket geometry scales to the visible axis span and to the font size, so the
    # annotation stack stays legible when the axis runs to 0 or the text is enlarged
    lo_data = min(d[y].min(), ref_line) if ref_line is not None else d[y].min()
    hi = d[y].max()
    lo_axis = bottom if bottom is not None else lo_data - 0.06 * (hi - lo_data)
    unit = 0.058 * (hi - lo_axis) * fontscale
    rows = []

    # inner brackets: WT vs mutant within a treatment group
    inner_top = {}
    for grp in groups:
        a, b = vals[(grp, genos[0])], vals[(grp, genos[1])]
        alt = 'two-sided' if grp == 'Control' else 'greater'   # treated: mutant > WT
        shared = a.index.intersection(b.index)
        p_unpaired = scipy.stats.ttest_ind(b, a, equal_var=False, alternative=alt)[1]
        if len(shared) >= MIN_PAIRS:
            p_paired = scipy.stats.ttest_rel(b.loc[shared], a.loc[shared], alternative=alt)[1]
            p_paired_two = scipy.stats.ttest_rel(b.loc[shared], a.loc[shared])[1]
        else:
            p_paired = p_paired_two = np.nan
        # paired test only when the design is actually complete in this group
        fully_paired = len(shared) >= MIN_PAIRS and len(shared) == len(a) == len(b)
        p_shown, shown_label = (p_paired, 'paired') if fully_paired else (p_unpaired, 'unpaired')

        rows.append(dict(Cohort=cohort, Comparison=f'{wt_label} vs {cohort} ({grp})', family='genotype',
                         test=alt, shown=shown_label, n_WT=len(a), n_mutant=len(b), n_pairs=len(shared),
                         mean_WT=a.mean(), mean_mutant=b.mean(), p_shown=p_shown,
                         p_paired=p_paired, p_welch_unpaired=p_unpaired, p_paired_two_sided=p_paired_two))

        top = max(a.max(), b.max()) + 0.45 * unit
        inner_top[grp] = top
        x0, x1 = pos[(grp, genos[0])], pos[(grp, genos[1])]
        ax.plot([x0, x0, x1, x1], [top, top + 0.22 * unit, top + 0.22 * unit, top], color='black', lw=1)
        ax.text((x0 + x1) / 2, top + 0.33 * unit, pval_text(p_shown) + f'\n({shown_label})',
                ha='center', va='bottom', fontsize=9.5 * fontscale)

    # outer brackets: Control vs treated within a genotype -> unpaired (different mice)
    base = max(inner_top.values()) + 2.5 * unit
    for tier, g in enumerate(genos):
        a, b = vals[(groups[0], g)], vals[(groups[1], g)]
        p = scipy.stats.ttest_ind(a, b, equal_var=False, alternative='greater')[1]  # treatment lowers value
        rows.append(dict(Cohort=cohort, Comparison=f'Control vs {drug} ({g})', family='treatment',
                         test='greater', shown='unpaired', n_WT=len(a), n_mutant=len(b), n_pairs=np.nan,
                         mean_WT=a.mean(), mean_mutant=b.mean(), p_shown=p, p_paired=np.nan,
                         p_welch_unpaired=p,
                         p_paired_two_sided=scipy.stats.ttest_ind(a, b, equal_var=False)[1]))
        yb = base + tier * 1.5 * unit
        x0, x1 = pos[(groups[0], g)], pos[(groups[1], g)]
        ax.plot([x0, x0, x1, x1], [yb, yb + 0.22 * unit, yb + 0.22 * unit, yb], color='black', lw=1)
        ax.text((x0 + x1) / 2, yb + 0.32 * unit, f'{g}:  ' + pval_text(p),
                ha='center', va='bottom', fontsize=9.5 * fontscale)

    ax.set_ylim(lo_axis, base + 1.5 * unit + 1.25 * unit)
    ax.set_xlim(-0.45, 2.45)
    ax.set_xticks([pos[(grp, g)] for grp in groups for g in genos])
    lbl = (lambda g: g.replace(' ', '\n')) if wrap_labels else (lambda g: g)
    ax.set_xticklabels([f'{lbl(g)}\n(n = {len(vals[(grp, g)])})' for grp in groups for g in genos],
                       fontsize=11.5 * fontscale)

    # treatment group labels under the genotype labels
    for grp in groups:
        xc = (pos[(grp, genos[0])] + pos[(grp, genos[1])]) / 2
        label_drop = -0.125 * fontscale - (0.055 if wrap_labels else 0)
        ax.annotate(grp, xy=(xc, label_drop), xycoords=('data', 'axes fraction'), ha='center', va='top',
                    fontsize=13 * fontscale, fontweight='bold',
                    color='#555555' if grp == 'Control' else DRUG_SHADES[drug][1], annotation_clip=False)
    ax.axvline(1.0, color='#DDDDDD', lw=1, zorder=0)

    ax.spines[['right', 'top']].set_visible(False)
    ax.tick_params(axis='y', which='major', labelsize=13 * fontscale)
    ax.tick_params(axis='x', length=0)
    ax.set_xlabel('')
    ax.set_ylabel(ylabel, fontsize=14 * fontscale)
    ax.set_title(title, fontsize=14 * fontscale, pad=10)
    fig.tight_layout()
    if outfile:
        fig.savefig(outfile, transparent=True, bbox_inches=_bbox(fig, force_square))
    return fig, pd.DataFrame(rows)


def _bbox(fig, force_square, pad=0.08):
    """Tight bounding box, optionally padded out to exactly 1:1.

    `bbox_inches='tight'` crops to the drawn content, so a square figure canvas does
    not give a square file — the width depends on how wide the y-axis label happens
    to be. When force_square is set, the tight box is expanded symmetrically on the
    shorter axis so the saved page is exactly square with the content centred.
    """
    if not force_square:
        return 'tight'
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer()).padded(pad)
    side = max(bb.width, bb.height)
    cx, cy = (bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2
    return Bbox.from_extents(cx - side / 2, cy - side / 2, cx + side / 2, cy + side / 2)
