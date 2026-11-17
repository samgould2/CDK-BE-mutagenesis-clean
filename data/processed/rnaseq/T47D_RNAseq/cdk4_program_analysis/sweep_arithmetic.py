"""
Derivation schedule: the NT sgRNA was transduced and sorted 2-3 weeks BEFORE the three
variant sgRNAs, which were transduced into the same CBE parental pool after it had been
carried in culture for those 2-3 weeks. So NT and the variants are two draws from the
same pool at two different times, each through its own sort bottleneck.

Question: can plain selection over 2-3 weeks take a 9p21-null / ER-null subclone from
minority to the ~96.5% implied by the residual parental signal? Or does the sort
bottleneck have to be doing the work?

Deterministic haploid selection:  odds(p1) = odds(p0) * (1+s)^g
"""
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
FIG = Path(__file__).resolve().parent / 'figures'

P1 = 0.965                       # final clone frequency (median residual parental ~3.5%)
TD_H = 36.0                      # T47D doubling time, hours (typical 32-40 h)
odds = lambda p: p / (1 - p)

p0 = np.linspace(0.01, 0.90, 400)
fig, ax = plt.subplots(figsize=(5.6, 3.7))
for days, c in [(14, '#4393c3'), (17.5, '#66bd63'), (21, '#d6604d')]:
    g = days * 24 / TD_H
    s = (odds(P1) / odds(p0)) ** (1 / g) - 1
    ax.plot(p0 * 100, s * 100, color=c, lw=1.6, label=f'{days:.0f} days  ({g:.0f} doublings)')
ax.axhspan(0, 10, color='#cfe8cf', alpha=.55, zorder=0)
ax.axhspan(10, 25, color='#fdf1c7', alpha=.6, zorder=0)
ax.axhspan(25, 100, color='#f8d7d3', alpha=.6, zorder=0)
for y, t in [(4, 'plausible for a fitness-neutral-ish passenger'),
             (16.5, 'large but possible for p16/ARF-null'),
             (40, 'implausible by selection alone')]:
    ax.text(89, y, t, fontsize=6.8, ha='right', va='center', color='#555555')
ax.set_xlabel('clone frequency in the parental pool when NT was sorted (%)', fontsize=8.5)
ax.set_ylabel('selective advantage required\nper doubling (%)', fontsize=8.5)
ax.set_ylim(0, 60); ax.set_xlim(1, 90)
ax.legend(fontsize=7.5, frameon=False, title=f'gap between NT and variant sorts', title_fontsize=7.5)
ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=8)
ax.set_title(f'Reaching {P1:.0%} by selection alone', fontsize=9)
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'sweep_arithmetic.{e}', dpi=300, bbox_inches='tight')
plt.close()

print(f'Final clone frequency assumed: {P1:.1%}  (from the ~3.5% residual parental signal)')
print(f'T47D doubling time assumed: {TD_H:.0f} h\n')
print(f"{'start freq':>11s} " + ' '.join(f'{d:>14s}' for d in ['14 d (9 dbl)', '17.5 d (12 dbl)', '21 d (14 dbl)']))
for q in [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90]:
    row = []
    for days in (14, 17.5, 21):
        g = days * 24 / TD_H
        row.append(f'{((odds(P1)/odds(q))**(1/g)-1)*100:13.0f}%')
    print(f'{q:10.0%}  ' + ' '.join(row))
print('\nNote: NT itself shows essentially no 9p21 loss (CDKN2A 21 TPM, 329 counts, 3/3 reps),')
print('so the clone must have been rare in the pool at the time of the NT sort.')
