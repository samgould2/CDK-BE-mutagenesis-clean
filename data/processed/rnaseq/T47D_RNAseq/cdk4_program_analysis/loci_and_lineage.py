"""
Two questions, three figures:

  1. genomewide_cnv_scan   - is the 9p21 loss the only segmental event, or one of many?
  2. locus_zoom            - does each affected region look like a DELETION (several
                             functionally unrelated neighbouring genes, sharp boundaries)
                             or like SILENCING (one transcription unit plus its own
                             antisense, immediate neighbours untouched)?
  3. lineage_programs      - what is the end state of the ER/luminal collapse: basal
                             switch, EMT, or loss of luminal identity without conversion?

Gene coordinates: Ensembl BioMart (GRCh38), cached in results/gene_coords.tsv.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

mpl.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})

HERE = Path(__file__).resolve().parent
BASE, OUT, FIG = HERE.parent, HERE / 'results', HERE / 'figures'

V = ['CDK4_D196N', 'CDK4_E56K', 'CDK4_W198']
LABEL = {'CDK4_D196N': 'D196N', 'CDK4_E56K': 'E56K', 'CDK4_W198': 'W198*'}
COLOR = {'CDK4_D196N': '#4393c3', 'CDK4_E56K': '#66bd63', 'CDK4_W198': '#d6604d'}

tpm = pd.read_csv(BASE / 'salmon.merged.gene_tpm.tsv', sep='\t', index_col=0)
gn = tpm.pop('gene_name')
meta = pd.Series(tpm.columns.str.extract(r'^T47D_(.+)_veh_\d+$')[0].values, index=tpm.columns)
gm = tpm.T.groupby(meta).mean().T
order = tpm.mean(axis=1).sort_values(ascending=False).index
s2i = pd.Series(order, index=gn.reindex(order).values); s2i = s2i[~s2i.index.duplicated()]

co = pd.read_csv(OUT / 'gene_coords.tsv', sep='\t')
co.columns = ['gene_id', 'symbol', 'chrom', 'start', 'end', 'biotype']
idx = pd.Series(tpm.index, index=[i.split('.')[0] for i in tpm.index])
co['mid'] = co.gene_id.map(idx)
co = co.dropna(subset=['mid'])
CH = [str(i) for i in range(1, 23)] + ['X']
co = co[co.chrom.isin(CH)]
for v in V:
    co[v] = np.log2((gm.loc[co.mid, v].values + 0.25) / (gm.loc[co.mid, 'CBE_NT'].values + 0.25))
co['nt'] = gm.loc[co.mid, 'CBE_NT'].values

# ══════════ 1. genome-wide scan ══════════════════════════════════════════════
e = co[co.nt > 2].copy()
e['chrom'] = pd.Categorical(e.chrom, categories=CH, ordered=True)
e = e.sort_values(['chrom', 'start']).reset_index(drop=True)
e['x'] = np.arange(len(e))
bounds = e.groupby('chrom', observed=True)['x'].agg(['min', 'max'])

fig, axes = plt.subplots(3, 1, figsize=(12, 6.2), sharex=True)
for ax, v in zip(axes, V):
    roll = e.groupby('chrom', observed=True)[v].transform(
        lambda s: s.rolling(7, center=True, min_periods=5).median())
    for k, (ch, b) in enumerate(bounds.iterrows()):
        if k % 2 == 0:
            ax.axvspan(b['min'], b['max'], color='#f2f2f2', zorder=0)
    ax.scatter(e.x, e[v].clip(-6, 4), s=1.6, color='#c8c8c8', lw=0, rasterized=True)
    ax.plot(e.x, roll, color=COLOR[v], lw=.9)
    ax.axhline(0, color='black', lw=.6)
    ax.set_ylim(-6.3, 4.3); ax.set_ylabel(f'{LABEL[v]}\nlog$_2$ vs NT', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=7)
    # mark the 9p21 block
    blk = e[(e.chrom == '9') & (e.start > 20.4e6) & (e.start < 22.1e6)]
    if len(blk):
        ax.annotate('9p21.3\nMTAP/CDKN2A', xy=(blk.x.mean(), -2.6), fontsize=6.5, ha='center', color='#a33',
                    xytext=(blk.x.mean(), -5.4), textcoords='data',
                    arrowprops=dict(arrowstyle='->', color='#a33', lw=.7))
    q = e[(e.chrom == '3') & (e.start > 178e6)]
    if len(q):
        ax.annotate('3q26-qter', xy=(q.x.mean(), -1.0), fontsize=6.5, ha='center', color='#36c',
                    xytext=(q.x.mean(), -4.2), textcoords='data',
                    arrowprops=dict(arrowstyle='->', color='#36c', lw=.7))
axes[-1].set_xticks(bounds.mean(axis=1).values)
axes[-1].set_xticklabels(bounds.index, fontsize=6.5)
axes[-1].set_xlabel('chromosome (expressed genes, NT TPM > 2, ordered by position)', fontsize=8.5)
axes[-1].set_xlim(0, len(e))
plt.tight_layout()
for ext in ('png', 'svg'):
    plt.savefig(FIG / f'genomewide_cnv_scan.{ext}', dpi=300, bbox_inches='tight')
plt.close()

# ══════════ 2. locus zoom: deletion vs silencing ═════════════════════════════
LOCI = [('9',  20.20, 22.35, '9p21.3  (MTAP / CDKN2A)',  (20.45, 21.99),
         ['MLLT3', 'FOCAD', 'HACD4', 'KLHL9', 'MTAP', 'CDKN2A', 'CDKN2B-AS1', 'TUSC1']),
        ('6',  151.15, 152.25, '6q25.1  (ESR1)',          (151.47, 151.70),
         ['AKAP12', 'ZBTB2', 'RMND1', 'DCPH1', 'CCDC170', 'ESR1', 'SYNE1']),
        ('10', 7.55, 8.35, '10p14  (GATA3)',              (8.010, 8.060),
         ['KIN', 'ATP5F1C', 'TAF3', 'GATA3-AS1', 'GATA3']),
        ('14', 37.05, 37.80, '14q21.1  (FOXA1)',          (37.45, 37.62),
         ['MIPOL1', 'FOXA1', 'TTC6'])]

fig, axes = plt.subplots(1, 4, figsize=(14.5, 3.8))
for ax, (ch, a, b, title, span, labels) in zip(axes, LOCI):
    d = co[(co.chrom == ch) & (co.start / 1e6 > a) & (co.start / 1e6 < b) & (co.nt > 0.5)].sort_values('start')
    ax.axvspan(span[0], span[1], color='#d6604d', alpha=.10, zorder=0)
    for v in V:
        ax.plot(d.start / 1e6, d[v], marker='o', ms=4.5, lw=1.0, color=COLOR[v],
                label=LABEL[v], alpha=.9, zorder=3)
    ax.axhline(0, color='black', lw=.7)
    ax.axhline(-1, color='#999999', lw=.7, ls='--')
    for r in d.itertuples():
        if r.symbol in labels:
            ax.annotate(r.symbol, (r.start / 1e6, min(getattr(r, v) for v in V)),
                        fontsize=6.5, style='italic', rotation=90, ha='center', va='top',
                        xytext=(0, -4), textcoords='offset points')
    ax.set_ylim(-7.6, 2.6)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(f'chr{ch} (Mb)', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=7.5)
    n_hit = int((d[V].max(axis=1) < -1).sum())
    ax.text(.03, .04, f'{n_hit} consecutive\nexpressed genes < −1', transform=ax.transAxes,
            fontsize=7, va='bottom', color='#a33')
axes[0].set_ylabel('log$_2$ (variant / CBE NT), mean TPM', fontsize=8.5)
axes[0].legend(fontsize=7.5, frameon=False, loc='upper left')
plt.tight_layout()
for ext in ('png', 'svg'):
    plt.savefig(FIG / f'locus_zoom_deletion_vs_silencing.{ext}', dpi=300, bbox_inches='tight')
plt.close()

# ══════════ 3. lineage programs ══════════════════════════════════════════════
PANELS = {
    'Luminal TF\nnetwork': ['ESR1', 'GATA3', 'FOXA1', 'SPDEF', 'AR', 'XBP1'],
    'ER target\ngenes': ['PGR', 'GREB1', 'TFF1', 'TFF3', 'STC2', 'CA12'],
    'Epithelial identity\n/ adhesion': ['CDH1', 'GRHL2', 'RAB25', 'ESRP1', 'CLDN4', 'EPCAM', 'KRT19'],
    'Basal /\nmyoepithelial': ['KRT5', 'KRT14', 'KRT17', 'TP63', 'SOX10', 'FOXC1', 'EGFR'],
    'Mesenchymal\n/ EMT': ['VIM', 'CDH2', 'ZEB1', 'SNAI2', 'FN1', 'AXL', 'LOXL2'],
}
genes, seps, ticks = [], [], []
for p, gs in PANELS.items():
    gs = [g for g in gs if g in s2i.index]
    ticks.append((len(genes) + len(gs) / 2 - .5, p))
    genes += gs; seps.append(len(genes) - .5)
M = pd.DataFrame({LABEL[v]: [np.log2((gm.loc[s2i[g], v] + .1) / (gm.loc[s2i[g], 'CBE_NT'] + .1)) for g in genes]
                  for v in V}, index=genes)
M.to_csv(OUT / 'lineage_program_log2FC.csv')

fig, ax = plt.subplots(figsize=(13, 2.5))
im = ax.imshow(M.T.values, cmap='RdBu_r', vmin=-6, vmax=6, aspect='auto')
ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=90, fontsize=7, style='italic')
ax.set_yticks(range(3)); ax.set_yticklabels(M.columns, fontsize=8.5)
for s in seps[:-1]:
    ax.axvline(s, color='black', lw=1.2)
for x, p in ticks:
    ax.text(x, -0.95, p, ha='center', va='bottom', fontsize=7.5)
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, pad=.012, aspect=10, shrink=.9)
cb.set_label('log$_2$ vs CBE NT', fontsize=8); cb.ax.tick_params(labelsize=7)
plt.tight_layout()
for ext in ('png', 'svg'):
    plt.savefig(FIG / f'lineage_programs.{ext}', dpi=300, bbox_inches='tight')
plt.close()

# ══════════ summary numbers ══════════════════════════════════════════════════
print('Consecutive expressed genes with log2 < -1 in all three lines, per locus:')
for ch, a, b, title, span, _ in LOCI:
    d = co[(co.chrom == ch) & (co.start / 1e6 > a) & (co.start / 1e6 < b) & (co.nt > 0.5)].sort_values('start')
    hit = d[(d[V] < -1).all(axis=1)]
    print(f'  {title:28s} {len(hit)} genes over {(hit.end.max()-hit.start.min())/1e3:7.0f} kb : '
          f'{", ".join(str(x) for x in hit.symbol)}')
print('\nFold change vs NT (median of the three lines):')
for p, gs in PANELS.items():
    vals = [2 ** M.loc[g].median() for g in gs if g in M.index]
    print(f'  {p.replace(chr(10), " "):32s} median {np.median(vals):6.2f}x   '
          f'range {min(vals):.2f}-{max(vals):.2f}x')
