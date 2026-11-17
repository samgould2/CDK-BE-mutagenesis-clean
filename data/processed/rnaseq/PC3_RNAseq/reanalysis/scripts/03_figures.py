"""Publication figures for PC3 CDK12 variant x CDK12/13 inhibitor RNA-seq.

Inputs: reanalysis/results (from 01_deseq2.R, 02_gsea.py)
Outputs: reanalysis/figures/*.pdf (vector, Illustrator-editable text) + *.png (300 dpi)
Elongation and targeted-GSEA figures live in 04_elongation.py.
Run:  cd reanalysis/scripts && python 03_figures.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import TwoSlopeNorm, ListedColormap, to_rgba
from sklearn.decomposition import PCA
from _style import (RES, FIG, GENOS, VARS, DRUGS, TREATS, GCOL, TLAB, TMARK, TFILL,
                    DIVERGING, INK, MUTED, GRID, BAR_EDGE, PADJ, LFC,
                    save, load, drug_con, meta as load_meta)

rng = np.random.default_rng(0)
meta = load_meta()
vst = pd.read_csv(os.path.join(RES, 'vst.csv'), index_col='gene_id')
gene_name = vst.pop('gene_name')
vst = vst[meta.index]


# ── Core NT drug-response signature (shared by both drugs) ───────────────────
def nt_sig(d):
    r = drug_con('NT', d)
    s = r[(r.padj < PADJ)]
    return set(s[s.log2FC_shrunk > LFC].index), set(s[s.log2FC_shrunk < -LFC].index)


(hu, hd), (zu, zd) = nt_sig('HQ461'), nt_sig('ZLY025')
CORE_UP, CORE_DN = sorted(hu & zu), sorted(hd & zd)
print(f'Core NT signature: {len(CORE_UP)} up, {len(CORE_DN)} down')


# ═════════════════════════════════════════════════════════════════════════════
# Fig 1 — PCA (all samples; DMSO only)
# ═════════════════════════════════════════════════════════════════════════════
def pca_panel(ax, samples, ntop=500):
    x = vst[samples]
    top = x.var(axis=1).sort_values(ascending=False).index[:ntop]
    p = PCA(n_components=2).fit(x.loc[top].T.values)
    pcs = p.transform(x.loc[top].T.values)
    for s, (pc1, pc2) in zip(samples, pcs):
        g, t = meta.loc[s, 'genotype'], meta.loc[s, 'treatment']
        ax.scatter(pc1, pc2, marker=TMARK[t], s=24,
                   facecolor='white' if t == 'Control' else GCOL[g],
                   edgecolor=GCOL[g], linewidth=0.9, zorder=3)
    ev = p.explained_variance_ratio_ * 100
    ax.set_xlabel(f'PC1 ({ev[0]:.0f}% variance)')
    ax.set_ylabel(f'PC2 ({ev[1]:.0f}% variance)')


fig, axes = plt.subplots(1, 2, figsize=(5.2, 2.3), gridspec_kw={'width_ratios': [1.35, 1], 'wspace': 0.35})
pca_panel(axes[0], list(meta.index))
axes[0].set_title('All samples (top 500 variable genes)')
pca_panel(axes[1], list(meta.index[meta.treatment == 'Control']))
axes[1].set_title('DMSO only')
h = [Patch(facecolor=GCOL[g], edgecolor='black', lw=0.6, label=g) for g in GENOS]
h += [Line2D([], [], ls='', marker=TMARK[t], ms=4.5, mfc='white' if t == 'Control' else MUTED,
             mec=INK, mew=0.8, label=TLAB[t]) for t in TREATS]
axes[1].legend(handles=h, loc='upper left', bbox_to_anchor=(1.02, 1.0), handlelength=1.2)
save(fig, 'Fig1_PCA')


# ═════════════════════════════════════════════════════════════════════════════
# Fig 2 — Number of drug-responsive genes per genotype
# ═════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(3.7, 2.2), sharey=True, gridspec_kw={'wspace': 0.12})
deg_rows = []
for ax, d in zip(axes, DRUGS):
    for i, g in enumerate(GENOS):
        r = drug_con(g, d)
        up = int(((r.padj < PADJ) & (r.log2FC_shrunk > LFC)).sum())
        dn = int(((r.padj < PADJ) & (r.log2FC_shrunk < -LFC)).sum())
        deg_rows.append({'drug': d, 'genotype': g, 'up': up, 'down': dn})
        ax.bar(i, up, width=0.72, color=GCOL[g], **BAR_EDGE)
        # face alpha only — passing alpha= would fade the black outline too
        ax.bar(i, -dn, width=0.72, color=to_rgba(GCOL[g], 0.45), **BAR_EDGE)
        ax.text(i, up + 90, f'{up:,}', ha='center', va='bottom', fontsize=5.5)
        ax.text(i, -dn - 90, f'{dn:,}', ha='center', va='top', fontsize=5.5)
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_xticks(range(4), GENOS, rotation=90)
    ax.set_title(f'{TLAB[d]} vs DMSO')
    ax.set_ylim(-5300, 4400)
    ax.set_yticks([-4000, -2000, 0, 2000, 4000], ['4,000', '2,000', '0', '2,000', '4,000'])
axes[0].set_ylabel('DEGs (padj < 0.05, |log$_2$FC| > 1)\n← down    up →')
axes[1].spines['left'].set_visible(False)
axes[1].tick_params(axis='y', length=0)
pd.DataFrame(deg_rows).to_csv(os.path.join(FIG, 'Fig2_DEG_counts.csv'), index=False)
save(fig, 'Fig2_DEG_counts')


# ── Fig 2b — same, plus baseline (variant vs NT in DMSO) as a third panel ────
# Same overall figure size, so panel widths are set proportional to bar count
# (4 + 4 + 3) and the bars stay the same physical width across all three.
# y-axis is shared throughout: baseline counts really are ~10x smaller than the
# drug effects, and that comparison is the point of the panel.
def deg_panel(ax, genos, counts, title):
    for i, g in enumerate(genos):
        up, dn = counts[g]
        ax.bar(i, up, width=0.72, color=GCOL[g], **BAR_EDGE)
        ax.bar(i, -dn, width=0.72, color=to_rgba(GCOL[g], 0.45), **BAR_EDGE)
        ax.text(i, up + 90, f'{up:,}', ha='center', va='bottom', fontsize=5.5)
        ax.text(i, -dn - 90, f'{dn:,}', ha='center', va='top', fontsize=5.5)
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_xticks(range(len(genos)), genos, rotation=90)
    ax.set_xlim(-0.65, len(genos) - 0.35)
    ax.set_title(title)
    ax.set_ylim(-5300, 4400)
    ax.set_yticks([-4000, -2000, 0, 2000, 4000], ['4,000', '2,000', '0', '2,000', '4,000'])


def deg_counts(res):
    return (int(((res.padj < PADJ) & (res.log2FC_shrunk > LFC)).sum()),
            int(((res.padj < PADJ) & (res.log2FC_shrunk < -LFC)).sum()))


fig, axes = plt.subplots(1, 3, figsize=(3.7, 2.2), sharey=True,
                         gridspec_kw={'wspace': 0.14, 'width_ratios': [4, 4, 3]})
for ax, d in zip(axes, DRUGS):
    deg_panel(ax, GENOS, {g: deg_counts(drug_con(g, d)) for g in GENOS}, f'{TLAB[d]} vs DMSO')
base_counts = {g: deg_counts(load(f'{g}_Control_vs_NT_Control')) for g in VARS}
deg_panel(axes[2], VARS, base_counts, 'DMSO, vs NT')
axes[0].set_ylabel('DEGs (padj < 0.05, |log$_2$FC| > 1)\n← down    up →')
for ax in axes[1:]:
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
pd.DataFrame(deg_rows + [{'drug': 'Control_vs_NT', 'genotype': g, 'up': u, 'down': dn}
                         for g, (u, dn) in base_counts.items()]
             ).to_csv(os.path.join(FIG, 'Fig2b_DEG_counts_with_baseline.csv'), index=False)
save(fig, 'Fig2b_DEG_counts_with_baseline')


# ═════════════════════════════════════════════════════════════════════════════
# Fig 3 — Fraction of the NT drug response retained by each variant
# score = mean VST(core UP) − mean VST(core DN), centred on NT DMSO
# ═════════════════════════════════════════════════════════════════════════════
score = vst.loc[CORE_UP].mean() - vst.loc[CORE_DN].mean()
score -= score[meta.condition == 'NT_Control'].mean()
sc = meta.assign(score=score)
cm = sc.groupby('condition').score.mean()

rows = []
for g in GENOS:
    for d in DRUGS:
        nt_delta = cm[f'NT_{d}'] - cm['NT_Control']
        for s, v in sc[(sc.genotype == g) & (sc.treatment == d)].score.items():
            rows.append({'genotype': g, 'drug': d, 'sample': s,
                         'pct_response': 100 * (v - cm[f'{g}_Control']) / nt_delta})
pct = pd.DataFrame(rows)
pct.to_csv(os.path.join(FIG, 'Fig3_pct_response_per_replicate.csv'), index=False)
sc[['genotype', 'treatment', 'score']].to_csv(os.path.join(FIG, 'Fig3_signature_score_per_sample.csv'))

fig, ax = plt.subplots(figsize=(2.3, 2.2))
for j, d in enumerate(DRUGS):
    for i, g in enumerate(GENOS):
        x = j * 5 + i
        v = pct[(pct.genotype == g) & (pct.drug == d)].pct_response
        ax.bar(x, v.mean(), width=0.78, color=GCOL[g], **BAR_EDGE)
        ax.errorbar(x, v.mean(), yerr=v.std(ddof=1), color=INK, lw=0.7, capsize=1.5, capthick=0.7)
        ax.scatter(x + rng.uniform(-0.15, 0.15, len(v)), v, s=6, color='white',
                   edgecolor=INK, lw=0.4, zorder=3)
        ax.text(x, max(v.max(), v.mean()) + 4, f'{v.mean():.0f}', ha='center', va='bottom', fontsize=5.5)
ax.set_xticks([j * 5 + i for j in range(2) for i in range(4)], GENOS * 2, rotation=90)
for j, d in enumerate(DRUGS):
    ax.text(j * 5 + 1.5, 1.02, TLAB[d], transform=ax.get_xaxis_transform(),
            ha='center', va='bottom', fontsize=7)
# 100% is the NT anchor, so ticks stop there; the limit extends a little further only to
# fit NT's own replicate scatter (97.4-102.2% of its mean) rather than clip it.
ax.set_ylim(min(-6, pct.pct_response.min() - 6), 112)
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.axhline(0, color=INK, lw=0.6)
ax.axhline(100, color=MUTED, lw=0.6, ls='--')
ax.set_ylabel('Drug response retained\n(% of NT signature shift)')
ax.set_title(f'Core NT signature\n{len(CORE_UP):,} up / {len(CORE_DN):,} down genes', fontsize=6.5, pad=14)
save(fig, 'Fig3_response_retained')


# ═════════════════════════════════════════════════════════════════════════════
# Fig 4 — Genome-wide LFC: variant drug response vs NT drug response
# slope through origin over NT DEGs = fraction of NT response retained
# ═════════════════════════════════════════════════════════════════════════════
LIM = 8
fig, axes = plt.subplots(2, 3, figsize=(5.2, 3.5), sharex=True, sharey=True,
                         gridspec_kw={'wspace': 0.15, 'hspace': 0.3})
slope_rows = []
for r_i, d in enumerate(DRUGS):
    nt = drug_con('NT', d)
    for c_i, g in enumerate(VARS):
        ax = axes[r_i, c_i]
        v = drug_con(g, d)
        j = nt[['log2FC_shrunk', 'padj']].join(v[['log2FC_shrunk']], rsuffix='_v').dropna(
            subset=['log2FC_shrunk', 'log2FC_shrunk_v'])
        sig = (j.padj < PADJ) & (j.log2FC_shrunk.abs() > LFC)
        x, y = j.log2FC_shrunk.clip(-LIM, LIM), j.log2FC_shrunk_v.clip(-LIM, LIM)
        ax.scatter(x[~sig], y[~sig], s=1, color='#CFCFCF', lw=0, rasterized=True)
        ax.scatter(x[sig], y[sig], s=1.5, color=GCOL[g], alpha=0.35, lw=0, rasterized=True)
        xs, ys = j.log2FC_shrunk[sig].values, j.log2FC_shrunk_v[sig].values
        b = (xs * ys).sum() / (xs * xs).sum()
        boots = [(xs[k] * ys[k]).sum() / (xs[k] ** 2).sum()
                 for k in (rng.integers(0, len(xs), len(xs)) for _ in range(1000))]
        lo, hi = np.percentile(boots, [2.5, 97.5])
        rr = np.corrcoef(xs, ys)[0, 1]
        slope_rows.append({'drug': d, 'genotype': g, 'n_NT_DEGs': len(xs), 'slope': b,
                           'ci_lo': lo, 'ci_hi': hi, 'pearson_r_NT_DEGs': rr,
                           'pearson_r_all_genes': np.corrcoef(j.log2FC_shrunk, j.log2FC_shrunk_v)[0, 1]})
        ax.plot([-LIM, LIM], [-LIM, LIM], ls='--', color=MUTED, lw=0.6)
        ax.plot([-LIM, LIM], [-LIM * b, LIM * b], color=INK, lw=1.0)
        ax.axhline(0, color='#C3C2B7', lw=0.4, zorder=0)
        ax.axvline(0, color='#C3C2B7', lw=0.4, zorder=0)
        ax.text(-LIM + 0.4, LIM - 0.4, f'slope = {b:.2f}\nr = {rr:.2f}', va='top', fontsize=6)
        ax.set_xlim(-LIM, LIM)
        ax.set_ylim(-LIM, LIM)
        if r_i == 0:
            ax.set_title(g, color=GCOL[g], fontweight='bold')
        if c_i == 0:
            ax.set_ylabel(f'{TLAB[d]}\nvariant log$_2$FC vs DMSO')
        if r_i == 1:
            ax.set_xlabel('NT log$_2$FC vs DMSO')
pd.DataFrame(slope_rows).to_csv(os.path.join(FIG, 'Fig4_response_slopes.csv'), index=False)
fig.text(0.92, 0.5, f'Colored: NT DEGs (padj < {PADJ}, |log$_2$FC| > {LFC:g})\n'
                    'Dashed: y = x (full response)\nSolid: fit through origin',
         fontsize=5.5, va='center', ha='left')
save(fig, 'Fig4_LFC_scatter_vs_NT')


# ═════════════════════════════════════════════════════════════════════════════
# Fig 5 — Heatmap: core signature genes, each genotype centred on its own DMSO
# ═════════════════════════════════════════════════════════════════════════════
order = list(meta.index)
genes = CORE_UP + CORE_DN
mat = vst.loc[genes, order].copy()
for g in GENOS:
    cols = meta.index[meta.genotype == g]
    ctrl = meta.index[(meta.genotype == g) & (meta.treatment == 'Control')]
    mat[cols] = mat[cols].sub(vst.loc[genes, ctrl].mean(axis=1), axis=0)
nt_mean = (drug_con('NT', 'HQ461').log2FC_shrunk + drug_con('NT', 'ZLY025').log2FC_shrunk) / 2
mat = mat.loc[sorted(CORE_UP, key=lambda k: -nt_mean[k]) + sorted(CORE_DN, key=lambda k: -nt_mean[k])]

fig = plt.figure(figsize=(3.6, 4.2))
gs = fig.add_gridspec(3, 2, height_ratios=[0.05, 0.05, 1], width_ratios=[1, 0.035], hspace=0.04, wspace=0.05)
ax_g, ax_t, ax_h, cax = (fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]),
                         fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1]))
ax_g.imshow([[GENOS.index(meta.loc[s, 'genotype']) for s in order]], aspect='auto',
            cmap=ListedColormap([GCOL[g] for g in GENOS]), vmin=-0.5, vmax=3.5, interpolation='nearest')
ax_t.imshow([[TREATS.index(meta.loc[s, 'treatment']) for s in order]], aspect='auto',
            cmap=ListedColormap([TFILL[t] for t in TREATS]), vmin=-0.5, vmax=2.5, interpolation='nearest')
for a, lab in [(ax_g, 'Genotype'), (ax_t, 'Treatment')]:
    a.set_xticks([])
    a.set_yticks([0], [lab], fontsize=6)
    a.tick_params(length=0)
    for sp in a.spines.values():
        sp.set_visible(False)
for i, g in enumerate(GENOS):
    ax_g.text(i * 9 + 4, -0.9, g, ha='center', va='bottom', fontsize=6.5, color=GCOL[g], fontweight='bold')
im = ax_h.imshow(mat.values.clip(-3, 3), aspect='auto', cmap=DIVERGING, norm=TwoSlopeNorm(0, -3, 3),
                 interpolation='nearest', rasterized=True)
for i in range(1, 4):
    for a in (ax_g, ax_t, ax_h):
        a.axvline(i * 9 - 0.5, color='white', lw=1.2)
ax_h.axhline(len(CORE_UP) - 0.5, color=INK, lw=0.6)
ax_h.set_xticks([])
ax_h.set_yticks([len(CORE_UP) / 2, len(CORE_UP) + len(CORE_DN) / 2],
                [f'NT up\n(n = {len(CORE_UP):,})', f'NT down\n(n = {len(CORE_DN):,})'], fontsize=6)
ax_h.tick_params(axis='y', length=0)
for sp in ax_h.spines.values():
    sp.set_visible(False)
cb = fig.colorbar(im, cax=cax, ticks=[-3, -1.5, 0, 1.5, 3])
cb.set_label('log$_2$ expression vs\nsame-genotype DMSO', fontsize=6)
cb.outline.set_linewidth(0.4)
cax.tick_params(labelsize=5.5, width=0.4)
ax_h.legend(handles=[Patch(facecolor=TFILL[t], edgecolor='#C3C2B7', lw=0.4, label=TLAB[t]) for t in TREATS],
            loc='upper center', bbox_to_anchor=(0.5, -0.015), ncol=3, fontsize=6,
            handlelength=1, columnspacing=1.4)
save(fig, 'Fig5_signature_heatmap')


# ═════════════════════════════════════════════════════════════════════════════
# Fig 8 — CDK12-dependent DDR gene panel
# ═════════════════════════════════════════════════════════════════════════════
DDR = ['BRCA1', 'BRCA2', 'ATR', 'ATM', 'FANCI', 'FANCD2', 'FANCA', 'PALB2', 'RAD51C', 'RAD51D',
       'CHEK1', 'NBN', 'MRE11', 'TOPBP1', 'ORC3', 'ERCC4']
gid = {n: i for i, n in gene_name.items()}
DDR = [g for g in DDR if g in gid]
drug_cons = [f'{g}_{d}_vs_{g}_Control' for d in DRUGS for g in GENOS]
cols = drug_cons + [f'{g}_Control_vs_NT_Control' for g in VARS]
lfc = pd.DataFrame({c: load(c).log2FC_shrunk.reindex([gid[g] for g in DDR]).values for c in cols}, index=DDR)
pad = pd.DataFrame({c: load(c).padj.reindex([gid[g] for g in DDR]).values for c in cols}, index=DDR)
nt_cons = [f'NT_{d}_vs_NT_Control' for d in DRUGS]
ordr = lfc[nt_cons].mean(axis=1).sort_values().index
lfc, pad = lfc.loc[ordr], pad.loc[ordr]
lfc.to_csv(os.path.join(FIG, 'Fig8_DDR_gene_LFC.csv'))

fig, ax = plt.subplots(figsize=(2.9, 2.6))
im = ax.imshow(lfc.values, aspect='auto', cmap=DIVERGING, norm=TwoSlopeNorm(0, -2, 2), interpolation='nearest')
for i in range(lfc.shape[0]):
    for j in range(lfc.shape[1]):
        if pad.iloc[i, j] < PADJ:
            ax.text(j, i, '•', ha='center', va='center', fontsize=5, color=INK)
labs = [g for d in DRUGS for g in GENOS] + VARS
ax.set_xticks(range(len(cols)), labs, rotation=90)
for t, g in zip(ax.get_xticklabels(), labs):
    t.set_color(GCOL[g])
ax.set_yticks(range(len(ordr)), ordr, fontstyle='italic')
for x in (3.5, 7.5):
    ax.axvline(x, color='white', lw=2)
for x, lab in [(1.5, TLAB['HQ461']), (5.5, TLAB['ZLY025']), (9, 'DMSO,\nvs NT')]:
    ax.text(x, -0.85, lab, ha='center', va='bottom', fontsize=6)
for sp in ax.spines.values():
    sp.set_visible(False)
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, shrink=0.4, anchor=(0, 1), pad=0.03, ticks=[-2, -1, 0, 1, 2])
cb.set_label('log$_2$FC', fontsize=6)
cb.outline.set_linewidth(0.4)
cb.ax.tick_params(labelsize=5.5, width=0.4)
ax.set_title('CDK12-dependent DDR genes\n(• padj < 0.05)', pad=22, fontsize=6.5)
save(fig, 'Fig8_DDR_genes')


# ═════════════════════════════════════════════════════════════════════════════
# Fig S — Hallmark GSEA dot plots (all 50 sets; the targeted version is Fig 6)
# ═════════════════════════════════════════════════════════════════════════════
gsea = pd.read_csv(os.path.join(RES, 'gsea', 'gsea_all_contrasts.csv'))
hm = gsea[gsea.library == 'hallmark'].copy()
hm['fdr_floor'] = hm.fdr.clip(lower=1e-3)
hm['term'] = hm.term.replace({'Pperoxisome': 'Peroxisome', 'heme Metabolism': 'Heme Metabolism'})


def dotplot(ax, contrasts, labels, terms, colors=None):
    sub = hm[hm.contrast.isin(contrasts)]
    nes = sub.pivot(index='term', columns='contrast', values='NES').reindex(index=terms, columns=contrasts)
    fdr = sub.pivot(index='term', columns='contrast', values='fdr_floor').reindex(index=terms, columns=contrasts)
    yy, xx = np.meshgrid(range(len(terms)), range(len(contrasts)), indexing='ij')
    sizes = (-np.log10(fdr.values)).clip(0, 3) / 3 * 45 + 1
    sig = fdr.values < 0.05
    s_ = ax.scatter(xx.ravel(), yy.ravel(), s=sizes.ravel(), c=nes.values.ravel(), cmap=DIVERGING,
                    norm=TwoSlopeNorm(0, -3, 3), edgecolor=np.where(sig.ravel(), INK, 'none'), linewidth=0.4)
    ax.set_xticks(range(len(contrasts)), labels, rotation=90)
    if colors:
        for t, c in zip(ax.get_xticklabels(), colors):
            t.set_color(c)
    ax.set_yticks(range(len(terms)), terms)
    ax.set_xlim(-0.6, len(contrasts) - 0.4)
    ax.set_ylim(len(terms) - 0.5, -0.5)
    ax.grid(color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    return s_


def size_legend(ax, loc):
    hs = [Line2D([], [], ls='', marker='o', ms=np.sqrt(v / 3 * 45 + 1), mfc=MUTED, mec='none', label=l)
          for v, l in [(-np.log10(0.25), '0.25'), (-np.log10(0.05), '0.05'), (2, '0.01'), (3, '≤0.001')]]
    ax.legend(handles=hs, title='FDR', loc='upper left', bbox_to_anchor=loc, labelspacing=0.8, title_fontsize=6)


nt_hm = hm[hm.contrast.isin(nt_cons)]
term_order = nt_hm[nt_hm.fdr < 0.05].groupby('term').NES.mean().sort_values(ascending=False).index.tolist()
fig, ax = plt.subplots(figsize=(2.7, 0.14 * len(term_order) + 0.9))
s_ = dotplot(ax, drug_cons, [g for d in DRUGS for g in GENOS], term_order,
             colors=[GCOL[g] for d in DRUGS for g in GENOS])
for j, d in enumerate(DRUGS):
    ax.text(j * 4 + 1.5, -1.2, TLAB[d], ha='center', va='bottom', fontsize=7)
ax.axvline(3.5, color=MUTED, lw=0.6)
cb = fig.colorbar(s_, ax=ax, shrink=0.25, anchor=(0, 1), pad=0.04)
cb.set_label('NES', fontsize=6)
cb.outline.set_linewidth(0.4)
cb.ax.tick_params(labelsize=5.5, width=0.4)
size_legend(ax, (1.02, 0.62))
ax.set_title('Hallmark GSEA: drug vs DMSO within genotype\n(sets with FDR < 0.05 in NT; outlined = FDR < 0.05)',
             pad=16, fontsize=6.5)
save(fig, 'FigS_hallmark_GSEA_drug_response')

base_cons = [f'{g}_Control_vs_NT_Control' for g in VARS]
b_hm = hm[hm.contrast.isin(base_cons)]
b_order = b_hm[b_hm.fdr < 0.05].groupby('term').NES.mean().sort_values(ascending=False).index.tolist()
fig, ax = plt.subplots(figsize=(1.9, 0.14 * len(b_order) + 0.9))
s_ = dotplot(ax, base_cons, VARS, b_order, colors=[GCOL[g] for g in VARS])
cb = fig.colorbar(s_, ax=ax, shrink=0.25, anchor=(0, 1), pad=0.08)
cb.set_label('NES', fontsize=6)
cb.outline.set_linewidth(0.4)
cb.ax.tick_params(labelsize=5.5, width=0.4)
size_legend(ax, (1.05, 0.62))
ax.set_title('Hallmark GSEA: variant vs NT (DMSO)', fontsize=6.5)
save(fig, 'FigS_hallmark_GSEA_baseline')


# ═════════════════════════════════════════════════════════════════════════════
# Fig S — volcano grids
# ═════════════════════════════════════════════════════════════════════════════
def volcano(ax, res, title, xlim=10, ylim=None):
    r = res.dropna(subset=['padj', 'log2FC_shrunk']).copy()
    y = -np.log10(r.padj.clip(lower=1e-300))
    up = (r.padj < PADJ) & (r.log2FC_shrunk > LFC)
    dn = (r.padj < PADJ) & (r.log2FC_shrunk < -LFC)
    x = r.log2FC_shrunk.clip(-xlim, xlim)
    ax.scatter(x[~(up | dn)], y[~(up | dn)], s=1, color='#CFCFCF', lw=0, rasterized=True)
    ax.scatter(x[up], y[up], s=1.5, color='#B2182B', lw=0, alpha=0.5, rasterized=True)
    ax.scatter(x[dn], y[dn], s=1.5, color='#2166AC', lw=0, alpha=0.5, rasterized=True)
    for v in (LFC, -LFC):
        ax.axvline(v, color=MUTED, lw=0.4, ls='--')
    ax.axhline(-np.log10(PADJ), color=MUTED, lw=0.4, ls='--')
    ax.set_xlim(-xlim, xlim)
    if ylim:
        ax.set_ylim(-3, ylim)
    ax.set_title(title, fontsize=6.5)
    ax.text(0.02, 0.98, f'{int(dn.sum()):,}', color='#2166AC', transform=ax.transAxes,
            va='top', ha='left', fontsize=6)
    ax.text(0.98, 0.98, f'{int(up.sum()):,}', color='#B2182B', transform=ax.transAxes,
            va='top', ha='right', fontsize=6)


fig, axes = plt.subplots(2, 4, figsize=(6.0, 3.1), sharex=True, sharey=True,
                         gridspec_kw={'wspace': 0.12, 'hspace': 0.32})
for i, d in enumerate(DRUGS):
    for j, g in enumerate(GENOS):
        volcano(axes[i, j], drug_con(g, d), g, ylim=160)
        axes[i, j].title.set_color(GCOL[g])
        if j == 0:
            axes[i, j].set_ylabel(f'{TLAB[d]} vs DMSO\n−log$_{{10}}$(padj)')
for j in range(4):
    axes[1, j].set_xlabel('log$_2$FC')
save(fig, 'FigS_volcano_drug_response')

fig, axes = plt.subplots(1, 3, figsize=(4.6, 1.8), sharex=True, sharey=True, gridspec_kw={'wspace': 0.12})
for j, g in enumerate(VARS):
    volcano(axes[j], load(f'{g}_Control_vs_NT_Control'), f'{g} vs NT (DMSO)', xlim=6)
    axes[j].title.set_color(GCOL[g])
    axes[j].set_xlabel('log$_2$FC')
axes[0].set_ylabel('−log$_{10}$(padj)')
save(fig, 'FigS_volcano_baseline')


# ═════════════════════════════════════════════════════════════════════════════
# Fig S — sample–sample correlation (QC)
# ═════════════════════════════════════════════════════════════════════════════
top = vst.var(axis=1).sort_values(ascending=False).index[:1000]
cc = vst.loc[top, order].corr(method='pearson')
fig, ax = plt.subplots(figsize=(4.2, 3.8))
im = ax.imshow(cc.values, cmap='Blues', vmin=cc.values.min(), vmax=1, interpolation='nearest')
lab = [f"{meta.loc[s, 'genotype']} {TLAB[meta.loc[s, 'treatment']]} r{meta.loc[s, 'rep']}" for s in order]
ax.set_xticks(range(len(order)), lab, rotation=90, fontsize=4.5)
ax.set_yticks(range(len(order)), lab, fontsize=4.5)
for t, s in zip(ax.get_xticklabels(), order):
    t.set_color(GCOL[meta.loc[s, 'genotype']])
for t, s in zip(ax.get_yticklabels(), order):
    t.set_color(GCOL[meta.loc[s, 'genotype']])
for i in range(1, 4):
    ax.axvline(i * 9 - 0.5, color='white', lw=1)
    ax.axhline(i * 9 - 0.5, color='white', lw=1)
cb = fig.colorbar(im, ax=ax, shrink=0.5)
cb.set_label('Pearson r (VST, top 1,000 variable genes)', fontsize=6)
cb.outline.set_linewidth(0.4)
ax.tick_params(length=0)
save(fig, 'FigS_sample_correlation')
print('All figures done.')
