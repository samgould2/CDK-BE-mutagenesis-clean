"""PC3 panel k: intrinsic (Control) NT vs CDK12 variant DE, corrected.
Same DESeq2 setup as pc3_deseq2_gsea.py (filter >=10 counts in >=3 samples, ~condition, Cook's refit),
but log2FC is the unshrunk MLE from each Wald contrast. The original lfc_shrink(coeff='condition[T.<test>]')
shrank against the design reference (G731K_Control), giving wrong log2FC for any contrast whose
reference is not G731K_Control. padj and stat are unaffected.
Usage: python pc3_panel_k_intrinsic.py <salmon.merged.gene_counts.tsv> <outdir>"""
import os, re, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Liberation Sans','Arial','Helvetica','DejaVu Sans'],
                     'pdf.fonttype':42,'svg.fonttype':'none'})
COUNTS, OUT = sys.argv[1:3]; os.makedirs(OUT, exist_ok=True)
cnt = pd.read_csv(COUNTS, sep='\t', index_col=0); names = cnt.pop('gene_name'); cnt = cnt.round().astype(int)
def parse(c):
    m = re.match(r'X?\d+_(\w+?)_(Control|HQ461|ZLY)Rep\d+', c); g, t = m.group(1), m.group(2)
    return g, ('ZLY025' if t == 'ZLY' else t)
meta = pd.DataFrame([parse(c) for c in cnt.columns], columns=['genotype','treatment'], index=cnt.columns)
meta['condition'] = meta.genotype + '_' + meta.treatment
cf = cnt.loc[(cnt >= 10).sum(axis=1) >= 3]
try:    dds = DeseqDataSet(counts=cf.T, metadata=meta[['condition']], design='~condition', refit_cooks=True, quiet=True)
except TypeError: dds = DeseqDataSet(counts=cf.T, metadata=meta[['condition']], design_factors='condition', refit_cooks=True)
dds.deseq2()
COMPS = [('NT_Control','I733V_Control','NT_vs_I733V_Control','Control: NT vs I733V'),
         ('NT_Control','G822E_Control','NT_vs_G822E_Control','Control: NT vs G822E'),
         ('NT_Control','G731K_Control','NT_vs_G731K_Control','Control: NT vs G731K')]
res = {}
for test, ref, name, title in COMPS:
    s = DeseqStats(dds, contrast=['condition', test, ref], alpha=0.05, quiet=True); s.summary()
    r = s.results_df.join(names).reset_index().rename(columns={'index':'gene_id'}).sort_values('padj')
    r.to_csv(f'{OUT}/{name}.csv', index=False); res[name] = (r, title)

LFC_T, PADJ_T, FLOOR = 1.0, 0.05, 1e-300
allx, ally = [], []
for r, _ in res.values():
    rc = r.dropna(subset=['padj','log2FoldChange']); allx += list(rc.log2FoldChange); ally += list(-np.log10(rc.padj.clip(lower=FLOOR)))
xl = np.ceil(max(abs(min(allx)), abs(max(allx)))) + 0.5; yl = np.ceil(max(ally)) + 2
fig, axes = plt.subplots(1, 3, figsize=(16, 5)); counts = []
for ax, (name, (r, title)) in zip(axes, res.items()):
    rc = r.dropna(subset=['padj','log2FoldChange']).copy(); rc['padj'] = rc.padj.clip(lower=FLOOR)
    up = (rc.padj < PADJ_T) & (rc.log2FoldChange > LFC_T); dn = (rc.padj < PADJ_T) & (rc.log2FoldChange < -LFC_T); ns = ~(up | dn)
    ax.scatter(rc.loc[ns,'log2FoldChange'], -np.log10(rc.loc[ns,'padj']), s=4, color='#AAAAAA', alpha=0.4, linewidths=0, rasterized=True)
    ax.scatter(rc.loc[up,'log2FoldChange'], -np.log10(rc.loc[up,'padj']), s=25, color='tab:red', alpha=0.7, linewidths=0.3, edgecolors='black', rasterized=True, label=f'{up.sum()} UP')
    ax.scatter(rc.loc[dn,'log2FoldChange'], -np.log10(rc.loc[dn,'padj']), s=25, color='tab:blue', alpha=0.7, linewidths=0.3, edgecolors='black', rasterized=True, label=f'{dn.sum()} DOWN')
    ax.set_xlim(-xl, xl); ax.set_ylim(0, yl); ax.legend(fontsize=10, frameon=True, framealpha=0.9, edgecolor='#CCCCCC')
    for v in (LFC_T, -LFC_T): ax.axvline(v, color='black', lw=0.5, ls='--', alpha=0.4)
    ax.axhline(-np.log10(PADJ_T), color='black', lw=0.5, ls='--', alpha=0.4)
    ax.spines[['top','right']].set_visible(False); ax.tick_params(labelsize=11)
    ax.set_xlabel('log$_2$ Fold Change', fontsize=12); ax.set_ylabel('-log$_{10}$(padj)', fontsize=12); ax.set_title(title, fontsize=12)
    counts.append((name, int(up.sum()), int(dn.sum()), len(rc)))
plt.tight_layout()
for ext in ('svg','pdf','png'): fig.savefig(f'{OUT}/PC3_panel_k_volcano_intrinsic.{ext}', dpi=300, bbox_inches='tight')
for c in counts: print(*c)
