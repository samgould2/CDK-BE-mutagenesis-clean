"""Panel f rebuilt with the official MSigDB Hallmark Androgen Response set.
1) Preranked GSEA (DESeq2 Wald stat, NT HQ461 vs NT Control) against all 50 Hallmarks.
2) Heatmap of the leading-edge genes: row z-score of log2(median-of-ratios norm counts + 1), NT samples.
Usage: python androgen_heatmap_official.py <analysis_dir> <hallmark.gmt> <outdir>"""
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns, gseapy as gp
matplotlib.rcParams.update({'font.family':['Avenir','Avenir Next','DejaVu Sans'],'font.weight':'medium',
                            'pdf.fonttype':42,'svg.fonttype':'none'})
BASE, GMT, OUT = sys.argv[1:4]; os.makedirs(OUT, exist_ok=True)
RES = os.path.join(BASE,'lncap_rnaseq_results'); TERM='HALLMARK_ANDROGEN_RESPONSE'
sets = {l.split('\t')[0]: [g.strip() for g in l.rstrip('\n').split('\t')[2:] if g.strip()] for l in open(GMT)}

de = (pd.read_csv(f'{RES}/deseq2_results/NT_HQ461_vs_NT_Control.csv')
        .dropna(subset=['stat','gene_name']).drop_duplicates('gene_name').set_index('gene_name'))
rnk = de['stat'].sort_values(ascending=False)
pre = gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=1000, seed=42,
                 min_size=5, max_size=500, verbose=False, outdir=None, threads=4)
r = pre.res2d.set_index('Term').loc[TERM]
lead = [g for g in str(r['Lead_genes']).split(';') if g]
print(f"{TERM}: NES={float(r['NES']):.2f} FDR={float(r['FDR q-val']):.3g} "
      f"NOM p={float(r['NOM p-val']):.3g} tested={r['Tag %']} leading edge n={len(lead)}")

cnt = pd.read_csv(os.path.join(BASE,'salmon.merged.gene_counts.tsv'), sep='\t', index_col=0)
names = cnt.pop('gene_name'); cnt = cnt.round()
lg = np.log(cnt[(cnt>0).all(axis=1)]); sf = np.exp(lg.sub(lg.mean(axis=1),axis=0).median())
norm = np.log2(cnt/sf+1); norm.index = names.values; norm = norm.groupby(level=0).mean()
samples = [c for g in ('NT_Control','NT_HQ461') for c in norm.columns if c.startswith(g)]
genes = sorted([g for g in lead if g in norm.index], key=lambda g: de.loc[g,'log2FoldChange'])
m = norm.loc[genes, samples]; z = m.sub(m.mean(axis=1),axis=0).div(m.std(axis=1),axis=0)

fig, ax = plt.subplots(figsize=(2.6, 0.17*len(genes)+0.9))
sns.heatmap(z.values, ax=ax, cmap='RdBu_r', center=0, vmin=-2, vmax=2, linewidths=0.8, linecolor='black',
            cbar=True, cbar_kws={'shrink':0.4,'pad':0.04}, yticklabels=genes)
ax.set_xticks([1.5,4.5]); ax.set_xticklabels(['Control','HQ461'], fontsize=8, rotation=0)
ax.set_yticklabels(genes, fontsize=8, rotation=0, style='italic'); ax.tick_params(length=0)
ax.set_title('Hallmark Androgen Response\n(leading edge, NT)', fontsize=8)
cb = ax.collections[0].colorbar; cb.set_label('Row z-score', fontsize=7); cb.ax.tick_params(labelsize=7)
plt.tight_layout(pad=0.4)
stem = f'{OUT}/LNCaP_NT_Hallmark_androgen_leadingedge_heatmap'
for ext in ('svg','pdf','png'): fig.savefig(f'{stem}.{ext}', bbox_inches='tight', dpi=300)

# source data
with open(f'{stem}_source_data.txt','w') as f:
    f.write(f"# Gene set: {TERM}, MSigDB Hallmark v7.0 ({len(sets[TERM])} genes)\n"
            f"# GSEA preranked (gseapy {gp.__version__}), rank = DESeq2 Wald stat, NT HQ461 vs NT Control; "
            f"all 50 Hallmarks, 1000 permutations, seed 42, size 5-500\n"
            f"# NES = {float(r['NES']):.3f}; FDR q = {float(r['FDR q-val']):.3g}; nominal p = {float(r['NOM p-val']):.3g}\n"
            f"# Heatmap rows = leading-edge genes, ordered by log2FC (most decreased first)\n"
            f"# Values = row z-score of log2(DESeq2 median-of-ratios normalized counts + 1), NT samples only\n")
    out = z.copy(); out.columns = [c+'_zscore' for c in samples]
    out = pd.concat([m.add_suffix('_log2norm'), out, de.loc[genes,['log2FoldChange','padj']]], axis=1)
    out.index.name = 'gene'; out.round(4).to_csv(f, sep='\t')
print(len(genes),'genes plotted:', ', '.join(genes))
