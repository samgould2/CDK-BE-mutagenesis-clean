"""Panel g: preranked GSEA on a single pre-specified set, GO:0006302 (double-strand break repair),
GO Biological Process 2021 (Enrichr library). Rank = DESeq2 Wald stat, HQ461: NT vs Cdk12 I733V.
Single set -> nominal p reported (no multiple-testing correction). Style matches plot_gsea_ddr_nt_vs_i733v.py.
Usage: python gsea_dsb_repair_official.py <analysis_dir> <GO_BP_2021.gmt> <outdir>"""
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt, gseapy as gs
matplotlib.rcParams.update({'font.family':['Avenir','Avenir Next','DejaVu Sans'],'font.weight':'medium',
                            'pdf.fonttype':42,'svg.fonttype':'none'})
BASE, GO, OUT = sys.argv[1:4]; os.makedirs(OUT, exist_ok=True)
TERM = 'double-strand break repair (GO:0006302)'
go = {l.split('\t')[0]: [g.strip() for g in l.rstrip('\n').split('\t')[2:] if g.strip()] for l in open(GO)}
r = (pd.read_csv(f'{BASE}/lncap_rnaseq_results/deseq2_results/NT_HQ461_vs_I733V_HQ461.csv')
       .dropna(subset=['stat','gene_name']).drop_duplicates('gene_name'))
rnk = r.set_index('gene_name')['stat'].sort_values(ascending=False)
pre = gs.prerank(rnk=rnk, gene_sets={TERM: go[TERM]}, permutation_num=1000, seed=42,
                 min_size=5, max_size=500, verbose=False, outdir=None)
t = pre.results[TERM]; es, hits, nes, p = np.array(t['RES']), np.array(t['hits']), t['nes'], t['pval']
x = np.arange(len(es))
fig, axes = plt.subplots(3, 1, figsize=(3.0, 2.8), sharex=True,
                         gridspec_kw={'height_ratios':[3,0.5,1.3], 'hspace':0.05})
ax = axes[0]; ax.plot(x, es, color='#2166ac', lw=1.5); ax.axhline(0, color='black', lw=0.6)
ax.set_ylabel('Enrichment score', fontsize=7)
ax.set_title('Double-strand break repair\nHQ461: NT vs Cdk12 I733V', fontsize=8)
ax.text(0.03, 0.05, f'NES = {nes:.2f}\np = {p:.3f}', transform=ax.transAxes, ha='left', va='bottom', fontsize=7)
ax.tick_params(labelsize=7); ax.spines[['top','right']].set_visible(False)
ax = axes[1]; ax.vlines(hits, 0, 1, color='black', lw=0.6)
ax.set_yticks([]); ax.spines[['top','right','left']].set_visible(False)
ax = axes[2]; ax.fill_between(x, rnk.values, 0, color='grey', lw=0); ax.axhline(0, color='black', lw=0.5)
ax.set_xlim(0, len(x)); ax.set_ylabel('Wald stat', fontsize=7); ax.set_xlabel('Gene rank', fontsize=7)
ax.tick_params(labelsize=7); ax.spines[['top','right']].set_visible(False)
ax.text(0.01, 0.95, 'higher in NT', transform=ax.transAxes, fontsize=6.5, va='top')
ax.text(0.99, 0.05, 'higher in I733V', transform=ax.transAxes, fontsize=6.5, ha='right')
stem = f'{OUT}/LNCaP_GSEA_GO_DSB_repair_HQ461_NT_vs_I733V'
for ext in ('svg','pdf','png'): fig.savefig(f'{stem}.{ext}', bbox_inches='tight', dpi=300)
with open(f'{stem}_source_data.txt','w') as f:
    f.write(f"# Gene set: {TERM}, GO Biological Process 2021 (Enrichr library); {len(go[TERM])} genes, {len(hits)} in ranked list\n"
            f"# Preranked GSEA, gseapy {gs.__version__}; rank = DESeq2 Wald stat, NT HQ461 vs Cdk12 I733V HQ461\n"
            f"# Single pre-specified set; 1000 permutations, seed 42\n"
            f"# NES = {nes:.3f}; nominal p = {p:.3g}\n# Leading edge: {t['lead_genes']}\n")
    d = pd.DataFrame({'rank':x+1,'gene':rnk.index,'wald_stat':rnk.values.round(4),'in_gene_set':0,'running_ES':es.round(5)})
    d.loc[hits,'in_gene_set'] = 1; d.to_csv(f, sep='\t', index=False)
print(f'NES={nes:.3f} p={p:.3g} hits={len(hits)}')
