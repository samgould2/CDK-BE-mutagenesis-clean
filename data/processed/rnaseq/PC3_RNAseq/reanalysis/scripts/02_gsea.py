"""Preranked GSEA (MSigDB Hallmark, all 50 sets + curated DDR/CDK12 sets) on every DESeq2 contrast.

Ranking metric: DESeq2 Wald statistic (unshrunk). Duplicate gene symbols collapsed to max |stat|.
Run:  cd reanalysis/scripts && python 02_gsea.py
"""
import os
import glob
import pandas as pd
import gseapy as gp

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
RES = os.path.join(BASE, 'results')
OUT = os.path.join(RES, 'gsea')
os.makedirs(OUT, exist_ok=True)
SIG_DIR = os.path.join(BASE, '..', 'pc3_rnaseq_results', 'gsea_signatures')


def read_gmt(path):
    sets = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            sets[p[0]] = [g for g in p[2:] if g]
    return sets


hallmark_path = os.path.join(OUT, 'MSigDB_Hallmark_2020.gmt')
if not os.path.exists(hallmark_path):
    hm = gp.get_library('MSigDB_Hallmark_2020', organism='Human')
    with open(hallmark_path, 'w') as f:
        for k, v in hm.items():
            f.write('\t'.join([k, 'MSigDB_Hallmark_2020'] + v) + '\n')
hallmark = read_gmt(hallmark_path)
curated = {**read_gmt(os.path.join(SIG_DIR, 'DDR_sets.gmt')),
           **read_gmt(os.path.join(SIG_DIR, 'CDK12_13_custom_sets.gmt'))}

files = sorted(glob.glob(os.path.join(RES, 'contrasts', '*.csv')) +
               glob.glob(os.path.join(RES, 'interaction', 'interaction_*.csv')))

rows = []
for fp in files:
    name = os.path.basename(fp).replace('.csv', '')
    df = pd.read_csv(fp).dropna(subset=['stat', 'gene_name'])
    df = df.loc[df['stat'].abs().groupby(df['gene_name']).idxmax()]
    rnk = df.set_index('gene_name')['stat'].sort_values(ascending=False)
    for lib_name, lib in [('hallmark', hallmark), ('curated', curated)]:
        pre = gp.prerank(rnk=rnk, gene_sets=lib, permutation_num=1000, seed=42,
                         min_size=10, max_size=1000, threads=6, outdir=None, verbose=False)
        r = pre.res2d.copy()
        r['contrast'] = name
        r['library'] = lib_name
        rows.append(r)
    print(name, 'done')

out = pd.concat(rows, ignore_index=True)
out = out.rename(columns={'Term': 'term', 'NOM p-val': 'nom_p', 'FDR q-val': 'fdr', 'FWER p-val': 'fwer'})
out[['contrast', 'library', 'term', 'ES', 'NES', 'nom_p', 'fdr', 'fwer', 'Tag %', 'Gene %', 'Lead_genes']].to_csv(
    os.path.join(OUT, 'gsea_all_contrasts.csv'), index=False)
print('Wrote', os.path.join(OUT, 'gsea_all_contrasts.csv'))
