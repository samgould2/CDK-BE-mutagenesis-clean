"""
Map the apparent 9p21.3 loss: expression of every gene on chr9:19-26 Mb in the
three CDK4 variant lines relative to CBE NT, ordered by genomic position, so the
boundaries of the affected interval can be read off directly.

Gene coordinates (GRCh38) come from the Ensembl REST API, cached in
results/chr9p21_genes.json.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
OUT  = HERE / 'results'
FIG  = HERE / 'figures'

CTRL     = 'CBE_NT'
VARIANTS = ['CDK4_D196N', 'CDK4_E56K', 'CDK4_W198']
LABEL = {'CDK4_D196N': 'D196N', 'CDK4_E56K': 'E56K', 'CDK4_W198': 'W198*'}
COLOR = {'CDK4_D196N': '#4393c3', 'CDK4_E56K': '#66bd63', 'CDK4_W198': '#d6604d'}

tpm = pd.read_csv(BASE / 'salmon.merged.gene_tpm.tsv', sep='\t', index_col=0)
counts = pd.read_csv(BASE / 'salmon.merged.gene_counts_scaled.tsv', sep='\t', index_col=0)
gene_name = tpm.pop('gene_name'); counts.pop('gene_name')
meta = pd.Series(tpm.columns.str.extract(r'^T47D_(.+)_veh_\d+$')[0].values, index=tpm.columns)
grp = lambda g: list(meta.index[meta == g])
gm = tpm.T.groupby(meta).mean().T

genes = pd.DataFrame(json.loads((OUT / 'chr9p21_genes.json').read_text()))
tpm_noversion = pd.Series(tpm.index, index=[i.split('.')[0] for i in tpm.index])
genes['matrix_id'] = genes.gene_id.map(tpm_noversion)
genes = genes.dropna(subset=['matrix_id'])

# keep genes with enough signal somewhere to be informative about copy number
expressed = genes[gm.loc[genes.matrix_id].max(axis=1).values > 1.0].copy()
expressed['ctrl_TPM'] = gm.loc[expressed.matrix_id, CTRL].values
for v in VARIANTS:
    expressed[v] = gm.loc[expressed.matrix_id, v].values
    expressed[v + '_log2ratio'] = np.log2((expressed[v] + 0.1) / (expressed.ctrl_TPM + 0.1))
expressed = expressed.sort_values('start')
expressed.to_csv(OUT / 'chr9p21_expression_map.csv', index=False)

print(f'{len(expressed)} genes on chr9:19-26 Mb with mean TPM > 1 in at least one genotype\n')
print(f"{'symbol':14s} {'Mb':>7s} {'biotype':16s} {'NT':>8s} " + ' '.join(f'{LABEL[v]:>8s}' for v in VARIANTS))
for r in expressed.itertuples():
    sym = r.symbol or r.gene_id
    print(f'{sym:14s} {r.start/1e6:7.3f} {r.biotype:16s} {r.ctrl_TPM:8.2f} ' +
          ' '.join(f'{getattr(r, v):8.2f}' for v in VARIANTS))

# raw counts for the core locus - is it zero, or just low?
print('\nRaw (scaled) counts per sample at the core locus:')
for sym in ['KLHL9', 'MTAP', 'CDKN2A', 'CDKN2B']:
    hits = gene_name[gene_name == sym].index
    for gid in hits:
        if gid in counts.index:
            print(f'  {sym:8s} ' + ' '.join(f'{counts.loc[gid, s]:8.1f}' for s in tpm.columns))
print('  samples: ' + ' '.join(f'{s.replace("T47D_","").replace("_veh",""):>8s}' for s in tpm.columns))

# ── figure: log2 ratio along the chromosome ─────────────────────────────────
fig, ax = plt.subplots(figsize=(9.5, 3.6))
x = expressed.start.values / 1e6
for v in VARIANTS:
    ax.plot(x, expressed[v + '_log2ratio'], marker='o', ms=4, lw=.9, color=COLOR[v], label=LABEL[v], alpha=.9)
ax.axhline(0, color='black', lw=.8)
ax.axhline(-1, color='#999999', lw=.8, ls='--')
lo, hi = 20.45, 22.00   # bounded by MLLT3 (retained) and CDKN2B-AS1 (retained)
ax.axvspan(lo, hi, color='#d6604d', alpha=.10, zorder=0)
ax.text((lo + hi) / 2, ax.get_ylim()[1], 'affected interval (~1.55 Mb)', ha='center', va='top', fontsize=8, color='#a33')
for sym in ['KLHL9', 'MTAP', 'CDKN2A', 'CDKN2B', 'IFNB1', 'TUSC1', 'ELAVL2', 'DMRTA1', 'MLLT3', 'FOCAD']:
    row = expressed[expressed.symbol == sym]
    if len(row):
        r = row.iloc[0]
        ax.annotate(sym, (r.start / 1e6, r[VARIANTS[2] + '_log2ratio']), fontsize=7, style='italic',
                    xytext=(0, -11), textcoords='offset points', ha='center')
ax.set_xlabel('chr9 position (Mb, GRCh38)', fontsize=9)
ax.set_ylabel('log$_2$ (variant / CBE NT)\nmean TPM', fontsize=9)
ax.legend(fontsize=8, frameon=False, loc='lower right')
ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=8)
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'chr9p21_expression_map.{e}', dpi=300, bbox_inches='tight')
plt.close()

# ── CDK4/6 pathway compensation panel ───────────────────────────────────────
PANEL = ['CDK4', 'CDK6', 'CCND1', 'CCND3', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CDKN2D',
         'CDKN1A', 'CDKN1B', 'RB1', 'CCNE1', 'CCNE2', 'CDK2', 'E2F1', 'MYC']
order = tpm.mean(axis=1).sort_values(ascending=False).index
sym2id = pd.Series(order, index=gene_name.reindex(order).values); sym2id = sym2id[~sym2id.index.duplicated()]
panel = pd.DataFrame({g: gm.loc[[sym2id[s] for s in PANEL if s in sym2id.index], g].values
                      for g in [CTRL] + VARIANTS},
                     index=[s for s in PANEL if s in sym2id.index])
panel.to_csv(OUT / 'cdk4_pathway_panel_TPM.csv')
print('\nCDK4/6-RB pathway (mean TPM):')
print(panel.round(1).to_string())
