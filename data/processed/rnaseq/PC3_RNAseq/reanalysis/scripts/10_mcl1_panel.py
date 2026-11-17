"""MCL1 expression panel, laid out to pair with a Western blot.

MCL1 is a canonical short-lived CDK9/CDK12-13 target. At 18.8 kb it sits mid-axis on the
gene-length distribution (decile 4 of 10), so its behaviour here is NOT a length effect --
it separates the two drugs directly.

One bar per condition in blot-lane order (genotype x treatment), per-replicate points shown,
significance from the DESeq2 drug-vs-same-genotype-DMSO contrast.

Run:  cd reanalysis/scripts && python 10_mcl1_panel.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import to_rgba
from _style import (RES, FIG, GENOS, DRUGS, TREATS, GCOL, TLAB, INK, MUTED, BAR_EDGE,
                    save, load, drug_con, meta as load_meta)

GENE, GENE_ID = 'MCL1', 'ENSG00000143384'
ALPHA = {'Control': 0.0, 'HQ461': 0.45, 'ZLY025': 1.0}   # fill depth by treatment

meta = load_meta()
nc = pd.read_csv(os.path.join(RES, 'normalized_counts.csv'), index_col='gene_id')
nc.pop('gene_name')
expr = nc.loc[GENE_ID, meta.index]


def stars(p):
    if not np.isfinite(p):
        return 'n.d.'
    return '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 0.05 else 'n.s.'


rows, rng = [], np.random.default_rng(0)
fig, ax = plt.subplots(figsize=(3.0, 2.3))
xt, xl = [], []
for i, g in enumerate(GENOS):
    base = expr[meta.index[(meta.genotype == g) & (meta.treatment == 'Control')]].mean()
    for j, t in enumerate(TREATS):
        x = i * 4 + j
        s = meta.index[(meta.genotype == g) & (meta.treatment == t)]
        v = expr[s]
        padj = drug_con(g, t).padj[GENE_ID] if t != 'Control' else np.nan
        rows.append({'genotype': g, 'treatment': t, 'mean_norm_count': v.mean(),
                     'sd': v.std(ddof=1), 'log2_vs_own_DMSO': np.log2(v.mean() / base),
                     'padj_vs_own_DMSO': padj,
                     **{f'rep{k + 1}': v.iloc[k] for k in range(len(v))}})
        ax.bar(x, v.mean(), width=0.8, color=to_rgba(GCOL[g], ALPHA[t]), **BAR_EDGE)
        ax.errorbar(x, v.mean(), yerr=v.std(ddof=1), color=INK, lw=0.7, capsize=1.6, capthick=0.7)
        ax.scatter(x + rng.uniform(-0.16, 0.16, len(v)), v, s=6, color='white',
                   edgecolor=INK, lw=0.4, zorder=3)
        if t != 'Control':
            ax.text(x, v.mean() + v.std(ddof=1) + 32, stars(padj), ha='center', va='bottom',
                    fontsize=5.5, color=INK)
        xt.append(x)
        xl.append(TLAB[t])
    ax.text(i * 4 + 1, 1035, g, ha='center', va='bottom', fontsize=7,
            color=GCOL[g], fontweight='bold')

ax.set_xticks(xt, xl, rotation=90, fontsize=5.5)
ax.set_xlim(-0.8, len(GENOS) * 4 - 1.2)
ax.set_ylim(0, 1120)
ax.set_ylabel(f'{GENE} expression\n(DESeq2 normalized counts)')
ax.axhline(0, color=INK, lw=0.6)
ax.legend(handles=[Patch(facecolor=to_rgba('#3D3D3D', ALPHA[t]), edgecolor='black',
                         lw=0.8, label=TLAB[t]) for t in TREATS],
          loc='upper left', bbox_to_anchor=(0.0, -0.30), ncol=3, fontsize=5.5,
          handlelength=1.1, columnspacing=1.2)
ax.set_title(f'{GENE} — stars: drug vs same-genotype DMSO', fontsize=6, pad=12)
save(fig, 'Fig12_MCL1_expression')

out = pd.DataFrame(rows)
out.to_csv(os.path.join(FIG, 'Fig12_MCL1_expression.csv'), index=False)
print(out.round(3).to_string(index=False))
