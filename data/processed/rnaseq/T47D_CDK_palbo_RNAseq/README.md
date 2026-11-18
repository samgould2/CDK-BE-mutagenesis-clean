# T47D CDK variant ± palbociclib RNA-seq

The experiment behind **figure 4m** and **supplementary 13e–13g**: six T47D
lines — the two non-targeting controls and four base-edited CDK variants —
grown in DMSO or palbociclib, three replicates each, 36 libraries
(nf-core/rnaseq, STAR + Salmon).

    ABE_NT   CBE_NT   CDK2_S261G   CDK4_I51V   CDK6_G53S   CDK6_R46W

The two CDK6 lines each carry a second edit the sample name does not record,
and every published panel prints the longer spelling: `CDK6_G53S` is
`CDK6 G53S_M54I`, `CDK6_R46W` is `CDK6 R44C_R46W`.

## What is here

| file | what it is |
|---|---|
| `salmon.merged.gene_tpm.csv.gz` | Salmon gene-level TPM, 78,932 genes × 36 libraries |
| `deseq2_results/<variant>_vs_<control>_<veh\|palbo>.csv.gz` | the eight DESeq2 contrasts, each variant against its matched NT arm in that treatment |
| `hallmark_gene_sets_cache.json` | the MSigDB Hallmark sets the volcano highlight and the GSEA dot plot use |
| `marker_gene_panel.ipynb`, `pca_plot.ipynb`, `volcano_plots.ipynb` | the notebooks as delivered |
| `tables/marker_gene_zscore_{DMSO,Palbociclib}.csv` | the z-scores `marker_gene_panel.ipynb` wrote |

`veh` is DMSO throughout, in file names and sample names alike.

## Running it

**The three notebooks do not run here.** They read literal `.csv` paths
relative to their own folder, and the tables beside them are gzipped —
`docs/DISCREPANCIES.md` item 17 is the general version of this. They are kept
unchanged as the record of what was run.

The version that runs is `cdkbe.source_data.rnaseq_t47d_palbo`, which is their
plotting code in the package, reading through `cdkbe.io`:

```python
from cdkbe.source_data import rnaseq_t47d_palbo as rp

fig, zscores = rp.panel_4m()      # the marker-gene heatmap
fig, pca = rp.panel_S13e()        # the PCA
fig, counts = rp.panel_S13g()     # the volcano grid
rp.gsea_ranks()                   # S13f's input; the dot plot is not redrawn
```

`tests/test_t47d_cdk_palbo.py` compares all three against the published SVGs
in `data/raw/collaborator_panels/T47D_MCF7/`, character for character.

## Two things to know before using it

* **This is not the T47D RNA-seq in `../T47D_RNAseq/`.** That one is figure 5:
  four CDK4 genotypes, vehicle only, twelve libraries, a different GENCODE
  release. The two share three sample *names* — `T47D_CBE_NT_veh_1`, `_2`,
  `_3` — and none of the libraries behind them. Joining on sample name pairs
  unrelated data without failing; `docs/DISCREPANCIES.md` item 36 has the
  numbers.
* **The marker-gene panel z-scores within each treatment arm**, not across the
  whole matrix, so the DMSO → palbociclib shift is deliberately not what it
  shows.

Delivered as three standalone notebook folders; `data/MANIFEST.csv` records
each file's path inside that delivery. The duplicate copy of the TPM matrix
and two superseded z-score tables were dropped —
`docs/DISCREPANCIES.md` item 35.
