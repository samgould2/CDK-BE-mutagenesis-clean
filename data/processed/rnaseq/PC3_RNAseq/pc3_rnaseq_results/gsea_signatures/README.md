# GSEA Signatures — PC3 CDK12 RNA-seq Analysis

Gene sets used across all GSEA analyses in this project.
All sets applied via `gseapy.prerank()` (1000 permutations, seed=42, min_size=5).

## Files

### `DDR_sets.gmt`
4 curated DDR/HR gene sets used for the DDR GSEA heatmaps (04_gsea_cdk12/).
Hardcoded inline in `scripts/regen_ddr_heatmaps.py`.

| Set | Source | n genes |
|-----|--------|---------|
| GOBP_DSB_REPAIR_VIA_HR | GO Biological Process (MSigDB) | 40 |
| KEGG_FANCONI_ANEMIA_PATHWAY | KEGG (MSigDB) | 37 |
| GOBP_DOUBLE_STRAND_BREAK_REPAIR | GO Biological Process (MSigDB) | 33 |
| GOBP_DNA_REPAIR | GO Biological Process (MSigDB) | 72 |

### `CDK12_13_custom_sets.gmt`
Curated CDK12/13 mechanistic signatures. Hardcoded in `scripts/pc3_deseq2_gsea.py`
and `scripts/make_gsea_svgs.py`. Used for GSEA_CDK12_13_heatmap and ssGSEA scoring.

| Set | Description | n genes |
|-----|-------------|---------|
| CDK12_LOSS_DN | Genes down-regulated upon CDK12 loss (DDR/replication) | 57 |
| CDK12_LOSS_UP | Genes up-regulated upon CDK12 loss (cGAS-STING/IFN) | 32 |
| CDK13_LOSS_DN | Genes down-regulated upon CDK13 loss (splicing/cell cycle) | 23 |
| CDK12_CDK13_SHARED_DN | Genes shared between CDK12 and CDK13 loss | 20 |
| HQ461_CDK12_CDK13_DN | Union of CDK12_LOSS_DN + CDK13_LOSS_DN + SHARED (sorted) | ~80 |
| CDK12_LONG_GENES | CDK12-sensitive long genes (>50 kb) used for elongation index | 24 |
| SHORT_HK_GENES_CTRL | Short housekeeping genes used as elongation index denominator | 13 |

### `hallmark_sets.gmt`
16 MSigDB Hallmark v2023.2 gene sets used in the Hallmark GSEA dot plots
(05_hallmark_gsea/). Hardcoded inline in `scripts/pc3_deseq2_gsea.py` and
`scripts/make_gsea_svgs.py` (subset; full MSigDB Hallmark library has 50 sets).

### `NT_drug_signatures/` (runtime-generated, not static GMT)
The NT drug sensitivity signatures are computed at runtime from the DESeq2 results:
- `NT_HQ461_DN`: genes with padj<0.05, LFC<-1 in NT_HQ461_vs_NT_Control
- `NT_HQ461_UP`: genes with padj<0.05, LFC>+1 in NT_HQ461_vs_NT_Control
- `NT_ZLY025_DN`: genes with padj<0.05, LFC<-1 in NT_ZLY025_vs_NT_Control
- `NT_ZLY025_UP`: genes with padj<0.05, LFC>+1 in NT_ZLY025_vs_NT_Control

The NT UP/DN intersection across both drugs (~1806 UP, ~2775 DN) defines the
core CDK inhibition signature used in NT_GSEA_enrichment and GSEA_preservation_score.

### Reactome (runtime-downloaded)
Reactome_2022 gene sets are downloaded at runtime via `gseapy.get_library('Reactome_2022')`.
The full library (~1600 gene sets) is not stored here; download requires internet access.

## Usage

All GMT files can be loaded directly into gseapy:
```python
import gseapy as gs
# From file:
result = gs.prerank(rnk=ranking, gene_sets='gsea_signatures/DDR_sets.gmt', ...)
# Or as dict (already in scripts):
result = gs.prerank(rnk=ranking, gene_sets=DDR_SETS, ...)
```
