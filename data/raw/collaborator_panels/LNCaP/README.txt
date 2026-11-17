qPCR_heatmap/
  LNCaP_CDK12_qPCR_DNArepair_heatmap_source_data.txt
    Per-replicate and summary values (mean +/- SEM) for the qPCR DNA repair heatmap.
    Stars = t-test on dCq vs 0 hr, Bonferroni x3.
  make_qpcr_heatmap.py
    Draws the heatmap from the qPCR pipeline's intermediate files. Run after
    scripts/run_full_analysis.py in qPCR_DNA_repair_analysis:
      python make_qpcr_heatmap.py <qPCR_DNA_repair_analysis_dir> [outdir]
  LNCaP_CDK12_qPCR_DNArepair_heatmap.svg/.pdf   Heatmap as produced by the script.

RNAseq_panel_d/
  deseq2.pca.vals.txt   PC1/PC2 per sample (source data for panel d). Made by
    nf-core/rnaseq 3.25.0's DESeq2 QC step (variance-stabilised counts, top 500
    most variable genes), not by pyDESeq2. PC1 31%, PC2 15% variance.
  PCA.pdf               Plot, drawn by lncap_deseq2_gsea.py from the file above.

RNAseq_panel_e/
  4 DESeq2 result tables (source data for panel e; numerator first):
    NT_HQ461_vs_NT_Control.csv          780 DE (187 up / 593 down)
    I733V_HQ461_vs_I733V_Control.csv     39 DE (12 up / 27 down)
    NT_HQ461_vs_I733V_HQ461.csv         286 DE (122 up / 164 down)
    NT_Control_vs_I733V_Control.csv       1 DE (0 up / 1 down)
    DE = padj < 0.05 and |log2FC| > 1. log2FC = unshrunk MLE (shrinkage bug fixed).
    Columns: gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, gene_name.
  volcano_all_comparisons.pdf   Plot.
  lncap_deseq2_gsea.py          Revised main script (DESeq2 + volcano + PCA plot).
    Run from the analysis folder holding salmon.merged.gene_counts.tsv.

RNAseq_panel_f/
  LNCaP_NT_Hallmark_androgen_leadingedge_heatmap.svg/.pdf/.png
    Replacement for panel f. Leading-edge genes of MSigDB HALLMARK_ANDROGEN_RESPONSE
    (NT HQ461 vs Control), row z-score of log2 normalized counts.
  ..._source_data.txt   Values plotted, plus log2FC/padj and GSEA settings.
  androgen_heatmap_official.py   Script. Usage:
    python androgen_heatmap_official.py <analysis_dir> <hallmark.gmt> <outdir>
  h.all.v7.0.symbols.gmt   Hallmark gene sets used (MSigDB v7.0).
    Before publication, rerun with the current MSigDB Hallmark .gmt and cite that version.

RNAseq_panel_g/
  LNCaP_GSEA_GO_DSB_repair_HQ461_NT_vs_I733V.svg/.pdf/.png
    Panel g. Preranked GSEA on a single pre-specified set, GO:0006302 (double-strand break
    repair; GO BP 2021, Enrichr library; includes BRCA1, BRCA2, ERCC4).
    HQ461: NT vs Cdk12 I733V. NES -1.30, nominal p = 0.045.
  ..._source_data.txt   GSEA summary + ranked list with set membership and running ES.
  gsea_dsb_repair_official.py   Script. Usage:
    python gsea_dsb_repair_official.py <analysis_dir> <GO_BP_2021.gmt> <outdir>
  enrichr.GO_Biological_Process_2021.gmt   Gene set library used.
