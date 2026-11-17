PC3 panel k - volcano plots, Control: NT vs CDK12 I733V / G822E / G731K

Source data: 3 DESeq2 tables (numerator = NT, reference = variant, all untreated)
  NT_vs_I733V_Control.csv   141 up / 52 down
  NT_vs_G822E_Control.csv   764 up / 429 down
  NT_vs_G731K_Control.csv   696 up / 363 down
  DE = padj < 0.05 and |log2FC| > 1. Columns: gene_id, baseMean, log2FoldChange, lfcSE,
  stat, pvalue, padj, gene_name. log2FC = unshrunk MLE from the Wald contrast.

Method: pyDESeq2 0.6.0 on salmon.merged.gene_counts.tsv (36 samples; 4 genotypes x
3 treatments x 3 reps), genes with >=10 counts in >=3 samples, ~condition (genotype_treatment),
Cook's refit, Wald test.

Script: python pc3_panel_k_intrinsic.py <salmon.merged.gene_counts.tsv> <outdir>

Correction vs the original analysis (pc3_deseq2_gsea.py): it called
lfc_shrink(coeff='condition[T.<test>]'), which shrinks against the design reference
(G731K_Control), not the stated reference. log2FC in NT_vs_I733V_Control and NT_vs_G822E_Control
were therefore effectively NT vs G731K (r = 0.999 with that file). padj/stat were unaffected and
are reproduced here (significance calls agree for >99.9% of genes). Old counts:
I733V 94/56, G822E 436/273, G731K 527/304 (G731K changes only because shrinkage is now off for all three).
