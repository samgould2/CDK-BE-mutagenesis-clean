# LNCaP CDK12 I733V x HQ461 — DESeq2 re-analysis
# 2 genotypes (NT, I733V) x 2 treatments (Control, HQ461) x 3 reps = 12 samples
# Library: 3' DGE (UMI-deduplicated, 3'-end reads, length correction disabled upstream)
#
# Mirrors the PC3 pipeline so the two cell lines are analysed identically:
#   1. ~ condition (4 levels)    -> pairwise contrasts
#   2. ~ genotype * treatment    -> interaction: does the HQ461 response differ in I733V?
#
# Run:  cd reanalysis/scripts && Rscript 01_deseq2.R

suppressPackageStartupMessages({
  library(DESeq2)
})

BASE <- Sys.getenv("LNCAP_BASE", normalizePath(file.path(getwd(), "..", "..")))
OUT  <- file.path(BASE, "reanalysis", "results")
dir.create(file.path(OUT, "contrasts"), recursive = TRUE, showWarnings = FALSE)

cts <- read.delim(file.path(BASE, "lncap_rnaseq_results", "salmon.merged.gene_counts.tsv"),
                  check.names = FALSE)
gene_names <- setNames(cts$gene_name, cts$gene_id)
mat <- round(as.matrix(cts[, -(1:2)]))
rownames(mat) <- cts$gene_id

# Columns: Cdk12I733V_Control_1 ... NT_HQ461_3
parse <- regmatches(colnames(mat), regexec("^(Cdk12I733V|NT)_(Control|HQ461)_(\\d)$", colnames(mat)))
stopifnot(all(lengths(parse) == 4))
meta <- data.frame(
  sample    = colnames(mat),
  genotype  = sub("^Cdk12", "", sapply(parse, `[`, 2)),
  treatment = sapply(parse, `[`, 3),
  rep       = as.integer(sapply(parse, `[`, 4)),
  stringsAsFactors = FALSE
)
meta$genotype  <- factor(meta$genotype,  levels = c("NT", "I733V"))
meta$treatment <- factor(meta$treatment, levels = c("Control", "HQ461"))
meta$condition <- factor(paste(meta$genotype, meta$treatment, sep = "_"),
                         levels = c("NT_Control", "NT_HQ461", "I733V_Control", "I733V_HQ461"))
rownames(meta) <- meta$sample
meta <- meta[order(meta$condition, meta$rep), ]
mat  <- mat[, meta$sample]
write.csv(meta, file.path(OUT, "sample_metadata.csv"), row.names = FALSE)

keep <- rowSums(mat >= 10) >= 3
mat  <- mat[keep, ]
cat("Genes retained:", nrow(mat), "of", nrow(cts), "\n")

dds <- DESeqDataSetFromMatrix(mat, meta, design = ~ condition)
dds <- DESeq(dds, quiet = TRUE)

vsd <- vst(dds, blind = FALSE)
write.csv(data.frame(gene_id = rownames(vsd), gene_name = gene_names[rownames(vsd)],
                     assay(vsd), check.names = FALSE),
          file.path(OUT, "vst.csv"), row.names = FALSE)
nc <- counts(dds, normalized = TRUE)
write.csv(data.frame(gene_id = rownames(nc), gene_name = gene_names[rownames(nc)], nc, check.names = FALSE),
          file.path(OUT, "normalized_counts.csv"), row.names = FALSE)

contrasts <- list(
  NT_HQ461_vs_NT_Control       = c("NT_HQ461", "NT_Control"),
  I733V_HQ461_vs_I733V_Control = c("I733V_HQ461", "I733V_Control"),
  I733V_Control_vs_NT_Control  = c("I733V_Control", "NT_Control"),
  I733V_HQ461_vs_NT_HQ461      = c("I733V_HQ461", "NT_HQ461")
)

summary_rows <- list()
for (nm in names(contrasts)) {
  cc <- contrasts[[nm]]
  r  <- results(dds, contrast = c("condition", cc[1], cc[2]), alpha = 0.05)
  rs <- lfcShrink(dds, contrast = c("condition", cc[1], cc[2]), res = r, type = "ashr", quiet = TRUE)
  df <- data.frame(gene_id = rownames(r), gene_name = gene_names[rownames(r)],
                   baseMean = r$baseMean, log2FC = r$log2FoldChange, lfcSE = r$lfcSE,
                   log2FC_shrunk = rs$log2FoldChange, stat = r$stat,
                   pvalue = r$pvalue, padj = r$padj)
  write.csv(df[order(df$pvalue), ], file.path(OUT, "contrasts", paste0(nm, ".csv")), row.names = FALSE)
  sig <- !is.na(df$padj) & df$padj < 0.05
  summary_rows[[nm]] <- data.frame(contrast = nm, numerator = cc[1], denominator = cc[2],
    up_padj05 = sum(sig & df$log2FC > 0), down_padj05 = sum(sig & df$log2FC < 0),
    up_lfc1 = sum(sig & df$log2FC_shrunk > 1), down_lfc1 = sum(sig & df$log2FC_shrunk < -1))
  cat(nm, "done\n")
}
write.csv(do.call(rbind, summary_rows), file.path(OUT, "DEG_summary.csv"), row.names = FALSE)

# Interaction: (I733V_HQ461 - I733V_Control) - (NT_HQ461 - NT_Control)
dds_int <- DESeqDataSetFromMatrix(mat, meta, design = ~ genotype + treatment + genotype:treatment)
dds_int <- DESeq(dds_int, quiet = TRUE)
coef <- "genotypeI733V.treatmentHQ461"
r  <- results(dds_int, name = coef, alpha = 0.05)
rs <- lfcShrink(dds_int, coef = coef, res = r, type = "ashr", quiet = TRUE)
df <- data.frame(gene_id = rownames(r), gene_name = gene_names[rownames(r)],
                 baseMean = r$baseMean, log2FC = r$log2FoldChange, lfcSE = r$lfcSE,
                 log2FC_shrunk = rs$log2FoldChange, stat = r$stat,
                 pvalue = r$pvalue, padj = r$padj)
write.csv(df[order(df$pvalue), ], file.path(OUT, "contrasts", "interaction_I733V_x_HQ461.csv"), row.names = FALSE)
cat("interaction padj<0.05:", sum(df$padj < 0.05, na.rm = TRUE), "\n")

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_deseq2.txt"))
cat("Done.\n")
