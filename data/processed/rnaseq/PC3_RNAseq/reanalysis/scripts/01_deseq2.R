# PC3 CDK12 variants x CDK12/13 degraders — DESeq2 re-analysis
# 4 genotypes (NT, I733V, G822E, G731K) x 3 treatments (Control, HQ461, ZLY025) x 3 reps
#
# Two models:
#   1. ~ condition (12 levels)          -> all pairwise contrasts (drug response, baseline, under-drug)
#   2. ~ genotype * treatment           -> interaction terms: does the drug response differ in variant vs NT?
# LFC shrinkage: ashr (apeglm not installed); unshrunk Wald stats kept for GSEA ranking.
#
# Run:  cd reanalysis/scripts && Rscript 01_deseq2.R

suppressPackageStartupMessages({
  library(DESeq2)
})

# BASE = "Comprehensive PC3 RNA-seq Analysis" folder; run from reanalysis/scripts or set PC3_BASE
BASE <- Sys.getenv("PC3_BASE", normalizePath(file.path(getwd(), "..", "..")))
OUT  <- file.path(BASE, "reanalysis", "results")
dir.create(file.path(OUT, "contrasts"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "interaction"), recursive = TRUE, showWarnings = FALSE)

# ── Counts + metadata ─────────────────────────────────────────────────────────
cts <- read.delim(file.path(BASE, "data", "salmon.merged.gene_counts.tsv"), check.names = FALSE)
gene_names <- setNames(cts$gene_name, cts$gene_id)
mat <- round(as.matrix(cts[, -(1:2)]))
rownames(mat) <- cts$gene_id

parse <- regmatches(colnames(mat), regexec("^X?(\\d+)_(\\w+?)_(Control|HQ461|ZLY)Rep(\\d)$", colnames(mat)))
meta <- data.frame(
  sample    = colnames(mat),
  num       = as.integer(sapply(parse, `[`, 2)),
  genotype  = sapply(parse, `[`, 3),
  treatment = sub("^ZLY$", "ZLY025", sapply(parse, `[`, 4)),
  rep       = as.integer(sapply(parse, `[`, 5)),
  stringsAsFactors = FALSE
)
meta$genotype  <- factor(meta$genotype,  levels = c("NT", "I733V", "G822E", "G731K"))
meta$treatment <- factor(meta$treatment, levels = c("Control", "HQ461", "ZLY025"))
meta$condition <- factor(paste(meta$genotype, meta$treatment, sep = "_"),
                         levels = as.vector(t(outer(levels(meta$genotype), levels(meta$treatment), paste, sep = "_"))))
rownames(meta) <- meta$sample
meta <- meta[order(meta$condition, meta$rep), ]
mat  <- mat[, meta$sample]
write.csv(meta, file.path(OUT, "sample_metadata.csv"), row.names = FALSE)

keep <- rowSums(mat >= 10) >= 3
mat  <- mat[keep, ]
cat("Genes retained:", nrow(mat), "\n")

# ── Model 1: condition ────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(mat, meta, design = ~ condition)
dds <- DESeq(dds, quiet = TRUE)
saveRDS(dds, file.path(OUT, "dds_condition.rds"))

vsd <- vst(dds, blind = FALSE)
vst_df <- data.frame(gene_id = rownames(vsd), gene_name = gene_names[rownames(vsd)], assay(vsd), check.names = FALSE)
write.csv(vst_df, file.path(OUT, "vst.csv"), row.names = FALSE)
nc <- counts(dds, normalized = TRUE)
write.csv(data.frame(gene_id = rownames(nc), gene_name = gene_names[rownames(nc)], nc, check.names = FALSE),
          file.path(OUT, "normalized_counts.csv"), row.names = FALSE)
write.csv(data.frame(sample = colnames(dds), size_factor = sizeFactors(dds)),
          file.path(OUT, "size_factors.csv"), row.names = FALSE)

genos <- c("NT", "I733V", "G822E", "G731K")
drugs <- c("HQ461", "ZLY025")
contrasts <- list()
for (g in genos) for (d in drugs)               contrasts[[paste0(g, "_", d, "_vs_", g, "_Control")]] <- c(paste0(g, "_", d), paste0(g, "_Control"))
for (g in genos[-1])                             contrasts[[paste0(g, "_Control_vs_NT_Control")]]     <- c(paste0(g, "_Control"), "NT_Control")
for (g in genos[-1]) for (d in drugs)            contrasts[[paste0(g, "_", d, "_vs_NT_", d)]]         <- c(paste0(g, "_", d), paste0("NT_", d))
contrasts[["NT_ZLY025_vs_NT_HQ461"]] <- c("NT_ZLY025", "NT_HQ461")

run_contrast <- function(num, den) {
  r  <- results(dds, contrast = c("condition", num, den), alpha = 0.05)
  rs <- lfcShrink(dds, contrast = c("condition", num, den), res = r, type = "ashr", quiet = TRUE)
  data.frame(gene_id = rownames(r), gene_name = gene_names[rownames(r)],
             baseMean = r$baseMean, log2FC = r$log2FoldChange, lfcSE = r$lfcSE,
             log2FC_shrunk = rs$log2FoldChange, stat = r$stat, pvalue = r$pvalue, padj = r$padj)
}

summary_rows <- list()
for (nm in names(contrasts)) {
  cc <- contrasts[[nm]]
  df <- run_contrast(cc[1], cc[2])
  write.csv(df[order(df$pvalue), ], file.path(OUT, "contrasts", paste0(nm, ".csv")), row.names = FALSE)
  sig <- !is.na(df$padj) & df$padj < 0.05
  summary_rows[[nm]] <- data.frame(contrast = nm, numerator = cc[1], denominator = cc[2],
    up_padj05 = sum(sig & df$log2FC > 0), down_padj05 = sum(sig & df$log2FC < 0),
    up_lfc1 = sum(sig & df$log2FC_shrunk > 1), down_lfc1 = sum(sig & df$log2FC_shrunk < -1))
  cat(nm, "done\n")
}
write.csv(do.call(rbind, summary_rows), file.path(OUT, "DEG_summary.csv"), row.names = FALSE)

# ── Model 2: genotype x treatment interaction ─────────────────────────────────
# Interaction coef genotypeX.treatmentD = (X_D - X_Control) - (NT_D - NT_Control)
# i.e. how much the drug response in variant X differs from the NT drug response.
dds_int <- DESeqDataSetFromMatrix(mat, meta, design = ~ genotype + treatment + genotype:treatment)
dds_int <- DESeq(dds_int, quiet = TRUE)
saveRDS(dds_int, file.path(OUT, "dds_interaction.rds"))
cat("Interaction coefs:", paste(resultsNames(dds_int), collapse = ", "), "\n")

int_rows <- list()
for (g in genos[-1]) for (d in drugs) {
  coef <- paste0("genotype", g, ".treatment", d)
  r  <- results(dds_int, name = coef, alpha = 0.05)
  rs <- lfcShrink(dds_int, coef = coef, res = r, type = "ashr", quiet = TRUE)
  df <- data.frame(gene_id = rownames(r), gene_name = gene_names[rownames(r)],
                   baseMean = r$baseMean, log2FC = r$log2FoldChange, lfcSE = r$lfcSE,
                   log2FC_shrunk = rs$log2FoldChange, stat = r$stat, pvalue = r$pvalue, padj = r$padj)
  nm <- paste0("interaction_", g, "_x_", d)
  write.csv(df[order(df$pvalue), ], file.path(OUT, "interaction", paste0(nm, ".csv")), row.names = FALSE)
  sig <- !is.na(df$padj) & df$padj < 0.05
  int_rows[[nm]] <- data.frame(contrast = nm, n_padj05 = sum(sig), up = sum(sig & df$log2FC > 0), down = sum(sig & df$log2FC < 0))
}

# Global LRT: does genotype modify drug response at all (any interaction term)?
dds_lrt <- DESeq(dds_int, test = "LRT", reduced = ~ genotype + treatment, quiet = TRUE)
r <- results(dds_lrt, alpha = 0.05)
lrt <- data.frame(gene_id = rownames(r), gene_name = gene_names[rownames(r)], baseMean = r$baseMean,
                  stat = r$stat, pvalue = r$pvalue, padj = r$padj)
write.csv(lrt[order(lrt$pvalue), ], file.path(OUT, "interaction", "LRT_any_genotype_x_treatment.csv"), row.names = FALSE)
int_rows[["LRT"]] <- data.frame(contrast = "LRT_any_interaction", n_padj05 = sum(r$padj < 0.05, na.rm = TRUE), up = NA, down = NA)
write.csv(do.call(rbind, int_rows), file.path(OUT, "interaction_summary.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_deseq2.txt"))
cat("Done.\n")
