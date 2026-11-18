# 05 — bulk RNA-seq

Nothing runs here. The DESeq2 results, TPM matrices and GSEA tables for the
PC3, LNCaP and the two T47D experiments ship under `data/processed/rnaseq/`,
together with the collaborators' own analysis notes, and
`cdkbe.source_data.rnaseq`, `rnaseq_t47d`, `rnaseq_t47d_palbo` and
`rnaseq_volcanoes` draw the panels from them.

The two T47D trees are different experiments and share three sample names
without sharing the libraries behind them — `T47D_RNAseq/` is figure 5's CDK4
series, `T47D_CDK_palbo_RNAseq/` is figure 4m's six lines ± palbociclib, and
that one carries its own README. `docs/DISCREPANCIES.md` item 36 has the
numbers.

The DESeq2 and GSEA runs themselves are upstream of this repository.
`docs/METHODS.md` §5 says what each contrast is.
