# 01 — library design

The guide library itself, and how it was designed.

Nothing runs here. The library ships as data: `data/reference/library/` has
the final gRNA library, `data/reference/library_generation/` the CDK
sequences and predicted edits it was built from, and `cdkbe.counts.library()`
is how every other stage reads it. Figure 1e and 1f are drawn by
`cdkbe.source_data.library` and `library_figures` from those tables.

The design code itself is upstream of this repository and is not shipped.
