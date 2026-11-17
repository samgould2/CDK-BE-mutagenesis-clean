# LNCaP CDK12 I733V RNA-seq — re-analysis (2026-10-31)

Re-analysis of the LNCaP dataset using the same DESeq2 + gene-length pipeline as the PC3
analysis (`../../Comprehensive PC3 RNA-seq Analysis/reanalysis/`), so the two cell lines are
processed identically.

## Design — note it is NOT the same as PC3

| | PC3 | LNCaP |
|---|---|---|
| Genotypes | NT, I733V, G822E, G731K | **NT, I733V only** |
| Treatments | DMSO, HQ461, ZLY-025 | **DMSO, HQ461 only** |
| Samples | 36 (4 × 3 × 3) | 12 (2 × 2 × 3) |
| Library | full-length bulk RNA-seq | **3′ DGE, UMI-deduplicated, length correction disabled** |

The PC3 `Fig10a` has three panels (HQ461 / ZLY-025 / baseline). LNCaP has no ZLY-025 arm, so the
recreated figure has two: HQ461 vs DMSO, and baseline I733V vs NT.

## Pipeline

```bash
cd reanalysis/scripts
Rscript 01_deseq2.R                 # DESeq2 contrasts + genotype×treatment interaction
python  02_elongation_shared_y.py   # the Fig10a recreation
```

## The resistance phenotype reproduces cleanly

| Contrast | up (padj<0.05) | down | up \|LFC\|>1 | down \|LFC\|>1 |
|---|---|---|---|---|
| NT HQ461 vs DMSO | 2,166 | 2,114 | 108 | 275 |
| I733V HQ461 vs DMSO | 116 | 281 | 6 | 3 |
| I733V vs NT (DMSO) | 280 | 150 | 0 | 0 |

The drug clearly works and I733V clearly blocks it. Canonical responses in NT, all abolished in
I733V: **CDKN1A +1.50** (padj 9×10⁻⁶³), **KLK3 −1.20** (3×10⁻⁶¹), **BRCA1 −0.95** (9×10⁻⁷),
AR −0.23 (8×10⁻³). In I733V the same genes are flat (BRCA1 −0.08, CDKN1A +0.12, both n.s.).
946 genes have a significant genotype×treatment interaction.

Baseline I733V vs NT has **zero** genes at |LFC| > 1, consistent with PC3 where I733V was the
mild variant with no baseline elongation defect.

---

## ⚠️ The gene-length result is not interpretable in this dataset

`Fig10a_elongation_shared_y_LNCaP` is flat. **This should not be reported as "LNCaP has no
elongation defect."** The library type differs from PC3 in a way that bears directly on the
measurement, and the comparison is confounded.

Top-decile response gap (same metric as PC3):

| Contrast | LNCaP gap (95% CI) | PC3 equivalent |
|---|---|---|
| NT HQ461 vs DMSO | **−0.024** (−0.045 to −0.010) | **−0.412** |
| I733V HQ461 vs DMSO | +0.002 (−0.001 to +0.006) | −0.044 |
| I733V vs NT (DMSO) | −0.003 (−0.007 to +0.001) | −0.004 |

The effect is **absent, not merely smaller**. Normalising for overall response magnitude does not
rescue it — the slope of log₂FC on log₁₀(span) divided by SD(log₂FC) is **0.018** in LNCaP versus
**−0.588** in PC3, a ~33× difference after scaling. The zoomed panel
(`Fig10a_elongation_zoom_LNCaP`, ~10× y-scale) shows the profile is a shallow hump, not a damped
monotonic decline, with replicate error bars comparable to the signal.

### What was ruled out

- **Failed experiment** — no. 4,280 genes change at padj < 0.05 and the canonical targets above
  all move in NT and are blocked in I733V.
- **Annotation mismatch** — no. Unversioned Ensembl IDs, same Ensembl 113 annotation as PC3,
  56,711 of 63,967 genes annotated; 9,914 expressed protein-coding genes pass the same filters
  (PC3: 8,547).
- **Baseline length–abundance coupling differing between libraries** — no. corr(log₁₀ baseMean,
  log₁₀ span) is −0.239 in LNCaP vs −0.260 in PC3; corr with mRNA length −0.303 vs −0.357.
  The two datasets have a similar standing length–abundance relationship.

### Candidate explanations, not distinguished by this data

1. **3′ DGE may be blind to the mechanism.** The elongation readout assumes a shortened
   transcript loses gene-level signal. In a 3′-end assay, one transcript gives one read at its
   poly(A) site; premature cleavage and polyadenylation at an intronic site still produces a
   polyadenylated 3′ end. Depending on whether those reads map to annotated transcripts, the
   transcript may still be counted — in which case the assay cannot see the event that
   full-length RNA-seq sees as lost 3′ coverage.
2. **Much weaker overall response.** SD(log₂FC) is 0.296 in LNCaP vs 1.188 in PC3 — a 4× smaller
   dynamic range. Dose and timepoint for the two experiments were not documented in the material
   I read; if they differ, the comparison is not like-for-like.
3. **Genuine cell-line biology** — different CDK12/13 dependence in LNCaP vs PC3.

I cannot separate these with the data in hand. The honest conclusion is that **this dataset
cannot test the elongation hypothesis**, not that the hypothesis fails here.

### What would resolve it

A 5′→3′ coverage metagene is not available for a 3′-end library. The decisive options are
(a) confirm the dose and timepoint match the PC3 experiment, and (b) if the elongation question
matters for LNCaP, generate full-length RNA-seq — 3′ DGE is the wrong assay for it regardless of
how the analysis is done.

Everything else in this dataset (DEG counts, canonical targets, the interaction model, the
resistance phenotype) is unaffected by this caveat and is directly comparable to PC3.

## Outputs

```
reanalysis/
├── results/
│   ├── contrasts/      4 pairwise contrasts + interaction_I733V_x_HQ461 (CSV)
│   ├── annotation/     Ensembl 113 gene features (copied from the PC3 analysis)
│   ├── vst.csv, normalized_counts.csv, sample_metadata.csv
│   ├── DEG_summary.csv, elongation_length_gaps.csv
└── figures/            Fig10a_elongation_shared_y_LNCaP, Fig10a_elongation_zoom_LNCaP (+ CSV)
```

PDFs are transparent-background, Type-42 text; PNGs are 300 dpi on white. Palette matches PC3
(NT dark grey, I733V blue).
