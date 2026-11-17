# Methods — gene-length ("elongation") quantification

Covers `Fig7_elongation_length`, `Fig9_baseline_elongation`, `Fig10*`,
`FigS_elongation_length_and_introns`, and the `LONG_GENES` / `INTRON_RICH` sets in
`Fig6_targeted_GSEA`.

Scripts: `04_elongation.py`, `05_elongation_length.py`, `06_elongation_combined.py`,
`07_elongation_shared_y.py`, `08_elongation_controls.py`; shared helpers in `_style.py`.

---

## 1. What is being measured, and what it is a proxy for

CDK12/13 phosphorylate the RNA Pol II CTD and promote processive elongation. When they are
inhibited, Pol II processivity drops and transcription is more likely to terminate prematurely —
canonically through use of intronic poly(A) sites. The predicted consequence at steady-state
mRNA level is a **length-dependent response**: transcripts that require the most elongation to
complete (long genes, intron-rich genes) are depleted, while short, intron-poor transcripts are
relatively spared or increased.

The quantification here measures exactly that prediction: **how the measured fold-change depends
on gene length**. It does *not* directly observe premature termination. See §7.

---

## 2. Gene annotation

Source: Ensembl release 113 GTF (`Homo_sapiens.GRCh38.113.gtf.gz`), which matches the GRCh38
Ensembl annotation the upstream nf-core/rnaseq Salmon quantification was built against.

Per gene, from the GTF `gene` feature:

- **`span_bp` = end − start + 1** — the genomic footprint of the locus, introns included.
  This is the primary axis. It is the right quantity for an elongation question: it is the
  distance Pol II must actually traverse. Mature transcript length is *not* the relevant axis
  and is not used as the primary one (it is recorded as `mrna_len`).

From the transcript tagged `Ensembl_canonical` (one per gene; present for all 78,932 genes):

- **`n_exons`**, **`n_introns` = n_exons − 1**, **`mrna_len`** (summed exon length),
  **`intronic_len` = span_bp − mrna_len**.

Using the single canonical transcript rather than a union or expression-weighted model is a
simplification: it makes intron counts well defined and reproducible, but it ignores
isoform usage in PC3 specifically.

Stored at `results/annotation/gene_features.tsv`.

---

## 3. Gene universe

Three filters, applied in order (`_style.expressed_pc`):

1. **`baseMean > 50`** — DESeq2 `baseMean` is the mean of normalized counts across **all 36
   samples**, not just the two groups of a contrast, so this filter is not biased toward any
   genotype or treatment. It removes genes where a decile median would be dominated by
   counting noise.
2. **`biotype == protein_coding`** — lncRNAs and small RNAs have very different length and
   processing behaviour and would confound a length axis.
3. Must have a length annotation.

Result: **8,547 genes** (of 16,290 in the DESeq2 output; 9,417 pass the `baseMean` filter alone).

---

## 4. Binning

Genes are split into **10 equal-count deciles** (855/854 genes each) by the binning variable:

```python
b = pd.qcut(values.rank(method='first'), 10, labels=False)
```

Ranking before `qcut` guarantees exactly equal-sized bins even with tied values. Each decile is
plotted at the **median** of its binning variable, so the x-position is a real data value rather
than a bin edge.

Length deciles span medians of **3.7 kb → 226.9 kb** (full range 207 bp → 1.44 Mb).

> **Tie caveat, intron axis only.** Intron count takes only 93 distinct values across 8,547
> genes, so `rank(method='first')` splits tied counts across adjacent decile boundaries
> arbitrarily (e.g. decile 0 spans 0–2 introns and decile 1 spans 2–4 — genes with exactly 2
> introns fall in both). Equal bin sizes are preserved but boundaries are not clean. Gene length
> has 8,204 distinct values of 8,547, so this does not apply to the primary axis.

---

## 5. Response metric — two variants

Both are used, deliberately, and they are not interchangeable.

### 5a. Contrast-level shrunken LFC (Fig7, Fig9, FigS, Fig6 GSEA)

The per-gene `log2FC_shrunk` from DESeq2 (`ashr` shrinkage), summarised as the **median within
each decile**. Median rather than mean because per-gene LFC distributions are heavy-tailed.

Uncertainty originally shown as the interquartile range across genes within a decile. **This was
dropped** — the IQR is gene-to-gene biological spread, not uncertainty in the decile estimate,
and it is enormous relative to the signal, which made the panels look far noisier than the
estimate actually is.

### 5b. Replicate-level VST ratios (Fig10, Fig10a, Fig10b)

A DESeq2 shrunken LFC is a single contrast-level estimate with no per-replicate decomposition,
so replicate error bars require a different quantity. For each replicate *i* of the test group:

```
log2ratio_i(gene) = VST(gene, replicate_i) − mean(VST(gene, reference replicates))
```

then the decile median is taken **per replicate**, giving three values per decile, reported as
**mean ± SD (n = 3)**. Reference group is the same genotype's DMSO samples for drug contrasts,
and NT DMSO for baseline contrasts.

**These two metrics are not numerically interchangeable.** VST shrinks low-count genes toward the
mean, so magnitudes are compressed: NT's shortest decile under HQ461 reads **+0.95** as a
shrunken LFC but **+0.71** as a VST ratio (~25% smaller). Shape, sign and genotype ordering are
identical. Read Fig10 within itself; do not quote its values against Fig7/Fig9.

Replicate SDs of the top-decile median are **0.012–0.081 log₂** — i.e. these profiles are highly
reproducible, far more so than the gene-level IQR implied.

---

## 6. Summary statistics and uncertainty

Two scalar summaries per condition, both in `results/elongation_slopes.csv` and
`results/elongation_length_gaps.csv`:

| Metric | Definition | Notes |
|---|---|---|
| **Elongation slope** | OLS slope of per-gene log2FC on log10(span), all genes | interpretable as Δlog₂FC per 10× length; assumes monotonicity |
| **Top-decile response gap** (used in figures) | median(log2FC \| top decile) − median(log2FC \| all genes) | no monotonicity assumption |

The **gap is the primary metric**. The slope is reported but is a poor summary where the profile
is non-monotonic — NT's intron profile under ZLY-025 dips at single-intron genes and then
declines, which a straight line fits badly.

The gap is referenced to the **overall median rather than the bottom decile** because
single-intron genes behave atypically under ZLY-025 and would otherwise dominate the contrast.

A useful consequence: because DESeq2 median-of-ratios normalization pins roughly the median gene
near zero (observed global median log2FC: −0.086 HQ461, −0.024 ZLY-025), the gap is **invariant to
any global compositional shift**. Only the *relative* difference across deciles is
normalization-independent; the absolute height of a curve is not. This matters because thousands
of genes change under these drugs, which strains the median-of-ratios assumption that most genes
do not. Interpret curve *shape*, not curve *offset*.

**Confidence intervals**: nonparametric bootstrap over genes, 2,000 resamples for the gap and
1,000 for the slope, percentile method (2.5/97.5), seeded `np.random.default_rng(0)`. These
intervals reflect gene-sampling uncertainty, **not** biological replicate uncertainty — for that,
see the n = 3 error bars in §5b. The two answer different questions and neither subsumes the other.

---

## 7. Controls (`08_elongation_controls.py`)

Gene length is not independent of other gene properties, so the trend was tested against the two
most plausible confounders.

**Length is correlated with expression** (*r* = −0.260 for log₁₀ span vs log₁₀ baseMean: longer
genes are less expressed) **and with intron number** (*r* = 0.569).

### Expression stratification

Length tertiles were recomputed **inside** each expression tertile, so only genes of similar
abundance are compared. Long-minus-short tertile median log₂FC:

| Condition | Global | Low expr | Mid expr | High expr |
|---|---|---|---|---|
| NT HQ461 | −1.021 | −0.570 | −0.655 | −0.855 |
| NT ZLY-025 | −0.780 | −0.236 | −0.341 | −0.607 |
| I733V HQ461 | −0.002 | +0.092 | +0.097 | −0.021 |
| G731K ZLY-025 | −0.026 | −0.000 | −0.013 | −0.038 |

The length effect **survives in every expression stratum** in NT, and the genotype ordering is
preserved. But it is **attenuated within strata** (NT HQ461: −0.57 to −0.86 vs −1.02 globally)
and is **clearly stronger in highly expressed genes**. So part of the global effect size is
carried by the length–expression correlation, and the effect is not uniform across the
expression range. Both are honest qualifications of the headline number.

> A methodological trap worth recording: the top-decile-vs-median **gap cannot be used within a
> stratum**. Because length and expression are anticorrelated, the median gene in the
> low-expression stratum is already long, so the gap collapses toward zero for a reason that has
> nothing to do with biology. Within-stratum top-vs-bottom tertile is the correct comparison.

### Multiple regression

Coefficient on log10(span) for NT as covariates are added:

| Model | β span | β expr | β introns |
|---|---|---|---|
| ~ span | −0.699 | | |
| ~ span + expr | −0.477 | +0.856 | |
| ~ span + introns | −0.723 | | +0.070 |
| ~ span + introns + expr | −0.522 | +0.858 | +0.129 |

Length retains ~**75%** of its coefficient when expression is controlled (−0.70 → −0.52) and is
essentially unchanged by intron number. Notably, once length is in the model the intron
coefficient is small and **positive**, so on these data intron number adds little independent
explanatory power beyond length — the two axes in `FigS_elongation_length_and_introns` are
largely telling the same story, not two independent ones. Expression level itself is a strong
positive predictor (β ≈ +0.86), which is why the stratified control above matters.

---

## 8. Limitations

1. **This is an inference, not a direct observation of premature termination.** A
   length-dependent decrease in steady-state mRNA is *consistent with* loss of processivity, but
   gene-level counts cannot distinguish it from length-correlated differences in mRNA stability,
   co-transcriptional processing, or decay. Nothing here localises where Pol II stops.
2. **Gene-level quantification only.** The direct readouts would be within-gene 5′→3′ coverage
   decay, intronic poly(A) site usage, or a last-exon/first-exon ratio. The Salmon output here is
   gene-level, so none of those were computable. Transcript-level or exon-level quantification
   from the same BAMs would permit them.
3. **Compositional normalization.** As in §6, absolute curve height is not identifiable; only
   shape is. The apparent *up*-regulation of short genes is partly a normalization consequence of
   large-scale downregulation elsewhere, and should not be reported as an absolute increase
   without spike-in normalization.
4. **Expression confound is real and only partly controlled** (§7).
5. **Canonical-transcript annotation** ignores PC3-specific isoform usage.
6. **Intron-axis ties** (§4).
7. **n = 3** per condition. Replicate SDs are small for the decile medians, but the baseline
   effects are small in absolute terms and individual deciles are noisy (see Fig10's baseline
   panel error bars).

## 9. What would make this direct

In rough order of effort: (a) 5′→3′ coverage metagene and last/first exon ratios from the
existing BAMs, stratified by gene length — cheapest, uses data already in hand; (b) 3′-end
sequencing or analysis of intronic poly(A) site usage to directly quantify premature cleavage;
(c) nascent transcription (PRO-seq/TT-seq) to measure Pol II processivity rather than infer it;
(d) spike-in normalization to make absolute direction interpretable.
