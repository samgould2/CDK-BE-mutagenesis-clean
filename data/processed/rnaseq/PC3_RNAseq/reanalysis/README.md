# PC3 CDK12 RNA-seq — re-analysis (2026-09-30)

Independent re-analysis of the 36-sample PC3 dataset described in `../README.md`.
Rebuilt from the raw Salmon count matrix using **R/DESeq2** rather than pydeseq2,
with a different set of summary statistics and a fresh figure set.

## Design

4 genotypes (NT, I733V, G822E, G731K) × 3 treatments (DMSO, HQ461, ZLY-025) × 3 replicates.
I733V, G822E and G731K are CDK12 kinase-domain variants recovered as resistance hits in the
base-editing screens (`notebooks/11_ChimeraX/cdk12_13_ddb1_hit_residues.csv`);
G822E is HQ461-specific in the screen data, G731K and I733V score as general resistance.

## Pipeline

```bash
cd reanalysis/scripts
Rscript 01_deseq2.R      # DESeq2: contrasts + genotype×treatment interaction    ~6 min
python  02_gsea.py       # preranked GSEA, all 50 Hallmark sets + curated DDR    ~5 min
python  03_figures.py    # main figures (PDF + 300-dpi PNG) + figure source data ~2 min
python  04_elongation.py # targeted GSEA + length & intron elongation analysis     ~3 min
python  05_elongation_length.py # length-only + baseline elongation figures        ~1 min
python  06_elongation_combined.py # combined drug + baseline, per-replicate error   ~1 min
python  07_elongation_shared_y.py # shared-y variants, with and without inset        ~1 min
python  08_elongation_controls.py # expression/intron confounder controls             ~1 min
python  09_length_heatmap_scatter.py # decile heatmap + per-gene density scatter       ~1 min
python  10_mcl1_panel.py # MCL1 per-condition panel, laid out to pair with a Western     ~1 min
```

`_style.py` holds the shared palette and figure defaults; edit it to restyle everything at once.

### What differs from the original analysis

| | Original | Here |
|---|---|---|
| DE engine | pydeseq2 | R/DESeq2 4.4 |
| Shrinkage | apeglm | ashr (apeglm not installed); unshrunk Wald stat kept for GSEA ranking |
| Contrasts | 17 | 18 pairwise + 6 interaction coefficients + a global LRT |
| Interaction model | none | `~ genotype * treatment` — directly tests whether the drug response differs by genotype |
| Hallmark sets | 16 hand-picked | all 50 (MSigDB Hallmark 2020 via Enrichr) |
| Resistance metric | Pearson *r*, GSEA "preservation score" | **response slope** (regression through origin over NT DEGs, bootstrap 95% CI) and **% signature shift retained** (per-replicate, so it carries error bars) |
| Gene length | 24-gene curated "elongation index" | genome-wide: log₂FC by gene-length **and intron-number** decile over 8,547 expressed protein-coding genes (Ensembl 113 canonical transcripts) |
| GSEA framing | 50-set Hallmark sweep as the main panel | targeted sets chosen from the expected mechanism (long/short, intron-rich/-poor, HR-DDR, cGAS-STING, immediate-early); the Hallmark sweep moves to supplemental |

Why the slope instead of Pearson *r*: *r* measures whether the variant moves the same genes
in the same rank order, and stays moderately high even when the response is strongly damped.
The slope measures *how much* of the NT response magnitude is reproduced, which is the actual
resistance quantity. Reporting both is useful — I733V under HQ461 has r = 0.45 but slope = 0.13,
i.e. it responds to the right genes, barely.

## Results

### Baseline (DMSO, variant vs NT)

- **G822E and G731K are near-identical at baseline** (LFC correlation *r* = 0.90 between the two
  variant-vs-NT contrasts) — 2,502/2,349 and 2,228/2,126 DEGs respectively. **I733V is much
  milder** (479/484 DEGs; *r* ≈ 0.3 with the other two). It shares the *direction* of the
  Hallmark program but at roughly a fifth of the effect size — and on the gene-length axis
  (Fig9) it shows no baseline defect at all, so it should not be called a hypomorph.
- All three are **up** for E2F Targets, G2-M Checkpoint, MYC Targets, DNA Repair and
  **down** for Hypoxia, Glycolysis, KRAS Signaling, Androgen Response, EMT, TNFα/NF-κB
  (FigS_hallmark_GSEA_baseline). Consistent with chronic replication stress plus loss of
  CDK12-dependent elongation at long, pause-prone genes.
- Interferon Alpha Response is up in G822E and G731K only — consistent with the cGAS-STING
  arm of CDK12 loss, and absent in the mild variant.

### Drug response

Fraction of the NT drug response retained (Fig3, Fig4; full numbers in
`figures/Fig4_response_slopes.csv` and `figures/Fig3_pct_response_per_replicate.csv`):

| Variant | HQ461 slope | HQ461 % retained | ZLY-025 slope | ZLY-025 % retained |
|---|---|---|---|---|
| I733V | 0.13 | 15% | 0.48 | 61% |
| G822E | 0.20 | 25% | 0.34 | 49% |
| G731K | 0.04 | 10% | 0.03 | 5% |

- **G731K is resistant to both drugs** and is the only variant that is. Slope ≈ 0.03–0.04 and
  44–47 up / 115–134 down DEGs vs ~2,600–3,100 / ~3,800–4,200 in NT (Fig2). Effectively a
  transcriptionally null response.
- **I733V and G822E are strongly resistant to HQ461 but substantially sensitive to ZLY-025**
  (~50–60% of the NT response retained). This drug-selective pattern is the most interesting
  result in the dataset and matches the screen classification of G822E as an HQ461-specific hit.
  It argues the two compounds engage CDK12 differently — expected if HQ461 acts as a molecular
  glue requiring a specific CDK12 surface, which these substitutions disrupt, whereas ZLY-025
  does not depend on that surface.
- The interaction model agrees: 6,000–9,700 genes per variant×drug have a significantly
  different drug response than NT (padj < 0.05), and 12,334 genes are significant in the global
  LRT (`results/interaction_summary.csv`).

### The elongation defect (Fig7) — the cleanest mechanistic panel

The expected consequence of CDK12/13 inhibition is loss of processive Pol II elongation, which
causes premature cleavage and polyadenylation at intronic poly(A) sites. Transcripts needing the
most elongation — long genes, intron-rich genes — should be lost preferentially, while short
intron-poor transcripts rise. That is exactly what NT does, on both axes:

- **Gene length.** In NT the response is monotonically length-dependent: median log₂FC ≈ +0.95 in
  the shortest decile under HQ461, falling to ≈ −0.6 in the longest. **All three variants flatten
  this to near zero under HQ461.** Under ZLY-025 I733V and G822E retain a partial gradient
  (long-gene gap −0.23 and −0.17 vs NT's −0.37) while G731K stays flat (−0.07).
- **Intron number.** The same pattern holds on an axis that is only partly collinear with length
  (log-scale *r* = 0.57): NT's intron-rich gap is −0.34 under HQ461, abolished in I733V (−0.04)
  and strongly reduced in G822E and G731K.

### Baseline elongation defect (Fig9) — direct evidence of hypomorphism

The same gene-length axis applied to the **no-drug** contrasts (variant vs NT in DMSO) separates
the variants on mechanism rather than on drug response:

| Variant | Baseline long-gene gap (95% CI) |
|---|---|
| I733V | −0.004 (−0.010 to +0.002) — **no defect** |
| G822E | −0.116 (−0.146 to −0.085) |
| G731K | −0.082 (−0.106 to −0.061) |

G822E and G731K show a monotonic decline across length deciles without any drug — median log₂FC
of roughly +0.10 for the shortest genes falling to −0.03 to −0.07 for the longest. That is a
chronic, partial elongation defect, i.e. hypomorphism measured directly on CDK12's own mechanistic
axis, and it agrees with the western-blot calls in the original analysis. The effect is ~4×
smaller than acute HQ461 in NT (−0.41), which is the expected size for a viable hypomorph.

**I733V is flat at baseline** (CI spans zero). Combined with its mild baseline DEG count (479/484
vs ~2,400 for the other two), this revises the original README's claim that I733V "clusters with
the confirmed hypomorphs": on the elongation axis it does not. Its resistance is therefore better
explained by altered drug engagement than by reduced basal CDK12 activity — which also fits its
being the *least* resistant variant to ZLY-025 (61% response retained) while strongly resistant
to HQ461 (15%).

Summary metric is the **top-decile response gap** (median log₂FC of the top decile minus the
median of all genes, bootstrap 95% CI). It is used in preference to an OLS slope because NT's
ZLY-025 intron profile is non-monotonic — single-intron genes behave atypically — where a linear
slope is a poor summary. Both metrics are in `results/elongation_slopes.csv`.

This reproduces CDK12's known elongation phenotype genome-wide rather than through a 24-gene
index, and discriminates the variants on the same axis as the resistance metrics.

Full methodology, controls and limitations for this analysis: **[METHODS_elongation.md](METHODS_elongation.md)**.
In brief, the length effect survives stratification by expression level (NT HQ461: −0.57 to −0.86
within expression tertiles vs −1.02 globally) and retains ~75% of its regression coefficient when
expression is controlled, but it is stronger in highly expressed genes and part of the global
effect size reflects the length–expression correlation (*r* = −0.26). Intron number adds little
independent signal once length is in the model.

### Targeted GSEA (Fig6)

Replacing the 50-set Hallmark sweep with sets chosen from the expected mechanism makes the result
legible in one panel: **LONG_GENES and INTRON_RICH are strongly negative, SHORT_GENES and
INTRON_POOR strongly positive in NT under both drugs**, with the enrichment collapsing in the
variants in the same rank order as every other resistance metric. Immediate-early genes (short,
intron-poor, stress-induced) are the most robustly up class. The full 50-set sweep is retained as
`FigS_hallmark_GSEA_drug_response`.

### MCL1 separates the two drugs (Fig12)

MCL1 is a canonical short-lived CDK9/CDK12-13 target, and at 18.8 kb it sits mid-distribution
on gene length (decile 4 of 10) — so its behaviour here is not a length effect.

| Genotype | HQ461 | ZLY-025 |
|---|---|---|
| NT | **+0.38** (padj 1×10⁻³) — *up* | **−1.26** (padj 2×10⁻²⁷) — 2.4-fold down |
| I733V | −0.06 (n.s.) | −0.24 (n.s.) |
| G822E | −0.23 (n.s.) | −0.35 (padj 6×10⁻³) |
| G731K | −0.04 (n.s.) | +0.02 (n.s.) |

The two drugs move MCL1 in **opposite directions in the same cells**, and all three variants
block the ZLY-025 effect. MCL1 is unchanged at baseline across genotypes (all n.s.), so this is
a drug response rather than a standing difference. Independent support for the HQ461-vs-ZLY-025
mechanistic difference, in a gene with direct functional consequence — MCL1 loss is the
apoptotic arm of CDK9/CDK12-13 inhibition, and its ~30 min protein half-life means mRNA loss
should track to protein quickly.

### DDR genes (Fig8)

Canonical CDK12-dependent DDR genes (*BRCA1/2*, *ATR*, *TOPBP1*, *PALB2*, *RAD51C/D*, *ERCC4*,
*NBN*, *CHEK1*) are suppressed in NT by both drugs and protected in G731K under both.
*FANCA*, *FANCI*, *FANCD2* and *ORC3* move **up**, so "DDR down" is not uniform — worth flagging
before that claim goes in a figure legend. At baseline G822E and G731K have modestly elevated
*CHEK1*, *RAD51D* and *ORC3*, the compensatory direction.

## Caveats

1. **`NT_ZLY025_vs_NT_HQ461` shows the two drugs are not equivalent in NT** (1,468 up / 2,072
   down at padj < 0.05, 1,008 genes with LFC < −1). NT LFC profiles correlate *r* = 0.86 between
   drugs, so they overlap heavily but ZLY-025 has real additional effects. Since dose was not
   normalized to equipotency, part of the I733V/G822E drug-selectivity could be a potency
   difference rather than a mechanistic one. A dose-response or a washout experiment would
   separate these. I'd treat the HQ461-vs-ZLY-025 contrast as the key follow-up.
2. **Library sizes vary ~2.5×** (4.2M–10.7M; the three G822E ZLY-025 replicates are the
   smallest). DESeq2 normalization handles this, but G822E ZLY-025 has the least power of any
   group.
3. **`I733V_ZLY025` replicates are less reproducible than other groups**
   (FigS_sample_correlation) — one replicate sits apart on PC1.
4. Gene "length" here is genomic span of the Ensembl 113 gene locus, not a transcript-weighted
   length. Fine for the decile trend; do not quote absolute values.
5. The original README's claim that a G822E/ZLY-025 GSEA NES of −1.0 with p = 1.0 is a
   "genuine biological absence" is not supported — that is gseapy's sentinel for a test that
   did not produce a usable enrichment, and it should not be interpreted. The slope and
   %-retained metrics here show G822E does respond to ZLY-025 (49% retained), which is the
   opposite conclusion.

## Figures

| File | Content |
|---|---|
| `Fig1_PCA` | PCA, all samples + DMSO only |
| `Fig2_DEG_counts` | DEGs per genotype × drug, up/down |
| `Fig2b_DEG_counts_with_baseline` | As Fig2 plus a third panel: variant vs NT in DMSO (shared y-axis) |
| `Fig3_response_retained` | % of NT response retained, per replicate (mean ± SD) |
| `Fig4_LFC_scatter_vs_NT` | Variant vs NT drug log₂FC, slope through origin + bootstrap CI |
| `Fig5_signature_heatmap` | Core NT signature genes, each genotype centred on its own DMSO |
| `Fig6_targeted_GSEA` | Targeted GSEA: elongation-dependent + CDK12-pathway sets |
| `Fig7_elongation_length` | log₂FC by gene-length decile + long-gene response gap (drug response) |
| `Fig8_DDR_genes` | CDK12-dependent DDR gene panel |
| `Fig9_baseline_elongation` | Gene-length profile with no drug (variant vs NT) — baseline hypomorphism |
| `Fig10_elongation_combined` | Fig7 + Fig9 in one 3-panel figure, replicate-level error bars, no bar panels |
| `Fig10a_elongation_shared_y` | As Fig10 but one y-scale throughout — baseline shown at true relative size |
| `Fig10b_elongation_shared_y_inset` | As Fig10a plus an inset zooming the variant-vs-NT baseline curves |
| `Fig12_MCL1_expression` | MCL1 normalized counts across all 12 conditions, blot-lane order |
| `Fig11_length_decile_heatmap` | Length decile × condition, colour = median log₂FC (all 11 contrasts) |
| `FigS_elongation_length_and_introns` | Both axes together: gene length and intron number |
| `FigS_length_scatter` | Per-gene density (hexbin) of log₂FC vs gene length + decile median |
| `FigS_hallmark_GSEA_drug_response` | Full 50-set Hallmark GSEA, drug vs DMSO within genotype |
| `FigS_hallmark_GSEA_baseline` | Hallmark GSEA, variant vs NT in DMSO |
| `FigS_volcano_drug_response`, `FigS_volcano_baseline` | Volcano grids |
| `FigS_sample_correlation` | Sample–sample Pearson correlation (QC) |

**Palette.** NT is neutral dark grey (`#3D3D3D`) as the reference genotype; variants are blue
`#1B6CA8` (I733V), gold `#E1A730` (G822E), plum `#8E4585` (G731K). The three variant hues were
checked for colourblind separation (deuteran/protan/tritan ΔE above the usable floor on all
pairs). Gold falls below 3:1 contrast on white, which the black bar outlines compensate for.

Each figure ships `.pdf` (vector, Type-42 text, Illustrator-editable, **transparent background**)
+ `.png` (300 dpi, white background for quick viewing).

In `Fig3_response_retained` the y-axis ticks stop at 100% because that is the NT anchor, not a
measured quantity; the axis limit runs slightly past it only so NT's own replicate scatter
(97.4–102.2% of its mean — the ~±2% noise floor of the metric) is shown rather than clipped.
Panels with derived statistics also write their source data as `.csv` next to the figure.

**Error bars in `Fig10`.** Mean ± SD across the three biological replicates, not the
gene-to-gene IQR shown in Fig7/Fig9 — the IQR is biological spread between genes, not
uncertainty in the decile estimate. Replicate-level log₂ ratios are built from the VST matrix
(each replicate against the mean of its reference group), because a DESeq2 shrunken LFC is a
single contrast-level estimate with no per-replicate decomposition. The two agree in shape and
sign but VST compresses magnitudes (it shrinks low-count genes toward the mean), so Fig10's
values run ~25% smaller than Fig7/Fig9 — read within the panel, not against those numbers.
The baseline panel carries its own expanded y-axis (the effect is ~4x smaller than acute drug);
compare its shape to the drug panels, not its amplitude.
Replicate SDs of the top-decile median are 0.01–0.08 log₂, i.e. the decile profiles are far more
reproducible than the gene-to-gene IQR implied. Per-replicate values are in
`figures/Fig10_replicate_decile_profiles.csv`.

Core NT signature used in Fig3/Fig5 = intersection of HQ461 and ZLY-025 NT DEGs
(padj < 0.05, |log₂FC| > 1): **2,098 up / 3,068 down**.

## Outputs

```
reanalysis/
├── results/
│   ├── contrasts/           18 pairwise DESeq2 contrasts (CSV)
│   ├── interaction/         6 genotype×treatment coefficients + global LRT
│   ├── gsea/                gsea_all_contrasts.csv, gsea_targeted.csv, *.gmt
│   ├── annotation/          Ensembl 113 gene spans + exon/intron structure
│   ├── elongation_slopes.csv, elongation_length_gaps.csv
│   ├── vst.csv, normalized_counts.csv, size_factors.csv, sample_metadata.csv
│   └── DEG_summary.csv, interaction_summary.csv
└── figures/                 PDF + PNG + per-figure source CSVs
```

`01_deseq2.R` also writes `dds_condition.rds` / `dds_interaction.rds` (~44 MB); these are
deleted after the run to keep the directory small and are regenerated by re-running the script.

Contrast CSV columns: `gene_id, gene_name, baseMean, log2FC` (unshrunk), `lfcSE`,
`log2FC_shrunk` (ashr — use for plotting/ranking by effect size), `stat` (Wald, used for GSEA
ranking), `pvalue`, `padj`.
