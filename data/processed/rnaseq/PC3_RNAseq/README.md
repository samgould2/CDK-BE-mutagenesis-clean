# PC3 CDK12 RNA-Seq Analysis

**Project:** CDK12 kinase domain variant functional characterization in PC3 prostate cancer cells  
**Analyst:** Jonuelle Acosta, Sánchez-Rivera Lab, MIT  
**Pipeline version:** September 2026

---

## Experimental Design

### Genotypes
| Genotype | Description |
|----------|-------------|
| NT | Non-targeting control (wild-type CDK12) |
| I733V | CDK12 kinase domain point mutant |
| G822E | CDK12 kinase domain point mutant — **confirmed hypomorphic** (reduced (P)-Ser2, slower growth) |
| G731K | CDK12 kinase domain point mutant — **confirmed hypomorphic** (reduced (P)-Ser2, slower growth) |

### Treatments
| Treatment | Description |
|-----------|-------------|
| Control | DMSO vehicle |
| HQ461 | CDK12/13 inhibitor |
| ZLY-025 | CDK12/13 inhibitor |

3 biological replicates per genotype × treatment = **36 total samples**

### Sample naming convention
`{num}_{genotype}_{treatment}Rep{rep}` — e.g., `5_G822E_HQ461Rep2` or `X12_NT_ControlRep1` (leading `X` stripped automatically)

---

## Folder Structure

```
Comprehensive PC3 RNA-seq Analysis/
├── README.md
├── data/
│   ├── salmon.merged.gene_counts.tsv     ← nf-core/rnaseq Salmon output
│   └── deseq2.pca.vals.txt               ← PCA coordinates
├── scripts/
│   ├── pc3_deseq2_gsea.py                ← main DESeq2 + GSEA pipeline (all 17 comparisons)
│   ├── nt_resistance_score.py            ← Pearson r genome-wide drug-response correlation
│   ├── nt_gsea_preservation_score.py     ← GSEA Signature Preservation Score (SPS) bar charts
│   ├── nt_gsea_enrichment_plots.py       ← GSEA enrichment curves + NES fold-change bars
│   ├── nt_lfc_scatter.py                 ← LFC scatter: variant vs NT drug response
│   ├── nt_deg_counts.py                  ← total DEG count bar chart per genotype × drug
│   ├── regen_ddr_heatmaps.py             ← regenerate DDR/CDK12 ssGSEA heatmaps
│   ├── make_concordance_figure.py        ← LFC concordance across conditions
│   ├── make_hypomorph_program.py         ← hypomorphic program ORA figure
│   └── make_upset_figure.py              ← DEG overlap UpSet plot
└── pc3_rnaseq_results/
    ├── deseq2_results/
    │   ├── drug_response/                ← within-genotype drug-effect DESeq2 CSVs (8)
    │   ├── resistance/                   ← NT vs mutant under drug CSVs (6)
    │   ├── intrinsic/                    ← NT vs mutant no-drug CSVs (3)
    │   └── supplemental/                 ← elongation index, ssGSEA scores, shared hypomorphic DEGs
    └── figures/
        ├── 01_qc/                        ← PCA
        ├── 02_volcano/                   ← volcano plots (drug effects, resistance, intrinsic)
        ├── 03_hypomorphic/               ← hypomorphic characterization figures
        ├── 04_gsea_cdk12/                ← CDK12 drug-response GSEA, SPS, LFC scatter
        ├── 05_hallmark_gsea/             ← MSigDB Hallmark GSEA dot plots
        ├── 06_reactome/                  ← Reactome GSEA dot plots
        ├── 07_dna_damage_stress/         ← DDR/stress ssGSEA heatmap
        └── 08_prostate_lineage/          ← prostate lineage ssGSEA heatmap
```

---

## Scripts

### `scripts/pc3_deseq2_gsea.py` — Main pipeline

Runs all 17 DESeq2 comparisons, volcano plots, PCA, Hallmark GSEA, Reactome GSEA, ssGSEA heatmaps, hypomorphic characterization figures, and elongation index.

```bash
cd "Comprehensive PC3 RNA-seq Analysis"
python scripts/pc3_deseq2_gsea.py
```

Runtime ~15–25 min. All outputs land in `pc3_rnaseq_results/`.

---

### `scripts/nt_resistance_score.py` — Pearson r resistance score

Computes genome-wide Pearson correlation between each variant's drug-response LFC and NT's drug-response LFC. r ≈ 1 = transcriptionally sensitive; r ≈ 0 = transcriptionally resistant.

```bash
python scripts/nt_resistance_score.py
```

**Output:** `figures/04_gsea_cdk12/NT_resistance_score.{pdf,png,svg}`

---

### `scripts/nt_lfc_scatter.py` — Genome-wide LFC scatter

Six panels (3 variants × 2 drugs), 3×2 grid. X-axis = NT drug LFC; Y-axis = variant drug LFC. NT DEGs (padj < 0.05, |LFC| > 1) highlighted; recovered DEGs (same-direction DE in variant) filled, lost DEGs open circles. Annotated with Pearson r, p-value, and DEG recovery %.

```bash
python scripts/nt_lfc_scatter.py
```

**Output:** `figures/04_gsea_cdk12/NT_LFC_scatter.{pdf,png,svg}`

---

### `scripts/nt_deg_counts.py` — DEG count bar chart

Total DEG counts (padj < 0.05, |LFC| > 1) per genotype for each drug. Two panels side by side (HQ461, ZLY-025). Bars ordered NT → G731K → I733V → G822E, making the resistance gradient immediately visible.

```bash
python scripts/nt_deg_counts.py
```

**Output:** `figures/04_gsea_cdk12/NT_DEG_counts.{pdf,png,svg}`

---

### `scripts/nt_gsea_preservation_score.py` — Signature Preservation Score

Runs preranked GSEA of the **core NT CDK inhibition signature** (intersection of HQ461 and ZLY-025 NT DEGs: ~1806 UP, ~2775 DN) against each variant's drug-response ranked list. Computes the Signature Preservation Score:

> **SPS = mean(NES_UP / NT_NES_UP, |NES_DN| / NT_NES_DN)** — normalized so NT self-score = 1.0

Only components with NOM p < 0.05 contribute to SPS; asymmetric cases (one component n.s.) are flagged with `†`. Bars annotated with FDR stars (*** < 0.001, ** < 0.01, * < 0.05). Sorted ascending by SPS (most resistant left).

```bash
python scripts/nt_gsea_preservation_score.py
```

**Output:** `figures/04_gsea_cdk12/GSEA_preservation_score.{pdf,png,svg}`  
**Raw scores:** `figures/04_gsea_cdk12/GSEA_preservation_score_raw.csv`

---

### `scripts/nt_gsea_enrichment_plots.py` — GSEA enrichment curves

Shows individual enrichment curves for the core NT UP signature in each variant condition, with NT reference curve overlaid in grey. Annotated with NES, ΔNES (variant − NT), and FDR. A separate NES fold-change bar plot summarizes enrichment attenuation across all conditions.

```bash
python scripts/nt_gsea_enrichment_plots.py
```

**Outputs:**
- `figures/04_gsea_cdk12/NT_GSEA_enrichment.{pdf,png,svg}` — 3×2 enrichment curve panels
- `figures/04_gsea_cdk12/NT_GSEA_NES_foldchange.{pdf,png,svg}` — NES fold-change bar plot

---

### `scripts/regen_ddr_heatmaps.py` — DDR heatmaps

Regenerates the DDR and CDK12/13 ssGSEA heatmaps split by drug (HQ461, ZLY-025) and variant effect, without re-running the full pipeline.

```bash
python scripts/regen_ddr_heatmaps.py
```

**Outputs:** `figures/04_gsea_cdk12/GSEA_DDR_heatmap_drug_HQ461.{pdf,png}`, `_drug_ZLY025`, `_variant_effect`

---

## Comparisons (17 total)

### Drug response — within each genotype (8, in `drug_response/`)
| Key | Interpretation |
|-----|----------------|
| `NT_HQ461_vs_NT_Control` | NT: HQ461 transcriptional effect |
| `NT_ZLY025_vs_NT_Control` | NT: ZLY-025 transcriptional effect |
| `I733V_HQ461_vs_I733V_Control` | I733V: HQ461 effect |
| `I733V_ZLY025_vs_I733V_Control` | I733V: ZLY-025 effect |
| `G822E_HQ461_vs_G822E_Control` | G822E: HQ461 effect |
| `G822E_ZLY025_vs_G822E_Control` | G822E: ZLY-025 effect |
| `G731K_HQ461_vs_G731K_Control` | G731K: HQ461 effect |
| `G731K_ZLY025_vs_G731K_Control` | G731K: ZLY-025 effect |

### Drug resistance — NT vs mutant under drug (6, in `resistance/`)
NT is the numerator. Positive LFC = higher in NT (sensitive); negative LFC = higher in mutant (resistant).

| Key | Interpretation |
|-----|----------------|
| `NT_vs_I733V_HQ461` | HQ461: NT vs I733V resistance |
| `NT_vs_G822E_HQ461` | HQ461: NT vs G822E resistance |
| `NT_vs_G731K_HQ461` | HQ461: NT vs G731K resistance |
| `NT_vs_I733V_ZLY025` | ZLY-025: NT vs I733V resistance |
| `NT_vs_G822E_ZLY025` | ZLY-025: NT vs G822E resistance |
| `NT_vs_G731K_ZLY025` | ZLY-025: NT vs G731K resistance |

### Intrinsic differences — no drug (3, in `intrinsic/`)
NT is the numerator. Positive LFC = higher in NT; negative LFC = higher in mutant at baseline.

| Key | Interpretation |
|-----|----------------|
| `NT_vs_I733V_Control` | Baseline: NT vs I733V |
| `NT_vs_G822E_Control` | Baseline: NT vs G822E |
| `NT_vs_G731K_Control` | Baseline: NT vs G731K |

---

## Output Figures

### `01_qc/PCA.{pdf,png}`
PCA of all 36 samples (VST-normalized). PC1 vs PC2, colored by condition, shaped by genotype.

### `02_volcano/volcano_drug_effects.{pdf,png}`
2×4 volcano grid (rows = HQ461, ZLY-025; cols = NT, I733V, G822E, G731K). Drug effect within each genotype. Coral = up, steel blue = down (padj < 0.05, |LFC| > 0.5).

### `02_volcano/volcano_resistance.{pdf,png}`
2×3 volcano grid (rows = HQ461, ZLY-025; cols = I733V, G822E, G731K). NT vs mutant under drug. Genes with negative LFC are relatively higher in the mutant (drug-resistant direction).

### `02_volcano/volcano_nodrug.{pdf,png}`
1×3 volcano grid. Baseline DEGs: NT vs each mutant under vehicle.

### `03_hypomorphic/` — Hypomorphic characterization

| Figure | Description |
|--------|-------------|
| `hypomorphic_ssgsea` | ssGSEA scores for curated DDR/HR gene sets across all 36 samples |
| `hypomorphic_empirical` | ssGSEA using empirical drug-response signatures on control samples ⚠️ |
| `hypomorphic_attenuation` | Drug-response transcriptional attenuation per variant |
| `hypomorphic_lfc_scatter` | NT LFC vs variant LFC scatter (2×3); includes Pearson r and OLS slope |
| `elongation_index` | Mean(long DDR genes) / mean(short HK genes) ratio by condition |
| `lfc_concordance` | LFC correlation across all pairwise comparisons |
| `lfc_concordance_de_upset` | UpSet of DEG overlaps with concordance |
| `de_gene_upset` | UpSet plot of DEG overlap across genotypes |
| `hypomorph_program_ora` | ORA of shared hypomorphic DEG program |

> ⚠️ `hypomorphic_empirical`: Acute drug-response signatures are inappropriate for detecting chronic hypomorphism. Chronic CDK12 hypomorphism upregulates DDR compensatory programs; acute inhibition suppresses them. Retained for reference only.

### `04_gsea_cdk12/` — CDK12 drug-response analysis

| Figure | Description |
|--------|-------------|
| `NT_DEG_counts` | Total DEG counts per genotype × drug; bars show resistance gradient at a glance |
| `NT_resistance_score` | Pearson r between variant and NT drug-response LFC profiles; r=1 = sensitive, r≈0 = resistant |
| `NT_LFC_scatter` | Genome-wide LFC scatter (3×2); NT DEGs highlighted by recovery status |
| `GSEA_preservation_score` | SPS bar chart; normalized to NT=1.0; † = asymmetric (one direction n.s.) |
| `NT_GSEA_enrichment` | Enrichment curves for core NT UP signature in each variant, with NT grey reference |
| `NT_GSEA_NES_foldchange` | NES fold-change (NES_variant / NES_NT) bars; NT anchor = 1.0 |
| `GSEA_CDK12_13_heatmap` | ssGSEA NES heatmap for CDK12/13-relevant gene sets |
| `GSEA_DDR_heatmap_drug_HQ461` | DDR ssGSEA heatmap, HQ461-treated conditions |
| `GSEA_DDR_heatmap_drug_ZLY025` | DDR ssGSEA heatmap, ZLY-025-treated conditions |
| `GSEA_DDR_heatmap_variant_effect` | DDR ssGSEA heatmap, variant vs NT contrasts |

### `05_hallmark_gsea/` — MSigDB Hallmark GSEA
Preranked GSEA (Wald stat) dot plots for NT drug response and baseline mutant vs NT contrasts. Dot size = −log₁₀(FDR); color = NES direction.

### `06_reactome/` — Reactome GSEA
Reactome preranked GSEA dot plots for baseline NT vs mutant comparisons (G731K, G822E, I733V).

### `07_dna_damage_stress/ddr_stress_ssgsea_heatmap.{pdf,png}`
ssGSEA NES heatmap for curated DDR and stress-response gene sets across all 12 conditions.

### `08_prostate_lineage/prostate_lineage_ssgsea_heatmap.{pdf,png}`
ssGSEA NES heatmap for prostate lineage gene sets across all 12 conditions.

---

## Core NT CDK Inhibition Signature

The shared CDK12 drug-response signature used for GSEA preservation analysis is defined as the **intersection** of HQ461 and ZLY-025 NT DEGs at padj < 0.05 and |LFC| > 1:

- **Core UP:** ~1,806 genes upregulated by both drugs in NT
- **Core DN:** ~2,775 genes downregulated by both drugs in NT

Using the intersection (rather than drug-specific signatures) ensures any enrichment in variant conditions reflects the canonical CDK12 transcriptional program and not drug-idiosyncratic effects.

---

## Key Biological Findings

### 1. Variant resistance phenotype (HQ461 and ZLY-025)

Ranked by transcriptional resistance (most resistant → most sensitive):

| Variant | Drug | Phenotype |
|---------|------|-----------|
| G731K | HQ461 | Near-complete resistance — lowest Pearson r and SPS |
| G731K | ZLY-025 | Near-complete resistance |
| G822E | HQ461 | Partial resistance |
| G822E | ZLY-025 | Asymmetric: UP program partially preserved, DN program absent (SPS †) |
| I733V | HQ461 | Partial resistance |
| I733V | ZLY-025 | Least resistant of the three variants |

G822E ZLY-025 DN program: NES ≈ −1.0 (gseapy sentinel value), NOM p = 1.0 — genuine biological absence of the downregulation program, not a numerical artifact.

### 2. Hypomorphic mutants show a chronic replication stress signature at baseline

G822E and G731K vs NT (no drug):
- **Upregulated in mutants:** E2F Targets, G2M Checkpoint, DNA Repair (FDR < 0.01) — compensatory ATR/CHK1 activation from chronic replication stress
- **Downregulated in mutants:** Hypoxia, Androgen Response, Inflammatory Response — consistent with CDK12's elongation function at long, pause-prone transcription units

### 3. G731K shows near-complete transcriptional drug resistance

Elongation index drops ~20% in NT under both drugs. G731K maintains near-baseline elongation and has the lowest Pearson r and SPS across both drugs.

### 4. G822E has asymmetric program preservation under ZLY-025

The UP transcriptional program is partially preserved in G822E under ZLY-025, but the DN program is completely absent (SPS flagged with †). This asymmetry may reflect mechanistic differences between HQ461 and ZLY-025.

### 5. I733V clusters with confirmed hypomorphs at baseline

Despite unconfirmed hypomorphism by western blot, I733V Hallmark GSEA clusters with G822E and G731K, suggesting mild CDK12 activity reduction below western blot sensitivity or altered substrate selectivity with convergent transcriptional consequences.

---

## Statistical Methods

| Step | Method |
|------|--------|
| Normalization | DESeq2 median-of-ratios (pydeseq2) |
| Differential expression | pydeseq2 Wald test + LFC shrinkage (apeglm) |
| DEG threshold | padj < 0.05, \|LFC\| > 0.5 (volcano); padj < 0.05, \|LFC\| > 1 (core signature) |
| padj floor | 1 × 10⁻³⁰⁰ (prevents Inf on −log₁₀ axis) |
| Pre-filtering | ≥ 10 counts in ≥ 3 samples |
| GSEA | gseapy prerank, Wald stat ranking, 1000 permutations, seed=42, max_size=5000 |
| Gene sets | MSigDB Hallmark v2023.2 (50 sets); Reactome; curated DDR/CDK12 sets |
| Core signature | NT HQ461 ∩ NT ZLY-025 DEGs (padj < 0.05, \|LFC\| > 1) |
| SPS | mean(NES_variant / NES_NT) for NOM p < 0.05 components; normalized to NT = 1.0 |
| Pearson r | Genome-wide drug-response LFC correlation (all shared genes with complete data) |
| PCA | sklearn on VST-normalized log counts |

---

## Figure Aesthetics

- Font: Liberation Sans (PDF-embedded as Type 42 — text editable in Illustrator)
- All figures: `.pdf` (vector, Illustrator-ready) + `.png` (300 dpi) + `.svg` where applicable

---

## Dependencies

```
pydeseq2    pandas    numpy    scipy    sklearn
matplotlib  gseapy    upsetplot
```

---

## Contact

Jonuelle Acosta — jonuelle@mit.edu  
Sánchez-Rivera Lab, MIT
