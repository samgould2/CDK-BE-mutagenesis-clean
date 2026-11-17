"""
03_build_logo_plots.py

Builds amino-acid sequence logo plots from the codon pileup data
produced by 02_align_and_pileup.py, in the style used in Sam Gould's
own analysis notebook (notebooks/09_analysis_and_viz/logo_plots.ipynb
in samgould2/CDK-BE-mutagenesis):
  - contiguous stretch of residues (not just the named mutation sites)
  - y-axis on a 0-100 "% reads" scale
  - non-WT glyphs highlighted in red

Usage:
    python 03_build_logo_plots.py

Requires: logomaker, matplotlib, pandas
Output: ../output/CDK12_aa_logo_plots.png
"""
import os
import pickle
from collections import Counter

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import logomaker

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "..", "output")

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

# WT amino acid at each residue plotted (used to decide which glyphs to
# highlight red as "non-WT")
WT_SEQ = {731: "G", 732: "I", 733: "I",
          816: "M", 817: "D", 818: "H", 819: "D", 820: "L", 821: "M", 822: "G"}

SAMPLE_LABELS = {
    "X3295C_1_G731K_.fastq": "G731K gRNA",
    "X3295C_2_I733V.fastq": "I733V gRNA",
    "X3295C_4_NT_51.fastq": "NT control",
    "X3295C_3_G822E.fastq": "G822E/K gRNA",
    "X3295C_5_NT_94.fastq": "NT control",
}

# plotting order: (fastq_name, panel_label)
SAMPLE_ORDER = [
    "X3295C_4_NT_51.fastq",
    "X3295C_1_G731K_.fastq",
    "X3295C_2_I733V.fastq",
    "X3295C_5_NT_94.fastq",
    "X3295C_3_G822E.fastq",
]

MIN_GLYPH_HIGHLIGHT_PCT = 1.0  # only highlight alt glyphs above this % as red


def build_freq_table(counters):
    """counters: {residue: Counter({codon: count})} -> DataFrame[residue, amino_acid] in %"""
    rows = {}
    for res, c in counters.items():
        total = sum(c.values())
        aa_counts = Counter()
        for codon, cnt in c.items():
            aa_counts[CODON_TABLE.get(codon, 'X')] += cnt
        rows[res] = {aa: 100.0 * cnt / total for aa, cnt in aa_counts.items()}
    return pd.DataFrame(rows).T.fillna(0.0).sort_index()


def main():
    pileup_path = os.path.join(OUT_DIR, "codon_pileups.pkl")
    with open(pileup_path, "rb") as f:
        all_results = pickle.load(f)

    try:
        plt.rc('font', family='Helvetica')
    except Exception:
        pass  # falls back to matplotlib default if Helvetica isn't installed

    fig, axes = plt.subplots(len(SAMPLE_ORDER), 1, figsize=(6, 2.3 * len(SAMPLE_ORDER)))

    for ax, fq_name in zip(axes, SAMPLE_ORDER):
        df = build_freq_table(all_results[fq_name])
        logo = logomaker.Logo(df, vpad=0.02, ax=ax)
        logo.style_xticks(anchor=0, spacing=1, rotation=0, fontsize=11)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_yticklabels([0, 25, 50, 75, 100], fontsize=10)
        ax.set_ylabel("% reads", fontsize=11)
        ax.set_xlabel("CDK12 residue", fontsize=11)
        ax.set_title(SAMPLE_LABELS[fq_name], fontsize=12, loc='left')

        for res in df.index:
            wt_aa = WT_SEQ.get(res)
            for aa in df.columns:
                val = df.loc[res, aa]
                if aa != wt_aa and val >= MIN_GLYPH_HIGHLIGHT_PCT:
                    logo.style_single_glyph(c=aa, p=res, color='tab:red')

        logo.style_spines(spines=['right', 'top'], visible=False)

    plt.tight_layout()
    out_path = os.path.join(OUT_DIR, "CDK12_aa_logo_plots.png")
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
