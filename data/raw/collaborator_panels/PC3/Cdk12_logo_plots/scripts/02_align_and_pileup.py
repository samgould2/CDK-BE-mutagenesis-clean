"""
02_align_and_pileup.py

Aligns each Plasmidsaurus Nanopore genotyping FASTQ to its correct
amplicon reference (minimap2 / mappy, map-ont preset), then pulls the
base calls at each target codon from every read individually (a true
per-read pileup, not relying on Plasmidsaurus's own allele-deconvolution
step, which can be unreliable for heterogeneous/pooled base-edited
populations).

Key correctness details (found the hard way -- see README):
  - ~40-45% of Nanopore reads map to the reverse strand. Reads must be
    reverse-complemented (sequence AND quality string) before walking
    the CIGAR, or you silently corrupt ~half the data.
  - Only positions with Phred quality >= MIN_QUAL are counted; a codon
    is only counted if ALL 3 of its bases pass this threshold.
  - Reads with a deletion inside the codon window are dropped for that
    residue (ambiguous/uncallable).

Usage:
    python 02_align_and_pileup.py

Requires: mappy
Output: writes codon_pileups.pkl (raw codon counts) and
        aa_frequency_table.csv (translated, %-scale) to ../output/
"""
import os
import pickle
from collections import Counter

import mappy as mp

from cds_mapping import codon_local_range

HERE = os.path.dirname(os.path.abspath(__file__))
REF_DIR = os.path.join(HERE, "..", "data", "reference")
FASTQ_DIR = os.path.join(HERE, "..", "data", "fastq")
OUT_DIR = os.path.join(HERE, "..", "output")
os.makedirs(OUT_DIR, exist_ok=True)

MIN_QUAL = 20  # Phred quality threshold per base

# residue -> which residues to pileup for each sample's amplicon.
# amp1 = CDK12_3351_2232 (703-750 region); amp2 = CDK12_3394_3395_3396 (807-870 region)
RESIDUES = {
    "amp1": list(range(731, 734)),   # 731, 732, 733
    "amp2": list(range(816, 823)),   # 816..822
}

# fastq filename -> (amplicon reference, human-readable sample label)
SAMPLES = {
    "X3295C_1_G731K_.fastq": ("amp1", "G731K gRNA"),
    "X3295C_2_I733V.fastq":  ("amp1", "I733V gRNA"),
    "X3295C_4_NT_51.fastq":  ("amp1", "NT control (amp1)"),
    "X3295C_3_G822E.fastq":  ("amp2", "G822E/K gRNA"),
    "X3295C_5_NT_94.fastq":  ("amp2", "NT control (amp2)"),
}

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

_COMP = str.maketrans("ACGTN", "TGCAN")
def revcomp(s):
    return s.translate(_COMP)[::-1]


def get_query_codon(hit, seq, qual, ref_start, ref_end, min_q=MIN_QUAL):
    """
    Walk a mappy alignment's CIGAR and return the query bases aligned to
    ref[ref_start:ref_end) (0-based half-open), or None if any base in
    that window fails the quality threshold, has a deletion, or the
    alignment doesn't fully span the window.

    Correctly handles reverse-strand hits by reverse-complementing the
    read (and reversing its quality string) before walking the CIGAR,
    since mappy's cigar/query-position bookkeeping for strand=-1 hits is
    relative to the reverse complement, not the raw read as given.
    """
    if hit.strand == -1:
        seq = revcomp(seq)
        qual = qual[::-1]
        q_st = len(seq) - hit.q_en
    else:
        q_st = hit.q_st

    r_pos, q_pos = hit.r_st, q_st
    mapping = {}
    for length, op in hit.cigar:
        if op in (0, 7, 8):  # match/mismatch (M, =, X)
            for i in range(length):
                rp = r_pos + i
                if ref_start <= rp < ref_end:
                    qbase = seq[q_pos + i]
                    qscore = ord(qual[q_pos + i]) - 33
                    mapping[rp] = qbase if qscore >= min_q else None
            r_pos += length
            q_pos += length
        elif op == 1:  # insertion in query
            q_pos += length
        elif op == 2:  # deletion in query
            for i in range(length):
                rp = r_pos + i
                if ref_start <= rp < ref_end:
                    mapping[rp] = '-'
            r_pos += length
        elif op == 3:  # N (skip)
            r_pos += length
        elif op == 4:  # soft clip
            q_pos += length
        # op 5 (hard clip) consumes neither ref nor query bases already in `seq`

    bases = [mapping.get(p) for p in range(ref_start, ref_end)]
    if any(b is None for b in bases):
        return None
    codon = ''.join(bases)
    return None if '-' in codon else codon


def pileup_sample(fastq_path, ref_fasta_path, residue_local_ranges):
    """Returns {residue_num: Counter({codon: count})} for one FASTQ file."""
    aligner = mp.Aligner(ref_fasta_path, preset="map-ont")
    counters = {res: Counter() for res in residue_local_ranges}
    for name, seq, qual in mp.fastx_read(fastq_path):
        hits = [h for h in aligner.map(seq) if h.is_primary]
        if not hits:
            continue
        hit = hits[0]
        for res, (s, e) in residue_local_ranges.items():
            if hit.r_st <= s and hit.r_en >= e:
                codon = get_query_codon(hit, seq, qual, s, e)
                if codon:
                    counters[res][codon] += 1
    return counters


def main():
    all_results = {}
    for fq_name, (amplicon, label) in SAMPLES.items():
        fastq_path = os.path.join(FASTQ_DIR, fq_name)
        ref_fasta_path = os.path.join(REF_DIR, f"ref_{amplicon}.fasta")
        residue_local_ranges = {
            r: codon_local_range(r, amplicon) for r in RESIDUES[amplicon]
        }
        print(f"Aligning {fq_name} -> {ref_fasta_path} ...")
        counters = pileup_sample(fastq_path, ref_fasta_path, residue_local_ranges)
        all_results[fq_name] = counters
        for res, c in counters.items():
            total = sum(c.values())
            top = c.most_common(3)
            print(f"  residue {res}: n={total}  top={top}")

    pickle_path = os.path.join(OUT_DIR, "codon_pileups.pkl")
    with open(pickle_path, "wb") as f:
        pickle.dump(all_results, f)
    print(f"\nSaved raw codon pileups to {pickle_path}")

    # also write a flat, translated, %-scale CSV for convenience/inspection
    import pandas as pd
    rows = {}
    for fq_name, counters in all_results.items():
        label = SAMPLES[fq_name][1]
        for res, c in counters.items():
            total = sum(c.values())
            aa_counts = Counter()
            for codon, cnt in c.items():
                aa_counts[CODON_TABLE.get(codon, 'X')] += cnt
            rows[(label, res)] = {aa: 100.0 * cnt / total for aa, cnt in aa_counts.items()}
    df = pd.DataFrame(rows).T.fillna(0.0)
    df.index.names = ["sample", "residue"]
    csv_path = os.path.join(OUT_DIR, "aa_frequency_table.csv")
    df.to_csv(csv_path)
    print(f"Saved translated frequency table to {csv_path}")


if __name__ == "__main__":
    main()
