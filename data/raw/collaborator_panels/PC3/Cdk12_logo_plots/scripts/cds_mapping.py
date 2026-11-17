"""
cds_mapping.py

Maps CDK12 (ENST00000447079 / CDK12-202) residue numbers to genomic
nucleotide coordinates, using the CDS exon-segment structure pulled from
the original SnapGene mapped-primers file.

Coordinates are 1-based, inclusive, on the same coordinate system as the
genomic sequence in `data/reference/*.gb` (i.e. positions within the
73,059 bp mapped region, matching what SnapGene reports).
"""

# CDS exon segments (genomic_start, genomic_end), 1-based inclusive,
# in transcript/coding order. Extracted from the SnapGene CDS feature
# block of CDK12-202_mapped_primers.dna.
CDS_SEGMENTS = [
    (587, 1632), (9394, 10278), (29072, 29248), (31266, 31405), (33039, 33209),
    (39765, 39954), (48220, 48276), (50044, 50145), (54246, 54323), (55955, 56071),
    (58471, 58602), (63189, 63400), (64379, 64831), (69119, 69831),
]

# Genomic start coordinate of each amplicon (matches primer_bind F start
# in the annotated GenBank files), used to convert genomic -> local
# amplicon coordinates.
AMPLICON_STARTS = {
    "amp1": 31215,  # CDK12_3351_2232 amplicon (covers residues ~703-750)
    "amp2": 39658,  # CDK12_3394_3395_3396 amplicon (covers residues ~807-870)
}


def codon_genomic_range(residue_num):
    """
    Return (genomic_start, genomic_end), 1-based inclusive, for the 3 nt
    of the codon encoding `residue_num`. Only resolves correctly for
    codons fully contained within a single CDS exon segment (true for
    all residues used in this analysis: 731-733 and 816-822).
    """
    cds_start = 3 * (residue_num - 1) + 1
    cds_end = 3 * residue_num
    cum = 0
    for seg_start, seg_end in CDS_SEGMENTS:
        seg_len = seg_end - seg_start + 1
        if cum < cds_end <= cum + seg_len and cum < cds_start <= cum + seg_len:
            off_start = cds_start - cum
            off_end = cds_end - cum
            return seg_start + off_start - 1, seg_start + off_end - 1
        cum += seg_len
    raise ValueError(
        f"Residue {residue_num} codon spans a splice junction or is "
        "outside the mapped CDS segments; not supported by this helper."
    )


def codon_local_range(residue_num, amplicon):
    """
    Return (local_start, local_end) 0-based half-open coordinates of the
    codon for `residue_num` within the given amplicon's own sequence
    (i.e. coordinates you can slice the amplicon FASTA/GenBank sequence
    with, or use directly against minimap2 alignments to that reference).
    """
    gs, ge = codon_genomic_range(residue_num)
    amp_start = AMPLICON_STARTS[amplicon]
    return gs - amp_start, ge - amp_start + 1


if __name__ == "__main__":
    # sanity check against the known validated mutation sites
    for r in [731, 732, 733, 816, 817, 818, 819, 820, 821, 822]:
        print(r, codon_genomic_range(r))
