"""
01_extract_reference_fasta.py

Extracts plain FASTA reference sequences from the annotated amplicon
GenBank files, for use as minimap2 alignment references.

Usage:
    python 01_extract_reference_fasta.py

Requires: biopython
"""
import os
from Bio import SeqIO

HERE = os.path.dirname(os.path.abspath(__file__))
REF_DIR = os.path.join(HERE, "..", "data", "reference")

GENBANK_FILES = {
    "amp1": "CDK12_3351_2232_amplicon_annotated.gb",
    "amp2": "CDK12_3394_3395_3396_amplicon_annotated.gb",
}

def main():
    for amp_name, gb_file in GENBANK_FILES.items():
        gb_path = os.path.join(REF_DIR, gb_file)
        record = SeqIO.read(gb_path, "genbank")
        fasta_path = os.path.join(REF_DIR, f"ref_{amp_name}.fasta")
        with open(fasta_path, "w") as f:
            f.write(f">{amp_name}\n{str(record.seq)}\n")
        print(f"Wrote {fasta_path} ({len(record.seq)} bp)")

if __name__ == "__main__":
    main()
