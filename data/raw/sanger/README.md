# Sanger sequencing

Every Sanger trace this repository ships, with one table over all of them.

`INDEX.csv` has one row per `.ab1`: where it lives, the sample label and
sequencing primer the vendor recorded, the gene, variant, treatment, selection
state and read direction parsed out of that label, and the run date,
instrument, well, order number, dye set, read length and median Phred score
read out of the ABIF file itself. `scripts/make_sanger_index.py` writes it and
documents every parsing rule; nothing here is typed by hand.

`basecalls.fasta.gz` carries the vendor basecall for each trace under the same
path, so the set can be searched or aligned without trace-viewing software.
The chromatograms are the primary record -- a base-edited locus is a mixed
peak, which a basecall flattens -- so quantification should go back to the
`.ab1`.

## What is here

| directory | what it is |
|---|---|
| `endogenous_editing/` (187 traces) | the edited locus amplified out of treated cells: the validation and post-selection genotyping behind figures 2, 3 and 4. Also holds the sequencing primer tables and the SnapGene transcript maps |
| `library_colonies/` (45 traces) | colonies from the library cloning, read at the U6 promoter -- the protospacer/barcode pairing QC. `colony_seq_identification.csv` is the call per colony |
| `guide_plasmids/` (56 traces) | the individual validation guide vectors, read at the U6 promoter. These ship under `data/reference/validation_guides/`, beside the designs they confirm, and are indexed here so that all three sets appear in one table |

## Reading them

    pip install biopython
    from Bio import SeqIO
    record = SeqIO.read("data/raw/sanger/<path from INDEX.csv>", "abi")

`record.annotations["abif_raw"]` holds the full ABIF header, including the
four dye traces, for quantifying a mixed peak.

## Provenance

The traces are vendor deliveries, copied unchanged -- `data/MANIFEST.csv`
records each one's source path and SHA-256, and `scripts/verify_data.py`
checks them. The sequencing primers are in
`endogenous_editing/cdk_seq_primers.xlsx` and
`endogenous_editing/cdk246_seq_primers.xlsx`, and the SnapGene maps beside
them show where each primer sits in the transcript.
