"""
Align the nine CDK kinase domains together with the ten kinases that carry
clinically observed resistance mutations (OncoKB R1/R2).

The CDK-only alignment used everywhere else
(`source_data/16_KLIFS/clustalo-*.aln-clustal_num`) has no coordinate for EGFR,
ABL1 and the rest, so clinical resistance can only be placed on it through the
85 KLIFS pocket indices.  This script builds a second, wider alignment in which
every one of those kinases has its own row, so a clinical mutation anywhere in
the kinase domain -- not just at a numbered pocket position -- lands on a
column that the CDK resistance profile can be read against.

Sequences and domain boundaries both come from the KLIFS json files already in
the repository (`uniprot.canonical_seq`, `pfam.start/end`), so nothing is
fetched from the network.  The domain is padded on both sides because several
CDK resistance sites sit just outside the Pfam kinase domain.

    python3 build_kinase_domain_msa.py

Output: source_data/16_KLIFS/kinase_domain_with_clinical.aln.fa
        source_data/16_KLIFS/kinase_domain_with_clinical_offsets.csv
"""
import json, os, subprocess, sys

import pandas as pd
from Bio import AlignIO
from .. import paths

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
# paths resolved from the repository root -- see _paths.py
from ._paths import (LIBRARY, KLIFS_DIR, DISORDER, ALPHAMISSENSE, DBNSFP, LFC_DIR, EDITING, FIG_DIR)  # noqa: E402
from ._paths import KLIFS_DIR, LIBRARY_GENERATION  # noqa: E402  repo-root paths

CDKS = ["CDK7", "CDK8", "CDK19", "CDK9", "CDK12", "CDK13", "CDK2", "CDK4", "CDK6"]
# the kinases with OncoKB R1/R2 resistance mutations, one json each
CLINICAL = ["ABL1", "ALK", "BTK", "EGFR", "KIT", "MET", "NTRK1", "NTRK3", "PDGFRA", "ROS1"]

# CDK resistance reaches past the Pfam domain on both sides (the CDK12/13 C-lobe
# sites in particular), so the slice is widened rather than cut at the domain edge
PAD = 60


def load_json(gene):
    for sub in ("", "oncokb_klifs"):
        p = os.path.join(KLIFS_DIR, sub, f"{gene}.json")
        if os.path.exists(p):
            return json.load(open(p))
    sys.exit(f"no KLIFS json for {gene}")


def domain_slice(gene):
    """(sequence of the padded kinase domain, 1-based residue number of its first residue)."""
    d = load_json(gene)
    seq = d["uniprot"]["canonical_seq"]
    lo = max(1, d["pfam"]["start"] - PAD)
    hi = min(len(seq), d["pfam"]["end"] + PAD)
    return seq[lo - 1:hi], lo, seq


def main():
    cdk_info = pd.read_csv(
        os.path.join(LIBRARY_GENERATION,
                     "cdks_info_updated.csv")).set_index("Gene")

    records, rows = [], []
    for gene in CDKS + CLINICAL:
        dom, lo, full = domain_slice(gene)
        if gene in CDKS:
            # the library's codon numbers are what the resistance calls use, so the
            # json sequence has to be the same protein the library was built on
            ref = cdk_info.loc[gene, "Protein"]
            assert full == ref, f"{gene}: KLIFS json sequence != library reference protein"
        records.append(f">{gene}\n{dom}")
        rows.append(dict(gene=gene, kind="CDK" if gene in CDKS else "clinical",
                         domain_start=lo, domain_end=lo + len(dom) - 1,
                         domain_len=len(dom), protein_len=len(full)))
        print(f"  {gene:<7} {rows[-1]['kind']:<8} residues {lo}-{lo + len(dom) - 1} "
              f"of {len(full)}")

    raw = os.path.join(KLIFS_DIR, "kinase_domain_with_clinical.input.fa")
    out = os.path.join(KLIFS_DIR, "kinase_domain_with_clinical.aln.fa")
    open(raw, "w").write("\n".join(records) + "\n")

    print(f"\nrunning mafft on {len(records)} kinase domains ...")
    with open(out, "w") as fh:
        subprocess.run(["mafft", "--auto", "--quiet", raw], stdout=fh, check=True)
    os.remove(raw)

    aln = AlignIO.read(out, "fasta")
    pd.DataFrame(rows).to_csv(
        os.path.join(KLIFS_DIR, "kinase_domain_with_clinical_offsets.csv"), index=False)

    print(f"alignment length {aln.get_alignment_length()}, {len(aln)} sequences")
    print(f"wrote {paths.rel(out)}")


if __name__ == "__main__":
    main()
