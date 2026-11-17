"""
Per-residue evolutionary conservation for the nine screened CDKs, from a deep
homolog search across all domains of life.

Steps: pool the taxon-stratified BLAST hits -> keep clear homologs -> one
representative per organism -> retrieve sequences + lineage from UniProt ->
MAFFT alignment -> Henikoff-weighted per-column conservation projected onto
human residue numbering -> join to the KLIFS pocket positions.
"""
import os, glob, json, subprocess, io, time
import numpy as np
import pandas as pd
import requests

from ._paths import (ALN, BLAST, CONS, EDITING, FIG, KLIFS, LFC, ORTH,  # noqa: F401
                      OUT, ROOT)


CDKS = {"CDK7": "P50613", "CDK8": "P49336", "CDK19": "Q9BWU1", "CDK9": "P50750",
        "CDK12": "Q9NYV4", "CDK13": "Q14004", "CDK2": "P24941", "CDK4": "P11802",
        "CDK6": "Q00534"}
AAS = list("ACDEFGHIKLMNPQRSTVWY")

# Clear-homolog thresholds. Deliberately permissive on identity so that the
# deep branches (fungi, plants, protists) survive, but strict on E-value and on
# subject length so that fragments and spurious kinase-domain hits are removed.
MIN_IDENTITY  = 25.0
MAX_EVALUE    = 1e-15
MIN_LENGTH    = 200
MAX_PER_BIN   = 220         # stops mammals from swamping the alignment
# Bit score per human kinase-domain residue. True CDK orthologs score ~2-3.5;
# generic kinase-fold proteins (archaeal/bacterial Ser-Thr kinases, unrelated
# human kinases) score <1. This is length-independent, so it behaves the same
# for the 346-residue CDK7 and the 1490-residue CDK12.
MIN_BITS_PER_RES = 1.3
MIN_DOMAIN_COV   = 0.60     # applied after alignment, over the Pfam domain

TAX_BINS = [("Mammalia", "Mammalia"), ("Aves", "Sauropsida"), ("Amphibia", "Amphibia"),
            ("Actinopterygii", "Actinopterygii"), ("Chordata", "Chordata"),
            ("Arthropoda", "Arthropoda"), ("Nematoda", "Nematoda"),
            ("Other Metazoa", "Metazoa"), ("Fungi", "Fungi"),
            ("Viridiplantae", "Viridiplantae"), ("Other Eukaryota", "Eukaryota"),
            ("Archaea", "Archaea"), ("Bacteria", "Bacteria")]


def kinase_domain(gene):
    """Pfam protein-kinase domain boundaries in human numbering."""
    j = json.load(open(os.path.join(KLIFS, f"{gene}.json")))
    return int(j["pfam"]["start"]), int(j["pfam"]["end"])


def tax_bin(lineage):
    for label, token in TAX_BINS:
        if token in lineage:
            return label
    return "Unclassified"


# ------------------------------------------------------------- hit pooling --
def pool_hits(gene):
    files = sorted(glob.glob(f"{BLAST}/{gene}_hits.tsv") +
                   glob.glob(f"{BLAST}/{gene}__*_hits.tsv"))
    frames = []
    for f in files:
        d = pd.read_csv(f, sep="\t")
        d["source_db"] = os.path.basename(f).replace(f"{gene}__", "").replace(
            f"{gene}_", "").replace("_hits.tsv", "")
        frames.append(d)
    if not frames:
        raise FileNotFoundError(f"no BLAST results for {gene}")
    d = pd.concat(frames, ignore_index=True)
    d = d.rename(columns={"Identities(%)": "pct_id", "Score(Bits)": "bits",
                          "E()": "evalue", "Length": "subj_len"})
    d = (d.sort_values("bits", ascending=False)
           .drop_duplicates("Accession", keep="first"))
    dstart, dend = kinase_domain(gene)
    d["bits_per_res"] = d["bits"] / (dend - dstart + 1)
    keep = ((d.pct_id >= MIN_IDENTITY) & (d.evalue <= MAX_EVALUE)
            & (d.subj_len >= MIN_LENGTH) & (d.bits_per_res >= MIN_BITS_PER_RES))
    return d[keep].copy(), d


# ------------------------------------------------------ sequence retrieval --
def fetch_uniprot(accs, batch=90):
    rows = []
    for i in range(0, len(accs), batch):
        chunk = accs[i:i + batch]
        for attempt in range(4):
            try:
                r = requests.get(
                    "https://rest.uniprot.org/uniprotkb/accessions",
                    params={"accessions": ",".join(chunk),
                            "fields": "accession,organism_name,organism_id,lineage,sequence",
                            "format": "tsv"}, timeout=180)
                r.raise_for_status()
                rows.append(pd.read_csv(io.StringIO(r.text), sep="\t"))
                break
            except Exception as e:
                if attempt == 3:
                    print(f"    fetch failed for batch {i}: {e}")
                time.sleep(4 * (attempt + 1))
    if not rows:
        return pd.DataFrame()
    d = pd.concat(rows, ignore_index=True)
    return d.rename(columns={"Entry": "Accession", "Organism": "organism",
                             "Organism (ID)": "taxid",
                             "Taxonomic lineage": "lineage", "Sequence": "sequence"})


def build_ortholog_set(gene, uniprot_id):
    hits, raw = pool_hits(gene)
    n_raw = len(raw)
    # record what was searched but rejected as too distant to be a clear homolog
    rej = raw[raw.bits_per_res < MIN_BITS_PER_RES]
    rej.to_csv(f"{ORTH}/{gene}_rejected_distant_hits.csv", index=False)
    seqs = fetch_uniprot(hits["Accession"].tolist())
    d = hits.merge(seqs, on="Accession", how="inner")
    d["lineage"] = d["lineage"].fillna("")
    d["tax_bin"] = d["lineage"].apply(tax_bin)

    # one representative (best bit score) per organism
    d = (d.sort_values("bits", ascending=False)
           .drop_duplicates("taxid", keep="first"))

    # stratified cap so deep clades are not drowned out by mammals
    d = (d.groupby("tax_bin", group_keys=False)
           .apply(lambda g: g.nlargest(min(len(g), MAX_PER_BIN), "bits")))

    # always include the human query itself, first
    human = fetch_uniprot([uniprot_id])
    human["tax_bin"] = "Mammalia"; human["pct_id"] = 100.0
    human["Accession"] = uniprot_id
    d = d[d.Accession != uniprot_id]
    d = pd.concat([human, d], ignore_index=True)

    d.to_csv(f"{ORTH}/{gene}_orthologs.csv", index=False)
    with open(f"{ORTH}/{gene}_orthologs.fasta", "w") as fh:
        for r in d.itertuples():
            if isinstance(r.sequence, str) and len(r.sequence) > 0:
                fh.write(f">{r.Accession}|{str(r.organism).replace(' ','_')}|{r.tax_bin}\n"
                         f"{r.sequence}\n")
    print(f"  {gene}: {n_raw} raw hits -> {len(d)} clear homologs "
          f"across {d.taxid.nunique()} organisms, {d.tax_bin.nunique()} clades "
          f"({len(rej)} rejected as distant kinase-fold only)")
    return d


# ------------------------------------------------------------- alignment ----
def align(gene, threads=8):
    fa, out = f"{ORTH}/{gene}_orthologs.fasta", f"{ALN}/{gene}_aligned.fasta"
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return out
    with open(out, "w") as fh:
        subprocess.run(["mafft", "--auto", "--anysymbol", "--quiet",
                        "--thread", str(threads), fa], stdout=fh, check=True)
    return out


def read_fasta(path):
    names, seqs, name, buf = [], [], None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                names.append(name); seqs.append("".join(buf))
            name, buf = line[1:], []
        else:
            buf.append(line.upper())
    if name is not None:
        names.append(name); seqs.append("".join(buf))
    return names, seqs


# ----------------------------------------------------------- conservation ---
def henikoff_weights(arr):
    """Henikoff & Henikoff position-based sequence weights (redundancy control)."""
    n_seq, n_col = arr.shape
    w = np.zeros(n_seq)
    for j in range(n_col):
        col = arr[:, j]
        vals, counts = np.unique(col[col >= 0], return_counts=True)
        if len(vals) == 0:
            continue
        r = len(vals)
        cmap = dict(zip(vals, counts))
        contrib = np.array([1.0 / (r * cmap[c]) if c in cmap else 0.0 for c in col])
        w += contrib
    return w / w.sum() * n_seq


def conservation_from_alignment(gene, uniprot_id):
    names, seqs = read_fasta(align(gene))
    lut = {a: i for i, a in enumerate(AAS)}
    arr = np.array([[lut.get(c, -1) for c in s] for s in seqs], dtype=np.int16)

    hi = next(i for i, n in enumerate(names) if n.startswith(uniprot_id))
    human = seqs[hi]

    # drop sequences that only partially cover the human kinase domain
    dstart, dend = kinase_domain(gene)
    res_no, dom_cols = 0, []
    for j, hc in enumerate(human):
        if hc == "-":
            continue
        res_no += 1
        if dstart <= res_no <= dend:
            dom_cols.append(j)
    dom_cols = np.array(dom_cols)
    cov = (arr[:, dom_cols] >= 0).mean(axis=1)
    keep = cov >= MIN_DOMAIN_COV
    keep[hi] = True
    n_dropped = int((~keep).sum())
    arr = arr[keep]
    names = [n for n, k in zip(names, keep) if k]
    hi = next(i for i, n in enumerate(names) if n.startswith(uniprot_id))
    print(f"  dropped {n_dropped} partial/fragmentary sequences "
          f"(<{MIN_DOMAIN_COV:.0%} kinase-domain coverage); {len(names)} retained")
    bins = [n.split("|")[-1] for n in names]

    w = henikoff_weights(arr)
    recs, res_no = [], 0
    for j, hc in enumerate(human):
        if hc == "-":
            continue
        res_no += 1
        col = arr[:, j]
        mask = col >= 0
        if mask.sum() == 0:
            continue
        ww = w[mask]
        p = np.bincount(col[mask], weights=ww, minlength=20)
        p = p / p.sum()
        nz = p[p > 0]
        H = -(nz * np.log2(nz)).sum()
        hidx = lut.get(hc, -1)
        recs.append({
            "gene": gene, "uniprot": uniprot_id,
            "human_position": res_no, "human_residue": hc,
            "conservation": 1 - H / np.log2(20),
            "shannon_entropy": H,
            "pct_identical_to_human": 100 * (p[hidx] if hidx >= 0 else np.nan),
            "consensus_residue": AAS[int(np.argmax(p))],
            "consensus_freq": float(p.max()),
            "n_seqs_with_residue": int(mask.sum()),
            "coverage": mask.sum() / len(names),
            "n_clades_covered": len({b for b, m in zip(bins, mask) if m}),
            "in_kinase_domain": dstart <= res_no <= dend,
        })
    d = pd.DataFrame(recs)
    d["n_seqs_total"] = len(names)
    d.to_csv(f"{CONS}/ortholog_conservation_{gene}.csv", index=False)
    return d


def main(genes=None):
    genes = genes or list(CDKS)
    out = []
    for g in genes:
        print(f"[{g}]")
        build_ortholog_set(g, CDKS[g])
        t0 = time.time()
        d = conservation_from_alignment(g, CDKS[g])
        print(f"  aligned + scored {len(d)} residues in {time.time()-t0:.0f}s")
        out.append(d)
    all_d = pd.concat(out, ignore_index=True)
    all_d.to_csv(f"{CONS}/ortholog_conservation_all_CDKs.csv", index=False)
    return all_d


if __name__ == "__main__":
    main()
