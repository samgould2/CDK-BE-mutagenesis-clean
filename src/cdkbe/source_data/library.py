"""Figure 1 panels: library coverage and the resistance variants it found.

1e — what fraction of each CDK's codons and single-amino-acid variants the
     library reaches, as a function of the sensor-editing cutoff.
     From editing_analysis_v2.ipynb cell 16, which looped 101 cutoffs x 9 genes
     and wrote into a pre-allocated wide frame. `plotting_df` was never saved,
     so this is a full re-derivation.

1f — every resistance variant, placed on its CDK's residue axis.
     From fig1_analysis.ipynb cells 4-13.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .. import counts, editing, io, paths

GENES = ["CDK7", "CDK8", "CDK19", "CDK9", "CDK12", "CDK13", "CDK2", "CDK4", "CDK6"]
AAS = list("ACDEFGHIKLMNPQRSTVWY") + ["*"]
POOL_OF = {"CDK7": "1", "CDK8": "1", "CDK19": "1", "CDK9": "1",
           "CDK12": "2", "CDK13": "2", "CDK2": "3", "CDK4": "3", "CDK6": "3"}

#: the drugs each CDK was screened against, for the "is this a resistance hit" call
DRUGS = {
    "CDK9": ["KB-0742_1500nM_DMSO", "KI-CDK9d-32N_1250nM_DMSO", "KI-CDK9d-32N_5000nM_DMSO",
             "KI-CDK9d-32_1000nM_DMSO", "KI-CDK9d-32_100nM_DMSO"],
    "CDK7": ["SY-5609_100nM_DMSO", "SY-5609_10nM_DMSO"],
    "CDK8": ["SEL120_4000nM_DMSO", "Senexin B_15000nM_DMSO", "Senexin B_2000nM_DMSO"],
    "CDK19": ["SEL120_4000nM_DMSO", "Senexin B_15000nM_DMSO", "Senexin B_2000nM_DMSO"],
    "CDK12": ["BSJ-4-116_DMSO", "CDK12-IN-2_DMSO", "HQ461_DMSO"],
    "CDK13": ["BSJ-4-116_DMSO", "CDK12-IN-2_DMSO", "HQ461_DMSO"],
    "CDK2": ["INX_DMSO", "TAGTO_DMSO"],
    "CDK4": ["ABEMA_DMSO", "ATIRMO_DMSO", "PALBO_DMSO", "RIBO_DMSO"],
    "CDK6": ["ABEMA_DMSO", "ATIRMO_DMSO", "PALBO_DMSO", "RIBO_DMSO"],
}

# ------------------------------------------------------------------ 1e


def _possible(cdk_info: pd.DataFrame):
    """Every codon and every single-aa substitution each protein could have."""
    codons, variants = {}, {}
    for g in GENES:
        seq = cdk_info.loc[cdk_info["Gene"] == g, "Protein"].values[0]
        codons[g] = [i + 1 for i in range(len(seq))]
        variants[g] = [f"{aa}{i + 1}{alt}" for i, aa in enumerate(seq)
                       for alt in AAS if alt != aa]
    return codons, variants


def _savs(sample: str, library: pd.DataFrame) -> pd.DataFrame:
    d = io.read(paths.EDITING / "SAVs_canonical" / f"{sample}_canonical.csv")
    return pd.merge(d, library[["gRNA_id", "Gene"]], on="gRNA_id")


def _codon_rates(savs: pd.DataFrame, library: pd.DataFrame) -> pd.DataFrame:
    return editing.codon_rates(savs, library)


def _venn(observed_cbe, observed_abe, possible) -> tuple:
    """(CBE-only, ABE-only, shared, missing) against what the protein could have."""
    s1, s2, s3 = set(observed_cbe), set(observed_abe), set(possible)
    return (len(s2 & s3 - s1),            # ABE only, and real
            len(s1 & s3 - s2),            # CBE only, and real
            len(s1 & s2 & s3),            # both
            len(s3 - s1 - s2))            # never reached


def figure_1e(cutoffs=range(0, 101)) -> pd.DataFrame:
    """Coverage of codons and variants vs the sensor-editing cutoff, per CDK."""
    library = counts.library()
    cdk_info = io.read(paths.REFERENCE / "library_generation" / "cdks_info_updated.csv")
    cdk_info = cdk_info[cdk_info["Gene"].isin(GENES)].reset_index(drop=True)
    possible_codons, possible_vars = _possible(cdk_info)

    savs = {f"{e}_subpool{p}": _savs(f"{e}_subpool{p}", library)
            for e in ("ABE", "CBE") for p in ("1", "2", "3")}
    codons = {k: _codon_rates(v, library) for k, v in savs.items()}

    rows = []
    for cutoff in cutoffs:
        for g in GENES:
            p = POOL_OF[g]
            A, C = savs[f"ABE_subpool{p}"], savs[f"CBE_subpool{p}"]
            cA, cC = codons[f"ABE_subpool{p}"], codons[f"CBE_subpool{p}"]
            row = dict(editing_cutoff=cutoff, Gene=g)
            for tag, (obs_c, obs_a, pos) in {
                    "vars": (C[(C["Gene"] == g) & (C["%Reads"] >= cutoff)]["HGVSp"],
                             A[(A["Gene"] == g) & (A["%Reads"] >= cutoff)]["HGVSp"],
                             possible_vars[g]),
                    "codons": (cC[(cC["Gene"] == g) & (cC["%Reads"] >= cutoff)]["Codon"],
                               cA[(cA["Gene"] == g) & (cA["%Reads"] >= cutoff)]["Codon"],
                               possible_codons[g])}.items():
                abe_only, cbe_only, shared, missing = _venn(obs_c, obs_a, pos)
                total = abe_only + cbe_only + shared + missing
                row[f"total_{tag}"] = 100 * (abe_only + cbe_only + shared) / total
                row[f"ABE_{tag}"] = 100 * (abe_only + shared) / total
                row[f"CBE_{tag}"] = 100 * (cbe_only + shared) / total
            rows.append(row)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ 1f

FDR_CUTOFF = 0.01
LFC_CUTOFF = 1
MIN_EDIT_PERC = 10      # a guide below this has no usable edit
MIN_SENSOR_READS = 100


def figure_1f() -> pd.DataFrame:
    """Resistance variants by residue: max LFC across that CDK's own drugs."""
    from .. import master_table as mt

    samples = ["CBE_subpool1", "ABE_subpool1", "CBE_CDK12_13", "ABE_CDK12_13",
               "CBE_CDK2_4_6", "ABE_CDK2_4_6"]
    edits = {s: editing.top_hgvsp(io.read(paths.EDITING / f"{s}_HGVSp_canonical.csv"))
             for s in samples}

    lfc = {e: io.read(mt.path(e)) for e in ("ABE", "CBE")}
    for d in lfc.values():
        d["gene"] = d["gRNA_id"].str.split("_").str[1]

    keep = ["gRNA_id", "Screen_ID", "Editor", "Base_Value", "gene", "max_LFC", "min_FDR"]
    parts = []
    for gene in GENES:
        conds = DRUGS[gene]
        lfc_cols = [f"LFC_median_{c}" for c in conds]
        fdr_cols = [f"FDR_{c}" for c in conds]
        for editor in ("CBE", "ABE"):
            d = lfc[editor]
            s = d[(d["gene"] == gene) & (d["classification"] == "targeting")].reset_index(drop=True)
            s = s[(s[fdr_cols] < FDR_CUTOFF).any(axis=1)].reset_index(drop=True)
            s = s[(s[lfc_cols] > LFC_CUTOFF).any(axis=1)].reset_index(drop=True)
            s = s.dropna(axis=1, how="all").reset_index(drop=True)
            if not len(s):
                continue
            s["max_LFC"] = s[lfc_cols].max(axis=1)
            s["min_FDR"] = s[fdr_cols].min(axis=1)
            parts.append(s[keep])
    comb = pd.concat(parts)

    pools = {"CBE": ["CBE_subpool1", "CBE_CDK12_13", "CBE_CDK2_4_6"],
             "ABE": ["ABE_subpool1", "ABE_CDK12_13", "ABE_CDK2_4_6"]}
    out = []
    for editor in ("CBE", "ABE"):
        sub = comb[comb["Editor"] == editor]
        rel = pd.concat([edits[s][edits[s]["gRNA_id"].isin(sub["gRNA_id"])]
                         for s in pools[editor]]).reset_index(drop=True)
        # a guide that barely edits, or has almost no sensor reads, cannot be placed
        rel = rel[(rel["%Reads"] >= MIN_EDIT_PERC) & (rel["#Reads"] >= MIN_SENSOR_READS)]
        merged = pd.merge(sub, rel, on="gRNA_id")
        merged = merged[merged["Codon"] != "WT"].copy()
        merged["Editor"] = editor
        out.append(merged[["gRNA_id", "Editor", "gene", "Codon", "HGVSp", "max_LFC",
                           "min_FDR", "%Reads", "#Reads", "WT_perc"]])
    d = pd.concat(out)
    d = d[d["gene"].isin(GENES)].copy()
    # the panel prints each gene's canonical length beside its track, and the
    # tracks are scaled by it, so the length belongs in the source table too
    from .library_figures import protein_lengths

    d["protein_length"] = d["gene"].map(protein_lengths())
    return d
