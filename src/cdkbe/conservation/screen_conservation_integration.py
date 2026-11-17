"""
Does evolutionary conservation predict what a base edit does in the screen?

Two outcomes are tested per gRNA:
  resistance  — enriched under the inhibitor that targets its CDK (drug vs DMSO)
  fitness     — depleted in the absence of drug (T0 vs plasmid, DMSO vs T0)

Edit identity comes from the canonical HGVSp editing tables and gRNA editing
efficiency is carried through the whole analysis, both as a filter and as a
covariate, because a guide that does not edit cannot score in either direction.

Conventions (thresholds, edit assignment, drug dictionary) follow
notebooks/09_analysis_and_viz/fig1_analysis.ipynb.
"""
import os, re, json
import numpy as np
import pandas as pd

from .. import editing as _editing
from ._paths import (ALN, BLAST, CONS, EDITING, FIG, KLIFS, LFC, ORTH,  # noqa: F401
                      OUT, ROOT)


GENES = ["CDK7", "CDK8", "CDK19", "CDK9", "CDK12", "CDK13", "CDK2", "CDK4", "CDK6"]
SAMPLES = ["CBE_subpool1", "ABE_subpool1", "CBE_CDK12_13", "ABE_CDK12_13",
           "CBE_CDK2_4_6", "ABE_CDK2_4_6"]

# CDK -> the inhibitors that target it (drug vs DMSO contrasts)
DRUG_DICT = {
    "CDK9":  ["KB-0742_1500nM_DMSO", "KI-CDK9d-32N_1250nM_DMSO", "KI-CDK9d-32N_5000nM_DMSO",
              "KI-CDK9d-32_1000nM_DMSO", "KI-CDK9d-32_100nM_DMSO"],
    "CDK7":  ["SY-5609_100nM_DMSO", "SY-5609_10nM_DMSO"],
    "CDK8":  ["SEL120_4000nM_DMSO", "Senexin B_15000nM_DMSO", "Senexin B_2000nM_DMSO"],
    "CDK19": ["SEL120_4000nM_DMSO", "Senexin B_15000nM_DMSO", "Senexin B_2000nM_DMSO"],
    "CDK12": ["BSJ-4-116_DMSO", "CDK12-IN-2_DMSO", "HQ461_DMSO"],
    "CDK13": ["BSJ-4-116_DMSO", "CDK12-IN-2_DMSO", "HQ461_DMSO"],
    "CDK2":  ["INX_DMSO", "TAGTO_DMSO"],
    "CDK4":  ["ABEMA_DMSO", "ATIRMO_DMSO", "PALBO_DMSO", "RIBO_DMSO"],
    "CDK6":  ["ABEMA_DMSO", "ATIRMO_DMSO", "PALBO_DMSO", "RIBO_DMSO"],
}
FITNESS_CONTRASTS = ["T0_Plasmid", "DMSO_T0"]

FDR_CUTOFF   = 0.01
LFC_CUTOFF   = 1.0     # enrichment for resistance / depletion for fitness
MIN_PCT_EDIT = 10.0    # gRNA editing efficiency floor (benchmark)
MIN_READS    = 100     # sensor read-depth floor (benchmark)


# --------------------------------------------------------------- edits ------
def extract_numbers(text):
    """First number in an HGVSp string. Shared implementation."""
    return _editing.first_number(text)


def top_HGVSp_compiler(edit):
    """Top non-WT edit per gRNA plus that gRNA's WT fraction.

    This was a sixth copy of the same reduction; it delegates to the shared
    implementation, which uses a stable sort so guides whose top two variants
    tie on %Reads keep the call the original `nlargest` gave them.
    """
    return _editing.top_hgvsp(edit)


def build_edits():
    frames = []
    for s in SAMPLES:
        e = pd.read_csv(f"{EDITING}/{s}_HGVSp_canonical.zip")
        t = top_HGVSp_compiler(e)
        t["sample"] = s
        t["Editor"] = s.split("_")[0]
        frames.append(t)
    d = pd.concat(frames, ignore_index=True)
    d = d[d.HGVSp != "WT"].copy()

    # editing efficiency thresholds (benchmark)
    n_before = len(d)
    d = d[(d["%Reads"] >= MIN_PCT_EDIT) & (d["#Reads"] >= MIN_READS)].copy()

    # unambiguous single substitutions only, so a residue can be attributed
    d["is_compound"] = d["HGVSp"].str.contains("_")
    d = d[~d.is_compound].copy()
    m = d["HGVSp"].str.extract(r"^([A-Z])(\d+)([A-Z*])$")
    d["ref_aa"], d["position"], d["alt_aa"] = m[0], m[1].astype(int), m[2]
    d["consequence"] = np.where(d.alt_aa == "*", "nonsense",
                        np.where(d.alt_aa == d.ref_aa, "silent", "missense"))
    d["gene"] = d["gRNA_id"].str.split("_").str[1]
    d = d.rename(columns={"%Reads": "edit_efficiency", "#Reads": "sensor_reads"})
    print(f"edits: {n_before} guides with an edit -> {len(d)} single-substitution guides "
          f"passing >={MIN_PCT_EDIT}% editing and >={MIN_READS} reads")
    return d[["gRNA_id", "gene", "Editor", "sample", "HGVSp", "ref_aa", "position",
              "alt_aa", "consequence", "edit_efficiency", "sensor_reads", "WT_perc"]]


# ------------------------------------------------------------- outcomes -----
def build_outcomes():
    """
    One row per guide. The master LFC tables are stacked: a guide appears once per
    (Screen_ID, reference) block, and each block populates only its own columns, so
    a guide's drug and fitness contrasts live on different rows and sometimes come
    from different screens. Everything is therefore aggregated back to the guide:
    a guide is a hit if any of its rows clears the thresholds.
    """
    frames = []
    for editor in ["ABE", "CBE"]:
        d = pd.read_csv(f"{LFC}/{editor}_master_LFC_table.csv")
        d["gene"] = d["gRNA_id"].str.split("_").str[1]
        d = d[(d.classification == "targeting") & d.gene.isin(GENES)].copy()

        # --- resistance: enriched under an inhibitor that targets this guide's CDK
        d["_res_lfc"] = np.nan
        d["_res_hit"] = False
        for gene, conds in DRUG_DICT.items():
            m = d.gene == gene
            if not m.any():
                continue
            lcols = [f"LFC_median_{c}" for c in conds if f"LFC_median_{c}" in d.columns]
            fcols = [f"FDR_{c}" for c in conds if f"FDR_{c}" in d.columns]
            if not lcols:
                continue
            L = d.loc[m, lcols].astype(float)
            F = d.loc[m, fcols].astype(float)
            d.loc[m, "_res_lfc"] = L.max(axis=1)
            d.loc[m, "_res_hit"] = ((F.values < FDR_CUTOFF) &
                                    (L.values > LFC_CUTOFF)).any(axis=1)

        # --- fitness: depleted in the absence of drug
        lf = d[[f"LFC_median_{c}" for c in FITNESS_CONTRASTS]].astype(float)
        ff = d[[f"FDR_{c}" for c in FITNESS_CONTRASTS]].astype(float)
        d["_fit_lfc"] = lf.min(axis=1)
        d["_fit_hit"] = ((ff.values < FDR_CUTOFF) & (lf.values < -LFC_CUTOFF)).any(axis=1)
        for c in FITNESS_CONTRASTS:
            d[f"LFC_{c}"] = d[f"LFC_median_{c}"].astype(float)
            d[f"FDRv_{c}"] = d[f"FDR_{c}"].astype(float)

        d["Editor"] = editor
        agg = {"_res_lfc": "max", "_res_hit": "any", "_fit_lfc": "min", "_fit_hit": "any"}
        for c in FITNESS_CONTRASTS:
            agg[f"LFC_{c}"] = "min"
            agg[f"FDRv_{c}"] = "min"
        g = (d.groupby(["gRNA_id", "gene", "Editor"], as_index=False)
               .agg(agg)
               .rename(columns={"_res_lfc": "max_drug_LFC", "_res_hit": "resistance_hit",
                                "_fit_lfc": "min_fitness_LFC", "_fit_hit": "fitness_hit"}))
        print(f"  {editor}: {len(d)} table rows -> {len(g)} unique targeting guides")
        frames.append(g)

    out = pd.concat(frames, ignore_index=True)
    print(f"outcomes: {len(out)} guides; "
          f"{out.resistance_hit.sum()} resistance, {out.fitness_hit.sum()} fitness hits")
    return out


# --------------------------------------------------------- conservation -----
def build_conservation():
    cons = pd.read_csv(f"{CONS}/ortholog_conservation_all_CDKs.csv")
    cons = cons.rename(columns={"human_position": "position",
                                "conservation": "ortholog_conservation"})
    klifs = pd.read_csv(os.path.join(KLIFS, "KLIFS_compiled.csv"))
    klifs["KLIFS_position"] = klifs["KLIFS_idx"].str.split(":").str[1].astype(int)
    kin = pd.read_csv(f"{CONS}/TABLE_kinome_KLIFS_conservation.csv")

    rows = []
    for _, r in klifs.iterrows():
        for g in GENES:
            if pd.notna(r[g]):
                rows.append({"gene": g, "position": int(r[g]),
                             "KLIFS_position": int(r["KLIFS_position"])})
    kmap = pd.DataFrame(rows).merge(
        kin[["KLIFS_position", "KLIFS_label", "region", "conservation", "mutability_class"]]
           .rename(columns={"conservation": "kinome_conservation"}),
        on="KLIFS_position", how="left")

    c = cons[["gene", "position", "human_residue", "ortholog_conservation",
              "pct_identical_to_human", "in_kinase_domain", "n_seqs_total"]] \
        .merge(kmap, on=["gene", "position"], how="left")
    c["in_KLIFS_pocket"] = c.KLIFS_position.notna()
    return c


def site_class(r):
    """
    Two-axis call. Only defined inside the KLIFS pocket, because the kinome axis
    needs a position that is alignable across kinases; elsewhere there is no
    cross-kinase comparison to make.
    """
    o, k = r["ortholog_conservation"], r["kinome_conservation"]
    if pd.isna(o):
        return "unmapped"
    if pd.isna(k):
        return "outside KLIFS pocket"
    if o >= 0.75 and k >= 0.50:
        return "Universally immutable"
    if o >= 0.75:
        return "CDK-specific, fixed"
    if o < 0.45:
        return "Tolerant in evolution"
    return "Intermediate"


def main():
    edits = build_edits()
    outcomes = build_outcomes()
    cons = build_conservation()

    d = (edits.merge(outcomes, on=["gRNA_id", "gene", "Editor"], how="inner")
              .merge(cons, on=["gene", "position"], how="left"))
    # sanity: the edited reference residue should match the UniProt sequence
    ok = d.human_residue.notna()
    agree = (d.loc[ok, "ref_aa"] == d.loc[ok, "human_residue"]).mean()
    print(f"merged: {len(d)} guides; reference residue agrees with UniProt for "
          f"{agree:.1%} of mapped guides")
    d = d[d.human_residue.isna() | (d.ref_aa == d.human_residue)].copy()

    d["site_class"] = d.apply(site_class, axis=1)
    d["edit_eff_tertile"] = pd.qcut(d.edit_efficiency, 3,
                                    labels=["low", "mid", "high"], duplicates="drop")
    d.to_csv(f"{CONS}/screen_conservation_merged.csv", index=False)
    print(f"wrote {CONS}/screen_conservation_merged.csv  ({len(d)} guides)")
    return d


if __name__ == "__main__":
    main()
