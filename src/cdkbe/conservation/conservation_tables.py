"""Build the reader-facing conservation tables (CSV + multi-sheet XLSX)."""
import os, json
import numpy as np
import pandas as pd
from . import kinome_conservation as kc

from ._paths import (ALN, BLAST, CONS, EDITING, FIG, KLIFS, LFC, ORTH,  # noqa: F401
                      OUT, ROOT)
CDKS = ["CDK7", "CDK8", "CDK19", "CDK9", "CDK12", "CDK13", "CDK2", "CDK4", "CDK6"]
GROUPS = ["AGC", "CAMK", "CK1", "CMGC", "STE", "TK", "TKL", "Other", "Atypical"]


def tolerated(freqs, cutoff=0.05):
    """Residues seen at >= cutoff across the kinome, as 'A 41% / V 22%'."""
    items = sorted(freqs.items(), key=lambda kv: -kv[1])
    return " / ".join(f"{a} {100*f:.0f}%" for a, f in items if f >= cutoff)


def main():
    master = pd.read_csv(os.path.join(OUT, "klifs_kinome_conservation.csv"))
    spec   = pd.read_csv(os.path.join(OUT, "klifs_group_specificity.csv"))
    freqs  = json.load(open(os.path.join(OUT, "klifs_position_aa_frequencies.json")))

    t = master.copy()
    t["tolerated_residues_>=5%"] = t["KLIFS_position"].map(
        lambda p: tolerated(freqs[str(p)]))
    t["n_residue_types_seen"] = t["KLIFS_position"].map(lambda p: len(freqs[str(p)]))

    # CDK-family excess constraint (size-matched null)
    cdk = spec[spec.group == "CDK family"].set_index("KLIFS_position")
    t["CDKfam_excess_conservation"] = t["KLIFS_position"].map(cdk["excess_conservation"])
    t["CDKfam_z_conservation"]      = t["KLIFS_position"].map(cdk["z_conservation"])
    t["CDKfam_z_KL_specificity"]    = t["KLIFS_position"].map(cdk["z_KL"])

    cols = (["KLIFS_position", "KLIFS_label", "region", "consensus_aa", "consensus_freq",
             "conservation", "shannon_entropy", "mutability_class",
             "tolerated_residues_>=5%", "n_residue_types_seen", "n_aa_90pct",
             "class_conservation", "majority_class", "frac_gap", "n_scored"]
            + [f"cons_{g}" for g in GROUPS]
            + ["cons_CDKfamily", "consensus_CDKfamily",
               "CDKfam_excess_conservation", "CDKfam_z_conservation",
               "CDKfam_z_KL_specificity"]
            + [f"{c}_residue" for c in CDKS])
    t = t[cols].round(4)
    t.to_csv(os.path.join(OUT, "TABLE_kinome_KLIFS_conservation.csv"), index=False)

    # ---- derived views
    immutable = t[t.mutability_class.isin(["Immutable (invariant)", "Strongly constrained"])] \
                  .sort_values("conservation", ascending=False)
    mutable   = t[t.mutability_class.isin(["Tolerant", "Hypervariable"])] \
                  .sort_values("conservation")

    by_group = spec.pivot(index="KLIFS_position", columns="group",
                          values="conservation").round(3)
    z_group  = spec.pivot(index="KLIFS_position", columns="group",
                          values="z_conservation").round(2)
    kl_group = spec.pivot(index="KLIFS_position", columns="group",
                          values="z_KL").round(2)
    for d in (by_group, z_group, kl_group):
        d.insert(0, "region", t.set_index("KLIFS_position")["region"])
        d.insert(1, "kinome_consensus", t.set_index("KLIFS_position")["consensus_aa"])

    cdk_map = t[["KLIFS_position", "KLIFS_label", "region", "consensus_aa",
                 "conservation", "mutability_class", "consensus_CDKfamily"]
                + [f"{c}_residue" for c in CDKS]]

    xlsx = os.path.join(OUT, "TABLE_kinome_KLIFS_conservation.xlsx")
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        t.to_excel(w, "Kinome_85_positions", index=False)
        immutable.to_excel(w, "Cannot_be_mutated", index=False)
        mutable.to_excel(w, "Can_be_mutated", index=False)
        by_group.to_excel(w, "Conservation_by_group")
        z_group.to_excel(w, "Group_excess_constraint_z")
        kl_group.to_excel(w, "Group_specificity_KL_z")
        cdk_map.to_excel(w, "Screened_CDK_residue_map", index=False)
    print("wrote", xlsx)

    print(f"\nIMMUTABLE / strongly constrained: {len(immutable)}/85 positions")
    print(immutable[["KLIFS_label", "consensus_aa", "consensus_freq", "conservation",
                     "tolerated_residues_>=5%"]].head(20).to_string(index=False))
    print(f"\nMOST MUTABLE: {len(mutable)}/85 positions")
    print(mutable[["KLIFS_label", "consensus_aa", "consensus_freq", "conservation",
                   "n_residue_types_seen", "tolerated_residues_>=5%"]].head(15).to_string(index=False))
    return t


if __name__ == "__main__":
    main()
