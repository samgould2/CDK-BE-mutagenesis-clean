"""
Kinome-wide evolutionary/structural conservation at the 85 KLIFS pocket positions.

Pulls the KLIFS alignment for every human kinase, computes per-position conservation
for the whole kinome and within each kinase group, and classifies positions as
mutable vs immutable. Writes tables to data/reference/conservation and figures to figures/.
"""
import json, os, requests
import numpy as np
import pandas as pd

from ._paths import (ALN, BLAST, CONS, EDITING, FIG, KLIFS, LFC, ORTH,  # noqa: F401
                      OUT, ROOT)


AAS = list("ACDEFGHIKLMNPQRSTVWY")

# KLIFS 85-position region layout (KLIFS numbering convention)
REGIONS = ([("I", 3), ("g.l", 6), ("II", 4), ("III", 6), ("aC", 11), ("b.l", 7),
            ("IV", 4), ("V", 3), ("GK", 1), ("hinge", 3), ("linker", 4), ("aD", 7),
            ("aE", 5), ("VI", 3), ("c.l", 8), ("VII", 3), ("VIII", 1),
            ("xDFG", 4), ("a.l", 2)])
REGION_OF_POS = []
for name, n in REGIONS:
    REGION_OF_POS += [name] * n
assert len(REGION_OF_POS) == 85

# Physicochemical classes used for "conservative substitution" tolerance
AA_CLASS = {}
for cls, members in {"hydrophobic": "AVLIMC", "aromatic": "FWY", "polar": "STNQ",
                     "positive": "KRH", "negative": "DE", "special": "GP"}.items():
    for a in members:
        AA_CLASS[a] = cls


# ------------------------------------------------------------------ data ----
def fetch_klifs(cache=None):
    """All human kinases in KLIFS with their 85-residue pocket alignment."""
    cache = cache or os.path.join(OUT, "klifs_human_kinases.json")
    if os.path.exists(cache):
        return json.load(open(cache))
    r = requests.get("https://klifs.net/api/kinase_information",
                     params={"species": "HUMAN"},
                     headers={"accept": "application/json"}, timeout=120)
    r.raise_for_status()
    json.dump(r.json(), open(cache, "w"))
    return r.json()


def build_matrix(kinases):
    """DataFrame: rows = kinases, columns = KLIFS positions 1..85."""
    rows, meta = [], []
    for k in kinases:
        p = k.get("pocket")
        if not isinstance(p, str) or len(p) != 85:
            continue
        rows.append(list(p))
        meta.append({"name": k["name"], "HGNC": k["HGNC"], "uniprot": k["uniprot"],
                     "group": k["group"], "family": k["family"]})
    mat = pd.DataFrame(rows, columns=range(1, 86))
    info = pd.DataFrame(meta)
    return mat, info


# ------------------------------------------------------- conservation -------
def position_stats(col):
    """Conservation metrics for one aligned column (a pandas Series of chars)."""
    vals = col.values
    n_tot = len(vals)
    aa = np.array([v for v in vals if v in AA_CLASS])          # drop gaps/'X'
    n = len(aa)
    out = {"n_seqs": n_tot, "n_scored": n, "frac_gap": 1 - n / n_tot}
    if n == 0:
        return {**out, "shannon_entropy": np.nan, "conservation": np.nan,
                "consensus_aa": "-", "consensus_freq": np.nan,
                "n_aa_90pct": np.nan, "class_conservation": np.nan}

    counts = pd.Series(aa).value_counts()
    freqs = counts / n
    H = float(-(freqs * np.log2(freqs)).sum())
    out["shannon_entropy"] = H
    out["conservation"] = 1 - H / np.log2(20)                  # 1 = invariant
    out["consensus_aa"] = counts.index[0]
    out["consensus_freq"] = float(freqs.iloc[0])
    # how many residue types are needed to cover 90% of the kinome at this site
    out["n_aa_90pct"] = int(np.searchsorted(np.cumsum(freqs.values), 0.90) + 1)
    # conservation of physicochemical class (tolerance to conservative swaps)
    cls = pd.Series([AA_CLASS[a] for a in aa]).value_counts() / n
    out["class_conservation"] = float(cls.iloc[0])
    out["majority_class"] = cls.index[0]
    out["aa_freqs"] = {a: round(float(freqs.get(a, 0.0)), 4) for a in AAS
                       if freqs.get(a, 0.0) > 0}
    return out


def conservation_table(mat, label):
    recs = []
    for pos in range(1, 86):
        s = position_stats(mat[pos])
        s["KLIFS_position"] = pos
        s["region"] = REGION_OF_POS[pos - 1]
        s["set"] = label
        recs.append(s)
    return pd.DataFrame(recs)


def classify(df, col="conservation"):
    """Mutability call from kinome-wide conservation."""
    def _c(v):
        if v >= 0.70:  return "Immutable (invariant)"
        if v >= 0.50:  return "Strongly constrained"
        if v >= 0.30:  return "Moderately constrained"
        if v >= 0.15:  return "Tolerant"
        return "Hypervariable"
    return df[col].apply(_c)


# ----------------------------------------------------------------- main -----
def main():
    kinases = fetch_klifs()
    mat, info = build_matrix(kinases)
    print(f"{len(mat)} human kinases with complete 85-residue KLIFS pockets")

    # ---- kinome-wide
    kinome = conservation_table(mat, "Kinome (all)")
    kinome["mutability_class"] = classify(kinome)

    # ---- per group (groups with enough members to be meaningful)
    per_group = []
    for grp, idx in info.groupby("group").groups.items():
        sub = mat.loc[idx]
        if len(sub) < 10:
            continue
        t = conservation_table(sub, grp)
        t["group"] = grp
        per_group.append(t)
    per_group = pd.concat(per_group, ignore_index=True)

    # ---- CDK family specifically (the family the screen targets)
    cdk_idx = info.index[info["family"].str.startswith("CDK", na=False)]
    cdk_tab = conservation_table(mat.loc[cdk_idx], "CDK family")
    cdk_tab["group"] = "CDK family"

    # ---- residue numbers of the 9 screened CDKs at each KLIFS position
    screened = pd.read_csv(os.path.join(KLIFS,
                                        "KLIFS_compiled.csv"))
    screened["KLIFS_position"] = screened["KLIFS_idx"].str.split(":").str[1].astype(int)

    # ---- master table
    master = kinome.drop(columns=["set"]).copy()
    wide = per_group.pivot(index="KLIFS_position", columns="group",
                           values="conservation")
    wide.columns = [f"cons_{c}" for c in wide.columns]
    master = master.merge(wide, on="KLIFS_position", how="left")
    ct = cdk_tab.set_index("KLIFS_position")
    master["cons_CDKfamily"] = master["KLIFS_position"].map(ct["conservation"])
    master["consensus_CDKfamily"] = master["KLIFS_position"].map(ct["consensus_aa"])

    for cdk in ["CDK7", "CDK8", "CDK19", "CDK9", "CDK12", "CDK13", "CDK2", "CDK4", "CDK6"]:
        m = screened.set_index("KLIFS_position")
        master[f"{cdk}_residue"] = master["KLIFS_position"].map(
            lambda p: f"{m.loc[p, cdk+'_seq']}{int(m.loc[p, cdk])}"
            if p in m.index and pd.notna(m.loc[p, cdk]) else "")

    master.insert(1, "KLIFS_label", master["region"] + ":" + master["KLIFS_position"].astype(str))

    aa_freq_json = master.pop("aa_freqs")
    master.to_csv(os.path.join(OUT, "klifs_kinome_conservation.csv"), index=False)
    per_group.drop(columns=["aa_freqs"]).to_csv(
        os.path.join(OUT, "klifs_conservation_by_group.csv"), index=False)
    json.dump({int(r.KLIFS_position): f for r, f in zip(master.itertuples(), aa_freq_json)},
              open(os.path.join(OUT, "klifs_position_aa_frequencies.json"), "w"), indent=1)
    mat.assign(**info[["name", "group", "family", "uniprot"]]).to_csv(
        os.path.join(OUT, "klifs_alignment_matrix_human.csv"), index=False)

    print(master[["KLIFS_label", "consensus_aa", "consensus_freq",
                  "conservation", "mutability_class"]].to_string())
    return master, per_group, cdk_tab, mat, info


if __name__ == "__main__":
    main()


# ------------------------------------- size-matched group specificity -------
def _encode(mat):
    """(n_kinases, 85) int array; 0-19 = amino acid, -1 = gap/unknown."""
    lut = {a: i for i, a in enumerate(AAS)}
    return np.array([[lut.get(c, -1) for c in row] for row in mat.values], dtype=np.int8)


def _col_stats(enc_rows):
    """Conservation and AA probability matrix for a set of rows."""
    n_pos = enc_rows.shape[1]
    cons = np.zeros(n_pos)
    probs = np.zeros((n_pos, 20))
    for j in range(n_pos):
        col = enc_rows[:, j]
        col = col[col >= 0]
        if len(col) == 0:
            cons[j] = np.nan
            continue
        cnt = np.bincount(col, minlength=20).astype(float)
        p = cnt / cnt.sum()
        nz = p[p > 0]
        cons[j] = 1 - (-(nz * np.log2(nz)).sum()) / np.log2(20)
        probs[j] = p
    return cons, probs


def specificity_table(mat, info, groups, n_perm=1000, seed=0, pseudo=0.05):
    """
    Group-specific constraint, corrected for the fact that smaller and more
    closely related subsets always look more conserved.

    For each group of size n we draw `n_perm` random size-matched subsets of the
    kinome and use them as the null. Reported per position:
      excess_conservation : observed − null mean
      z_conservation      : (observed − null mean) / null sd
      KL_vs_kinome        : KL(group ‖ kinome) in bits (residue-preference shift)
      z_KL                : same, standardised against the size-matched null
    """
    rng = np.random.default_rng(seed)
    enc = _encode(mat)
    n_kin = enc.shape[0]
    kin_cons, kin_probs = _col_stats(enc)
    kp = (kin_probs + pseudo) / (kin_probs + pseudo).sum(axis=1, keepdims=True)

    def kl(probs):
        q = (probs + pseudo) / (probs + pseudo).sum(axis=1, keepdims=True)
        return (q * np.log2(q / kp)).sum(axis=1)

    recs = []
    for grp, rows in groups.items():
        idx = np.asarray(rows)
        n = len(idx)
        obs_cons, obs_probs = _col_stats(enc[idx])
        obs_kl = kl(obs_probs)

        null_cons = np.empty((n_perm, 85))
        null_kl = np.empty((n_perm, 85))
        for b in range(n_perm):
            samp = rng.choice(n_kin, size=n, replace=False)
            c, pr = _col_stats(enc[samp])
            null_cons[b] = c
            null_kl[b] = kl(pr)

        mc, sc = np.nanmean(null_cons, 0), np.nanstd(null_cons, 0)
        mk, sk = np.nanmean(null_kl, 0), np.nanstd(null_kl, 0)
        for j in range(85):
            recs.append({
                "group": grp, "KLIFS_position": j + 1, "region": REGION_OF_POS[j],
                "n_kinases": n,
                "conservation": obs_cons[j],
                "null_conservation": mc[j],
                "excess_conservation": obs_cons[j] - mc[j],
                "z_conservation": (obs_cons[j] - mc[j]) / sc[j] if sc[j] > 0 else np.nan,
                "KL_vs_kinome": obs_kl[j],
                "z_KL": (obs_kl[j] - mk[j]) / sk[j] if sk[j] > 0 else np.nan,
            })
    return pd.DataFrame(recs)
