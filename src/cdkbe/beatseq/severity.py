"""Severity map and the pre-registered tests.

Everything here is a WITHIN-TARGETING contrast or is calibrated against a
pseudo-guide null built from NT cells, so none of it depends on the leave-one-out
artifact in the delivered pseudobulk.

Tests, in the order the plan specifies them:
  1. LOF (nonsense + start-loss) vs missense on the cross-fitted statistic
  2. severity calibrated to a pseudo-guide null -> the severity ladder
  3. severity x tier interaction against the bulk untreated-fitness axis
     (pre-registered: rho > 0 for untreated-fitness alleles, ~0 for drug-conditional)
  4. KLIFS kinase-pocket residues vs the rest
  5. CDK9 L156F, one-sided, plus the four shared-protospacer ABE/CBE contrasts
  6. cross-cell-line replication of the effect vector
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from . import norm, nullcal, paths  # noqa: E402

ARMS = paths.ARMS
CL_OF = {a: a.split("/")[0] for a in ARMS}


def build_null(genes, m: pd.DataFrame, n_draws: int = 1000):
    """Pseudo-guide null for D, drawn at sizes that SPAN each arm's guide cell counts.

    CORRECTED. The previous version hardcoded [40,80,150,300,500] / [25,40,70,110,160] and
    justified a single pooled null per arm with "D is n-unbiased". D is unbiased in its
    MEAN, which is what the shipped check (spearmanr(n_cells, D)) tested. Its sampling SD is
    strongly n-dependent -- about eleven-fold across the sizes drawn -- so the null must be
    size-matched, and the sizes must cover the guides being scored. Under the old sizes,
    54%, 69%, 24% and 0% of guides in the four arms lay outside the drawn range.
    """
    out = []
    for arm in ARMS:
        sizes = nullcal.sizes_for_guides(m.loc[m.arm == arm, "n_cells"])
        print(f"  {arm:11s} null sizes {sizes} x {n_draws}", flush=True)
        out.append(norm.pseudo_guide_null(CL_OF[arm], arm, genes, sizes, n_draws=n_draws))
    return pd.concat(out, ignore_index=True)


def add_severity(m: pd.DataFrame, null: pd.DataFrame) -> pd.DataFrame:
    """Severity = D standardised against the null AT THAT GUIDE'S CELL COUNT.

    The superseded pooled calibration is retained as severity_pooled / p_pooled / q_pooled
    so the two can be compared; `severity`, `p` and `q` are the corrected values.
    """
    m = m.copy()
    # --- superseded pooled calibration, kept for comparison ------------------
    stats = null.groupby("arm")["D"].agg(["mean", "std"])
    m["null_mean_pooled"] = m.arm.map(stats["mean"])
    m["null_sd_pooled"] = m.arm.map(stats["std"])
    m["severity_pooled"] = (m.D - m.null_mean_pooled) / m.null_sd_pooled
    m["p_pooled"] = [
        float((null.loc[null.arm == r.arm, "D"].to_numpy() >= r.D).mean())
        for _, r in m.iterrows()]
    m["q_pooled"] = np.nan
    for arm, idx in m.groupby("arm").groups.items():
        m.loc[idx, "q_pooled"] = nullcal.bh(m.loc[idx, "p_pooled"].to_numpy())

    # --- corrected size-matched calibration ---------------------------------
    models, std_draws = nullcal.fit_size_model(null)
    m = nullcal.standardise(m, models, std_draws, d_col="D", prefix="")
    m["p_emp"] = m["p"]                     # kept for backward compatibility
    m.attrs["models"] = models
    return m


def report(m: pd.DataFrame, null: pd.DataFrame) -> None:
    tgt = m[~m.is_control]
    print("\n" + "=" * 78)
    print("1. LOF vs MISSENSE on the cross-fitted statistic")
    print("=" * 78)
    for arm, g in tgt.groupby("arm"):
        lof, mis = g[g.is_LOF].severity, g[g.variant_class == "missense"].severity
        if len(lof) < 3:
            continue
        p = mannwhitneyu(lof, mis, alternative="greater")[1]
        print(f"  {arm:11s} LOF n={len(lof):2d} med={lof.median():+.2f} | "
              f"missense n={len(mis):3d} med={mis.median():+.2f}   p={p:.4f}")
    lof, mis = tgt[tgt.is_LOF], tgt[tgt.variant_class == "missense"]
    # rank within arm so arms with different scales can be pooled
    tgt = tgt.copy()
    tgt["r"] = tgt.groupby("arm")["severity"].rank(pct=True)
    p = mannwhitneyu(tgt[tgt.is_LOF].r, tgt[tgt.variant_class == "missense"].r,
                     alternative="greater")[1]
    print(f"  POOLED (within-arm ranks): LOF n={len(lof)} med={tgt[tgt.is_LOF].r.median():.3f} | "
          f"missense n={len(mis)} med={tgt[tgt.variant_class=='missense'].r.median():.3f}  p={p:.4f}")

    print("\n" + "=" * 78)
    print("2. CONTROL CALIBRATION — do the 9 real NT guides sit inside the null?")
    print("=" * 78)
    for arm, g in m[m.is_control].groupby("arm"):
        nd = null[null.arm == arm]["D"]
        print(f"  {arm:11s} NT guides n={len(g)} median severity={g.severity.median():+.2f} "
              f"(range {g.severity.min():+.2f} to {g.severity.max():+.2f})   "
              f"null sd={nd.std():.2e}")

    print("\n" + "=" * 78)
    print("3. SEVERITY x TIER — the pre-registered interaction")
    print("=" * 78)
    print("   prediction: rho > 0 for untreated-fitness alleles, ~0 for drug-conditional")
    for tier, g in tgt.groupby("tier"):
        g = g.dropna(subset=["LFC_median_T0_Plasmid", "severity"])
        if len(g) < 10:
            continue
        r, p = spearmanr(g.LFC_median_T0_Plasmid.abs(), g.severity)
        rs, ps = spearmanr(g.LFC_median_T0_Plasmid, g.severity)
        print(f"  {tier:20s} n={len(g):3d}  |LFC_T0| vs severity rho={r:+.3f} p={p:.3f}   "
              f"signed rho={rs:+.3f} p={ps:.3f}")

    print("\n" + "=" * 78)
    print("4. KLIFS kinase-pocket residues")
    print("=" * 78)
    k = tgt.dropna(subset=["KLIFS_residue"])
    for arm, g in k.groupby("arm"):
        a, b = g[g.KLIFS_residue == 1].severity, g[g.KLIFS_residue == 0].severity
        if len(a) < 5:
            continue
        p = mannwhitneyu(a, b, alternative="greater")[1]
        print(f"  {arm:11s} pocket n={len(a):2d} med={a.median():+.2f} | "
              f"non-pocket n={len(b):2d} med={b.median():+.2f}  p={p:.3f}")
    kk = k.copy(); kk["r"] = kk.groupby("arm")["severity"].rank(pct=True)
    p = mannwhitneyu(kk[kk.KLIFS_residue == 1].r, kk[kk.KLIFS_residue == 0].r,
                     alternative="greater")[1]
    print(f"  POOLED (within-arm ranks) p={p:.4f}")

    print("\n" + "=" * 78)
    print("5. CDK9 L156F and the shared-protospacer contrasts")
    print("=" * 78)
    print("   Same protospacer, same integration biology, different installed allele --")
    print("   any difference cannot be a guide or off-target effect.")
    for gid, alleles in paths.SHARED_PROTOSPACER_PAIRS.items():
        rows = tgt[tgt.gRNA_id == gid]
        if rows.empty:
            continue
        print(f"\n  {gid}   ABE {alleles['ABE']}  vs  CBE {alleles['CBE']}")
        for _, r in rows.sort_values("arm").iterrows():
            print(f"    {r.arm:11s} {str(r.Mutation):16s} n={int(r.n_cells):4d} "
                  f"severity={r.severity:+6.2f}  p_emp={r.p_emp:.3f}  "
                  f"LFC_T0={r.LFC_median_T0_Plasmid:+.2f}  edit={r.any_edit:.0f}%")

    print("\n" + "=" * 78)
    print("6. CROSS-CELL-LINE REPLICATION of the effect vector")
    print("=" * 78)
    for ed in ["ABE", "CBE"]:
        pa = norm.load_profile(f"A549/{ed}")
        pt = norm.load_profile(f"T47D/{ed}")
        shared = [g for g in pa.index if g in set(pt.index)]
        ctrl_ids = set(m[m.is_control].guide_identity)
        tg = [g for g in shared if g not in ctrl_ids]
        A = pa.loc[tg].to_numpy(); B = pt.loc[tg].to_numpy()
        A = (A - A.mean(0)) / (A.std(0) + 1e-9)
        B = (B - B.mean(0)) / (B.std(0) + 1e-9)
        C = np.corrcoef(A, B)[:len(tg), len(tg):]
        diag = np.diag(C)
        off = C[~np.eye(len(tg), dtype=bool)]
        rng = np.random.default_rng(0)
        perm = np.array([np.diag(C[:, rng.permutation(len(tg))]).mean() for _ in range(2000)])
        p = float((perm >= diag.mean()).mean())
        print(f"  {ed}: {len(tg)} shared targeting guides   same-guide median r={np.median(diag):+.4f} | "
              f"different-guide median r={np.median(off):+.4f}")
        print(f"      permutation p={p:.4f}  (z={(diag.mean()-perm.mean())/perm.std():+.2f})")
        top = pd.Series(diag, index=tg).sort_values(ascending=False).head(8)
        print("      best replicating: " + ", ".join(
            f"{g.split('|')[0]} {v:.3f}" for g, v in top.items()))


if __name__ == "__main__":
    paths.ensure_dirs()
    genes = np.load(paths.RESULTS / "gene_universe.npy", allow_pickle=True)
    m = pd.read_csv(paths.RESULTS / "master_guide_table.csv")
    null_path = paths.RESULTS / "pseudo_guide_null.csv"
    rebuild = "--rebuild-null" in sys.argv or not null_path.exists()
    if rebuild:
        null = build_null(genes, m)
        null.to_csv(null_path, index=False)
    else:
        null = pd.read_csv(null_path)
    print(f"pseudo-guide null: {len(null)} draws across {null.arm.nunique()} arms")

    # D is n-free in its MEAN (this was the only check originally run) ...
    for arm, g in null.groupby("arm"):
        r = spearmanr(g.n_cells, g.D)[0]
        print(f"  {arm:11s} null rho(D, n)={r:+.3f}  mean={g.D.mean():+.2e} sd={g.D.std():.2e}")
    # ... but its SD is NOT, which is why the null must be size-matched.
    m = add_severity(m, null)
    rep = nullcal.calibration_report(m, null, m.attrs["models"], prefix="",
                                     old_sev="severity_pooled", old_p="p_pooled")
    rep.to_csv(paths.RESULTS / "severity_calibration.csv", index=False)
    m.to_csv(paths.RESULTS / "severity.csv", index=False)

    t = m[~m.is_control]
    print(f"\n  HIT COUNTS   q<0.05: {int((t.q_pooled < 0.05).sum())} (pooled) -> "
          f"{int((t.q < 0.05).sum())} (size-matched)")
    print(f"               q<0.10: {int((t.q_pooled < 0.10).sum())} -> {int((t.q < 0.10).sum())}")
    print(f"               raw p<0.05: {int((t.p_pooled < 0.05).sum())} -> "
          f"{int((t.p < 0.05).sum())}   (chance {0.05*len(t):.1f})")
    report(m, null)
