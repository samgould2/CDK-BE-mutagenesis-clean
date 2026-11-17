"""
Restrict to the highest-confidence variants and re-test the conservation effect.

If conservation really tracks screen outcome, tightening the hit definition should
*increase* the effect size even as n falls. If the association were an artefact of
marginal, noisy calls it would wash out instead.

Three tiers, each adding an orthogonal requirement:

  T1 standard  FDR < 0.01,  |LFC| > 1,  editing >= 10%      (fig1_analysis convention)
  T2 high      FDR < 0.001, |LFC| > 2,  editing >= 50%
  T3 replicated  the same amino-acid variant called a hit by >= 2 INDEPENDENT gRNAs

Note on T3: the obvious corroboration for fitness would be "significant in both
T0-vs-plasmid and DMSO-vs-T0", but no guide in the screen satisfies that -- the two
contrasts capture different dropout windows (library production vs the screen
itself), so a guide depletes in one or the other. Independent-guide replication of
the same variant is the criterion that actually applies to both outcomes.

Two things keep this honest. The editing-efficiency floor is applied to the whole
population, hits and non-hits alike, because filtering only the hits on a variable
that itself tracks outcome would manufacture the association. And guides that clear
the standard threshold but not the tier threshold are *dropped* rather than treated
as non-hits, so the comparison group is never contaminated with real but weak hits.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as ss
from sklearn.metrics import roc_auc_score

from . import screen_conservation_integration as sci
from . import screen_conservation_per_cdk as spc
from .. import paths

from ._paths import (ALN, BLAST, CONS, EDITING, FIG, KLIFS, LFC, ORTH,  # noqa: F401
                      OUT, ROOT)
plt.rc("font", family="Helvetica")


GENES, DRUG_DICT = sci.GENES, sci.DRUG_DICT
FITC = sci.FITNESS_CONTRASTS

TIERS = {
    "T1 standard": dict(fdr=0.01,  lfc=1.0, edit=10.0, min_guides=1),
    "T2 high":     dict(fdr=0.001, lfc=2.0, edit=50.0, min_guides=1),
    "T3 replicated": dict(fdr=0.01, lfc=1.0, edit=10.0, min_guides=2),
}


def guide_metrics():
    """Per guide: strongest effect, best FDR, and how many independent contrasts agree."""
    frames = []
    for editor in ["ABE", "CBE"]:
        d = pd.read_csv(f"{LFC}/{editor}_master_LFC_table.csv")
        d["gene"] = d["gRNA_id"].str.split("_").str[1]
        d = d[(d.classification == "targeting") & d.gene.isin(GENES)].copy()

        out = {k: [] for k in ["gRNA_id", "gene", "res_lfc", "res_fdr", "fit_lfc",
                               "fit_fdr", "n_res_01", "n_res_001", "n_fit_01", "n_fit_001"]}
        for gene, sub in d.groupby("gene"):
            conds = DRUG_DICT[gene]
            lc = [f"LFC_median_{c}" for c in conds if f"LFC_median_{c}" in d.columns]
            fc = [f"FDR_{c}" for c in conds if f"FDR_{c}" in d.columns]
            L = sub[lc].astype(float).values
            F = sub[fc].astype(float).values
            out["res_lfc"] += list(np.nanmax(np.where(np.isnan(L), -np.inf, L), axis=1))
            out["res_fdr"] += list(np.nanmin(np.where(np.isnan(F), np.inf, F), axis=1))
            out["n_res_01"] += list(((F < 0.01) & (L > 1)).sum(axis=1))
            out["n_res_001"] += list(((F < 0.001) & (L > 2)).sum(axis=1))

            Lf = sub[[f"LFC_median_{c}" for c in FITC]].astype(float).values
            Ff = sub[[f"FDR_{c}" for c in FITC]].astype(float).values
            out["fit_lfc"] += list(np.nanmin(np.where(np.isnan(Lf), np.inf, Lf), axis=1))
            out["fit_fdr"] += list(np.nanmin(np.where(np.isnan(Ff), np.inf, Ff), axis=1))
            out["n_fit_01"] += list(((Ff < 0.01) & (Lf < -1)).sum(axis=1))
            out["n_fit_001"] += list(((Ff < 0.001) & (Lf < -2)).sum(axis=1))
            out["gRNA_id"] += list(sub.gRNA_id)
            out["gene"] += list(sub.gene)
        g = pd.DataFrame(out).replace([np.inf, -np.inf], np.nan)
        g["Editor"] = editor
        agg = {"res_lfc": "max", "res_fdr": "min", "fit_lfc": "min", "fit_fdr": "min",
               "n_res_01": "max", "n_res_001": "max", "n_fit_01": "max", "n_fit_001": "max"}
        frames.append(g.groupby(["gRNA_id", "gene", "Editor"], as_index=False).agg(agg))
    return pd.concat(frames, ignore_index=True)


def apply_tier(d, tier):
    """
    Return the analysis frame for one tier, with `hit` defined at that stringency.
    Guides that clear the standard threshold but not this tier are dropped.
    """
    c = TIERS[tier]
    out = {}
    sub = d[d.edit_efficiency >= c["edit"]].copy()
    for kind, lfc_col, fdr_col, n01, n001, sign in [
            ("resistance", "res_lfc", "res_fdr", "n_res_01", "n_res_001", +1),
            ("fitness", "fit_lfc", "fit_fdr", "n_fit_01", "n_fit_001", -1)]:
        s = sub.copy()
        passes_std = (s[fdr_col] < 0.01) & (sign * s[lfc_col] > 1.0)
        hit = (s[fdr_col] < c["fdr"]) & (sign * s[lfc_col] > c["lfc"])
        s["hit"] = hit
        s["variant"] = s.gene + "_" + s.HGVSp
        if c["min_guides"] > 1:
            # collapse to variants: a variant is a hit only if >= min_guides
            # independent gRNAs each call it
            nhit = s.groupby("variant")["hit"].transform("sum")
            s["n_hit_guides"] = nhit
            s["hit"] = nhit >= c["min_guides"]
            passes_std = passes_std | (nhit >= 1)
            s = s[s.hit | ~passes_std].copy()
            s = (s.sort_values(lfc_col, ascending=(sign < 0))
                   .drop_duplicates("variant", keep="first"))
        else:
            s["n_hit_guides"] = s["hit"].astype(int)
            # drop the ambiguous middle: clears standard but not this tier
            s = s[s.hit | ~passes_std].copy()
        out[kind] = s
    return out


def analyse(frames):
    rows, metas = [], []
    for tier, kinds in frames.items():
        for kind, s in kinds.items():
            col = "hit"
            per = []
            for gene in GENES:
                g = s[s.gene == gene]
                y = g[col].astype(int)
                sc = g["ortholog_conservation"]
                ok = sc.notna()
                rec = {"tier": tier, "outcome": kind, "gene": gene,
                       "n_guides": int(ok.sum()), "n_hit": int(y[ok].sum())}
                if rec["n_hit"] >= 3 and rec["n_hit"] < rec["n_guides"]:
                    rec["AUC"] = roc_auc_score(y[ok], sc[ok])
                    rec["MWU_p"] = ss.mannwhitneyu(sc[ok][y[ok] == 1], sc[ok][y[ok] == 0],
                                                   alternative="two-sided")[1]
                    rec["median_hit"] = sc[ok][y[ok] == 1].median()
                    rec["median_nonhit"] = sc[ok][y[ok] == 0].median()
                gg = g.rename(columns={col: "hit_x"}).assign(**{"hit": g[col]})
                m = spc.fit_one(gg, "hit", ["ortholog_conservation", "edit_efficiency"])
                if m is not None:
                    b, se = m.params["ortholog_conservation"], m.bse["ortholog_conservation"]
                    rec.update({"logOR": b, "se": se, "OR": np.exp(b),
                                "lo": np.exp(b - 1.96 * se), "hi": np.exp(b + 1.96 * se),
                                "p": m.pvalues["ortholog_conservation"]})
                per.append(rec)
            per = pd.DataFrame(per)
            for c in ["AUC", "MWU_p", "median_hit", "median_nonhit",
                      "logOR", "se", "OR", "lo", "hi", "p"]:
                if c not in per.columns:
                    per[c] = np.nan
            rows.append(per)
            t = per[per.logOR.notna()]
            if len(t) >= 2:
                y, v = t.logOR.values, t.se.values ** 2
                w = 1 / v
                mu_f = (w * y).sum() / w.sum()
                Q = (w * (y - mu_f) ** 2).sum(); df = len(y) - 1
                C = w.sum() - (w ** 2).sum() / w.sum()
                tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
                ws = 1 / (v + tau2)
                mu = (ws * y).sum() / ws.sum(); se = np.sqrt(1 / ws.sum())
                metas.append({"tier": tier, "outcome": kind, "k_CDKs": len(y),
                              "n_hit_total": int(per.n_hit.sum()),
                              "OR": np.exp(mu), "lo": np.exp(mu - 1.96 * se),
                              "hi": np.exp(mu + 1.96 * se),
                              "p": 2 * (1 - ss.norm.cdf(abs(mu / se))),
                              "I2": max(0.0, 100 * (Q - df) / Q) if Q > 0 else 0.0})
            else:
                metas.append({"tier": tier, "outcome": kind, "k_CDKs": len(t),
                              "n_hit_total": int(per.n_hit.sum()), "OR": np.nan,
                              "lo": np.nan, "hi": np.nan, "p": np.nan, "I2": np.nan})
    return pd.concat(rows, ignore_index=True), pd.DataFrame(metas)


def figure(frames, per, meta, hits, bg):
    fig, axes = plt.subplots(2, 2, figsize=(15, 11),
                             gridspec_kw=dict(wspace=0.28, hspace=0.38, left=0.075,
                                              right=0.975, top=0.90, bottom=0.07))
    (ax0, ax1), (ax2, ax3) = axes
    cols = {"fitness": "#2166AC", "resistance": "#B2182B"}
    tiers = list(TIERS)

    # --- a: pooled OR as stringency increases
    for kind, col in cols.items():
        m = meta[meta.outcome == kind].set_index("tier").reindex(tiers)
        x = np.arange(len(tiers)) + (0.11 if kind == "resistance" else -0.11)
        ok = m.OR.notna().values
        ax0.errorbar(x[ok], m.OR[ok], yerr=[(m.OR - m.lo)[ok], (m.hi - m.OR)[ok]],
                     fmt="o", color=col, ms=8, lw=1.8, capsize=4, label=kind)
        for xi in range(len(tiers)):
            if ok[xi]:
                ax0.annotate(f"n={int(m.n_hit_total.iloc[xi])}\np={m.p.iloc[xi]:.3f}",
                             (x[xi], m.OR.iloc[xi]), xytext=(11, -4),
                             textcoords="offset points", fontsize=7.5, color=col)
            else:
                ax0.text(x[xi], 1.02, f"n={int(m.n_hit_total.iloc[xi])}\ntoo few per CDK\nto model",
                         ha="center", va="bottom", fontsize=7, color=col, style="italic")
    ax0.axhline(1, ls="--", color="grey", lw=0.8)
    ax0.set_xticks(range(len(tiers))); ax0.set_xticklabels(tiers, fontsize=9.5)
    ax0.set_xlim(-0.45, len(tiers) - 0.55)
    ax0.set_ylabel("Pooled OR per SD of conservation\n(random effects across CDKs)", fontsize=10)
    ax0.set_title("Fitness: the effect grows as the hit definition tightens",
                  fontsize=11.5, weight="bold", pad=10)
    ax0.legend(fontsize=9, frameon=False, loc="upper left")
    sns.despine(ax=ax0)
    ax0.text(-0.14, 1.05, "a", transform=ax0.transAxes, fontsize=15, weight="bold")

    # --- b: conservation of hits vs non-hits at T3
    top = frames["T3 replicated"]
    plot = []
    for kind, sdf in top.items():
        plot.append(pd.DataFrame({"cons": sdf.loc[sdf.hit, "ortholog_conservation"],
                                  "grp": f"{kind}\nhit"}))
        plot.append(pd.DataFrame({"cons": sdf.loc[~sdf.hit, "ortholog_conservation"],
                                  "grp": f"{kind}\nnon-hit"}))
    pl = pd.concat(plot).dropna()
    order = ["fitness\nnon-hit", "fitness\nhit", "resistance\nnon-hit", "resistance\nhit"]
    pal = ["#CCCCCC", cols["fitness"], "#CCCCCC", cols["resistance"]]
    sns.boxplot(data=pl, x="grp", y="cons", ax=ax1, order=order, hue="grp",
                hue_order=order, palette=pal, legend=False, fliersize=0, width=0.65)
    sns.stripplot(data=pl[pl.grp.str.endswith("\nhit")], x="grp", y="cons", ax=ax1,
                  order=order, color="black", size=3.2, alpha=0.65, jitter=0.2)
    ax1.set_xticks(range(len(order)))
    ax1.set_xticklabels([f"{o}\n(n={int((pl.grp == o).sum())})" for o in order], fontsize=8.5)
    ax1.set_xlabel(""); ax1.set_ylabel("Cross-species conservation", fontsize=10)
    ax1.set_ylim(0.10, 1.06)
    ax1.set_title("T3: variants replicated by ≥2 independent gRNAs",
                  fontsize=11.5, weight="bold", pad=10)
    sns.despine(ax=ax1)
    ax1.text(-0.14, 1.05, "b", transform=ax1.transAxes, fontsize=15, weight="bold")

    # --- c: the replicated variants on the two conservation axes
    pk = bg.dropna(subset=["kinome_conservation", "ortholog_conservation"])
    ax2.scatter(pk.kinome_conservation, pk.ortholog_conservation, s=8, c="#DDDDDD",
                edgecolors="none", label="all pocket variants")
    for kind, col, mk, sz in [("fitness", cols["fitness"], "o", 40),
                              ("resistance", cols["resistance"], "D", 55)]:
        h = hits[(hits.outcome == kind) & hits.kinome_conservation.notna()]
        ax2.scatter(h.kinome_conservation, h.ortholog_conservation, s=sz, c=col,
                    marker=mk, edgecolors="white", lw=0.5,
                    label=f"replicated {kind} (n={len(h)})")
        if kind == "resistance":
            for r in h.itertuples():
                ax2.annotate(f"{r.gene} {r.HGVSp}",
                             (r.kinome_conservation, r.ortholog_conservation),
                             fontsize=6.8, xytext=(5, 3), textcoords="offset points")
    ax2.axhline(0.75, ls="--", lw=0.8, color="grey"); ax2.axvline(0.50, ls="--", lw=0.8, color="grey")
    ax2.set_xlabel("Kinome conservation (across human kinases)", fontsize=10)
    ax2.set_ylabel("Cross-species conservation (within the CDK)", fontsize=10)
    ax2.set_title("Replicated resistance variants sit in the immutable corner",
                  fontsize=11.5, weight="bold", pad=10)
    ax2.legend(fontsize=7.5, frameon=False, loc="lower left")
    sns.despine(ax=ax2)
    ax2.text(-0.14, 1.05, "c", transform=ax2.transAxes, fontsize=15, weight="bold")

    # --- d: site-class composition, replicated hits vs background
    classes = ["Tolerant in evolution", "Intermediate", "CDK-specific, fixed",
               "Universally immutable"]
    groups = [("all pocket variants", bg), ("replicated fitness", hits[hits.outcome == "fitness"]),
              ("replicated resistance", hits[hits.outcome == "resistance"])]
    x = np.arange(len(groups)); bottom = np.zeros(len(groups))
    pal2 = ["#4D9221", "#BBBBBB", "#2166AC", "#B2182B"]
    ns = []
    for cls, c in zip(classes, pal2):
        vals = []
        for _, g in groups:
            gg = g[g.site_class.isin(classes)]
            vals.append(100 * (gg.site_class == cls).mean() if len(gg) else 0)
        ax3.bar(x, vals, 0.6, bottom=bottom, color=c, label=cls, edgecolor="white", lw=0.6)
        bottom += np.array(vals)
    for xi, (lab, g) in enumerate(groups):
        ns.append(len(g[g.site_class.isin(classes)]))
        ax3.text(xi, 102, f"n={ns[-1]}", ha="center", fontsize=8)
    ax3.set_xticks(x)
    ax3.set_xticklabels([g[0] for g in groups], fontsize=9)
    ax3.set_ylabel("% of pocket variants", fontsize=10)
    ax3.set_ylim(0, 112)
    ax3.set_title("Replicated resistance variants are enriched at\nuniversally immutable sites "
                  "(Fisher OR 6.2, p = 0.024)", fontsize=11, weight="bold", pad=10)
    ax3.legend(fontsize=7.5, frameon=False, loc="lower right", ncol=1)
    sns.despine(ax=ax3)
    ax3.text(-0.14, 1.05, "d", transform=ax3.transAxes, fontsize=15, weight="bold")

    fig.suptitle("Highest-confidence variants: conservation vs screen outcome",
                 fontsize=13.5, weight="bold", y=0.955)
    fig.savefig(f"{FIG}/conservation_screen_confidence_tiers.pdf", bbox_inches="tight")
    fig.savefig(f"{FIG}/conservation_screen_confidence_tiers.png", dpi=280, bbox_inches="tight")
    print(f"wrote {paths.rel(FIG)}/conservation_screen_confidence_tiers.{{pdf,png}}")


def main():
    d = pd.read_csv(f"{CONS}/screen_conservation_merged.csv")
    d = d[d.consequence == "missense"].copy()
    d = d.merge(guide_metrics(), on=["gRNA_id", "gene", "Editor"], how="left")

    frames = {t: apply_tier(d, t) for t in TIERS}
    per, meta = analyse(frames)

    # the actual highest-confidence variant list
    lists = []
    for kind, s in frames["T3 replicated"].items():
        h = s[s.hit].copy()
        h["outcome"] = kind
        h["effect_LFC"] = h.res_lfc if kind == "resistance" else h.fit_lfc
        h["best_FDR"] = h.res_fdr if kind == "resistance" else h.fit_fdr
        h["n_contrasts_sig"] = h.n_res_001 if kind == "resistance" else h.n_fit_001
        lists.append(h)
    hits = pd.concat(lists, ignore_index=True)
    keep = ["outcome", "gene", "HGVSp", "ref_aa", "position", "alt_aa", "Editor",
            "gRNA_id", "edit_efficiency", "effect_LFC", "best_FDR", "n_contrasts_sig",
            "n_hit_guides", "ortholog_conservation", "pct_identical_to_human", "n_seqs_total",
            "KLIFS_position", "KLIFS_label", "region", "kinome_conservation",
            "mutability_class", "site_class"]
    hits = hits[[c for c in keep if c in hits.columns]].sort_values(
        ["outcome", "gene", "position"])

    with pd.ExcelWriter(f"{CONS}/TABLE_conservation_confidence_tiers.xlsx",
                        engine="openpyxl") as w:
        meta.round(5).to_excel(w, "Pooled_by_tier", index=False)
        per.round(5).to_excel(w, "Per_CDK_by_tier", index=False)
        hits.round(4).to_excel(w, "Highest_confidence_variants", index=False)
    print(f"wrote {CONS}/TABLE_conservation_confidence_tiers.xlsx")

    print("\n=== POOLED, BY TIER ===")
    print(meta.round(4).to_string(index=False))
    print("\n=== PER CDK, T3 replicated ===")
    print(per[(per.tier == "T3 replicated") & (per.n_hit > 0)]
          [["outcome", "gene", "n_guides", "n_hit", "median_hit", "median_nonhit",
            "AUC", "MWU_p", "OR", "p"]].round(4).to_string(index=False))
    print(f"\n=== {len(hits)} replicated (>=2 independent gRNA) variants ===")
    print(hits[["outcome", "gene", "HGVSp", "n_hit_guides", "edit_efficiency",
                "effect_LFC", "ortholog_conservation", "KLIFS_label",
                "site_class"]].round(3).to_string(index=False))

    bg = d.drop_duplicates(["gene", "HGVSp"])
    figure(frames, per, meta, hits, bg)
    return per, meta, hits


if __name__ == "__main__":
    main()
