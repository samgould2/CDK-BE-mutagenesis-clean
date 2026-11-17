"""Statistics + figures: does conservation predict resistance / fitness in the screen?"""
import os, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as ss
from statsmodels.discrete.discrete_model import Logit
from sklearn.metrics import roc_curve, roc_auc_score

from . import screen_conservation_integration as sci
from .. import paths

from ._paths import (ALN, BLAST, CONS, EDITING, FIG, KLIFS, LFC, ORTH,  # noqa: F401
                      OUT, ROOT)
plt.rc("font", family="Helvetica")


OUTCOMES = {"fitness_hit": ("Fitness (depleted, no drug)", "#2166AC"),
            "resistance_hit": ("Resistance (enriched under inhibitor)", "#B2182B")}


def wilson(k, n, z=1.96):
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return p, max(c - h, 0), min(c + h, 1)


# ------------------------------------------------------------- statistics ---
def univariate(d):
    rows = []
    for outcome, (label, _) in OUTCOMES.items():
        for var, vname in [("ortholog_conservation", "Cross-species conservation"),
                           ("kinome_conservation", "Kinome conservation"),
                           ("edit_efficiency", "gRNA editing efficiency")]:
            s = d[[outcome, var]].dropna()
            a = s.loc[s[outcome], var]
            b = s.loc[~s[outcome].astype(bool), var]
            if len(a) < 5 or len(b) < 5:
                continue
            u, p = ss.mannwhitneyu(a, b, alternative="two-sided")
            auc = roc_auc_score(s[outcome].astype(int), s[var])
            rows.append({"outcome": label, "variable": vname,
                         "n_hit": len(a), "n_nonhit": len(b),
                         "median_hit": a.median(), "median_nonhit": b.median(),
                         "AUC": auc, "MWU_p": p})
    return pd.DataFrame(rows)


def logistic(d, outcome, predictors, cluster="gene"):
    """Logistic regression with gene-clustered standard errors; ORs per SD."""
    s = d[[outcome, cluster] + predictors].dropna().copy()
    X = s[predictors].astype(float).copy()
    scales = {}
    for c in predictors:
        if X[c].nunique() > 2:
            mu, sd = X[c].mean(), X[c].std()
            X[c] = (X[c] - mu) / sd
            scales[c] = "per SD"
        else:
            scales[c] = "yes vs no"
    X.insert(0, "const", 1.0)
    y = s[outcome].astype(int)
    # gene-clustered SEs need >= 2 genes; fall back to heteroscedasticity-robust
    if s[cluster].nunique() >= 2:
        m = Logit(y, X).fit(disp=0, cov_type="cluster",
                            cov_kwds={"groups": s[cluster]})
        cov_used = f"cluster({cluster})"
    else:
        m = Logit(y, X).fit(disp=0, cov_type="HC1")
        cov_used = "HC1"
    out = pd.DataFrame({"term": X.columns, "coef": m.params.values,
                        "se": m.bse.values, "p": m.pvalues.values})
    out["OR"] = np.exp(out.coef)
    out["lo"] = np.exp(out.coef - 1.96 * out.se)
    out["hi"] = np.exp(out.coef + 1.96 * out.se)
    out["scale"] = out.term.map(scales).fillna("")
    out["outcome"] = outcome
    out["cov_type"] = cov_used
    out["n"] = len(s)
    out["n_hit"] = int(y.sum())
    return out[out.term != "const"], m


def hit_rate_table(d, outcome, var="ortholog_conservation", nbins=5):
    s = d[[outcome, var, "edit_eff_tertile"]].dropna().copy()
    s["bin"] = pd.qcut(s[var], nbins, duplicates="drop")
    g = (s.groupby(["bin", "edit_eff_tertile"], observed=True)[outcome]
           .agg(["mean", "sum", "size"]).reset_index()
           .rename(columns={"mean": "hit_rate", "sum": "n_hit", "size": "n"}))
    g["bin_mid"] = g["bin"].apply(lambda b: b.mid)
    g["outcome"] = outcome
    return g


# ---------------------------------------------------------------- figures ---
def figure5(d, uni, fits, rates):
    fig = plt.figure(figsize=(16.5, 10.5))
    gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.30,
                          left=0.06, right=0.98, top=0.91, bottom=0.08)

    # a/b: conservation by outcome
    for i, (outcome, (label, col)) in enumerate(OUTCOMES.items()):
        ax = fig.add_subplot(gs[0, i])
        s = d[[outcome, "ortholog_conservation"]].dropna()
        s["grp"] = np.where(s[outcome], "hit", "not a hit")
        sns.violinplot(data=s, x="grp", y="ortholog_conservation", ax=ax,
                       order=["not a hit", "hit"], hue="grp",
                       hue_order=["not a hit", "hit"], palette=["#CCCCCC", col],
                       legend=False, cut=0, inner="quartile", linewidth=1)
        r = uni[(uni.outcome == label) & (uni.variable == "Cross-species conservation")]
        if len(r):
            p = r.MWU_p.iloc[0]
            ax.set_title(f"{label}\nMWU p = {p:.2e}   AUC = {r.AUC.iloc[0]:.3f}",
                         fontsize=10.5, weight="bold")
        n_hit = int(s[outcome].sum())
        ax.set_xticklabels([f"not a hit\n(n={len(s)-n_hit})", f"hit\n(n={n_hit})"])
        ax.set_xlabel(""); ax.set_ylabel("Cross-species conservation", fontsize=10)
        sns.despine(ax=ax)
        ax.text(-0.14, 1.06, "ab"[i], transform=ax.transAxes, fontsize=15, weight="bold")

    # c: ROC
    axc = fig.add_subplot(gs[0, 2])
    for outcome, (label, col) in OUTCOMES.items():
        for var, ls in [("ortholog_conservation", "-"), ("edit_efficiency", ":")]:
            s = d[[outcome, var]].dropna()
            if s[outcome].sum() < 5:
                continue
            fpr, tpr, _ = roc_curve(s[outcome].astype(int), s[var])
            auc = roc_auc_score(s[outcome].astype(int), s[var])
            nm = "conservation" if var == "ortholog_conservation" else "editing efficiency"
            axc.plot(fpr, tpr, ls, color=col, lw=1.8,
                     label=f"{label.split(' (')[0]} — {nm} ({auc:.2f})")
    axc.plot([0, 1], [0, 1], "--", color="grey", lw=0.8)
    axc.set_xlabel("False positive rate", fontsize=10)
    axc.set_ylabel("True positive rate", fontsize=10)
    axc.set_title("Conservation vs editing efficiency\nas a predictor", fontsize=10.5, weight="bold")
    axc.legend(fontsize=7.2, frameon=False, loc="lower right")
    sns.despine(ax=axc)
    axc.text(-0.14, 1.06, "c", transform=axc.transAxes, fontsize=15, weight="bold")

    # d/e: hit rate vs conservation, stratified by editing efficiency
    for i, (outcome, (label, col)) in enumerate(OUTCOMES.items()):
        ax = fig.add_subplot(gs[1, i])
        g = rates[rates.outcome == outcome]
        for tert, mk in zip(["low", "mid", "high"], ["o", "s", "^"]):
            gg = g[g.edit_eff_tertile == tert].sort_values("bin_mid")
            ax.plot(gg.bin_mid, 100 * gg.hit_rate, marker=mk, color=col,
                    alpha={"low": 0.35, "mid": 0.65, "high": 1.0}[tert],
                    lw=1.6, ms=5, label=f"{tert} editing")
        ax.set_xlabel("Cross-species conservation (quintile midpoint)", fontsize=9.5)
        ax.set_ylabel(f"{label.split(' (')[0]} hit rate (%)", fontsize=10)
        ax.set_title(f"{label.split(' (')[0]} hit rate vs conservation,\n"
                     "stratified by gRNA editing efficiency", fontsize=10.5, weight="bold")
        ax.legend(fontsize=7.5, frameon=False)
        sns.despine(ax=ax)
        ax.text(-0.14, 1.06, "de"[i], transform=ax.transAxes, fontsize=15, weight="bold")

    # f: odds-ratio forest
    axf = fig.add_subplot(gs[1, 2])
    lbl = {"ortholog_conservation": "Cross-species conservation",
           "kinome_conservation": "Kinome conservation",
           "edit_efficiency": "gRNA editing efficiency",
           "is_CBE": "CBE (vs ABE)", "in_pocket": "In KLIFS pocket"}
    y, ticks, tickpos = 0.0, [], []
    for outcome, (label, col) in OUTCOMES.items():
        f = fits[fits.outcome == outcome]
        for r in f.itertuples():
            axf.plot([r.lo, r.hi], [y, y], color=col, lw=1.8)
            axf.plot(r.OR, y, "o", color=col, ms=6)
            ticks.append(lbl.get(r.term, r.term)); tickpos.append(y)
            if r.p < 0.05:
                axf.text(r.hi * 1.06, y, "*" * (1 + (r.p < 0.01) + (r.p < 0.001)),
                         va="center", fontsize=9, color=col)
            y += 1
        y += 0.7
    axf.axvline(1, ls="--", color="grey", lw=0.8)
    axf.set_yticks(tickpos); axf.set_yticklabels(ticks, fontsize=8.5)
    axf.invert_yaxis(); axf.set_xscale("log")
    axf.set_xlabel("Odds ratio (95% CI, gene-clustered)", fontsize=10)
    axf.set_title("Multivariable logistic regression\n"
                  "blue = fitness, red = resistance", fontsize=10.5, weight="bold")
    sns.despine(ax=axf)
    axf.text(-0.14, 1.06, "f", transform=axf.transAxes, fontsize=15, weight="bold")

    fig.suptitle("Evolutionary conservation as a predictor of base-editing screen outcome",
                 fontsize=14, weight="bold", y=0.975)
    fig.savefig(f"{FIG}/conservation_vs_screen_outcome.pdf", bbox_inches="tight")
    fig.savefig(f"{FIG}/conservation_vs_screen_outcome.png", dpi=280, bbox_inches="tight")
    print(f"wrote {paths.rel(FIG)}/conservation_vs_screen_outcome.{{pdf,png}}")


def figure6(d, d_all):
    """Two-axis view. `d` is the missense set; `d_all` keeps nonsense for the control."""
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2),
                             gridspec_kw=dict(wspace=0.30, left=0.06, right=0.98,
                                              top=0.85, bottom=0.14))
    ax0, ax1, ax2 = axes

    p = d.dropna(subset=["ortholog_conservation", "kinome_conservation"])
    ax0.scatter(p.kinome_conservation, p.ortholog_conservation, s=9, c="#DDDDDD",
                label="all pocket edits", edgecolors="none")
    f = p[p.fitness_hit]
    r = p[p.resistance_hit]
    ax0.scatter(f.kinome_conservation, f.ortholog_conservation, s=22, c="#2166AC",
                label=f"fitness hit (n={len(f)})", edgecolors="white", lw=0.3)
    ax0.scatter(r.kinome_conservation, r.ortholog_conservation, s=34, c="#B2182B",
                marker="D", label=f"resistance hit (n={len(r)})", edgecolors="white", lw=0.4)
    ax0.axhline(0.75, ls="--", lw=0.8, color="grey"); ax0.axvline(0.50, ls="--", lw=0.8, color="grey")
    ax0.set_xlabel("Kinome conservation (across human kinases)", fontsize=10)
    ax0.set_ylabel("Cross-species conservation (within the CDK)", fontsize=10)
    ax0.set_title("Where the two hit classes sit", fontsize=11.5, weight="bold")
    ax0.legend(fontsize=8, frameon=False, loc="lower left")
    sns.despine(ax=ax0)
    ax0.text(-0.13, 1.05, "a", transform=ax0.transAxes, fontsize=15, weight="bold")

    # hit rate by site class
    order = ["Tolerant in evolution", "Intermediate", "CDK-specific, fixed",
             "Universally immutable"]
    sub = d[d.site_class.isin(order)]
    n = sub.groupby("site_class").size().reindex(order)
    x = np.arange(len(order))
    for off, col, key, lab in [(-0.2, "#2166AC", "fitness_hit", "fitness"),
                               (0.2, "#B2182B", "resistance_hit", "resistance")]:
        k = sub.groupby("site_class")[key].sum().reindex(order)
        rate, lo, hi = [], [], []
        for kk, nn_ in zip(k, n):
            p_, l_, h_ = wilson(kk, nn_)
            rate.append(100 * p_); lo.append(100 * (p_ - l_)); hi.append(100 * (h_ - p_))
        ax1.bar(x + off, rate, 0.4, color=col, label=lab)
        ax1.errorbar(x + off, rate, yerr=[lo, hi], fmt="none", ecolor="black",
                     elinewidth=0.9, capsize=2.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"{o}\n(n={n[o]})" for o in order], fontsize=8, rotation=20, ha="right")
    ax1.set_ylabel("Hit rate (%)", fontsize=10)
    ax1.set_title("Hit rate by conservation class\n(Wilson 95% CI)", fontsize=11, weight="bold")
    ax1.legend(fontsize=8.5, frameon=False)
    sns.despine(ax=ax1)
    ax1.text(-0.13, 1.05, "b", transform=ax1.transAxes, fontsize=15, weight="bold")

    # consequence control: nonsense truncates and should deplete regardless of
    # how conserved the specific residue is -- taken from the full set
    cc = d_all[d_all.consequence.isin(["missense", "nonsense", "silent"])]
    order_c = [c for c in ["silent", "missense", "nonsense"]
               if c in set(cc.consequence)]
    sns.boxplot(data=cc, x="consequence", y="min_fitness_LFC", ax=ax2,
                order=order_c, hue="consequence", hue_order=order_c, legend=False,
                palette=["#CCCCCC", "#7FAFD4", "#2166AC"][:len(order_c)], fliersize=1.5)
    ax2.axhline(0, ls="--", color="black", lw=0.7)
    ax2.set_xlabel(""); ax2.set_ylabel("Fitness LFC (most depleted contrast)", fontsize=10)
    nn = cc.consequence.value_counts()
    ax2.set_xticklabels([f"{c}\n(n={nn.get(c,0)})" for c in order_c], fontsize=9)
    ax2.set_title("Positive control: nonsense edits deplete\n(why they are excluded from the primary test)",
                  fontsize=11, weight="bold")
    sns.despine(ax=ax2)
    ax2.text(-0.13, 1.05, "c", transform=ax2.transAxes, fontsize=15, weight="bold")

    fig.savefig(f"{FIG}/conservation_screen_two_axis.pdf", bbox_inches="tight")
    fig.savefig(f"{FIG}/conservation_screen_two_axis.png", dpi=280, bbox_inches="tight")
    print(f"wrote {paths.rel(FIG)}/conservation_screen_two_axis.{{pdf,png}}")


def main():
    d = sci.main()
    # primary analysis: missense only (nonsense deplete for reasons unrelated to
    # the conservation of the specific residue)
    mis = d[d.consequence == "missense"].copy()
    mis["is_CBE"] = (mis.Editor == "CBE").astype(int)
    mis["in_pocket"] = mis.in_KLIFS_pocket.fillna(False).astype(int)
    print(f"\nprimary set: {len(mis)} missense guides "
          f"({mis.fitness_hit.sum()} fitness, {mis.resistance_hit.sum()} resistance hits)")

    uni = univariate(mis)
    fits, models = [], {}
    for outcome in OUTCOMES:
        f, m = logistic(mis, outcome,
                        ["ortholog_conservation", "edit_efficiency", "is_CBE", "in_pocket"])
        fits.append(f); models[outcome] = m
    fits = pd.concat(fits, ignore_index=True)

    # pocket-restricted model, where kinome conservation is defined
    pocket = mis[mis.in_pocket == 1]
    pk_fits = []
    for outcome in OUTCOMES:
        if pocket[outcome].sum() >= 10:
            f, _ = logistic(pocket, outcome,
                            ["ortholog_conservation", "kinome_conservation",
                             "edit_efficiency", "is_CBE"])
            f["model"] = "KLIFS pocket only"
            pk_fits.append(f)
    pk_fits = pd.concat(pk_fits, ignore_index=True) if pk_fits else pd.DataFrame()

    rates = pd.concat([hit_rate_table(mis, o) for o in OUTCOMES], ignore_index=True)

    with pd.ExcelWriter(f"{CONS}/TABLE_conservation_vs_screen.xlsx", engine="openpyxl") as w:
        uni.round(5).to_excel(w, "Univariate", index=False)
        fits.round(5).to_excel(w, "Logistic_all_missense", index=False)
        if len(pk_fits):
            pk_fits.round(5).to_excel(w, "Logistic_pocket_only", index=False)
        rates.astype({"bin": str}).round(5).to_excel(w, "Hit_rate_by_bin", index=False)
        mis.round(4).to_excel(w, "Merged_guide_level", index=False)
    print(f"wrote {CONS}/TABLE_conservation_vs_screen.xlsx")

    print("\n=== UNIVARIATE ===");  print(uni.round(4).to_string(index=False))
    print("\n=== LOGISTIC (all missense) ===")
    print(fits[["outcome", "term", "OR", "lo", "hi", "p", "n", "n_hit"]].round(4).to_string(index=False))
    if len(pk_fits):
        print("\n=== LOGISTIC (KLIFS pocket only) ===")
        print(pk_fits[["outcome", "term", "OR", "lo", "hi", "p", "n", "n_hit"]].round(4).to_string(index=False))

    figure5(mis, uni, fits, rates)
    figure6(mis, d)
    return d, mis, uni, fits, rates


if __name__ == "__main__":
    main()
