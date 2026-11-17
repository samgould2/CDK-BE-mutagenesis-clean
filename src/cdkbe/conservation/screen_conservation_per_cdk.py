"""
Per-CDK test: does conservation predict screen outcome *within* each CDK?

Each CDK has its own inhibitors, its own essentiality and its own conservation
profile, so a pooled analysis can be driven entirely by differences between CDKs
rather than by conservation. Everything here is therefore estimated within a
single CDK and only then combined, using a random-effects meta-analysis that
keeps the between-CDK heterogeneity visible.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as ss
from statsmodels.discrete.discrete_model import Logit
from sklearn.metrics import roc_auc_score

from . import screen_conservation_integration as sci
from . import screen_conservation_stats as scs
from .. import paths

from ._paths import (ALN, BLAST, CONS, EDITING, FIG, KLIFS, LFC, ORTH,  # noqa: F401
                      OUT, ROOT)
plt.rc("font", family="Helvetica")


GENES = sci.GENES
OUTCOMES = {"fitness_hit": ("Fitness", "#2166AC"),
            "resistance_hit": ("Resistance", "#B2182B")}
MIN_HITS = 8          # below this a within-CDK model is not interpretable


def fit_one(sub, outcome, predictors):
    s = sub[[outcome] + predictors].dropna().copy()
    y = s[outcome].astype(int)
    if y.sum() < MIN_HITS or y.sum() == len(y):
        return None
    X = s[predictors].astype(float).copy()
    for c in predictors:
        sd = X[c].std()
        if sd == 0 or np.isnan(sd):
            return None
        X[c] = (X[c] - X[c].mean()) / sd
    X.insert(0, "const", 1.0)
    try:
        m = Logit(y, X).fit(disp=0, cov_type="HC1")
    except Exception:
        return None
    if not np.all(np.isfinite(m.bse.values)):
        return None
    return m


def per_cdk(d, predictors=("ortholog_conservation", "edit_efficiency")):
    predictors = list(predictors)
    rows = []
    for gene in GENES:
        sub = d[d.gene == gene]
        for outcome, (label, _) in OUTCOMES.items():
            s = sub[[outcome, "ortholog_conservation"]].dropna()
            y = s[outcome].astype(int)
            rec = {"gene": gene, "outcome": outcome, "outcome_label": label,
                   "n_guides": len(s), "n_hit": int(y.sum())}
            if y.sum() >= 3 and y.sum() < len(y):
                rec["AUC"] = roc_auc_score(y, s.ortholog_conservation)
                rec["MWU_p"] = ss.mannwhitneyu(
                    s.loc[y == 1, "ortholog_conservation"],
                    s.loc[y == 0, "ortholog_conservation"],
                    alternative="two-sided")[1]
                rec["median_hit"] = s.loc[y == 1, "ortholog_conservation"].median()
                rec["median_nonhit"] = s.loc[y == 0, "ortholog_conservation"].median()
            m = fit_one(sub, outcome, predictors)
            if m is not None:
                b = m.params["ortholog_conservation"]
                se = m.bse["ortholog_conservation"]
                rec.update({"logOR": b, "se": se, "OR": np.exp(b),
                            "lo": np.exp(b - 1.96 * se), "hi": np.exp(b + 1.96 * se),
                            "p": m.pvalues["ortholog_conservation"],
                            "OR_edit_eff": np.exp(m.params["edit_efficiency"]),
                            "p_edit_eff": m.pvalues["edit_efficiency"]})
            rows.append(rec)
    return pd.DataFrame(rows)


def meta(tab, outcome):
    """DerSimonian-Laird random-effects pooling of the within-CDK log odds ratios."""
    t = tab[(tab.outcome == outcome) & tab.logOR.notna()]
    if len(t) < 2:
        return None
    y, v = t.logOR.values, t.se.values ** 2
    w = 1 / v
    mu_fixed = (w * y).sum() / w.sum()
    Q = (w * (y - mu_fixed) ** 2).sum()
    df = len(y) - 1
    C = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    ws = 1 / (v + tau2)
    mu = (ws * y).sum() / ws.sum()
    se = np.sqrt(1 / ws.sum())
    return {"outcome": outcome, "k_CDKs": len(y), "OR": np.exp(mu),
            "lo": np.exp(mu - 1.96 * se), "hi": np.exp(mu + 1.96 * se),
            "p": 2 * (1 - ss.norm.cdf(abs(mu / se))),
            "Q": Q, "Q_p": 1 - ss.chi2.cdf(Q, df) if df > 0 else np.nan,
            "I2": max(0.0, 100 * (Q - df) / Q) if Q > 0 else 0.0, "tau2": tau2}


def figure(d, tab, metas):
    fig = plt.figure(figsize=(16.5, 11.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.25], hspace=0.34, wspace=0.28,
                          left=0.075, right=0.98, top=0.92, bottom=0.06)

    # --- a/b: per-CDK forest plots
    for i, (outcome, (label, col)) in enumerate(OUTCOMES.items()):
        ax = fig.add_subplot(gs[0, i])
        t = tab[(tab.outcome == outcome)].set_index("gene").reindex(GENES)
        ypos, ticks = [], []
        for j, gene in enumerate(GENES):
            r = t.loc[gene]
            y = j
            ypos.append(y)
            if pd.notna(r.get("OR", np.nan)):
                ax.plot([r["lo"], r["hi"]], [y, y], color=col, lw=1.7)
                ax.plot(r["OR"], y, "o", color=col, ms=6)
                ticks.append(f"{gene}  ({int(r['n_hit'])}/{int(r['n_guides'])})")
            else:
                ax.text(1.0, y, "too few hits", fontsize=7, color="grey",
                        ha="center", va="center", style="italic")
                ticks.append(f"{gene}  ({int(r['n_hit'])}/{int(r['n_guides'])})")
        m = metas.get(outcome)
        if m:
            y = len(GENES) + 0.6
            ax.plot([m["lo"], m["hi"]], [y, y], color="black", lw=2.4)
            ax.plot(m["OR"], y, "D", color="black", ms=8)
            ypos.append(y)
            ticks.append(f"Random-effects ({m['k_CDKs']} CDKs)")
            ax.text(m["hi"] * 1.08, y, f"OR {m['OR']:.2f}  p={m['p']:.3f}",
                    va="center", fontsize=7.5)
        ax.axvline(1, ls="--", color="grey", lw=0.8)
        ax.set_yticks(ypos); ax.set_yticklabels(ticks, fontsize=8.5)
        ax.invert_yaxis(); ax.set_xscale("log")
        ax.set_xlabel("Odds ratio per SD of conservation\n(adjusted for editing efficiency)",
                      fontsize=9.5)
        ax.set_title(f"{label}: conservation effect within each CDK",
                     fontsize=11.5, weight="bold")
        sns.despine(ax=ax)
        ax.text(-0.30, 1.05, "ab"[i], transform=ax.transAxes, fontsize=15, weight="bold")

    # --- c: per-CDK AUC
    axc = fig.add_subplot(gs[0, 2])
    w = 0.38
    x = np.arange(len(GENES))
    for off, (outcome, (label, col)) in zip([-w/2, w/2], OUTCOMES.items()):
        t = tab[tab.outcome == outcome].set_index("gene").reindex(GENES)
        axc.bar(x + off, t["AUC"].fillna(0), w, color=col, label=label)
        for xi, (a, nh) in enumerate(zip(t["AUC"], t["n_hit"])):
            if pd.notna(a):
                axc.text(xi + off, a + 0.015, f"{int(nh)}", ha="center", fontsize=6.5)
    axc.axhline(0.5, ls="--", color="grey", lw=0.8)
    axc.set_xticks(x); axc.set_xticklabels(GENES, rotation=45, fontsize=8.5, ha="right")
    axc.set_ylabel("AUC (conservation alone)", fontsize=10)
    axc.set_ylim(0, 1.0)
    axc.set_title("Discrimination per CDK\n(number above bar = hits)",
                  fontsize=11.5, weight="bold")
    axc.legend(fontsize=8, frameon=False, loc="upper left")
    sns.despine(ax=axc)
    axc.text(-0.16, 1.05, "c", transform=axc.transAxes, fontsize=15, weight="bold")

    # --- d: per-CDK conservation distributions, hit vs not
    gsl = gs[1, :].subgridspec(1, len(GENES), wspace=0.12)
    for j, gene in enumerate(GENES):
        ax = fig.add_subplot(gsl[0, j])
        sub = d[d.gene == gene]
        plot = []
        for outcome, (label, col) in OUTCOMES.items():
            s = sub[[outcome, "ortholog_conservation"]].dropna()
            plot.append(pd.DataFrame({"cons": s.loc[s[outcome], "ortholog_conservation"],
                                      "grp": label}))
        plot.append(pd.DataFrame({"cons": sub["ortholog_conservation"].dropna(),
                                  "grp": "all"}))
        pl = pd.concat(plot)
        sns.boxplot(data=pl, x="grp", y="cons", ax=ax, order=["all", "Fitness", "Resistance"],
                    hue="grp", hue_order=["all", "Fitness", "Resistance"], legend=False,
                    palette=["#DDDDDD", "#2166AC", "#B2182B"], fliersize=0, width=0.7)
        sns.stripplot(data=pl[pl.grp != "all"], x="grp", y="cons", ax=ax,
                      order=["all", "Fitness", "Resistance"], hue="grp",
                      hue_order=["Fitness", "Resistance"], legend=False,
                      palette=["#0B3D75", "#7A0E1E"], size=2.2, alpha=0.6, jitter=0.22)
        ax.set_title(gene, fontsize=10, weight="bold")
        ax.set_xlabel(""); ax.set_ylim(0.1, 1.05)
        ax.set_xticks(range(3))
        ax.set_xticklabels(["all", "fit", "res"], fontsize=7.5)
        if j == 0:
            ax.set_ylabel("Cross-species conservation", fontsize=10)
        else:
            ax.set_ylabel(""); ax.set_yticklabels([])
        sns.despine(ax=ax)
        if j == 0:
            ax.text(-0.55, 1.10, "d", transform=ax.transAxes, fontsize=15, weight="bold")

    fig.suptitle("Conservation and screen outcome, estimated separately within each CDK",
                 fontsize=14, weight="bold", y=0.965)
    fig.savefig(f"{FIG}/conservation_vs_screen_per_CDK.pdf", bbox_inches="tight")
    fig.savefig(f"{FIG}/conservation_vs_screen_per_CDK.png", dpi=280, bbox_inches="tight")
    print(f"wrote {paths.rel(FIG)}/conservation_vs_screen_per_CDK.{{pdf,png}}")


def main():
    d = pd.read_csv(f"{CONS}/screen_conservation_merged.csv")
    mis = d[d.consequence == "missense"].copy()
    tab = per_cdk(mis)
    metas = {o: meta(tab, o) for o in OUTCOMES}
    metas = {k: v for k, v in metas.items() if v}
    mt = pd.DataFrame(metas.values())

    with pd.ExcelWriter(f"{CONS}/TABLE_conservation_vs_screen_per_CDK.xlsx",
                        engine="openpyxl") as w:
        tab.round(5).to_excel(w, "Per_CDK", index=False)
        mt.round(5).to_excel(w, "Random_effects_meta", index=False)
    print(f"wrote {CONS}/TABLE_conservation_vs_screen_per_CDK.xlsx")

    cols = ["gene", "outcome_label", "n_guides", "n_hit", "median_hit",
            "median_nonhit", "AUC", "MWU_p", "OR", "lo", "hi", "p"]
    for o, (label, _) in OUTCOMES.items():
        print(f"\n=== {label} (per CDK) ===")
        print(tab[tab.outcome == o][cols].round(4).to_string(index=False))
    print("\n=== Random-effects meta-analysis ===")
    print(mt.round(4).to_string(index=False))

    figure(mis, tab, metas)
    return tab, mt


if __name__ == "__main__":
    main()
