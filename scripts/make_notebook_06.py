#!/usr/bin/env python3
"""Generate notebooks/06_structure_and_conservation/conservation.ipynb."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _lib.nbwrite import md, code, write, HEADER

ROOT = Path(__file__).resolve().parents[1]

cells = [
md("""# Evolutionary conservation of the kinase pocket

Two complementary views of which residues evolution is allowed to change.

1. **Across the human kinome** — the KLIFS 85-residue pocket alignment for every
   human kinase, scored overall and within each kinase group. Answers: *at this
   pocket position, do different kinases use different residues?*
2. **Across species, within each screened CDK** — a deep homolog search over all
   domains of life, aligned with MAFFT. Answers: *has evolution ever changed
   this residue in this particular CDK?*

A residue fixed on **both** axes is catalytic core and cannot be mutated. One
fixed across species but variable across the kinome is a CDK-specific
functional site — the interesting class for selective inhibitors and for
interpreting the base-editing screen.

**In** `data/reference/conservation/` (cached BLAST hits, ortholog alignments,
KLIFS tables). **Out** tables back into that directory, figures into
`figures/`.

The modules this drives live in `cdkbe.conservation`. They were already
standalone scripts; what changed is that their paths resolve from the
repository root rather than from `__file__/../..`."""),

code(HEADER + """
from cdkbe import paths
from cdkbe.conservation import (
    kinome_conservation as kc,
    kinome_conservation_figures as kcf,
    conservation_tables as ct,
    cdk_ortholog_conservation as co,
    cdk_ortholog_figures as cof,
)

pd.set_option("display.width", 200)
paths.ensure(paths.FIGURES)"""),

md("""## 1. Kinome-wide KLIFS conservation

Scores all 85 pocket positions across every human kinase. The kinase list is
cached in `data/reference/conservation/klifs_human_kinases.json`, so this does
not need the KLIFS API unless the cache is removed."""),

code("""master, per_group, cdk_tab, mat, info = kc.main()
print(mat.shape, "kinases x KLIFS positions")
info.group.value_counts()"""),

code("""# size-matched null: smaller groups always look more conserved, so correct for it
groups = {g: list(i) for g, i in info.groupby("group").groups.items() if len(i) >= 10}
groups["CDK family"] = list(info.index[info["family"].str.startswith("CDK", na=False)])

spec = kc.specificity_table(mat, info, groups, n_perm=1000)
spec.to_csv(Path(kc.OUT) / "klifs_group_specificity.csv", index=False)

spec[spec.group == "CDK family"].nlargest(10, "z_KL")[
    ["KLIFS_position", "region", "KL_vs_kinome", "z_KL", "z_conservation"]]"""),

code("""kcf.figure1(master, per_group, cdk_tab, mat)
kcf.figure2(master, spec)"""),

code("""tab = ct.main()
tab.head()"""),

md("""## 2. Deep ortholog conservation for the nine screened CDKs

The BLAST searches are run once by `cdkbe.conservation.blast_deep_search`,
taxon-stratified against UniProtKB because a single all-of-UniProt search
saturates its 1000-hit cap inside the vertebrates. This reuses the cached hit
tables in `data/reference/conservation/blast/`."""),

code("""cons = co.main()
cons.groupby("gene")["n_seqs_total"].first()"""),

code("""cons_all, merged, kin = cof.load_all()
out, comp = cof.build_tables(cons_all, merged, kin)
comp"""),

code("""cof.figure3(merged, comp, out)
cof.figure4(cons_all, merged)"""),

md("""## 3. Residues that cannot be mutated"""),

code("""out[out.site_class == "Universally immutable (catalytic core)"] \\
    .sort_values(["gene", "KLIFS_position"]).head(30)"""),

md("""## 4. CDK-specific sites — fixed across species, variable across the kinome"""),

code("""out[out.site_class == "CDK-specific, evolutionarily fixed"] \\
    .sort_values("ortholog_conservation", ascending=False).head(30)"""),

md("""## 5. Does conservation predict what the edit does in the screen?

Edit identity comes from the **canonical** HGVSp tables, filtered as figure 1
does (>=10% editing, >=100 sensor reads), restricted to unambiguous single
substitutions. Outcomes are per guide:

* **resistance** — enriched (FDR < 0.01, LFC > 1) under an inhibitor targeting
  that CDK
* **fitness** — depleted (FDR < 0.01, LFC < -1) with no drug, in T0 vs plasmid
  or DMSO vs T0

Two things are deliberate. Nonsense edits are held out of the primary test —
they truncate the protein, so they deplete regardless of how conserved the
edited residue is. And gRNA editing efficiency is carried as a covariate, since
a guide that does not edit cannot score in either direction."""),

code("""from cdkbe.conservation import screen_conservation_stats as scs

d, mis, uni, fits, rates = scs.main()
uni"""),

md("""### Per CDK

Each CDK has its own inhibitors, its own essentiality and its own conservation
profile, so a pooled estimate can be driven entirely by differences *between*
CDKs. The primary analysis runs within each CDK and is combined with a
random-effects meta-analysis, which keeps the between-CDK heterogeneity (I²)
visible."""),

code("""from cdkbe.conservation import screen_conservation_per_cdk as spc

per_cdk, meta = spc.main()
per_cdk"""),

code("""meta   # random-effects pooling of the within-CDK odds ratios"""),

md("""### Highest-confidence variants only

If the conservation effect is real it should *grow* as the hit definition
tightens, not wash out. Three tiers:

| | stats | editing | corroboration |
|---|---|---|---|
| T1 standard | FDR < 0.01, abs(LFC) > 1 | >=10% | — |
| T2 high | FDR < 0.001, abs(LFC) > 2 | >=50% | — |
| T3 replicated | FDR < 0.01, abs(LFC) > 1 | >=10% | same variant called by >=2 independent gRNAs |

Two design points. The editing floor is applied to the **whole population**,
hits and non-hits alike — filtering only the hits on a variable that itself
tracks outcome would manufacture the association. And guides that clear the
standard threshold but not the tier threshold are dropped rather than counted
as non-hits, so the comparison group is never contaminated with real-but-weak
hits.

T3 uses independent-gRNA replication rather than "significant in both fitness
contrasts", because **no guide in the screen is significant in both** T0-vs-plasmid
and DMSO-vs-T0 — those two contrasts capture different dropout windows."""),

code("""from cdkbe.conservation import screen_conservation_confidence as scc

per_tier, meta_tier, top_variants = scc.main()
meta_tier"""),

code("""top_variants[top_variants.outcome == "resistance"]   # the replicated resistance alleles"""),
]

p = write(ROOT / "notebooks" / "06_structure_and_conservation" / "conservation.ipynb", cells)
print(f"wrote {p.relative_to(ROOT)} ({len(cells)} cells)")
