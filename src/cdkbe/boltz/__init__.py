"""CDK9-ligand interaction fingerprints from Boltz-2 predicted complexes.

Supplementary figure 6's quantitative panels. Ported from
`manuscript_source_data/cdk9_boltz/cdk9_ifp.ipynb`, with its method written up
in `cdk9_idp_README.md` beside it.

    4 ligands  x  3 CDK9 variants (WT, L156F, I157V)
    x 3 random seeds  x 3 structural samples  =  108 predicted structures

Confidence scores are per structure; affinity outputs are per seed, so there
are 36 affinity records rather than 108. The three seeds are 635925368,
1821385089 and 2036692060.

WHAT IS IN THIS REPOSITORY. The delivered screen is 786 MB, almost all of it
MSA and embedding output the analysis never opens. The 108 `.cif` structures
and the 144 confidence/affinity JSONs it does open are copied into
`data/raw/cdk9_boltz/`, gzipped, at 8 MB -- so every panel rebuilds from the
predicted structures rather than from a shipped intermediate.

    paths   where the copied tree lives, and the design constants
    ifp     the loaders, the interaction fingerprint, and the weighted PCA

`cdkbe.source_data.cdk9_boltz` draws the panels.

These scores characterise atom-class proximity. They are not validated hydrogen
bonds, electrostatics, stacking geometries or refined energy terms, and the
binding probabilities are model output, not measurements.
"""
