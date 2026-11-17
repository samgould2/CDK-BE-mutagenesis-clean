# External resistance-variant databases

Used by `notebooks/09_analysis_and_viz/CDK_resistance_patterns.ipynb` to compare the
base-editing resistance calls against curated clinical resistance, on the shared
KLIFS + inter-element-gap coordinate.

## CIViC (present)

`nightly-ClinicalEvidenceSummaries.tsv` — downloaded 2026-10-13 from
<https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv>.
Open access, no registration. Refresh with:

```bash
curl -L -o source_data/24_resistance_databases/nightly-ClinicalEvidenceSummaries.tsv \
  https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv
```

The notebook keeps rows with `significance == "Resistance"` and
`evidence_status == "accepted"` whose molecular profile is a single amino-acid
substitution, then maps them onto KLIFS positions for the kinases that have a json
in `source_data/16_KLIFS/oncokb_klifs/`.

## COSMIC (not included — needs a login)

COSMIC's resistance-mutation export sits behind a registered-user download; it cannot be
fetched unauthenticated, so it is not committed here. Academic registration is free.

1. Sign in at <https://cancer.sanger.ac.uk/cosmic/download>
2. Download **Resistance Mutations** for GRCh38 (file `Cosmic_ResistanceMutations_Tsv_v*.tar`)
3. Extract it and place the `.tsv` here as `Cosmic_ResistanceMutations.tsv`

The notebook detects that file automatically and adds COSMIC to the comparison; if it is
absent the COSMIC panel is skipped with a message and everything else still runs.

Expected columns (v98+): `GENE_SYMBOL`, `MUTATION_AA`, `DRUG_NAME`, `MUTATION_SOMATIC_STATUS`.
