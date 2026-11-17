# Sequencing pipeline

How `data/raw/` was produced. These three steps run on a SLURM cluster, once
per sample, and are **upstream of everything else in this repository** — the
notebooks start from their output.

You do not need to run them to reproduce any figure.

The FASTQs they take as input are deposited at NCBI under BioProject
[PRJNA1390412](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1390412); they
are not in this repository.

```
FASTQ ──step1──▶ reads split by guide ──step2──▶ CRISPResso per guide ──step3──▶ per-guide HGVSp calls
                 + guide count tables                                            = data/raw/screen_editing/
                 = data/raw/screen_counts/
```

| step | script | what it does |
|---|---|---|
| 1 | `step1_sensor_extraction.py` | matches each read pair to a guide by protospacer and barcode, writes per-guide FASTQs and the guide count table |
| 2 | `step2_crispresso.py` | runs CRISPResso per guide against a WT and a dummy-HDR amplicon |
| 3 | `step3_hgvsp_quant.py` | reads each guide's allele table and writes one row per distinct edit with its protein consequence |

`barcode_counter.py` is the standalone barcode counter used for QC.

## Deduplication

The old `cluster_scripts/` tree carried **six copies of each step**, one per
screen. They were byte-identical except for a single output directory name
(`crispresso` vs `crispresso3`), so there is now one copy of each:

- step 2 takes the directory as an optional third argument
- step 3 reads it from `CRISPRESSO_DIR`, defaulting to `./crispresso`

That collapsed 59 MB to 1.7 MB, with the inputs gzipped.

## Per-screen job scripts

`jobs/<screen>/` holds the SLURM wrappers and the config file naming each
sample's FASTQ pair, as they were submitted:

```
jobs/ABE/                 subpool 1, ABE arm
jobs/CDK12_13/            subpool 2
jobs/CDK2_4_6/            subpool 3
jobs/KB_compound_mutant/  compound-mutant screen
jobs/SY-5609/             SY-5609 arm
```

Each is a SLURM array job indexed by the config file, e.g.

```bash
sbatch jobs/ABE/sensor_extraction_guide_counts.sh
```

**The paths and the conda environments in those scripts are site-specific** —
they point at the BMC cluster and at environments built there. Treat them as a
record of how the data were generated rather than something that will run
elsewhere unchanged.

## Inputs

`inputs/` holds what the steps read, gzipped and deduplicated:

| file | used by | what it is |
|---|---|---|
| `CDK_library_CLUSTER.csv` | step 1 | the guide library, with protospacer and barcode |
| `CDK_pool{1,2,3}_targeting_dummy_sensor_alt.csv` | step 2 | WT and edited sensor amplicons per guide |
| `master_df_HGVSp_quant_pool{1,2,3}.csv` | step 3 | codon table per guide for translating edits |
| `cdks_info_updated.csv` | step 3 | CDK protein sequences and exon coordinates |
| `crispresso_quant_blank.csv` | step 3 | the empty frame used when a guide has no reads |

The same library also ships, in its analysis form, as
`data/reference/library/CDK_library_final.csv`.
