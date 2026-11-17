"""Central path and design configuration for the CDK BEAT-seq re-analysis.

The delivered bundle ships outputs only -- no analysis code came with it -- so this
module reconstructs the experimental design from the READMEs and the object metadata.

Design: 164 guides (92 ABE + 72 CBE) across CDK2/4/6/7/8/9/12/13/19 + NT controls,
in A549 and T47D, untreated, 10x 3' SC3Pv4, CellRanger 10.0.0, GRCh38.
"""

from pathlib import Path

from .. import io as _io
from .. import paths as _repo

ROOT = _repo.ROOT
#: the derived tables, copied into this repository so figure 5 is self-contained
RESULTS = _repo.PROCESSED / "beatseq"
FIGURES = _repo.FIGURES


#: The delivered bundle's derived tables are copied into this repository, under
#: `outputs/` and `metadata/` mirroring the two delivered trees. Every figure-5
#: and single-cell supplementary panel reads from these; nothing reaches outside
#: the repository. Some were gzipped on the way in, so each constant goes through
#: `io.resolve`, which accepts either spelling.
OUT = RESULTS / "outputs"
META = RESULTS / "metadata"


def _t(root: Path, rel: str) -> Path:
    """A delivered table, by its path in the original bundle."""
    return _io.resolve(root / rel)


#: The single-cell objects themselves are NOT in this repository: ~34 GB of
#: .h5ad that nothing drawing a figure reads -- every panel comes from the
#: derived tables above. These paths document where each object came from;
#: opening one raises a normal FileNotFoundError naming the path.
#: See data/external/NOT_INCLUDED.csv.
_ABSENT = _repo.EXTERNAL / "CDK-BEAT-seq" / "CDK_Beatseq_outputs_20260813" / "cdk_lib" / "beatseq-cdk"

# --- analysis-ready objects -------------------------------------------------
# Only these are used. The norm_pop_* objects (~34 GB) and norm_pop_all.hdf (44 GB)
# are deliberately never touched: recomputing normalization from the singlets takes
# seconds, and the .hdf twins are redundant re-encodings for the perturbseq package.
SINGLETS = {
    "A549": _ABSENT / "05_primary_analysis/outputs/cdk_singlets_A549.h5ad",
    "T47D": _ABSENT / "05_primary_analysis/outputs/cdk_singlets_T47D.h5ad",
}
MEAN_POP = _ABSENT / "05_primary_analysis/outputs/mean_pop_all.h5ad"

# --- delivered tables we reuse rather than recompute -------------------------
T_GUIDE_DE = _t(OUT, "05_primary_analysis/outputs/tables/de/guide_de_summary.csv")
T_PSEUDOBULK = _t(OUT, "05_primary_analysis/outputs/tables/de/pseudobulk_control_norm.csv.gz")
T_PB_ANNOT = _t(OUT, "05_primary_analysis/outputs/tables/de/pseudobulk_guide_annotation.csv")
T_CELLCYCLE_PERCELL = _t(OUT, "05_primary_analysis/outputs/tables/cellcycle/cellcycle_percell.csv.gz")
T_CELLCYCLE_POWER = _t(OUT, "05_primary_analysis/outputs/tables/cellcycle/cellcycle_power_summary.csv")
T_CELLS_PER_GUIDE = _t(OUT, "plotting/figures/qc/cells_per_guide.csv")

# Guide-level annotation. NOTE: `tier`, `variant_class` and `LFC_median_T0_Plasmid`
# are NOT in the singlet obs -- they live here and in mean_pop_all.h5ad obs.
T_UNTREATED_EXPECT = _t(META, "02_library_qc/outputs/untreated_expectation.csv")
T_GUIDE_EDITING = _t(META, "02_library_qc/outputs/guide_editing_qc.csv")
T_SHARED_PROTO = _t(META, "02_library_qc/outputs/shared_protospacers.csv")
GUIDE_LIBRARY = {
    "ABE": _t(META, "01_library_design/outputs/cdk_ABE_guide_library.csv"),
    "CBE": _t(META, "01_library_design/outputs/cdk_CBE_guide_library.csv"),
}

# --- design -----------------------------------------------------------------
# gem_group is perfectly nested inside arm. A549 arms have 2 gem groups (true
# replicates); both T47D arms have exactly one, so T47D cannot support any
# statistic that relies on technical replication.
ARMS = ["A549/ABE", "A549/CBE", "T47D/ABE", "T47D/CBE"]
ARM_GEM_GROUPS = {
    "A549/ABE": [1, 2],
    "A549/CBE": [3, 4],
    "T47D/ABE": [5],
    "T47D/CBE": [6],
}
CELL_LINES = ["A549", "T47D"]
EDITORS = ["ABE", "CBE"]

CDK_GENES = ["CDK2", "CDK4", "CDK6", "CDK7", "CDK8", "CDK9", "CDK12", "CDK13", "CDK19"]

# Functional subpools from the parent paper (10.1101/2025.11.02.685764).
SUBPOOLS = {
    "transcription_initiation": ["CDK7", "CDK8", "CDK9", "CDK19"],
    "elongation_termination": ["CDK12", "CDK13"],
    "cell_cycle": ["CDK2", "CDK4", "CDK6"],
}

# --- known data hazards -----------------------------------------------------
# Four targeting protospacers appear in BOTH editor arms installing DIFFERENT
# alleles. Same spacer, same integration biology, different amino-acid change --
# these are the cleanest controlled contrasts in the experiment. Any cross-arm
# join must key on (gRNA_id, editor); joining on sgRNA or id alone silently
# collapses these pairs.
SHARED_PROTOSPACER_PAIRS = {
    "gRNA_CDK2_targ_7379": {"ABE": "V225A", "CBE": "E224K"},
    "gRNA_CDK6_targ_8160": {"ABE": "V76A_V77A", "CBE": "V77I"},
    "gRNA_CDK9_targ_1201": {"ABE": "I157V", "CBE": "L156F"},  # L156F = paper's headline variant
    "gRNA_CDK12_targ_2287": {"ABE": "M816V_D817G", "CBE": "H818Y"},
}

# Two NT controls share a protospacer (GGCTAGGTCCGGTAAGTGCG), so gRNA_nt_9243 is
# never recovered and the effective number of independent NT units is 9, not 10.
# That is the hard ceiling on any guide-level permutation null.
NT_DUPLICATE_PAIR = ("gRNA_nt_8877", "gRNA_nt_9243")
N_EFFECTIVE_NT = 9

# Genes excluded from every global/aggregate score. Ribosomal and mitochondrial
# genes are short, highly expressed and depth-sensitive -- they will single-handedly
# manufacture a spurious gene-length slope. IEGs are simultaneously the best
# short-half-life markers and the worst dissociation-stress artifact, so they are
# excluded from global fits and only analysed as a declared, caveated panel.
EXCLUDE_PREFIXES = ("MT-", "MT_", "RPL", "RPS", "MRPL", "MRPS")
IEG_GENES = ["FOS", "FOSB", "JUN", "JUNB", "EGR1", "IER2", "DUSP1", "NR4A1", "ZFP36"]


def arm_cell_line(arm: str) -> str:
    return arm.split("/")[0]


def arm_editor(arm: str) -> str:
    return arm.split("/")[1]


def ensure_dirs() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
