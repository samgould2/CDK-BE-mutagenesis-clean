"""Single-cell CROP-seq (BEAT-seq) analysis behind figure 5.

The delivered bundle shipped outputs only -- no analysis code came with it --
so `paths` reconstructs the design from the READMEs and the object metadata,
and the modules here are the re-analysis built on top of the delivered tables.

    paths        design, and where every delivered table lives
    singlet_io   the only module that opens a .h5ad; those are not vendored
    norm         held-out control normalisation and the n-unbiased statistic D
    nullcal      size-matched calibration of the pseudo-guide null
    severity     the severity map and the pre-registered tests
    fig_panels   all 46 source panels, one file each
    style        the manuscript's figure conventions

`cdkbe.source_data.figure5` maps these panels onto the published figure.
"""
