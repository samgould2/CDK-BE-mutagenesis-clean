"""Shared matplotlib style for the T47D cell-cycle notebooks in this folder.

Modeled on the visual conventions of the lab's CDK12/CDK13 base-editing-screen
figure: bold sans-serif type for categorical (gene/variant) labels, open axes
(no top/right spines, a slightly heavier axis and tick weight than the
matplotlib default), and a muted/dusty color palette rather than pastel
fills.

Usage, at the top of a notebook, before any figures are created:

    import figure_style
    figure_style.apply_style()
"""

import logging

import matplotlib.pyplot as plt

# Sans-serif fallback chain (Arial first, to match the reference figure;
# falls through to fonts that are actually available on this machine). Arial
# and Helvetica are frequently missing on Linux/Mac systems -- matplotlib
# still falls back correctly, it's just noisy about it (see apply_style).
FONT_FAMILY = ["Arial", "Helvetica", "Nimbus Sans", "DejaVu Sans"]


def apply_style():
    """Call once per notebook/script, before creating any figures."""
    # Silence matplotlib's "findfont: Font family 'Arial'/'Helvetica' not
    # found" warnings -- harmless (it still falls back to an available font
    # in FONT_FAMILY and renders/saves fine), just noisy on machines where
    # those two aren't installed.
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

    plt.rcParams.update({
        "font.family": FONT_FAMILY,
        "font.size": 11,
        "axes.linewidth": 1.3,
        "axes.edgecolor": "black",
        "xtick.major.width": 1.3,
        "ytick.major.width": 1.3,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.color": "black",
        "ytick.color": "black",
        "legend.frameon": False,
    })


# Cell-cycle phase fills (stacking order G0 -> G1 -> S -> G2/M), muted/dusty
# rather than pastel, in the same tonal family as the reference figure.
PHASE_COLORS = {
    "G0": "#D9A441",   # dusty gold
    "G1": "#5B84B1",   # steel blue
    "S": "#9A82B0",    # dusty mauve/purple
    "G2M": "#7FB49B",  # sage green
}

# Treatment marker/line colors, shared across the dot-plot and fold-change
# figures (only the keys a given script uses need to be present).
TREATMENT_COLORS = {
    "DMSO": "#8C8C86",   # neutral grey
    "Palbo": "#D9702E",  # burnt orange
    "ss": "#3E6FA8",     # steel navy ("serum free")
}
