"""Publication matplotlib style for the fekappa figures."""

from __future__ import annotations

import matplotlib as mpl
from cycler import cycler

# Colorblind-safe palette drawn from the matplotlib "tab" colors.
PUBLICATION_COLORS: tuple[str, ...] = (
    "tab:blue",
    "tab:red",
    "tab:green",
    "tab:orange",
)


def set_publication_style() -> None:
    """Apply the project-wide publication matplotlib rcParams.

    Sets a serif / Computer-Modern-math style suitable for the manuscript,
    a tighter axes line width, sensible screen vs. saved DPI, and a
    colorblind-safe color cycle. Mutates the global :data:`matplotlib.rcParams`.
    """
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 11,
            "axes.linewidth": 0.8,
            "mathtext.fontset": "cm",
            "figure.dpi": 130,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.prop_cycle": cycler(color=list(PUBLICATION_COLORS)),
        }
    )
