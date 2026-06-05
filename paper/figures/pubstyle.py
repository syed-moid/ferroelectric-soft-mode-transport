"""Thin helpers for the manuscript figure scripts.

NO restyling. Each figure reproduces its original notebook plotting cell verbatim
(colors, markers, sizes, line styles, fonts, ticks, figure size, scales/limits) via
``fekappa.plots.set_publication_style`` -- the exact style the originals were rendered
with. This module only provides: a path shim, a margin panel-label, and a dual
PDF + 600-dpi-PNG saver. Deterministic (Agg backend).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless + deterministic

OUT_DIR = Path(__file__).resolve().parent / "out"


def add_fekappa_to_path() -> Path:
    """Put the repo's ``src/`` on ``sys.path`` so ``fekappa`` imports here. Returns repo root."""
    here = Path(__file__).resolve()
    for p in here.parents:
        if (p / "pyproject.toml").exists():
            src = p / "src"
            if src.exists() and str(src) not in sys.path:
                sys.path.insert(0, str(src))
            return p
    return here.parent


def panel_label(ax, letter: str):
    """Place ``(a)``/``(b)`` in the OUTER top-left margin (just above the axes box),
    so it cannot overlap an interior legend or the data."""
    text = letter if letter.startswith("(") else f"({letter})"
    ax.annotate(text, xy=(0, 1), xycoords="axes fraction",
                xytext=(0, 6), textcoords="offset points",
                fontsize=13, fontweight="bold", ha="left", va="bottom",
                annotation_clip=False)


def save(fig, name: str):
    """Write ``out/<name>.pdf`` (vector) and ``out/<name>.png`` (600 dpi). Returns paths.
    Uses ``bbox_inches='tight'`` exactly as the original notebook saves did."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pdf = OUT_DIR / f"{name}.pdf"
    png = OUT_DIR / f"{name}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=600, bbox_inches="tight")
    return pdf, png
