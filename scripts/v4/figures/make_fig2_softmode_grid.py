#!/usr/bin/env python3
"""Fig. 2: SCPH soft-mode frequencies at 300 K against the q1 grid, with the archived 2015 fit added.

Same panels, series and sources as make_figure_drafts.fig_softmodes, plus the archived 2015 fit in the configuration
of the 2015 study (archived force constants, mixed-space non-analytic term, Born data of that study), on the grids
2, 4, 8, 12 (soft-mode tables written by scripts/v4/scph_softmodes.py from those SCPH runs). This figure needs the
files of the 2015 study, which are not in this release.
Colour: categorical slot 4 of the reference palette, with its own marker (triangle) and line style (dotted).
Usage (repository root): make_fig2_softmode_grid.py KPETA_CUM_A KPETA_CUM_B B3_SOFTMODES_CSV SOURCE_BANDS OUTDIR
                         ARCHIVED_SOFTMODES_CSV [ARCHIVED_SOFTMODES_CSV ...]
  The archived soft-mode tables hold runs named e_q2, e_q4, e_q8, e_q12. Writes OUTDIR/fig_softmode_vs_grid_v2.png.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_figure_drafts as M  # noqa: E402

plt = M.plt
ARCH_COL = "#eda100"


def archived_e(files):
    rows = []
    for f in files:
        d = pd.read_csv(f)
        d = d[(d.T_K == 300) & d.run.isin(("e_q2", "e_q4", "e_q8", "e_q12")) & d.converged]
        rows.append(d)
    d = pd.concat(rows)
    d["q1"] = d.run.str.split("_q").str[1].astype(int)
    return d.sort_values("q1").set_index("q1")


LABEL = "archived 2015 fit, configuration of the 2015 study"
OUTNAME = "fig_softmode_vs_grid_v2.png"


def main():
    M.ARGS.update(cum_a=sys.argv[1], cum_b=sys.argv[2], b3=sys.argv[3])
    M.SRC, M.OUT = Path(sys.argv[4]), Path(sys.argv[5])
    M.OUT.mkdir(parents=True, exist_ok=True)
    e = archived_e(sys.argv[6:])
    orig_savefig = plt.Figure.savefig
    # draw the unchanged figure, then add the archived series before saving under the new name
    def savefig(fig, path, *a, **k):
        axs = fig.axes
        for ax, pt in zip(axs, ("Gamma", "M", "R")):
            full = set(e.index) >= {2, 4, 8, 12}        # line only when no grid between points is missing
            ax.plot(e.index, e[pt], color=ARCH_COL, lw=2, ls=(0, (1, 1.5)) if full else "none", marker="^", ms=6,
                    mec="white", mew=1,
                    label=LABEL)
        for lg in fig.legends:
            lg.remove()
        h, lab = axs[2].get_legend_handles_labels()
        fig.legend(h, lab, loc="lower center", ncol=3, fontsize=7.5, bbox_to_anchor=(0.5, -0.01))
        fig.tight_layout(rect=(0, 0.2, 1, 1))
        return orig_savefig(fig, M.OUT / OUTNAME, *a, **k)
    plt.Figure.savefig = savefig
    try:
        M.fig_softmodes()
    finally:
        plt.Figure.savefig = orig_savefig
    print("written", OUTNAME + "; archived grids:", list(e.index))


if __name__ == "__main__":
    main()
