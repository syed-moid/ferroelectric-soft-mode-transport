#!/usr/bin/env python3
"""Fig. S4: kappa_P against temperature in the configuration of the 2015 study, next to the 2015 values.

Archived fit, mixed-space non-analytic term, Born data of the 2015 study, q1 8^3, 12^3 transport; our runs next to
the numerical values behind Fig. 7 of the 2015 paper. One linear axis; 200 K, the lowest converged temperature of
the SCPH series, is marked with an open symbol and a note. Reads the table written by scripts/v4/p5_kappa_T.py. This
figure shows values of the files of the 2015 study, which are not in this release.
Usage: make_kappa_T_2015.py KAPPA_T_CSV OUT_PNG
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_figure_drafts as M  # noqa: E402

plt = M.plt
SRC = Path(sys.argv[1]) if len(sys.argv) > 2 else None
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else None


def main():
    d = pd.read_csv(SRC, comment="#")
    ours = d.dropna(subset=["kappa_P_xx"])
    fig, ax = plt.subplots(figsize=(3.6, 2.8), dpi=200)
    ax.plot(d.T_K, d.kappa_2015, color=M.COL["2015"], lw=2, ls=(0, (4, 2)), marker="s", ms=5, mec="white", mew=1,
            label="2015 study (numerical data of its Fig. 7)")
    ax.plot(ours.T_K, ours.kappa_P_xx, color=M.COL["A"], lw=2, marker="o", ms=5, mec="white", mew=1,
            label="this work, configuration of the 2015 study")
    low = ours[ours.T_K == 200]
    if len(low):
        ax.plot(low.T_K, low.kappa_P_xx, ls="none", marker="o", ms=7, mfc="white", mec=M.COL["A"], mew=1.5)
        ref200 = d[d.T_K == 200]  # keep the 2015 square at 200 K visible above the open circle
        ax.plot(ref200.T_K, ref200.kappa_2015, ls="none", marker="s", ms=5, color=M.COL["2015"], mec="white", mew=1,
                zorder=5)
        ax.annotate("lowest converged\nSCPH temperature", (200, float(low.kappa_P_xx.iloc[0])), xytext=(235, 4.6),
                    fontsize=6.5, color=M.INK2, arrowprops=dict(arrowstyle="-", color=M.INK2, lw=0.6))
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"$\kappa_\mathrm{P}$ (W m$^{-1}$ K$^{-1}$)")
    ax.set_xlim(150, 1050)
    ax.set_ylim(0, None)
    ax.grid(axis="y", color=M.GRID, lw=0.6)
    ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT)
    print("written:", OUT)


if __name__ == "__main__":
    main()
