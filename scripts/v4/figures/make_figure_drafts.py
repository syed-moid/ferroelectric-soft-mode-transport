#!/usr/bin/env python3
"""Figure drafts from existing data. No T4 data.

fig_kp_eta_cumulative.png : cumulative share of kappa_P carried by modes with eta' <= eta, 300 K, 16^3, models A
                            and B; group-averaged widths. Sources: the kp_eta cumulative tables of A and B
fig_softmode_vs_grid.png  : SCPH frequencies at Gamma (TO1), M and R (lowest branch) at 300 K against the q1 grid,
                            for model A, model B (H1), its companion fit (H2) and the source study (Table I).
                            Sources: data/processed/v3/r_point_check.csv (A), b3_scph_softmodes.csv (H1, H2)
                           , data/literature/STO/tadano_2015/scph_convergence_300K.csv (2015, meV), and the
                            source study of model B's data: its one-shot SCPH series on 2³ and 12³ at 300 K
                            (SOURCE_BANDS/SCP-oneshot-q{222,121212}/band.300.dat of the model B dataset,
                            THz; same extraction rule as scripts/v4/b3_source_compare.py).
Usage (repository root): make_figure_drafts.py KPETA_CUM_A KPETA_CUM_B B3_SOFTMODES_CSV SOURCE_BANDS OUTDIR
  SOURCE_BANDS: data/SrTiO3/figures/bands of the unpacked model B dataset (Materials Cloud Archive 2026.91).
"""
import sys
from pathlib import Path

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OUT = Path(".")
ARGS = {}
MEV_CM = 8.065544
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
COL = {"A": "#2a78d6", "B": "#eb6834", "H2": "#1baf7a", "2015": "#7a7974"}   # reference palette slots 1-3; gray
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "legend.frameon": False})


def fig_cumulative():
    c = pd.concat([pd.read_csv(ARGS["cum_a"]), pd.read_csv(ARGS["cum_b"])], ignore_index=True)
    fig, ax = plt.subplots(figsize=(3.6, 2.8), dpi=200)
    for run, lab, key in (("a345", "model A", "A"), ("H1", "model B", "B")):
        g = c[(c.run == run) & (c.kind == "mode")].sort_values("eta")
        x = np.r_[1e-4, g.eta.to_numpy()]
        y = np.r_[0.0, g.cum_share_kP_le_eta.to_numpy()]
        ax.step(x, 100 * y, where="post", color=COL[key], lw=2, label=lab)
    for t in (0.02, 0.05, 0.1):
        ax.axvline(t, color=INK2, lw=0.8, ls=(0, (3, 3)))
        ax.text(t, 2, f" η = {t}", rotation=90, va="bottom", ha="right", color=INK2, fontsize=7)
    ax.set_xscale("log")
    ax.set_xlim(2e-3, 0.3)
    ax.set_xticks([0.003, 0.01, 0.03, 0.1, 0.3])
    ax.set_xticklabels(["0.003", "0.01", "0.03", "0.1", "0.3"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_ylim(0, 101)
    ax.set_xlabel("damping ratio η = Γ / Ω  (group-averaged)")
    ax.set_ylabel(r"share of $\kappa_\mathrm{P}$ with $\eta' \leq \eta$  (%)")
    ax.legend(loc="upper left")
    ax.set_title("300 K, 16³ mesh", fontsize=8, color=INK2, loc="left")
    fig.tight_layout()
    fig.savefig(OUT / "fig_kp_eta_cumulative.png")
    plt.close(fig)


SRC = Path(".")
THZ = 33.35641


def source_oneshot(path):
    """Gamma = 4th branch at Gamma, M and R = lowest branch; path G-X-M-G-R-M (as scripts/b3_source_compare.py)."""
    blocks, cur = [], []
    for line in open(path):
        if line.startswith("#"):
            continue
        if not line.strip():
            if cur:
                blocks.append(np.array(cur))
                cur = []
            continue
        cur.append([float(x) for x in line.split()])
    if cur:
        blocks.append(np.array(cur))
    x, b = blocks[0][:, 0], np.array([bl[:, 1] for bl in blocks]) * THZ
    at = lambda xx: b[:, np.argmin(abs(x - xx))]  # noqa: E731
    return {"Gamma": at(1.7071)[3], "R": at(2.5731)[0], "M": at(1.0)[0]}


def fig_softmodes():
    a = pd.read_csv("data/processed/v3/r_point_check.csv", comment="#")
    a = a[(a.PM == "on") & (a.T_K == 300)].set_index("q1")
    b = pd.read_csv(ARGS["b3"])
    b = b[b.T_K == 300].set_index("run")
    t = pd.read_csv("data/literature/STO/tadano_2015/scph_convergence_300K.csv")
    t = t[t.method == "SCPH"].copy()
    t["q1"] = t.branch.str.extract(r"_(\d+)x").astype(int)
    grids = [2, 4, 8, 12]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), dpi=200, sharex=True)
    for ax, pt, ttl in zip(axs, ("Gamma", "M", "R"), ("Γ (TO1)", "M (lowest)", "R (lowest)")):
        series = {
            "model A": ("A", [a.loc[q, pt] for q in grids], grids),
            "model B (H1)": ("B", [b.loc[f"H1_q{q}", pt] for q in grids], grids),
            "companion fit (H2)": ("H2", [b.loc[f"H2_q{q}", pt] for q in grids], grids),
        }
        tq = t[t.q_label == pt].sort_values("q1")
        for name, (key, y, x) in series.items():
            ax.plot(x, y, color=COL[key], lw=2, marker="o", ms=4.5, mec="white", mew=1, label=name,
                    ls="-" if key != "H2" else (0, (4, 2)))
        ax.plot(tq.q1, tq.omega_meV * MEV_CM, color=COL["2015"], lw=1.2, marker="s", ms=4.5, mfc="white",
                mew=1.2, label="2015 study (Table I)")
        os_ = [source_oneshot(SRC / f"SCP-oneshot-q{q}{q}{q}" / "band.300.dat")[pt] for q in (2, 12)]
        ax.plot([2, 12], os_, color=INK, ls="none", marker="D", ms=4.5, mfc="white", mew=1.2,
                label="source study of model B's data, one-shot")
        ax.set_title(ttl, fontsize=9, color=INK, loc="left")
        ax.set_xticks([2, 4, 6, 8, 10, 12])
        ax.set_xticklabels(["2³", "4³", "6³", "8³", "10³", "12³"])
        ax.set_xlabel("SCPH q₁ grid")
    axs[0].set_ylabel("frequency at 300 K (cm⁻¹)")
    h, l = axs[2].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=7.5, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.16, 1, 1))
    fig.savefig(OUT / "fig_softmode_vs_grid.png")
    plt.close(fig)


if __name__ == "__main__":
    ARGS.update(cum_a=sys.argv[1], cum_b=sys.argv[2], b3=sys.argv[3])
    SRC, OUT = Path(sys.argv[4]), Path(sys.argv[5])
    OUT.mkdir(parents=True, exist_ok=True)
    fig_cumulative()
    fig_softmodes()
    print("written:", *sorted(p.name for p in OUT.glob("*.png")))
