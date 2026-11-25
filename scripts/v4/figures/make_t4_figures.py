#!/usr/bin/env python3
"""T4 figure, one per model: who carries the heat, and who takes part in scattering the
carriers. 300 K, 16^3, frozen primary definition c_min 0.5, TRISYM = 1, group-averaged widths. Same axes and scales
for both models. Data only; no interpretation in the figure.

Panels
 a  kappa_P carried by each family (all modes)                    <- A/t4/t4_{A,B}_summary.json
 b  share of kappa_P carried by modes with eta above a threshold  <- A/t4/t4_{A,B}_modes.csv (all modes)
 c  per carrier: share of its linewidth from processes with >= 1 soft partner, against the carrier's frequency;
    colour = the carrier's own family; marker area ~ its kappa_P share  <- A/t4/t4_fig_{A,B}.csv
 d  carriers' linewidth by number of soft partners (kappa_P-weighted mean)   <- summary json
 e  carriers' linewidth by partner family pair (kappa_P-weighted mean)       <- summary json

Control figure, fig_t5_control.png: share of the carriers' linewidth per partner mode (% per 1000
modes), rot and non-rot (neither ∪ pol) partners, per frequency band, both models  <- A/t5/t5_{A,B}.json

T5b figure, fig_t5b_control.png: the same quantity in the bands below 80, 80-95, 95-105 and the merged
band above 105 cm^-1, both models  <- A/t5b/t5b_{A,B}.json

Usage: make_t4_figures.py DIR_A DIR_B OUTDIR [t5b]
  DIR_A, DIR_B: processed outputs of model A and model B (data/processed/v4/a, data/processed/v4/b); A/ above
  stands for DIR_A or DIR_B.
Writes OUTDIR/fig_t4_{A,B}.png and fig_t5_control.png; with "t5b" only fig_t5b_control.png
"""
import json
from pathlib import Path

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import sys  # noqa: E402

DIRS = {}
OUT = Path(".")
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#ffffff"
FAMCOL = {"neither": "#7a7974", "pol": "#2a78d6", "rot": "#1baf7a"}   # neutral gray; reference palette slots 1, 3
BAR = "#52514e"
DEF = "primary_c0.5"
PAIRS = ["neither/neither", "neither/rot", "neither/pol", "rot/rot", "pol/rot", "pol/pol"]
plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "legend.frameon": False,
                     "axes.titlesize": 8, "axes.titleweight": "bold", "axes.titlelocation": "left"})


def hbars(ax, labels, vals, colors, xmax=100):
    y = np.arange(len(labels))[::-1]
    ax.barh(y, vals, height=0.6, color=colors, edgecolor=SURF, linewidth=1)
    for yi, v in zip(y, vals):
        ax.text(v + 1.5, yi, f"{v:.1f} %" if v >= 0.05 else f"{v:.2f} %", va="center", fontsize=7, color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(0, xmax)
    ax.grid(axis="y", visible=False)
    ax.tick_params(axis="y", length=0)


def figure(model):
    s = json.load(open(DIRS[model] / f"t4/t4_{model}_summary.json"))
    d = s["definitions"][DEF]
    M = pd.read_csv(DIRS[model] / f"t4/t4_{model}_modes.csv")
    M = M[M.definition == DEF]
    F = pd.read_csv(DIRS[model] / f"t4/t4_fig_{model}.csv")

    fig = plt.figure(figsize=(7.2, 5.0), dpi=200)
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], hspace=0.55, wspace=0.75)
    fig.suptitle(f"Model {model}: 300 K, 16³, primary definition (c_min 0.5)", x=0.01, ha="left",
                 fontsize=9, fontweight="bold", color=INK)
    fig.text(0.01, 0.935, "Who carries the heat", fontsize=8, color=INK2)
    fig.text(0.01, 0.47, f"Who takes part in scattering the carriers ({s['carriers']} modes carrying "
             f"{100 * s['carrier_kappa_share']:.1f} % of " + r"$\kappa_\mathrm{P}$)", fontsize=8, color=INK2)

    ax = fig.add_subplot(gs[0, 0])
    k = d["kappa_P_share_by_family"]
    hbars(ax, ["neither", "pol", "rot"], [100 * k[f] for f in ("neither", "pol", "rot")],
          [FAMCOL[f] for f in ("neither", "pol", "rot")])
    ax.set_xlabel(r"share of $\kappa_\mathrm{P}$ (%)")
    ax.set_title(r"a  $\kappa_\mathrm{P}$ by family")

    ax = fig.add_subplot(gs[0, 1:])
    e = M.eta.fillna(0).to_numpy()
    w = M.kappa_mode_W_mK.to_numpy() / M.kappa_mode_W_mK.sum()
    o = np.argsort(e)
    xs, ys = e[o], 1 - np.cumsum(w[o])                      # share of kappa_P with eta' > x
    ax.step(np.r_[1e-4, xs], 100 * np.r_[1.0, ys], where="post", color=INK, lw=2)
    ax.axvline(0.1, color=INK2, lw=0.8, ls=(0, (3, 3)))
    ax.text(0.1, 60, f" η > 0.1: {100 * s['share_kP_eta_gt_0.1']:.4f} %", fontsize=7, color=INK2, ha="left")
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 0.3)
    ax.set_xticks([0.001, 0.003, 0.01, 0.03, 0.1, 0.3])
    ax.set_xticklabels(["0.001", "0.003", "0.01", "0.03", "0.1", "0.3"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_ylim(0, 101)
    ax.set_xlabel("damping ratio η = Γ / Ω (group-averaged)")
    ax.set_ylabel(r"share of $\kappa_\mathrm{P}$ with $\eta' > \eta$ (%)")
    ax.set_title(r"b  $\kappa_\mathrm{P}$ by damping ratio (all modes)")

    ax = fig.add_subplot(gs[1, 0])
    soft = 100 * (F.nsoft1_frac + F.nsoft2_frac)
    for fam in ("neither", "pol", "rot"):
        sel = F.carrier_family == fam
        if sel.any():
            ax.scatter(F.omega_cm1[sel], soft[sel], s=8 + 900 * F.kappa_share[sel], color=FAMCOL[fam],
                       edgecolor=SURF, linewidth=0.6, label=f"{fam} ({int(sel.sum())})", zorder=3)
    mean = 100 * (d["kappa_weighted_nsoft_frac"]["1"] + d["kappa_weighted_nsoft_frac"]["2"])
    ax.axhline(mean, color=INK2, lw=0.8, ls=(0, (3, 3)), label=f"κ-weighted mean {mean:.1f} %")
    ax.set_xlim(0, 450)
    ax.set_ylim(0, 100)
    ax.set_xlabel("carrier frequency Ω (cm⁻¹)")
    ax.set_ylabel("linewidth with ≥ 1 soft partner (%)")
    ax.set_title("c  per carrier")
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], ls="", marker="o", ms=4, color=FAMCOL[f], label=f"{f} ({int((F.carrier_family == f).sum())})")
          for f in ("neither", "pol", "rot") if (F.carrier_family == f).any()]
    hs.append(Line2D([], [], color=INK2, lw=0.8, ls=(0, (3, 3)), label=f"mean {mean:.1f} %"))
    ax.legend(handles=hs, title="carrier family", title_fontsize=6, fontsize=6, loc="upper right",
              handlelength=1.2, borderaxespad=0.2)

    ax = fig.add_subplot(gs[1, 1])
    n = d["kappa_weighted_nsoft_frac"]
    hbars(ax, ["no soft", "one soft", "both soft"], [100 * n[x] for x in "012"], [BAR] * 3)
    ax.set_xlabel("share of carriers' linewidth (%)")
    ax.set_title("d  by soft partners")

    ax = fig.add_subplot(gs[1, 2])
    p = d["kappa_weighted_pair_frac"]
    hbars(ax, PAIRS, [max(100 * p[x], 0.0) for x in PAIRS], [BAR] * len(PAIRS))
    ax.set_xlabel("share of carriers' linewidth (%)")
    ax.set_title("e  by partner families")

    fig.savefig(OUT / f"fig_t4_{model}.png", facecolor=SURF, bbox_inches="tight")
    plt.close(fig)


def fig_control(ctl="t5"):
    if ctl == "t5":
        bands = ["<105", "105-145", "145-160", ">=160"]
        blab = ["< 105", "105–145", "145–160", "≥ 160"]
        ymax = 33
    else:
        bands = ["<80", "80-95", "95-105", ">=105"]
        blab = ["< 80", "80–95", "95–105", "≥ 105"]
        ymax = 66
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), dpi=200, sharey=True)
    for ax, m in zip(axs, "AB"):
        d = json.load(open(DIRS[m] / f"{ctl}/{ctl}_{m}.json"))["bands"]
        x = np.arange(len(bands))
        for off, key, col, nm in ((-0.19, "rot", FAMCOL["rot"], "rot partners"),
                                  (0.19, "nonrot", FAMCOL["neither"], "non-rot partners (neither ∪ pol)")):
            v = [d[b][key]["per_mode"]["total"] for b in bands]
            n = [d[b][key]["N_modes"] for b in bands]
            y = [1e5 * a if a is not None else 0.0 for a in v]
            ax.bar(x + off, y, width=0.36, color=col, edgecolor=SURF, linewidth=1, label=nm)
            for xi, yi, a, ni in zip(x + off, y, v, n):
                ax.text(xi, yi + 0.4 * ymax / 33, "—" if a is None else f"{yi:.1f}", ha="center", va="bottom", fontsize=6.5,
                        color=INK)
        nr = [d[b]["rot"]["N_modes"] for b in bands]
        nn = [d[b]["nonrot"]["N_modes"] for b in bands]
        ax.set_xticks(x)
        ax.set_xticklabels([f"{bl}\nN {a} / {c}" for bl, a, c in zip(blab, nr, nn)], fontsize=7)
        ax.tick_params(axis="x", length=0)
        ax.set_axisbelow(True)
        ax.set_xlabel("partner frequency band (cm⁻¹)")
        ax.set_title(f"Model {m}")
        ax.grid(axis="x", visible=False)
        ax.set_ylim(0, ymax)
    axs[0].set_ylabel("share of carriers' linewidth\nper partner mode (% per 1000 modes)")
    fig.text(0.5, 0.09, "N: partner modes in the band on the full 16³ mesh, rot / non-rot; — : no rot modes in the band",
             ha="center", fontsize=6.5, color=INK2)
    h, lb = axs[0].get_legend_handles_labels()
    fig.legend(h, lb, loc="lower center", ncol=2, fontsize=7, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    fig.savefig(OUT / f"fig_{ctl}_control.png", facecolor=SURF)
    plt.close(fig)


if __name__ == "__main__":
    DIRS.update(A=Path(sys.argv[1]), B=Path(sys.argv[2]))
    OUT = Path(sys.argv[3])
    OUT.mkdir(parents=True, exist_ok=True)
    if sys.argv[4:] == ["t5b"]:
        fig_control("t5b")
        print(OUT / "fig_t5b_control.png")
        sys.exit()
    for m in "AB":
        figure(m)
    fig_control()
    print(OUT / "fig_t4_A.png", OUT / "fig_t4_B.png", OUT / "fig_t5_control.png")
