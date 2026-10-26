"""Figures 1-4 of Manuscript v8 (draft 3), from production tables only (own SCPH q1 = 12^3, PM).

Fig. 1  convergence with transport mesh: kappa_P (a) and kappa_C (b) vs mesh at 300, 500, 800 K
Fig. 2  damping ratio: eta vs Omega on 16^3 at 300, 650, 1000 K (a); kappa-weighted fractions eta > 0.05, > 0.1 vs T (b)
Fig. 3  particle conductivity: kappa_P(T) on 16^3 with the published SCPH values and experiment, all literature
        points read from figure (a); sector shares, eigenvector-character rule (b)
Fig. 4  coherent conductivity: kappa_C(T) per mesh with the 300 K range over >= 15^3 (a); hard/hard block vs the
        other five coherence blocks, 16^3, character rule (b)

Colour (dataviz skill; palette validated with a Python port of validate_palette.js, white surface):
  temperature and mesh are ordered -> one-hue blue ramp (steps 350 / 500 / 700; distinct markers, Fig. 1 also line styles);
  sectors -> categorical slots 1-3 (blue / orange / aqua; all-pairs CVD dE 9.2, normal-vision 24.0), hard = aqua
  in Figs. 3 and 4; calculation vs experiment in Fig. 3(a) -> neutral ink with distinct markers.
Inputs: data/processed/sto_production/{kappa_mesh_T,eta_validity_T,kappa_sectors_T,coherence_blocks_T}.csv,
mode_table_T.csv.gz; data/processed/revision/exp_datasets.csv. Output: out/v8_fig{1..4}.{pdf,png}.
"""
from __future__ import annotations

import pubstyle

ROOT = pubstyle.add_fekappa_to_path()

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from fekappa.plots import set_publication_style  # noqa: E402

set_publication_style()
PROD = ROOT / "data/processed/sto_production"
PJ = ROOT / "data/processed/revision"
rd = lambda f: pd.read_csv(f, comment="#")
km, eta, sec, coh = (rd(PROD / f) for f in ("kappa_mesh_T.csv", "eta_validity_T.csv", "kappa_sectors_T.csv",
                                             "coherence_blocks_T.csv"))
RAMP3 = ["#5598e7", "#256abf", "#0d366b"]          # blue steps 350 / 500 / 700 (ordinal; darker for print and grayscale)
SECTOR = {"acoustic": "#2a78d6", "soft": "#eb6834", "hard": "#1baf7a"}
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#9a9893"
MAIN = (km.T_K >= 300)


def style(ax):
    ax.grid(True, color="0.9", lw=0.6, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ------------------------------------------------------------------ Fig. 1
fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.3))
for (col, lab), ax, letter in zip((("kappa_P", r"$\kappa_\mathrm{P}$"), ("kappa_C", r"$\kappa_\mathrm{C}$")), axs, "ab"):
    for T, c, mk, ls in zip((300, 500, 800), RAMP3, ("o", "s", "^"), ("-", "--", ":")):
        x = km[km.T_K == T].sort_values("mesh")
        ax.plot(x.mesh, x[col], marker=mk, ms=6, lw=1.6, ls=ls, color=c, label=f"{T} K", zorder=3)
    ax.set_xlabel(r"transport mesh $N$ ($N^3$)")
    ax.set_ylabel(lab + r" (W m$^{-1}$ K$^{-1}$)")
    ax.set_xticks([8, 12, 16, 20, 24])
    style(ax); pubstyle.panel_label(ax, letter)
    ax.axvline(16, color=INK3, lw=0.9, ls=":", zorder=1)          # designated reference mesh
    ax.annotate(r"reference ($16^3$)", xy=(16, 0.0), xycoords=("data", "axes fraction"), xytext=(3, 3),
                textcoords="offset points", va="bottom", ha="left", fontsize=7.5, color=INK2)
fig.legend(*axs[0].get_legend_handles_labels(), frameon=False, ncol=3, loc="upper center",
           bbox_to_anchor=(0.5, 1.06))
pubstyle.save(fig, "v8_fig1"); plt.close(fig)

# ------------------------------------------------------------------ Fig. 2
mt = rd(PROD / "mode_table_T.csv.gz")
fig, (a, b) = plt.subplots(1, 2, figsize=(8.8, 3.4), gridspec_kw={"width_ratios": [1.45, 1]})
for T, c in zip((300, 650, 1000), RAMP3):
    x = mt[(mt.T_K == T) & (mt.mesh == 16) & (mt.omega_cm1 > 1) & (mt.eta > 0)]
    a.scatter(x.omega_cm1, x.eta, s=5, color=c, alpha=0.55, linewidths=0, rasterized=True, label=f"{T} K", zorder=3)
a.axhline(0.5, color=INK2, lw=1.0, ls="--", zorder=2)
a.annotate(r"$\eta = 0.5$", xy=(1, 0.5), xycoords=("axes fraction", "data"), xytext=(-4, 3),
           textcoords="offset points", ha="right", va="bottom", fontsize=8, color=INK2)
a.set_yscale("log"); a.set_ylim(1e-3, 1)
a.set_xlabel(r"$\Omega$ (cm$^{-1}$)"); a.set_ylabel(r"$\eta = \Gamma/\Omega$")
leg = a.legend(frameon=False, loc="lower right", markerscale=3, title=r"$16^3$")
for h in leg.legend_handles:
    h.set_alpha(1)
style(a); pubstyle.panel_label(a, "a")
e = eta[(eta.mesh == 16) & (eta.T_K >= 300)].sort_values("T_K")
b.plot(e.T_K, 100 * e.frac_eta_gt_0p05_kappa, color=INK, marker="o", ms=5, lw=1.6, label=r"$\eta > 0.05$")
b.plot(e.T_K, 100 * e.frac_eta_gt_0p1_kappa, color=INK3, marker="s", ms=5, lw=1.6, ls="--", label=r"$\eta > 0.1$")
b.set_xlabel(r"$T$ (K)"); b.set_ylabel(r"$\kappa$-weighted fraction (%)")
b.legend(frameon=False); style(b); pubstyle.panel_label(b, "b")
pubstyle.save(fig, "v8_fig2"); plt.close(fig)

# ------------------------------------------------------------------ Fig. 3
ex = rd(PJ / "exp_datasets.csv")
fig, (a, b) = plt.subplots(1, 2, figsize=(8.8, 3.4))
m16 = km[(km.mesh == 16) & MAIN].sort_values("T_K")
a.plot(m16.T_K, m16.kappa_P, color=INK, lw=1.8, marker="o", ms=4, label=r"this work, $\kappa_\mathrm{P}$ ($16^3$, reference mesh)", zorder=4)
t = ex[ex.dataset == "Tadano2015_calc"]
a.plot(t.T_K, t.kappa_W_mK, color=INK2, lw=1.4, ls="--", marker="^", ms=5, mfc="white",
       label="SCPH, Ref. [1] (read from figure)", zorder=3)
mu = ex[ex.dataset == "Muta2005"]
a.errorbar(mu.T_K, mu.kappa_W_mK, yerr=mu.unc_W_mK, fmt="D", ms=5, color=INK, mfc="white", lw=0.8, capsize=0,
           label="Muta et al. (read from figure)", zorder=5)
po = ex[ex.dataset == "Popuri2014"]
a.errorbar(po.T_K, po.kappa_W_mK, yerr=[po.kappa_W_mK - po.bar_lo, po.bar_hi - po.kappa_W_mK], fmt="o", ms=5,
           color=INK3, mfc="white", lw=0.8, capsize=2, label="Popuri et al. (read from figure)", zorder=5)
tb = ex[ex.dataset == "Tachibana2008"]
a.errorbar(tb.T_K, tb.kappa_W_mK, yerr=tb.unc_W_mK, fmt="s", ms=6, color=INK, lw=1.0, capsize=3,
           label="Tachibana et al.", zorder=6)
a.set_xlabel(r"$T$ (K)"); a.set_ylabel(r"$\kappa$ (W m$^{-1}$ K$^{-1}$)"); a.set_xlim(280, 1020)
a.legend(frameon=False, fontsize=7.5); style(a); pubstyle.panel_label(a, "a")
s = sec[(sec.mesh == 16) & (sec.labels == "char") & (sec.T_K >= 300)].sort_values("T_K")
for col, name, ls in (("kappa_P_A", "acoustic", "-"), ("kappa_P_soft", "soft", "--"), ("kappa_P_hard", "hard", "-")):
    y = 100 * s[col] / s.kappa_P
    b.plot(s.T_K, y, color=SECTOR[name], lw=2, ls=ls, zorder=3)
    b.annotate(f"{name} optical" if name != "acoustic" else name, xy=(s.T_K.iloc[-1], y.iloc[-1]),
               xytext=(4, 0), textcoords="offset points", va="center", fontsize=8, color=INK)
b.set_xlabel(r"$T$ (K)"); b.set_ylabel(r"share of $\kappa_\mathrm{P}$ (%)"); b.set_ylim(0, 100); b.set_xlim(280, 1150)
b.set_xticks([300, 500, 700, 900])
style(b); pubstyle.panel_label(b, "b")
pubstyle.save(fig, "v8_fig3"); plt.close(fig)

# ------------------------------------------------------------------ Fig. 4
fig, (a, b) = plt.subplots(1, 2, figsize=(8.8, 3.4))
# markers = calculated values; thin lines only join them (no interpolated or expected values are drawn)
for m, c, mk in zip((12, 16, 20), RAMP3, ("o", "s", "^")):
    x = km[(km.mesh == m) & MAIN].sort_values("T_K")
    a.plot(x.T_K, x.kappa_C, color=c, lw=0.8, ls=":" if m == 20 else "-", alpha=0.7, zorder=2)
    a.plot(x.T_K, x.kappa_C, color=c, marker=mk, ms=5, lw=0, label=f"${m}^3$", zorder=3)
r300 = km[(km.T_K == 300) & (km.mesh >= 15)].kappa_C
a.vlines(300, r300.min(), r300.max(), color=INK, lw=3, alpha=0.35, zorder=4)
a.plot([300] * len(r300), r300, ls="none", marker="_", ms=9, mew=1.4, color=INK, zorder=5)
a.annotate("observed mesh spread\n" + r"300 K, $15^3$–$24^3$", xy=(300, r300.min()), xytext=(8, -2),
           textcoords="offset points", fontsize=7.5, color=INK, va="top")
a.set_xlabel(r"$T$ (K)"); a.set_ylabel(r"$\kappa_\mathrm{C}$ (W m$^{-1}$ K$^{-1}$)"); a.set_ylim(0, None)
a.legend(frameon=False, loc="lower right", title="mesh"); style(a); pubstyle.panel_label(a, "a")
c16 = coh[(coh.mesh == 16) & (coh.labels == "char") & (coh.T_K >= 300)].sort_values("T_K")
other = c16.kappa_C_pairsum - c16.kC_hard_hard
b.plot(c16.T_K, c16.kappa_C_pairsum, color=INK, lw=1.8, label=r"$\kappa_\mathrm{C}$ (all pairs)", marker="o", ms=3, zorder=3)
b.plot(c16.T_K, c16.kC_hard_hard, color=SECTOR["hard"], lw=2.2, label="hard/hard", marker="o", ms=3, zorder=4)
b.plot(c16.T_K, other, color=INK3, lw=1.6, ls="--", label="other five contributions", marker="o", ms=3, zorder=3)
b.set_xlabel(r"$T$ (K)"); b.set_ylabel(r"$\kappa_\mathrm{C}$ (W m$^{-1}$ K$^{-1}$)")
b.set_ylim(0, 0.8)
b.legend(frameon=False, loc="upper right"); style(b); pubstyle.panel_label(b, "b")
pubstyle.save(fig, "v8_fig4"); plt.close(fig)
print("wrote out/v8_fig1..4 (.pdf, .png)")
