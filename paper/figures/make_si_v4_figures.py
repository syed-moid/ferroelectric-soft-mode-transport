"""Figures S1-S8 of Supplementary Information v4, from production, revision and SI tables only.

S1  SCPH soft-mode frequencies (Gamma-TO1, M, R) against T for q1 = 8^3, 12^3, 16^3 (PM; converged points only)
S2  benchmark matrix: kappa_P(300 K) against the SCPH grid q1, with and without polarization mixing, 8^3 and 12^3 meshes
S3  kappa_C convergence: change between meshes against r of the coarser mesh (a); top-20 pair share against mesh (b)
S4  hard/hard share of kappa_C against T per mesh, eigenvector-character rule (a) and frequency rule (b)
S5  largest damping ratio against T per mesh (a); kappa-weighted fraction with eta > 0.05 per mesh (b)
S6  four-phonon: Gamma_4 against epsilon at 4^3, unmodified and audit build (a); eta_4 of every sampled case (b)
S7  R-point mode: Omega_R(300 K) against q1 with Tadano & Tsuneyuki Table I (a); Omega_R(T) against experiment (b)
S8  damping ratio against mean free path, 16^3, 300 K, with the 2.8 A and d(Ti-O) lines

Colour as in make_v8_figures.py: ordered quantities (T, mesh, q1) on the one-hue blue ramp; sectors on categorical slots
1-3; calculation vs reference in neutral ink with distinct markers. Output: out/si_fig{1..8}.{pdf,png}.
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
SI = ROOT / "data/processed/si_v4"
rd = lambda f: pd.read_csv(f, comment="#")
RAMP3 = ["#5598e7", "#256abf", "#0d366b"]          # blue steps 350 / 500 / 700 (ordinal; darker for print and grayscale)
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#9a9893"
SECTOR = {"acoustic": "#2a78d6", "soft": "#eb6834", "hard": "#1baf7a"}
MEV = 8.0655


def style(ax):
    ax.grid(True, color="0.9", lw=0.6, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ------------------------------------------------------------------ S1
rpc = rd(PJ / "r_point_check.csv")
fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.2), sharex=True)
for ax, (col, lab), letter in zip(axs, (("Gamma", "Γ-TO1"), ("M", "M"), ("R", "R")), "abc"):
    for q1, c, mk in zip((8, 12, 16), RAMP3, ("o", "s", "^")):
        x = rpc[(rpc.PM == "on") & (rpc.q1 == q1) & rpc.converged].sort_values("T_K")
        ax.plot(x.T_K, x[col], color=c, marker=mk, ms=4, lw=1.4, label=f"$q_1 = {q1}^3$", zorder=3)
    ax.axvline(300, color=INK3, lw=0.8, ls=":", zorder=1)
    ax.set_xlabel(r"$T$ (K)"); ax.set_ylabel(fr"$\Omega$ ({lab}) (cm$^{{-1}}$)")
    style(ax); pubstyle.panel_label(ax, letter)
fig.legend(*axs[0].get_legend_handles_labels(), frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.07))
pubstyle.save(fig, "si_fig1"); plt.close(fig)

# ------------------------------------------------------------------ S2
bm = rd(PJ / "c1_benchmark_matrix.csv")
fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.2), sharey=True)
for ax, mesh, letter in zip(axs, (8, 12), "ab"):
    for pm, mk, mfc, lab in ((True, "o", INK, "with polarization mixing"), (False, "s", "white", "without")):
        x = bm[(bm.transport_mesh == mesh) & (bm.PM == pm)].sort_values("kmesh_scph")
        ax.plot(x.kmesh_scph, x.kappa_P, color=INK, marker=mk, mfc=mfc, ms=6, lw=1.2, label=lab, zorder=3)
    ax.axhline(8.7, color=INK2, lw=1.0, ls="--", zorder=2)
    ax.annotate("Tadano & Tsuneyuki (read from figure)", xy=(8, 8.7), xytext=(0, -4), textcoords="offset points",
                ha="right", va="top", fontsize=7.5, color=INK2)
    ax.set_xticks([2, 4, 8]); ax.set_xticklabels(["2³", "4³", "8³"])
    ax.set_xlabel(r"SCPH grid $q_1$"); ax.set_title(f"transport mesh {mesh}³", fontsize=9)
    style(ax); pubstyle.panel_label(ax, letter)
axs[0].set_ylabel(r"$\kappa_\mathrm{P}$(300 K) (W m$^{-1}$ K$^{-1}$)")
axs[1].legend(frameon=False, fontsize=8, loc="upper right")
pubstyle.save(fig, "si_fig2"); plt.close(fig)

# ------------------------------------------------------------------ S3
pr = rd(SI / "s3_pairs.csv"); d2 = rd(PJ / "d2_summary.csv")
fig, (a, b) = plt.subplots(1, 2, figsize=(8.6, 3.3))
a.axhspan(-2.5, 2.5, color="0.94", zorder=0)
a.axvline(0.88, color=INK3, lw=0.9, ls="--", zorder=1)
for _, r in pr.iterrows():
    reg = isinstance(r.registered, str) and r.registered
    a.plot(r.r_coarser, r.dkC_pct, marker="D" if reg else "o", ms=7 if reg else 6, color=INK, mfc="white" if reg else INK,
           zorder=3)
    dx, dy, ha = {300: (6, 7, "left"), 400: (6, -8, "left"), 650: (-6, -9, "right"), 1000: (-6, -11, "left")}.get(
        r.T_K if r.pair == "16³→20³" else 0, (6, 0, "left"))
    a.annotate(f"{r.T_K} K, {r.pair}" + (f"\n{r.registered}" if reg else ""), xy=(r.r_coarser, r.dkC_pct),
               xytext=(dx, dy), textcoords="offset points", fontsize=7, va="center", ha=ha, color=INK2)
a.set_xlabel(r"$r$ of the coarser mesh"); a.set_ylabel(r"$\Delta\kappa_\mathrm{C}$ (%)"); a.set_xlim(0.35, 1.45); a.set_ylim(-17, 26)
style(a); pubstyle.panel_label(a, "a")
for T, c, mk in zip((300, 500, 800), RAMP3, ("o", "s", "^")):
    x = d2[d2.T_K == T].sort_values("mesh")
    b.plot(x.mesh, 100 * x.top20_share, color=c, marker=mk, ms=5, lw=1.4, label=f"{T} K", zorder=3)
b.set_xlabel(r"transport mesh $N$ ($N^3$)"); b.set_ylabel(r"share of $\kappa_\mathrm{C}$ in the top-20 pairs (%)")
b.set_xticks([12, 16, 20, 24]); b.set_ylim(0, None); b.legend(frameon=False)
style(b); pubstyle.panel_label(b, "b")
pubstyle.save(fig, "si_fig3"); plt.close(fig)

# ------------------------------------------------------------------ S4
hh = rd(SI / "s5_hh_share.csv")
fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.2), sharey=True)
for ax, col, title, letter in zip(axs, ("hh_char_pct", "hh_freq_pct"), ("eigenvector-character rule", "frequency rule"), "ab"):
    for m, c, mk in zip((12, 16, 20), RAMP3, ("o", "s", "^")):
        x = hh[(hh.mesh == m)].dropna(subset=[col]).sort_values("T_K")
        ax.plot(x.T_K, x[col], color=c, marker=mk, ms=4, lw=1.4, label=f"{m}³", zorder=3)
    o = hh[hh.mesh.isin([15, 17, 24])].dropna(subset=[col])
    ax.plot(o.T_K, o[col], ls="none", marker="x", ms=6, color=INK, label="15³, 17³, 24³", zorder=4)
    ax.set_xlabel(r"$T$ (K)"); ax.set_title(title, fontsize=9); style(ax); pubstyle.panel_label(ax, letter)
axs[0].set_ylabel(r"hard/hard share of $\kappa_\mathrm{C}$ (%)"); axs[0].set_ylim(0, 100)
axs[1].legend(frameon=False, fontsize=8, loc="lower right", title="mesh")
pubstyle.save(fig, "si_fig4"); plt.close(fig)

# ------------------------------------------------------------------ S5
ev = rd(PROD / "eta_validity_T.csv"); ev = ev[ev.T_K >= 300]
fig, (a, b) = plt.subplots(1, 2, figsize=(8.4, 3.2))
for m, c, mk in zip((12, 16, 20), RAMP3, ("o", "s", "^")):
    x = ev[ev.mesh == m].sort_values("T_K")
    a.plot(x.T_K, x.eta_max, color=c, marker=mk, ms=4, lw=1.4, label=f"{m}³", zorder=3)
    b.plot(x.T_K, 100 * x.frac_eta_gt_0p05_kappa, color=c, marker=mk, ms=4, lw=1.4, label=f"{m}³", zorder=3)
a.set_ylabel(r"largest $\eta = \Gamma/\Omega$")
b.set_ylabel(r"$\kappa$-weighted fraction, $\eta > 0.05$ (%)"); b.set_ylim(0, None)
for ax, letter in ((a, "a"), (b, "b")):
    ax.set_xlabel(r"$T$ (K)"); style(ax); pubstyle.panel_label(ax, letter)
a.legend(frameon=False, title="mesh", loc="upper right", ncol=3, fontsize=8)
a.set_ylim(0, 0.27)
pubstyle.save(fig, "si_fig5"); plt.close(fig)

# ------------------------------------------------------------------ S6
g = rd(SI / "s7_gate_4cubed.csv"); cases = rd(PJ / "eta4_margin_cases.csv")
fig, (a, b) = plt.subplots(1, 2, figsize=(8.6, 3.3))
for (mode, lab), mk in zip((("Gamma-TO1", "Γ-TO1"), ("M", "M"), ("R", "R")), ("o", "s", "^")):
    x = g[g["mode"] == mode]
    a.plot(x.eps, x.G4_unpatched, color=INK3, marker=mk, ms=5, lw=1.2, ls="--", zorder=3)
    a.plot(x.eps, x.G4_patched, color=INK, marker=mk, ms=5, lw=1.4, zorder=4, label=lab)
a.set_yscale("log"); a.set_xlabel(r"$\varepsilon$ (cm$^{-1}$)"); a.set_ylabel(r"$\Gamma_4$ (cm$^{-1}$), $4^3$, 300 K")
a.set_xticks([2, 5, 10])
a.annotate("unmodified", xy=(10, 268), xytext=(-4, 6), textcoords="offset points", ha="right", fontsize=7.5, color=INK2)
a.annotate("identity set excluded", xy=(10, 2.5), xytext=(-4, -12), textcoords="offset points", ha="right", fontsize=7.5,
           color=INK)
a.legend(frameon=False, loc="center right", fontsize=8); style(a); pubstyle.panel_label(a, "a")
GRP = {"soft": ("soft modes (rule)", "o", INK), "hard optic (1000 K check)": ("hard optical, 1000 K (rule)", "s", INK),
       "hard optic (descriptive)": ("off-symmetry, 8³ (descriptive)", "D", INK),
       "heat carrier": ("heat carriers (300, 800 K)", "^", "white")}
for grp, (lab, mk, fc) in GRP.items():
    x = cases[cases.group.str.startswith(grp)]
    b.scatter(x.omega_cm1, x.eta4, marker=mk, s=22, facecolors=fc, edgecolors=INK, linewidths=0.8, label=lab, zorder=3)
b.set_xscale("log"); b.set_ylim(0, 0.1); b.set_xlim(20, 700)
b.set_xticks([30, 100, 300]); b.set_xticklabels(["30", "100", "300"]); b.minorticks_off()
b.set_xlabel(r"$\Omega$ (cm$^{-1}$)"); b.set_ylabel(r"$\eta_4 = \Gamma_4/\Omega$")
b.legend(frameon=False, fontsize=7.5, loc="upper right"); style(b); pubstyle.panel_label(b, "b")
pubstyle.save(fig, "si_fig6"); plt.close(fig)

# ------------------------------------------------------------------ S7
td = rd(SI / "s7b_tadano.csv"); ex = rd(SI / "s7b_experiment.csv")
fig, (a, b) = plt.subplots(1, 2, figsize=(8.4, 3.2))
x = td[td.PM == "on"].sort_values("q1")
a.plot(x.q1, x.R_ours, color=INK, marker="o", ms=6, lw=1.4, label="this work (PM)", zorder=3)
y = x.dropna(subset=["R_tadano"])
a.plot(y.q1, y.R_tadano, color=INK2, marker="^", mfc="white", ms=6, lw=1.2, ls="--", label="Tadano & Tsuneyuki, Table I", zorder=3)
a.set_xticks([2, 4, 8, 12, 16]); a.set_xlabel(r"SCPH grid $q_1$ ($q_1^3$)"); a.set_ylabel(r"$\Omega_\mathrm{R}$(300 K) (cm$^{-1}$)")
a.legend(frameon=False, loc="center right", fontsize=8); style(a); pubstyle.panel_label(a, "a")
r12 = rpc[(rpc.PM == "on") & (rpc.q1 == 12) & rpc.converged & (rpc.T_K >= 250)].sort_values("T_K")
b.plot(r12.T_K, r12.R, color=INK, lw=1.6, marker="o", ms=3.5, label=r"this work ($q_1 = 12^3$)", zorder=3)
b.errorbar(ex.T_K, ex.R_exp_cm_read, yerr=ex.R_exp_unc_cm, fmt="s", ms=6, mfc="white", color=INK2, capsize=3,
           label="Fauqué et al. (read from figure)", zorder=4)
b.set_xlabel(r"$T$ (K)"); b.set_ylabel(r"$\Omega_\mathrm{R}$ (cm$^{-1}$)"); b.set_ylim(0, None)
b.legend(frameon=False, loc="lower right", fontsize=8); style(b); pubstyle.panel_label(b, "b")
pubstyle.save(fig, "si_fig7"); plt.close(fig)

# ------------------------------------------------------------------ S8
mt = rd(PROD / "mode_table_T.csv.gz")
x = mt[(mt.T_K == 300) & (mt.mesh == 16) & (mt.omega_cm1 > 1e-3)].copy()
v = np.sqrt(x.v2_star_m2s2 / x.multiplicity)
x["mfp_A"] = v / (2 * 2 * np.pi * 2.99792458e10 * x.gamma_hwhm_cm1) * 1e10
fig, ax = plt.subplots(figsize=(4.8, 3.4))
for name, c in SECTOR.items():
    y = x[x.sector_char == {"acoustic": "A"}.get(name, name)]
    ax.scatter(y.mfp_A, y.eta, s=5, color=c, alpha=0.6, linewidths=0, rasterized=True, label=name, zorder=3)
for cut, lab, dx, ha in ((2.8, "2.8 Å", 3, "left"), (3.8963 / 2, "d(Ti–O)", -3, "right")):
    ax.axvline(cut, color=INK2, lw=0.9, ls="--", zorder=2)
    ax.annotate(lab, xy=(cut, 0.95), xycoords=("data", "axes fraction"), xytext=(dx, 0), textcoords="offset points",
                fontsize=7.5, color=INK2, ha=ha)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-2, 2e3)  # modes with vanishing |v| lie off the left edge
ax.set_xlabel(r"mean free path $|\mathbf{v}|\tau$ (Å)"); ax.set_ylabel(r"$\eta = \Gamma/\Omega$")
leg = ax.legend(frameon=False, markerscale=3, loc="upper right", fontsize=8, title="character rule", title_fontsize=8)
for h in leg.legend_handles:
    h.set_alpha(1)
style(ax)
pubstyle.save(fig, "si_fig8"); plt.close(fig)
print("wrote out/si_fig1..8 (.pdf, .png)")
