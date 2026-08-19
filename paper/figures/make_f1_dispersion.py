"""Figure F1 -- SCPH-renormalized dispersion vs bare harmonic; stabilization.

Panel (a): phonon dispersion along G-X-M-G-R from the committed path tables --
bare harmonic (STO_anharm.xml fc2; unstable branches plotted at -|omega|,
dashed gray) vs SCPH-coupled at 300 K (solid). Linewidth whiskers (+/- Gamma,
HWHM) at the on-path 8x8x8 mesh modes of the Gamma-X segment.

Panel (b): zone-center soft-mode hardening Omega_0^2(T) -- computed (this
dataset) vs the digitized Tadano 2015 Fig. 5 SCPH curves (with/without
polarization mixing, PM) and experimental references; annotation records
the stabilization count 11/525 -> 0/525 at 300 K.

Data: data/processed/sto_alamode/{dispersion_path_bare,dispersion_path_scph300,
mode_kappa_T,eta_validity_T}.csv + data/literature/STO/tadano_2015/scph_omega_T.csv.
"""

from __future__ import annotations

import pubstyle

ROOT = pubstyle.add_fekappa_to_path()

import matplotlib.pyplot as plt
import pandas as pd

from fekappa.plots import set_publication_style

set_publication_style()

TIDY = ROOT / "data" / "processed" / "sto_alamode"
CM1_TO_MEV = 1.0 / 8.06554

bare = pd.read_csv(TIDY / "dispersion_path_bare.csv", comment="#")
scph = pd.read_csv(TIDY / "dispersion_path_scph300.csv", comment="#")
# Path breakpoints from the CSV header (G, X, M, G, R in 1/Ang).
with open(TIDY / "dispersion_path_scph300.csv") as fh:
    breaks_line = next(line for line in fh if line.startswith("# path breakpoints"))
pairs = [tok.strip().split("=") for tok in breaks_line.split(":")[1].split(",")]
break_pos = [float(v) for _, v in pairs]
break_lab = [{"G": r"$\Gamma$"}.get(k, k) for k, _ in pairs]

modes = pd.read_csv(TIDY / "mode_kappa_T.csv", comment="#")
m300 = modes[modes["T_K"] == 300.0]
# On-path mesh modes of the Gamma-X segment: q_frac = (0, 0, xi), xi = 0..1/2
# (q_index 1-5 on the 8x8x8 list); by cubic symmetry these sit at
# x = 2*xi*|GX| on the plotted [100] segment.
onpath = m300[m300["q_index"].between(1, 5)].copy()
xi = (onpath["q_index"] - 1) / 8.0
onpath["x_path"] = 2.0 * xi * break_pos[1]

branch_cols = [c for c in scph.columns if c.startswith("branch_")]

fig, (ax_a, ax_b) = plt.subplots(
    1, 2, figsize=(9.6, 3.6), gridspec_kw={"width_ratios": [1.9, 1.0]}
)

for col in branch_cols:
    ax_a.plot(bare["q_ang1"], bare[col], color="0.65", lw=0.9, ls="--", zorder=1)
    ax_a.plot(scph["q_ang1"], scph[col], color="tab:blue", lw=1.1, zorder=2)
ax_a.errorbar(
    onpath["x_path"], onpath["omega_cm1"], yerr=onpath["gamma_hwhm_cm1"],
    fmt="o", ms=2.5, lw=0.9, color="tab:orange", zorder=3,
    label=r"mesh modes $\pm\,\Gamma$ (HWHM)",
)
ax_a.axhline(0.0, color="0.4", lw=0.6)
for x in break_pos[1:-1]:
    ax_a.axvline(x, color="0.85", lw=0.6, zorder=0)
ax_a.set_xticks(break_pos, break_lab)
ax_a.set_xlim(break_pos[0], break_pos[-1])
ax_a.set_ylabel(r"$\omega$, $\Omega$ (cm$^{-1}$)")
handles = [
    plt.Line2D([], [], color="0.65", lw=0.9, ls="--",
               label=r"bare harmonic ($-|\omega|$ if unstable)"),
    plt.Line2D([], [], color="tab:blue", lw=1.1, label="SCPH-coupled, 300 K"),
    plt.Line2D([], [], color="tab:orange", marker="o", ms=2.5, lw=0.9,
               label=r"mesh modes $\pm\,\Gamma$ (HWHM)"),
]
ax_a.legend(handles=handles, loc="upper left", fontsize=8, framealpha=0.9)
pubstyle.panel_label(ax_a, "a")

val = pd.read_csv(TIDY / "eta_validity_T.csv", comment="#")
omega0_sq = (val["Omega0_cm1"] * CM1_TO_MEV) ** 2
ax_b.plot(val["T_K"], omega0_sq, "o-", ms=3.5, color="tab:blue",
          label="computed (fixed-cell SCPH)")

tad = pd.read_csv(ROOT / "data" / "literature" / "STO" / "tadano_2015" / "scph_omega_T.csv")
styles = {
    "SCPH_without_PM": ("s--", "tab:purple", "Tadano SCPH(-PM)"),
    "SCPH_with_PM": ("^--", "tab:green", "Tadano SCPH(+PM)"),
    "exp_Servoin": ("v:", "0.35", "exp. (Servoin)"),
}
for method, (fmt, color, label) in styles.items():
    c = tad[(tad["mode"] == "TO1") & (tad["method"] == method)].sort_values("T_K")
    ax_b.plot(c["T_K"], c["Omega_sq_meV2"], fmt, ms=4, color=color, lw=1, label=label)
pubstyle.mark_phase_boundary(ax_b)
ax_b.set_xlabel("$T$ (K)")
ax_b.set_ylabel(r"$\Omega_0^2$ (meV$^2$)")
ax_b.set_xlim(0, 1000)
ax_b.legend(fontsize=7.5, loc="upper left")
ax_b.annotate(
    "unstable modes\nbare: 11/525\nSCPH: 0/525",
    xy=(0.97, 0.05), xycoords="axes fraction", ha="right", va="bottom", fontsize=8,
    bbox={"boxstyle": "round", "fc": "white", "ec": "0.7"},
)
pubstyle.panel_label(ax_b, "b")

ax_a.set_xlabel("wave vector")
fig.tight_layout()
paths = pubstyle.save(fig, "f1_dispersion")
print("wrote:", *paths)
n_neg = int((bare[branch_cols].to_numpy() < -1.0).any(axis=0).sum())
print(f"bare branches dipping below -1 cm^-1 on the path: {n_neg}")