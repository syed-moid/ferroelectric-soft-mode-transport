"""Figure F3 -- kappa(T): five-way decomposition, negative control, benchmarks.

Stacked areas: particle-like sectors (acoustic / soft-optic / hard-optic,
boundary 175 cm^-1) + the directly computed Wigner coherence channel kappa_C.
Dashed: bare-RTA kappa (negative control -- same anharmonic FCs, harmonic part
not renormalized). Markers at 300 K: Tadano 2015 first-principles 8.7 W/mK
(12x12x12) and Tachibana 2008 experiment ~11 W/mK.

Data: data/processed/sto_alamode/{kappa_five_way_T,bare_vs_scph_control}.csv +
data/literature/STO/{tadano_2015,tachibana_2008}/kappa_total_T.csv.
"""

from __future__ import annotations

import pubstyle

ROOT = pubstyle.add_fekappa_to_path()

import matplotlib.pyplot as plt
import pandas as pd

from fekappa.plots import set_publication_style

set_publication_style()

TIDY = ROOT / "data" / "processed" / "sto_alamode"
LIT = ROOT / "data" / "literature" / "STO"

five = pd.read_csv(TIDY / "kappa_five_way_T.csv", comment="#")
five = five[five["boundary_cm1"] == 175.0].sort_values("T_K")
control = pd.read_csv(TIDY / "bare_vs_scph_control.csv", comment="#")
bare = control[control["record"] == "kappa_T"].sort_values("T_K")
tadano = pd.read_csv(LIT / "tadano_2015" / "kappa_total_T.csv")
tachibana = pd.read_csv(LIT / "tachibana_2008" / "kappa_total_T.csv")

fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.stackplot(
    five["T_K"],
    five["kappa_acoustic_P"], five["kappa_soft_optic_P"],
    five["kappa_hard_optic_P"], five["kappa_C"],
    labels=[
        r"$\kappa_{\mathrm{ac}}^{P}$ (acoustic)",
        r"$\kappa_{\mathrm{soft}}^{P}$ (soft optic, $\Omega<175$ cm$^{-1}$)",
        r"$\kappa_{\mathrm{hard}}^{P}$ (hard optic)",
        r"$\kappa_{C}$ (Wigner coherence, direct)",
    ],
    colors=["#4c72b0", "#dd8452", "#55a868", "#c44e52"],
    alpha=0.85,
)
ax.plot(bare["T_K"], bare["bare"], "--", color="0.25", lw=1.4,
        label="bare-RTA (negative control)")
ax.errorbar(
    tadano["T_K"], tadano["kappa_W_per_m_K"], yerr=tadano["kappa_uncertainty"],
    fmt="D", ms=6, color="tab:purple", zorder=5, capsize=3,
    label=r"Tadano 2015 (FP, $12^3$)",
)
ax.errorbar(
    tachibana["T_K"], tachibana["kappa_W_per_m_K"], yerr=tachibana["kappa_uncertainty"],
    fmt="*", ms=11, color="k", zorder=5, capsize=3,
    label=r"Tachibana 2008 (exp., $\approx$)",
)
ax.set_xlim(0, 1000)
ax.set_ylim(0, 24)
ax.set_xlabel("$T$ (K)")
ax.set_ylabel(r"$\kappa$ (W m$^{-1}$ K$^{-1}$)")
ax.legend(fontsize=8.5, loc="upper right")

fig.tight_layout()
paths = pubstyle.save(fig, "f3_kappa_decomposition")
print("wrote:", *paths)