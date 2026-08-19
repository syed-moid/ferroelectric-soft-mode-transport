"""Figure 3 -- kappa(T): 3+6 decomposition, negative control, benchmarks.

Panel (a): stacked particle sectors (acoustic / soft-optic / hard-optic,
boundary 175 cm^-1) + the directly computed Wigner coherence channel kappa_C
(sum of the six pair-sector blocks); dashed bare-RTA negative control; 300 K
benchmark markers. Panel (b): the six coherence blocks kappa^C_{sigma sigma'}
versus temperature. The light vertical line marks the ~105 K cubic-tetragonal
boundary (fixed-structure cubic reference below it).

Data: data/processed/sto_alamode/{kappa_five_way_T,kappa_coherence_blocks_T,
bare_vs_scph_control}.csv + data/literature/STO/{tadano_2015,tachibana_2008}/
kappa_total_T.csv.
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

blocks = pd.read_csv(TIDY / "kappa_coherence_blocks_T.csv", comment="#")
blocks = blocks[blocks["T_K"] > 0].sort_values("T_K")

fig, (ax, ax_b) = plt.subplots(
    1, 2, figsize=(9.8, 4.0), gridspec_kw={"width_ratios": [1.35, 1.0]}
)
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
pubstyle.mark_phase_boundary(ax)
ax.set_xlim(0, 1000)
ax.set_ylim(0, 24)
ax.set_xlabel("$T$ (K)")
ax.set_ylabel(r"$\kappa$ (W m$^{-1}$ K$^{-1}$)")
ax.legend(fontsize=8, loc="upper right")
pubstyle.panel_label(ax, "a")

BLOCK_STYLE = [
    ("kappa_C_hard_hard", "hard/hard", "#55a868", "-"),
    ("kappa_C_A_soft", "ac/soft", "#dd8452", "-"),
    ("kappa_C_soft_soft", "soft/soft", "#c44e52", "-"),
    ("kappa_C_A_A", "ac/ac", "#4c72b0", "-"),
    ("kappa_C_A_hard", "ac/hard", "#8172b3", "--"),
    ("kappa_C_hard_soft", "hard/soft", "0.4", ":"),
]
for col, label, color, ls in BLOCK_STYLE:
    ax_b.plot(blocks["T_K"], blocks[col], ls, color=color, lw=1.3,
              marker="o", ms=2.5, label=label)
pubstyle.mark_phase_boundary(ax_b, label=True)
ax_b.set_yscale("log")
ax_b.set_ylim(1e-4, 1.0)
ax_b.set_xlim(0, 1000)
ax_b.set_xlabel("$T$ (K)")
ax_b.set_ylabel(r"$\kappa^{C}_{\sigma\sigma'}$ (W m$^{-1}$ K$^{-1}$)")
ax_b.legend(fontsize=7.5, loc="lower right", ncols=2)
pubstyle.panel_label(ax_b, "b")

fig.tight_layout()
paths = pubstyle.save(fig, "f3_kappa_decomposition")
print("wrote:", *paths)