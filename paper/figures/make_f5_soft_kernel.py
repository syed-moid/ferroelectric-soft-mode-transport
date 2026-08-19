"""Figure F5 -- reduced soft-mode kernel: computed vs literature-input version.

kappa_soft^P(T) from the recomputed kernel (computed Omega_0, Gamma_0,
Lambda^2; Debye cutoff with the Lambda^2 +/-30% band), the Vogt-anchored
labeled variant (literature Omega_0(T), computed Gamma_0 and Lambda^2), the
literature-only value (0.137 [0.102, 0.167] W/mK at ~300 K,
committed in section4_eta_kappa_triples.csv), the direct soft-optic sector
sum from the zone decomposition, and kappa_total for scale (log axis).

Data: data/processed/sto_alamode/{soft_kernel_T,kappa_five_way_T,kappa_rta_T}.csv
+ data/processed/section4_eta_kappa_triples.csv.
"""

from __future__ import annotations

import pubstyle

ROOT = pubstyle.add_fekappa_to_path()

import matplotlib.pyplot as plt
import pandas as pd

from fekappa.plots import set_publication_style

set_publication_style()

TIDY = ROOT / "data" / "processed" / "sto_alamode"
kernel = pd.read_csv(TIDY / "soft_kernel_T.csv", comment="#")
kernel = kernel[kernel["q_max"] == "debye"]
computed = kernel[kernel["variant"] == "computed"].sort_values("T_K")
anchored = kernel[kernel["variant"] == "vogt_anchored"].sort_values("T_K")
five = pd.read_csv(TIDY / "kappa_five_way_T.csv", comment="#")
five = five[(five["boundary_cm1"] == 175.0) & (five["T_K"] > 0)].sort_values("T_K")
kappa_tot = pd.read_csv(TIDY / "kappa_rta_T.csv", comment="#")
kappa_tot = kappa_tot[kappa_tot["T_K"] > 0]
triples = pd.read_csv(ROOT / "data" / "processed" / "section4_eta_kappa_triples.csv")
old = triples[triples["material"] == "STO"].iloc[0]

fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.fill_between(computed["T_K"], computed["kappa_soft_lo"], computed["kappa_soft_hi"],
                color="tab:blue", alpha=0.20, lw=0)
ax.plot(computed["T_K"], computed["kappa_soft_W_mK"], "o-", ms=4, color="tab:blue",
        label=r"$\kappa^{P}_{\rm soft}$ kernel, computed inputs ($\pm 30\%\ \Lambda^2$)")
ax.plot(anchored["T_K"], anchored["kappa_soft_W_mK"], "s--", ms=4, color="tab:green",
        label=r"kernel, Vogt-anchored $\Omega_0$ (labeled variant)")
ax.errorbar(
    [old["T_K"]], [old["kappa_soft_W_m_K"]],
    yerr=[[old["kappa_soft_W_m_K"] - old["kappa_soft_lo"]],
          [old["kappa_soft_hi"] - old["kappa_soft_W_m_K"]]],
    fmt="D", ms=7, color="tab:red", capsize=4, zorder=5,
    label="literature-only kernel",
)
ax.plot(five["T_K"], five["kappa_soft_optic_P"], "^:", ms=4, color="tab:orange",
        label=r"soft-optic sector, zone sum (this work)")
ax.plot(kappa_tot["T_K"], kappa_tot["kappa_W_mK"], "-", color="0.4", lw=1.2,
        label=r"$\kappa_{\rm total}$ (scale)")
pubstyle.mark_phase_boundary(ax)
ax.set_yscale("log")
ax.set_xlim(0, 1000)
ax.set_ylim(3e-2, 40)
ax.set_xlabel("$T$ (K)")
ax.set_ylabel(r"$\kappa$ (W m$^{-1}$ K$^{-1}$)")
ax.legend(fontsize=8, loc="upper right")

fig.tight_layout()
paths = pubstyle.save(fig, "f5_soft_kernel")
print("wrote:", *paths)