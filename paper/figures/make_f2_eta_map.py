"""Figure F2 -- zone-resolved damping map eta_qs(T) and validity summary.

Panel (a): eta = Gamma/Omega vs Omega for every mode at every temperature
(color = T, log-y), with the narrow-width (eta = 0.5) and overdamped (eta = 1)
lines. Panel (b): zone-wide max and median eta vs T (computed inputs), the
Vogt-rescaled soft-branch max (sensitivity, 150-300 K), and the zone-center
soft-mode eta both ways; the validity fractions are identically zero at every
temperature (annotated).

Data: data/processed/sto_alamode/{eta_map_T,eta_validity_T,
eta_validity_sensitivity}.csv.
"""

from __future__ import annotations

import pubstyle

ROOT = pubstyle.add_fekappa_to_path()

import matplotlib.pyplot as plt
import pandas as pd

from fekappa.plots import set_publication_style

set_publication_style()

TIDY = ROOT / "data" / "processed" / "sto_alamode"
eta_map = pd.read_csv(TIDY / "eta_map_T.csv", comment="#")
val = pd.read_csv(TIDY / "eta_validity_T.csv", comment="#")
sens = pd.read_csv(TIDY / "eta_validity_sensitivity.csv", comment="#")

fig, (ax_a, ax_b) = plt.subplots(
    1, 2, figsize=(9.2, 3.5), gridspec_kw={"width_ratios": [1.5, 1.0]}
)

nonzero = eta_map[eta_map["eta"] > 0]
sc = ax_a.scatter(
    nonzero["omega_cm1"], nonzero["eta"], c=nonzero["T_K"], cmap="viridis",
    s=4, alpha=0.5, linewidths=0, rasterized=True,
)
ax_a.axhline(0.5, ls="--", color="tab:red", lw=1)
ax_a.axhline(1.0, ls="-", color="tab:red", lw=1)
ax_a.annotate(r"$\eta = 0.5$ (narrow-width limit)", xy=(0.98, 0.5),
              xycoords=("axes fraction", "data"), xytext=(0, 3),
              textcoords="offset points", ha="right", fontsize=8, color="tab:red")
ax_a.annotate(r"$\eta = 1$ (overdamped)", xy=(0.98, 1.0),
              xycoords=("axes fraction", "data"), xytext=(0, 3),
              textcoords="offset points", ha="right", fontsize=8, color="tab:red")
ax_a.set_yscale("log")
ax_a.set_ylim(1e-5, 3.0)
ax_a.set_xlabel(r"$\Omega$ (cm$^{-1}$)")
ax_a.set_ylabel(r"$\eta = \Gamma/\Omega$")
fig.colorbar(sc, ax=ax_a, label="$T$ (K)", pad=0.02)
pubstyle.panel_label(ax_a, "a")

ax_b.plot(val["T_K"], val["eta_max"], "o-", ms=3.5, color="tab:blue",
          label=r"zone max $\eta$")
ax_b.plot(val["T_K"], val["eta_median"], "s-", ms=3.5, color="tab:cyan",
          label=r"zone median $\eta$")
ax_b.plot(val["T_K"], val["eta0_computed"], "^-", ms=3.5, color="tab:orange",
          label=r"$\eta_0$ soft mode (computed $\Omega_0$)")
ax_b.plot(val["T_K"], val["eta0_vogt_anchored"], "v--", ms=4, color="tab:red",
          label=r"$\eta_0$ soft mode (Vogt-anchored)")
ax_b.plot(sens["T_K"], sens["eta_max_rescaled"], "x:", ms=5, color="0.35",
          label=r"zone max, soft branch Vogt-rescaled")
ax_b.axhline(0.5, ls="--", color="tab:red", lw=1)
pubstyle.mark_phase_boundary(ax_b)
ax_b.set_yscale("log")
ax_b.set_ylim(1e-4, 1.0)
ax_b.set_xlabel("$T$ (K)")
ax_b.set_ylabel(r"$\eta$")
ax_b.legend(fontsize=7, loc="lower right")
ax_b.annotate(
    r"$\eta > 0.5$ fraction $= 0$" + "\n(count- and $\\kappa$-weighted,\nall $T$)",
    xy=(0.03, 0.97), xycoords="axes fraction", ha="left", va="top", fontsize=8,
    bbox={"boxstyle": "round", "fc": "white", "ec": "0.7"},
)
pubstyle.panel_label(ax_b, "b")

fig.tight_layout()
paths = pubstyle.save(fig, "f2_eta_map")
print("wrote:", *paths)