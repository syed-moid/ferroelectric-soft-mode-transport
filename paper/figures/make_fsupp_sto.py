"""Supplementary figures S1-S3 for the computed SrTiO3 dataset.

S1 -- mesh convergence: kappa_P and kappa_C at 300 K on 8x8x8 vs 12x12x12,
      plus the eta-map summary shift.
S2 -- small-q Lambda^2 fit: Omega^2 vs q^2 at 300 K (fit line + points from the
      committed path table's Gamma-X segment) and Lambda^2(T) with Omega_0(T).
S3 -- band-resolved tau(300 K) validation: bare vs SCPH-coupled band medians
      (p10-p90 whiskers) against the digitized Tadano 2015 Fig. 7 ranges.

Data: data/processed/sto_alamode/{mesh_convergence,lambda_sq_T,
dispersion_path_scph300}.csv + data/processed/alamode_crosscheck/*.csv.
"""

from __future__ import annotations

import re

import pubstyle

ROOT = pubstyle.add_fekappa_to_path()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fekappa.plots import set_publication_style

set_publication_style()

TIDY = ROOT / "data" / "processed" / "sto_alamode"
CROSS = ROOT / "data" / "processed" / "alamode_crosscheck"
CM1_TO_MEV = 1.0 / 8.06554

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12.6, 3.6))

# --- S1: mesh convergence ---------------------------------------------------
mesh = pd.read_csv(TIDY / "mesh_convergence.csv", comment="#")
x = np.arange(len(mesh))
w = 0.36
ax1.bar(x - w / 2, mesh["kappa_W_mK"], w, color="tab:blue", label=r"$\kappa_P$")
ax1.bar(x + w / 2, mesh["kappa_coherent_W_mK"], w, color="tab:red",
        label=r"$\kappa_C$")
for xi, r in zip(x, mesh.itertuples()):
    ax1.annotate(f"{r.kappa_W_mK:.2f}", (xi - w / 2, r.kappa_W_mK), ha="center",
                 va="bottom", fontsize=8)
    ax1.annotate(f"{r.kappa_coherent_W_mK:.2f}", (xi + w / 2, r.kappa_coherent_W_mK),
                 ha="center", va="bottom", fontsize=8)
ax1.set_xticks(x, mesh["mesh"])
ax1.set_ylabel(r"$\kappa$(300 K) (W m$^{-1}$ K$^{-1}$)")
ax1.set_ylim(0, 12.5)
med = mesh["median_eta"].tolist()
ax1.annotate(
    rf"median $\eta$: {med[0]:.4f} $\to$ {med[1]:.4f}",
    xy=(0.5, 0.97), xycoords="axes fraction", ha="center", va="top", fontsize=8,
)
ax1.legend(fontsize=8, loc="center right")
pubstyle.panel_label(ax1, "a")

# --- S2: Lambda^2 fit --------------------------------------------------------
lam = pd.read_csv(TIDY / "lambda_sq_T.csv", comment="#")
disp = pd.read_csv(TIDY / "dispersion_path_scph300.csv", comment="#")
r300 = lam[lam["T_K"] == 300.0].iloc[0]
seg = disp[disp["q_ang1"] <= 0.25]
omega_to1 = 0.5 * (seg["branch_4"] + seg["branch_5"]) * CM1_TO_MEV
ax2.plot(seg["q_ang1"] ** 2, omega_to1**2, "o", ms=4, color="tab:blue",
         label="TO1 doublet (SCPH, 300 K)")
qq = np.linspace(0, 0.25**2, 50)
ax2.plot(qq, r300["Omega0_meV"] ** 2 + r300["lambda_sq_meV2_ang2"] * qq, "-",
         color="tab:orange",
         label=rf"fit: $\Lambda^2$ = {r300['lambda_sq_meV2_ang2']:.0f} meV$^2\,$Å$^2$")
ax2.axvline(r300["qmax_fit_ang1"] ** 2, ls=":", color="0.5", lw=1)
ax2.annotate(r"$q_{\rm max}$ fit", xy=(r300["qmax_fit_ang1"] ** 2, 0.02),
             xycoords=("data", "axes fraction"), fontsize=8, rotation=90,
             ha="right", va="bottom", color="0.4")
ax2.set_xlabel(r"$q^2$ (Å$^{-2}$)")
ax2.set_ylabel(r"$\Omega^2$ (meV$^2$)")
ax2.legend(fontsize=8, loc="upper left")
ax2in = ax2.inset_axes([0.55, 0.14, 0.42, 0.38])
ax2in.plot(lam["T_K"], lam["lambda_sq_meV2_ang2"], "s-", ms=3, color="tab:purple")
ax2in.set_xlabel("$T$ (K)", fontsize=7)
ax2in.set_ylabel(r"$\Lambda^2$", fontsize=7)
ax2in.tick_params(labelsize=7)
pubstyle.panel_label(ax2, "b")

# --- S3: band-resolved tau validation ----------------------------------------
bare = pd.read_csv(CROSS / "alamode_sto_tau300K_bands_bare.csv", comment="#")
scph = pd.read_csv(CROSS / "alamode_sto_tau300K_bands_scph.csv", comment="#")


def paper_range(text: str):
    """Leading 'a-b' (ps) from the digitized-range annotation, else None."""
    m = re.match(r"\s*([\d.]+)-([\d.]+)\s*ps", str(text))
    return (float(m.group(1)), float(m.group(2))) if m else None


for i, (df, color, label) in enumerate(
    ((bare, "0.5", "bare-RTA"), (scph, "tab:blue", "SCPH-coupled"))
):
    ok = df[df["n_modes"] > 0]
    centers = 0.5 * (ok["band_lo_cm1"] + ok["band_hi_cm1"]) + (i - 0.5) * 12
    yerr = np.array([ok["tau_median_ps"] - ok["tau_p10_ps"],
                     ok["tau_p90_ps"] - ok["tau_median_ps"]])
    ax3.errorbar(centers, ok["tau_median_ps"], yerr=yerr, fmt="o", ms=4,
                 color=color, capsize=2.5, lw=1, label=label)
drawn = False
for _, r in scph.iterrows():
    rng = paper_range(r["paper_range_ps"])
    if rng is None:
        continue
    ax3.fill_between([r["band_lo_cm1"], r["band_hi_cm1"]], rng[0], rng[1],
                     color="tab:green", alpha=0.18, lw=0,
                     label=None if drawn else "Tadano 2015 Fig. 7 range")
    drawn = True
ax3.set_yscale("log")
ax3.set_xlabel(r"$\Omega$ band (cm$^{-1}$)")
ax3.set_ylabel(r"$\tau$(300 K) (ps)")
ax3.legend(fontsize=8, loc="upper right")
pubstyle.panel_label(ax3, "c")

fig.tight_layout()
paths = pubstyle.save(fig, "fsupp_sto")
print("wrote:", *paths)