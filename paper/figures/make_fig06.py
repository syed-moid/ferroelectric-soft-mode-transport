"""Figure 6 -- soft-mode transport kernel kappa_soft(T) (verbatim NB 08 s2 envelopes).

Verbatim NB 08 s2: STO (tab:blue) and KTaO3 (tab:green) kappa_soft(T) with their
Lambda^2 +/- frac-unc envelopes, same labels (eta / envelope %), recomputed via
fekappa.kappa_kernel.K_soft. Per the Fig-6 spec the BaTiO3 order-of-magnitude point
(tab:red, "BTO (723 K)" label) is included.

Surgical edits only: title removed; "soft-mode sector only" prose note removed
(it's in the caption). Legend kept.
"""

from __future__ import annotations

import warnings

import pubstyle

ROOT = pubstyle.add_fekappa_to_path()

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from fekappa.plots import set_publication_style  # noqa: E402
from fekappa import eta as eta_mod  # noqa: E402
from fekappa import kappa_kernel as kk  # noqa: E402
from fekappa import literature_io as lit  # noqa: E402

set_publication_style()

Ts, es = eta_mod.eta_sto_vogt(); Tk, ek = eta_mod.eta_kto_vogt()
eta_STO, eta_KTO = float(es[-1]), float(ek[-1])
has = lit.load_harada_axe_shirane_1971_mode_parameters()
bto_to = (has[(has["mode"] == "TO") & (has["q_label"] == 0.2)].drop_duplicates("T_K").sort_values("T_K"))
eta_BTO = float((bto_to["gamma_meV"] / bto_to["omega_meV"]).to_numpy().mean())

sto_o = lit.load_vogt1995_sto_omega_T(); sto_o = sto_o[sto_o["phase_mode"] == "cubic_F1u"].sort_values("T_K")
kto_o = lit.load_vogt1995_kto_omega_T(); kto_o = kto_o[kto_o["method"] == "exp_HyperRaman_Vogt1995"].sort_values("T_K")


def envelope(material, T, Om, eta_val):
    m = kk.MATERIALS[material]; f = m["lambda_sq_frac_unc"]; L2 = m["lambda_sq"]
    O2 = np.asarray(Om, float) ** 2
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if isinstance(L2, dict):
            return {d: tuple(kk.K_soft(O2, v * s, eta_val, T, m["V_cell_ang3"]) for s in (1.0, 1 - f, 1 + f))
                    for d, v in L2.items()}
        return tuple(kk.K_soft(O2, L2 * s, eta_val, T, m["V_cell_ang3"]) for s in (1.0, 1 - f, 1 + f))


Tsto = sto_o["T_K"].to_numpy(); Osto = sto_o["Omega_meV"].to_numpy()
Tkto = kto_o["T_K"].to_numpy(); Okto = kto_o["Omega_meV"].to_numpy()
Tbto = bto_to["T_K"].to_numpy(); Obto = bto_to["omega_meV"].to_numpy()
sto_ce, sto_lo, sto_hi = envelope("STO", Tsto, Osto, eta_STO)
kto_ce, kto_lo, kto_hi = envelope("KTO", Tkto, Okto, eta_KTO)
bto_env = envelope("BTO", Tbto, Obto, eta_BTO)

fig, ax = plt.subplots(figsize=(7, 4))
ax.fill_between(Tsto, sto_lo, sto_hi, color="tab:blue", alpha=0.25)
ax.plot(Tsto, sto_ce, "-o", color="tab:blue", ms=3, label=r"STO  ($\eta$=%.3f, $\Lambda^2$=889$\pm$30%%)" % eta_STO)
ax.fill_between(Tkto, kto_lo, kto_hi, color="tab:green", alpha=0.2)
ax.plot(Tkto, kto_ce, "-s", color="tab:green", ms=3, label=r"KTaO$_3$  ($\eta$=%.3f, proxy $\pm$50%%)" % eta_KTO)
# BaTiO3 order-of-magnitude point (Fig-6 spec): soft->stiff bracket at its cubic T.
b_lo = min(v[1][-1] for v in bto_env.values()); b_hi = max(v[2][-1] for v in bto_env.values())
Tb = float(Tbto[-1]); b_mid = float(np.sqrt(b_lo * b_hi))
ax.errorbar([Tb], [b_mid], yerr=[[b_mid - b_lo], [b_hi - b_mid]], fmt="D", color="tab:red",
            ms=6, capsize=3, label="BTO (overdamped, order-of-mag)")
ax.annotate("BTO (723 K)", (Tb, b_mid), xytext=(-6, 16), textcoords="offset points",
            ha="right", color="tab:red", fontsize=8)
ax.set_xlim(0, 800)
ax.set_xlabel("T (K)"); ax.set_ylabel(r"$\kappa_{\rm soft}$ (W m$^{-1}$ K$^{-1}$)")
# SURGICAL: title removed; "soft-mode sector only" prose note removed.
ax.legend(loc="upper right", fontsize=7)
fig.tight_layout()
pdf, png = pubstyle.save(fig, "fig06")

print("MANIFEST  fig06  (verbatim NB08 s2 envelopes + BTO order-of-mag point)")
print(f"  outputs : {pdf} | {png}")
print(f"  values  : STO {sto_ce[-1]:.3f} | KTaO3 {kto_ce[-1]:.3f} | BTO [{b_lo:.4f},{b_hi:.4f}] @{Tb:.0f}K")
print("  REMOVED:")
print('    title: "Soft-mode particle-like κ (K$_P$) -- underdamped QPE"')
print('    prose note: "soft-mode sector only (≪ total ~10 W m⁻¹K⁻¹)"')
print("  KEPT (verbatim): STO tab:blue / KTaO3 tab:green envelopes + eta/Lambda^2 legend labels.")
print("  PER Fig-6 SPEC: BTO order-of-mag point (tab:red) + 'BTO (723 K)' label kept;")
print("    legend at upper-right (empty region; clears the T<300 curves); xlim 0-800 to show BTO.")
