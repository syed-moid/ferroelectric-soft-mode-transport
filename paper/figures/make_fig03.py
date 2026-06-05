"""Figure 3 -- BTO hybrid negative result, side-by-side 2-panel (a = SCPH vs hybrid + delta, b = eta).

Verbatim reproduction of the NB 09 plotting (same scalar-SCPH calibration A=208.28
from the preserved snapshot, same unchanged hybrid CSVs, same colors: SCPH tab:red,
hybrid gray 0.4, delta tab:purple, HAS tab:red). NB 09's s2 (Omega^2 overlay + delta)
is assembled into panel (a) on a twin axis; s3 (eta) is panel (b).

Surgical edits only: titles removed; disclaimer dropped entirely (it lives in the
caption); margin "(a)"/"(b)". Data annotations the caption references (+2806%, the
eta=0.5 crossing) are kept.
"""

from __future__ import annotations

import glob
import sys

import pubstyle

ROOT = pubstyle.add_fekappa_to_path()
ARCHIVE = ROOT / "archive" / "scph_calibration"
if str(ARCHIVE) not in sys.path:
    sys.path.insert(0, str(ARCHIVE))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from fekappa.plots import set_publication_style  # noqa: E402
import scph_calibration_snapshot as scph_old  # noqa: E402

set_publication_style()

rows = []
for f in sorted(glob.glob(str(ROOT / "data" / "raw" / "BTO" / "phonon_modes_*K.csv"))):
    d = pd.read_csv(f)
    so = d[d["mode_type"] == "polar_soft"].iloc[0]
    rows.append(dict(T_K=float(so["temperature_K"]), omega_meV=float(so["omega_Q_meV"]),
                     gamma_hwhm_meV=float(so["Gamma_Q_HWHM_meV"])))
hyb = pd.DataFrame(rows).sort_values("T_K").reset_index(drop=True)
hyb["omega_sq_meV2"] = hyb["omega_meV"] ** 2
hyb["eta"] = hyb["gamma_hwhm_meV"] / hyb["omega_meV"]
E0_SQ, A = -424.0, scph_old.calibrate_A_from_ehsan(-424.0, 21.08, 500.0)
hyb["scph_sq_meV2"] = [scph_old.scph_iterate(E0_SQ, A, float(T))[0] for T in hyb["T_K"]]
hyb["delta_pct"] = (hyb["scph_sq_meV2"] - hyb["omega_sq_meV2"]) / hyb["omega_sq_meV2"] * 100.0
d400 = float(hyb.loc[hyb["T_K"] == 400, "delta_pct"].iloc[0])
ETA_BTO_HAS = 2.27

fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.2, 3.8))

# --- (a) Omega^2: SCPH vs hybrid, with delta on a twin axis (NB 09 s2) -----
Tg = np.linspace(100, 650, 60)
curve = scph_old.scph_curve(Tg, E0_SQ, A)
axA.plot(Tg, curve["E_sq_meV2"], "-", color="tab:red", label="scalar SCPH (Ehsan-calib.)")
axA.plot(hyb["T_K"], hyb["omega_sq_meV2"], "o", color="0.4", label="hybrid reconstruction")
axA.set_xlabel("T (K)"); axA.set_ylabel(r"$\Omega^2$ (meV$^2$)"); axA.legend(fontsize=7)
axd = axA.twinx()
axd.plot(hyb["T_K"], hyb["delta_pct"], "-s", color="tab:purple")
axd.annotate(f"+{d400:.0f}%", (400, d400), xytext=(8, -4), textcoords="offset points", color="tab:purple")
axd.set_ylabel(r"$\delta = (\Omega^2_{\rm SCPH}-\Omega^2_{\rm data})/\Omega^2_{\rm data}$  (%)",
               color="tab:purple")
axd.tick_params(axis="y", colors="tab:purple")
pubstyle.panel_label(axA, "a")

# --- (b) eta(T) crossing 0.5 (NB 09 s3) -----------------------------------
axB.axhspan(0.5, 3.0, color="0.92", label=r"overdamped ($\eta>0.5$)")
axB.axhline(0.5, ls="--", color="0.5", lw=1)
axB.plot(hyb["T_K"], hyb["eta"], "-o", color="0.4", label=r"hybrid $\eta=\Gamma/\Omega$")
axB.axhline(ETA_BTO_HAS, color="tab:red", lw=1.5, label=r"HAS cubic $\eta_{\rm BTO}\approx2.27$ (verified)")
axB.set_xlabel("T (K)"); axB.set_ylabel(r"$\eta = \Gamma_{\rm HWHM}/\Omega$")
axB.legend(fontsize=7)
pubstyle.panel_label(axB, "b")

fig.tight_layout()
pdf, png = pubstyle.save(fig, "fig03")

print("MANIFEST  fig03  (verbatim NB09 s2+s3, side-by-side 2-panel)")
print(f"  outputs : {pdf} | {png}")
print(f"  calibration: A={A:.2f} meV^2 (Ehsan); delta(400K)={d400:+.0f}%")
print("  REMOVED:")
print('    (a) title: "SCPH vs hybrid (incommensurate)" + (delta subpanel) "Residual blow-up [preserved negative result]"')
print('    (b) title: "BTO soft mode is overdamped (structural, not poorly-fit)"')
print("    disclaimer footer: DROPPED entirely (now in caption)")
print("  KEPT (verbatim): NB09 colors (SCPH tab:red, hybrid 0.4, delta tab:purple, HAS tab:red),")
print("    +2806% annotation, eta=0.5 line + overdamped span, all legend entries.")
print("  ADDED: margin (a)/(b); s2's two subpanels assembled onto one twin-axis panel (a).")
