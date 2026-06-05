"""Figure 5 -- cubic BTO directional Lambda^2 (verbatim NB 06 cells s1 + s4).

Surgical edits only: title removed; bottom-right repo text box removed. Everything
else verbatim -- tab:purple HAS data + digitized refit, tab:blue locked A(100)=972,
tab:red stiff [110]/[111]=4750, markers, axis limits, legend at upper-left.
"""

from __future__ import annotations

import pubstyle

ROOT = pubstyle.add_fekappa_to_path()

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from fekappa.plots import set_publication_style  # noqa: E402
from fekappa import literature_io as lit  # noqa: E402

set_publication_style()


def fit_lambda_sq(q_inv_ang, omega_meV):
    q = np.asarray(q_inv_ang, float)
    w = np.asarray(omega_meV, float)
    A = np.vstack([np.ones_like(q), q**2]).T
    (c0, c1), *_ = np.linalg.lstsq(A, w**2, rcond=None)
    return float(np.sqrt(max(c0, 0.0))), float(c1)


cd = lit.load_harada_axe_shirane_1971_cubic_dispersion()
to = cd[cd["branch"] == "TO"].sort_values("q_inv_angstrom")
q = to["q_inv_angstrom"].to_numpy(float)
w = to["omega_meV"].to_numpy(float)
Omega0_fit, Lambda_sq_fit = fit_lambda_sq(q, w)
ns = lit.load_harada_axe_shirane_1971_Lambda_sq()
Lambda_sq_locked = ns.central_meV2_A2

q_ref = np.linspace(0, 0.52, 100)
fig, ax = plt.subplots(figsize=(7.5, 5.2))
ax.scatter(q**2, w**2, s=80, color="tab:purple", zorder=5,
           label="HAS 1971 Fig. 5 soft-TO ([100], cubic 423 K)")
ax.plot(q_ref**2, Omega0_fit**2 + Lambda_sq_fit * q_ref**2, "-", color="tab:purple",
        label=rf"refit: $\Omega_0^2$+{Lambda_sq_fit:.0f}$q^2$ (digitized)")
ax.plot(q_ref**2, 3.0**2 + Lambda_sq_locked * q_ref**2, "--", color="tab:blue",
        label="locked A(100)=972 (HAS text)")
ax.plot(q_ref**2, 3.0**2 + 4750.0 * q_ref**2, ":", color="tab:red",
        label="stiff [110]/[111] A=4750 (HAS text)")
ax.set_xlabel(r"$q^2$  (Å$^{-2}$)")
ax.set_ylabel(r"$\Omega^2$  (meV$^2$)")
# SURGICAL: title removed; bottom-right text box removed.
ax.set_xlim(-0.005, 0.27)
ax.set_ylim(0, 260)
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout()
pdf, png = pubstyle.save(fig, "fig05")

print("MANIFEST  fig05  (verbatim NB06 s1+s4)")
print(f"  outputs : {pdf} | {png}")
print(f"  values  : digitized refit Lambda^2={Lambda_sq_fit:.0f}; locked={Lambda_sq_locked:.0f}; stiff=4750")
print("  REMOVED:")
print('    title: "Figure — cubic BaTiO$_3$ directional $\\Lambda^2$ (HAS 1971)"')
print('    text box: "soft [100] refit {N} vs text 972 ({±N}%)\\nstiff/soft anisotropy ~5x\\n'
      'soft TO OVERDAMPED (eta=Gamma/Omega~2.3)"')
print("  MOVED: legend loc 'upper left' -> 'lower right'")
print("  KEPT (verbatim): tab:purple data+refit, tab:blue locked, tab:red stiff; markers; limits.")
