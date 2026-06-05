"""Figure 2 -- STO soft-mode damping Gamma(T): YS vs Vogt (verbatim NB 07 cell).

Surgical edits only: title removed; bottom-right repo text box removed. Everything
else verbatim (three lines, tab: colors, markers, y-label, legend at upper-left).
"""

from __future__ import annotations

import pubstyle

ROOT = pubstyle.add_fekappa_to_path()

import matplotlib.pyplot as plt  # noqa: E402
from fekappa.plots import set_publication_style  # noqa: E402
from fekappa import literature_io as lit  # noqa: E402

set_publication_style()

ysg = lit.load_yamada_shirane_gamma_T().sort_values("T_K")
vgg = lit.load_vogt1995_sto_gamma_T(); vggc = vgg[vgg["phase_mode"] == "cubic_F1u"].sort_values("T_K")
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(ysg["T_K"], ysg["gamma_meV"], marker="o", color="tab:blue", ls="-", label="YS 1969 (neutron, HWHM)")
ax.plot(vggc["T_K"], vggc["gamma_meV"], marker="s", color="tab:red", ls="--", label="Vogt plotted γ (FWHM)")
ax.plot(vggc["T_K"], vggc["gamma_meV"] / 2, marker="D", color="tab:green", ls="-", label="Vogt Γ_HWHM = γ/2 (canonical)")
ax.set_xlabel("Temperature (K)"); ax.set_ylabel(r"$\Gamma$ (meV)")
# SURGICAL: title removed; bottom-right text box removed.
ax.legend(fontsize=8, loc="upper left")
fig.tight_layout()
pdf, png = pubstyle.save(fig, "fig02")

print("MANIFEST  fig02  (verbatim NB07 YS-vs-Vogt cell)")
print(f"  outputs : {pdf} | {png}")
print("  REMOVED:")
print('    title: "Figure — STO soft-mode Γ(T): YS 1969 vs Vogt 1995 (HWHM basis)"')
print('    text box: "@250K: raw γ {N}x; HWHM-basis {N}x\\nYS finite-q + central-peak '
      'contaminated\\nVogt true q=0, central peak (δ) separated"')
print("  MOVED: none (legend stays loc='upper left')")
print("  KEPT (verbatim): 3 lines (YS HWHM, Vogt FWHM, Vogt γ/2), tab: colors, markers, y-label Γ (meV)")
