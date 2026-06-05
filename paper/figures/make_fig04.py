"""Figure 4 (MAIN) -- damping-axis classification (verbatim NB 08 cluster cell).

Surgical edit only: title removed. Everything else verbatim -- tab:blue (QPE) /
tab:red (BTO), the log-eta strip, the empty-interior span + label, the eta=0.5 line,
the QPE/BTO legend, and the STO / KTaO3 / BTO point labels. eta values are pulled
from the committed sources (fekappa.eta + HAS Table I), not hardcoded.
"""

from __future__ import annotations

import pubstyle

ROOT = pubstyle.add_fekappa_to_path()

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from fekappa.plots import set_publication_style  # noqa: E402
from fekappa import eta as eta_mod  # noqa: E402
from fekappa import literature_io as lit  # noqa: E402

set_publication_style()

Ts, es = eta_mod.eta_sto_vogt()
Tk, ek = eta_mod.eta_kto_vogt()
eta_STO, eta_KTO = float(es[-1]), float(ek[-1])
has = lit.load_harada_axe_shirane_1971_mode_parameters()
bto_to = (has[(has["mode"] == "TO") & (has["q_label"] == 0.2)]
          .drop_duplicates("T_K").sort_values("T_K"))
eta_BTO = float((bto_to["gamma_meV"] / bto_to["omega_meV"]).to_numpy().mean())

THR = 0.5
INTERIOR = (0.5, 2.0)
fig, ax = plt.subplots(figsize=(7, 2.6))
ax.axvspan(*INTERIOR, color="0.9", label="empty interior (extrapolated)")
ax.axvline(THR, ls="--", color="0.5", lw=1)
ax.scatter([eta_STO, eta_KTO], [1, 1], s=90, color="tab:blue", zorder=3, label="QPE (underdamped)")
ax.scatter([eta_BTO], [1], s=90, color="tab:red", zorder=3, label="BTO (overdamped)")
ax.annotate("STO", (eta_STO, 1), xytext=(0, 10), textcoords="offset points", ha="center")
# KTaO3 label moved BELOW its marker so it no longer overlaps the STO label.
ax.annotate("KTaO$_3$", (eta_KTO, 1), xytext=(0, -10), textcoords="offset points", ha="center", va="top")
ax.annotate("BTO", (eta_BTO, 1), xytext=(0, 10), textcoords="offset points", ha="center")
ax.set_xscale("log"); ax.set_yticks([]); ax.set_ylim(0.95, 1.15)
ax.set_xlabel(r"$\eta = \Gamma_{\rm HWHM}/\Omega$")
# SURGICAL: title removed; legend relocated "lower right" -> "upper center" because the
# original's lower-right legend box overlapped the BTO data point (collision).
ax.legend(loc="upper center", fontsize=7)
fig.tight_layout()
pdf, png = pubstyle.save(fig, "fig04")

print("MANIFEST  fig04  (verbatim NB08 cluster cell; MAIN)")
print(f"  outputs : {pdf} | {png}")
print(f"  eta (committed sources): STO={eta_STO:.3f}  KTaO3={eta_KTO:.3f}  BTO={eta_BTO:.2f}")
print("  REMOVED: title \"Three-material damping cluster\"")
print("  MOVED: legend 'lower right' -> 'upper center' (lower-right box overlapped the BTO point)")
print("  LABEL FIX: KTaO3 point label moved BELOW its marker (STO stays above) to end the overlap.")
print("  KEPT (verbatim): tab:blue/tab:red, log-eta strip, empty-interior span+label,")
print("    eta=0.5 line, QPE/BTO legend, STO/KTaO3/BTO point labels.")
print("  NOTE: the original cell has NO separate 'underdamped/overdamped' text and NO 'eta>~1'")
print("    line (those were in the earlier restyled draft). Left out to stay verbatim --")
print("    say if you want them added back as region/threshold labels.")
