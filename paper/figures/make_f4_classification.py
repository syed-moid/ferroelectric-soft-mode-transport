"""Figure F4 -- damping-axis classification, updated with the computed STO point.

Based on the v5 Fig. 4 cluster strip (make_fig04.py). The literature-based STO
point remains PRIMARY (apples-to-apples with the KTO/BTO literature points);
the computed zone-center eta_0(300 K) from the SCPH-coupled dataset is added
alongside as an open marker with a connector. The offset follows the verified
chain: computed ~ Tadano fixed-cell SCPH > Tadano +expansion > experiment
(METHODS_AND_CAVEATS.md C8).

Data: fekappa.eta (Vogt literature values), HAS Table I (BTO), and
data/processed/sto_alamode/eta_validity_T.csv (computed point).
"""

from __future__ import annotations

import pubstyle

ROOT = pubstyle.add_fekappa_to_path()

import matplotlib.pyplot as plt
import pandas as pd

from fekappa import eta as eta_mod
from fekappa import literature_io as lit
from fekappa.plots import set_publication_style

set_publication_style()

_Ts, es = eta_mod.eta_sto_vogt()
_Tk, ek = eta_mod.eta_kto_vogt()
eta_STO, eta_KTO = float(es[-1]), float(ek[-1])
has = lit.load_harada_axe_shirane_1971_mode_parameters()
bto_to = (has[(has["mode"] == "TO") & (has["q_label"] == 0.2)]
          .drop_duplicates("T_K").sort_values("T_K"))
eta_BTO = float((bto_to["gamma_meV"] / bto_to["omega_meV"]).to_numpy().mean())

val = pd.read_csv(
    ROOT / "data" / "processed" / "sto_alamode" / "eta_validity_T.csv", comment="#"
)
eta_STO_comp = float(val.loc[val["T_K"] == 300.0, "eta0_computed"].item())

THR = 0.5
INTERIOR = (0.5, 2.0)
fig, ax = plt.subplots(figsize=(7, 2.6))
ax.axvspan(*INTERIOR, color="0.9", label="empty interior (extrapolated)")
ax.axvline(THR, ls="--", color="0.5", lw=1)
ax.scatter([eta_STO, eta_KTO], [1, 1], s=90, color="tab:blue", zorder=3,
           label="QPE (underdamped, literature)")
ax.scatter([eta_BTO], [1], s=90, color="tab:red", zorder=3,
           label="BTO (overdamped, literature)")
# Computed STO point (this work), connected to its literature benchmark.
ax.scatter([eta_STO_comp], [1], s=90, facecolors="none", edgecolors="tab:blue",
           linewidths=1.6, zorder=3, label="STO computed (this work)")
ax.annotate(
    "", xy=(eta_STO_comp, 1.045), xytext=(eta_STO, 1.045),
    arrowprops={"arrowstyle": "->", "color": "tab:blue", "lw": 1.0,
                "shrinkA": 2, "shrinkB": 2},
)
ax.annotate("STO", (eta_STO, 1), xytext=(0, 10), textcoords="offset points",
            ha="center")
ax.annotate("STO\n(computed)", (eta_STO_comp, 1), xytext=(0, -12),
            textcoords="offset points", ha="center", va="top", fontsize=8)
ax.annotate("KTaO$_3$", (eta_KTO, 1), xytext=(0, -10), textcoords="offset points",
            ha="center", va="top")
ax.annotate("BTO", (eta_BTO, 1), xytext=(0, 10), textcoords="offset points",
            ha="center")
ax.set_xscale("log")
ax.set_xlim(0.012, 3.5)
ax.set_yticks([])
ax.set_ylim(0.9, 1.15)
ax.set_xlabel(r"$\eta = \Gamma_{\rm HWHM}/\Omega$")
ax.legend(loc="upper center", fontsize=7, ncols=2)
fig.tight_layout()
pdf, png = pubstyle.save(fig, "f4_classification")

print("wrote:", pdf, png)
print(f"eta: STO(lit)={eta_STO:.3f}  STO(computed)={eta_STO_comp:.4f}  "
      f"KTaO3={eta_KTO:.3f}  BTO={eta_BTO:.2f}")