"""Figure 1 -- soft-mode omega(T) triangulation, stacked 2x1 (a = STO, b = KTaO3).

Verbatim reproduction of the original NB 07 (STO) and NB 07b (KTaO3) triangulation
cells -- same colors, markers, marker sizes, line styles, fonts (set_publication_style),
ticks, axis scales/limits. The ONLY changes are the three surgical edits:
  (a) STO  : title removed; legend moved upper-left -> lower-right (inside); margin "(a)".
  (b) KTaO3: title removed; repo text box removed; legend kept at original lower-right;
             margin "(b)".
The two original single panels are assembled into one 2x1 stack, each panel at its
original size/aspect.
"""

from __future__ import annotations

import pubstyle

ROOT = pubstyle.add_fekappa_to_path()

import matplotlib.pyplot as plt  # noqa: E402
from fekappa.plots import set_publication_style  # noqa: E402
from fekappa import literature_io as lit  # noqa: E402

set_publication_style()

fig, (axA, axB) = plt.subplots(2, 1, figsize=(9, 10.8))

# =========================== (a) SrTiO3 (NB 07 s1) ========================
ax = axA
ta = lit.load_tadano2015_scph_omega_T(); ta1 = ta[ta["mode"] == "TO1"]
vd = lit.load_verdi_sscha_omega_T()
fq = lit.load_fauque2022_omega_T(); fq_to = fq[(fq["mode"] == "TO") & (fq["method"] == "experimental")].sort_values("T_K")
sc = lit.load_schmidt2025_omega_T(); sc16 = sc[(sc["isotope"] == "ST16O") & (sc["method"] == "SSCHA")].sort_values("T_K")
ys = lit.load_yamada_shirane_omega_T(); ys0 = ys[ys["condition_direction"] == "q_0_TO_soft"].sort_values("T_K")
vg = lit.load_vogt1995_sto_omega_T(); vgf = vg[vg["phase_mode"] == "cubic_F1u"].sort_values("T_K")

sty = {"SCPH_with_PM": ("o", "tab:blue", "-", "Tadano SCPH(+PM)"),
       "SCPH_without_PM": ("s", "tab:cyan", "--", "Tadano SCPH(-PM)"),
       "TDEP": ("D", "tab:green", ":", "Tadano TDEP")}
for m, (mk, co, ls, lab) in sty.items():
    sub = ta1[ta1["method"] == m].sort_values("T_K")
    ax.plot(sub["T_K"], sub["Omega_meV"], marker=mk, color=co, ls=ls, ms=6, label=lab)
vm = vd[vd["method"] == "MLFF_SSCHA"].sort_values("T_K")
ax.plot(vm["T_K"], vm["Omega_meV"], marker="x", color="tab:orange", ls="-.", label="Verdi MLFF+SSCHA (shifted -4.93)")
ax.plot(fq_to["T_K"], fq_to["omega_meV"], marker="*", color="black", ls="-", ms=10, label="Fauqué Γ-TO exp (INS)")
ax.plot(sc16["T_K"], sc16["Omega_meV"], marker="P", color="tab:brown", ls="-", ms=6, label="Schmidt ST¹⁶O SSCHA")
ax.plot(ys0["T_K"], ys0["omega_meV"], marker="o", color="dimgray", ls="-", ms=5, mfc="none", label="Yamada-Shirane exp (Table II)")
ax.plot(vgf["T_K"], vgf["Omega_meV"], marker="D", color="tab:olive", ls="-", ms=5, label="Vogt 1995 HR exp (cubic F1u)")
ax.set_xlim(0, 1450); ax.set_xlabel("Temperature (K)"); ax.set_ylabel(r"$\hbar\Omega_{\rm TO1}$ (meV)")
ax.legend(fontsize=7.5, loc="lower right", ncol=2)   # SURGICAL: moved from loc="upper left"
pubstyle.panel_label(ax, "a")                          # SURGICAL: margin label

# =========================== (b) KTaO3 (NB 07b s1) =======================
ax = axB
mei = lit.load_meier2022_omega_T()
ran = lit.load_ranalli2023_omega_T()
vog = lit.load_vogt1995_kto_omega_T(); vog95 = vog[vog["method"] == "exp_HyperRaman_Vogt1995"].sort_values("T_K")

for m, st in {"Calculated": ("o", "tab:blue", "-", "Meier SSCHA+MLIP (calc)"),
              "Shirane_exp": ("s", "tab:red", "", "Meier: Shirane INS exp")}.items():
    s = mei[mei["method"] == m].sort_values("T_K"); ax.plot(s["T_K"], s["Omega_meV"], marker=st[0], color=st[1], ls=st[2], ms=6, label=st[3])
rm = ran[ran["method"] == "MLFF+SSCHA"].sort_values("T_K"); rx = ran[ran["method"] == "expt"].sort_values("T_K")
ax.plot(rm["T_K"], rm["Omega_meV"], marker="^", color="tab:green", ls="--", ms=6, label="Ranalli MLFF+SSCHA (shifted -4.30)")
ax.plot(rx["T_K"], rx["Omega_meV"], marker="v", color="tab:purple", ls="", ms=6, mfc="none", label="Ranalli INS exp")
ax.plot(vog95["T_K"], vog95["Omega_meV"], marker="D", color="tab:olive", ls="-", ms=5, label="Vogt 1995 HR exp (q=0)")
ax.set_xlim(0, 310); ax.set_xlabel("Temperature (K)"); ax.set_ylabel(r"$\hbar\Omega_{\rm TO}$ (meV)")
# SURGICAL: title removed; repo tie-breaker text box removed.
ax.legend(fontsize=7.5, loc="lower right")             # kept at original location
pubstyle.panel_label(ax, "b")                          # SURGICAL: margin label

fig.tight_layout()
pdf, png = pubstyle.save(fig, "fig01")

# ------------------------- removed / moved list ---------------------------
print("MANIFEST  fig01  (verbatim originals NB07 s1 + NB07b s1, stacked 2x1)")
print(f"  outputs : {pdf} | {png}")
print("  REMOVED:")
print('    (a) title: "Figure — STO TO1 Ω(T): multi-source triangulation (raw, no fits)"')
print('    (b) title: "Figure — KTaO$_3$ TO Ω(T): three-source triangulation + Vogt tie-breaker"')
print('    (b) text box: "Vogt exp (q=0) tie-breaks toward Meier-like high-T freq;\\n'
      'Ranalli -4.30 shift over-applied >150 K. Lambda^2_KTO: all q=0 -> bracket."')
print("  MOVED:")
print('    (a) legend: loc "upper left" (ncol=2) -> "lower right" (ncol=2, fontsize unchanged)')
print("  KEPT (verbatim): all colors/markers/sizes/linestyles, both legends' entries")
print('    incl. "(shifted -4.93)" / "(shifted -4.30)"; (b) legend stays at original "lower right";')
print("    all axis labels, limits, set_publication_style fonts/ticks.")
print("  ADDED: margin (a)/(b) labels; 2x1 assembly only.")
