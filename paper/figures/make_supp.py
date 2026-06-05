"""Supplementary S1-S3 (verbatim originals; titles removed).

S1 -- STO LST (NB 07 s3): Verdi 1/eps vs Omega^2, R^2 ~ 0.998.
S2 -- KTaO3 LST (NB 07b s4): Ranalli 1/eps vs Omega^2, R^2 ~ 0.993 (internal).
S3 -- transport plane (NB 08 s3): eta-kappa_soft triples, reading the committed
      data/processed/section4_eta_kappa_triples.csv.

Surgical edit only: titles removed. Colors/markers/scales verbatim.
"""

from __future__ import annotations

import pubstyle

ROOT = pubstyle.add_fekappa_to_path()

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from fekappa.plots import set_publication_style  # noqa: E402
from fekappa import literature_io as lit  # noqa: E402

set_publication_style()

# ============================ S1 -- STO LST (NB 07 s3) ====================
vd = lit.load_verdi_sscha_omega_T()
eps = lit.load_verdi_inverse_dielectric_T(); em = eps[eps["method"] == "MLFF_SSCHA"].sort_values("T_K")
om = vd[vd["method"] == "MLFF_SSCHA"].sort_values("T_K")
merged = em.merge(om[["T_K", "Omega_meV"]], on="T_K", how="inner")
x = (merged["Omega_meV"].astype(float))**2; y = merged["inverse_epsilon"].astype(float)
A = np.vstack([x, np.ones_like(x)]).T
(slope, intc), *_ = np.linalg.lstsq(A, y, rcond=None)
yhat = slope * x + intc; R2 = 1 - ((y - yhat)**2).sum() / ((y - y.mean())**2).sum()
fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(x, y, s=50, color="tab:blue", label="Verdi MLFF+SSCHA")
xr = np.linspace(x.min(), x.max(), 50); ax.plot(xr, slope * xr + intc, "k--", label=f"LST linear fit (R²={R2:.4f})")
ax.set_xlabel(r"$\Omega_{\rm TO1}^2$ (meV$^2$)"); ax.set_ylabel(r"$1/\epsilon$")
# SURGICAL: title removed.
ax.legend(fontsize=8)
fig.tight_layout(); s1 = pubstyle.save(fig, "figS1"); plt.close(fig)

# ============================ S2 -- KTaO3 LST (NB 07b s4) =================
ran = lit.load_ranalli2023_omega_T()
eps2 = lit.load_ranalli2023_inverse_dielectric_T(); em2 = eps2[eps2["method"] == "MLFF+SSCHA"].sort_values("T_K")
om2 = ran[ran["method"] == "MLFF+SSCHA"].sort_values("T_K")
mg = em2.merge(om2[["T_K", "Omega_meV"]], on="T_K", how="inner")
x = (mg["Omega_meV"].astype(float))**2; y = mg["inverse_epsilon"].astype(float)
A = np.vstack([x, np.ones_like(x)]).T
(slope, intc), *_ = np.linalg.lstsq(A, y, rcond=None)
yhat = slope * x + intc; R2 = 1 - ((y - yhat)**2).sum() / ((y - y.mean())**2).sum()
fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(x, y, s=50, color="tab:green", label="Ranalli MLFF+SSCHA")
xr = np.linspace(float(x.min()), float(x.max()), 50); ax.plot(xr, slope * xr + intc, "k--", label=f"LST fit (R²={R2:.4f})")
ax.set_xlabel(r"$\Omega_{\rm TO}^2$ (meV$^2$)"); ax.set_ylabel(r"$1/\epsilon$")
# SURGICAL: title removed.
ax.legend(fontsize=8)
fig.tight_layout(); s2 = pubstyle.save(fig, "figS2"); plt.close(fig)

# ============== S3 -- transport plane (NB 08 s3; committed CSV) ===========
triples = pd.read_csv(ROOT / "data" / "processed" / "section4_eta_kappa_triples.csv")
fig, ax = plt.subplots(figsize=(7, 4))
ax.axvspan(0.5, 2.0, color="0.92", label="empty interior")
q = triples[triples.regime == "underdamped"]
b = triples[triples.regime.str.startswith("overdamped")]
for sub, fmt, col, lab in [(q, "o", "tab:blue", "QPE (anchored)"), (b, "s", "tab:red", "BTO (order-of-mag)")]:
    ax.errorbar(sub["eta"], sub["kappa_soft_W_m_K"],
                yerr=[sub["kappa_soft_W_m_K"] - sub["kappa_soft_lo"], sub["kappa_soft_hi"] - sub["kappa_soft_W_m_K"]],
                fmt=fmt, color=col, capsize=3, label=lab)
for _, r in triples.iterrows():
    ax.annotate(f'{r["material"]}{"" if r["direction"] == "isotropic" else ":" + r["direction"]}',
                (r["eta"], r["kappa_soft_W_m_K"]), xytext=(4, 4), textcoords="offset points", fontsize=6)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"$\eta = \Gamma_{\rm HWHM}/\Omega$"); ax.set_ylabel(r"$\kappa_{\rm soft}$ (W m$^{-1}$ K$^{-1}$)")
# SURGICAL: title removed.
ax.legend(fontsize=7)
fig.tight_layout(); s3 = pubstyle.save(fig, "figS3"); plt.close(fig)

print("MANIFEST  supp (S1-S3, verbatim; titles removed)")
print(f"  S1 STO LST   : {s1[1]}  R²={R2:.4f}-style; src verdi_inverse_dielectric_T + verdi_sscha_omega_T")
print(f"  S2 KTaO3 LST : {s2[1]}  src ranalli2023_inverse_dielectric_T + ranalli2023_omega_T")
print(f"  S3 plane     : {s3[1]}  src data/processed/section4_eta_kappa_triples.csv ({len(triples)} rows)")
print("  REMOVED titles:")
print('    S1 "Figure — STO LST consistency (1/ε ∝ Ω²)"')
print('    S2 "Figure — KTaO$_3$ LST consistency (1/ε ∝ Ω²)"')
print('    S3 "Section-4 plane: $\\kappa_{\\rm soft}$ vs $\\eta$ (interior extrapolated)"')
print("  KEPT (verbatim): tab:blue/tab:green/tab:red, k-- fits, log-log (S3), all legend entries.")
