"""Short mean free path vs damping ratio on production data (claim (i) of Guo & Guo 2023).
Guo & Guo classify 75.6 % of modes at 300 K as "diffusion" phonons, whose MFP they state to be below 2.8 A (p.7). We
test the temporal criterion eta = Gamma/Omega. Here, per mode of the 16^3 production mesh (own SCPH q1 = 12^3, PM):
  |v| = sqrt(v2_star / multiplicity)  (v2_star is the star-summed |v|^2, src/fekappa/kappa_decomposition.py)
  tau = 1 / (2 Gamma), Gamma = 2 pi c * gamma_hwhm_cm1 (angular units)
  MFP = |v| tau
Reported at 300 and 900 K: the mesh-weighted fraction of modes and the kappa_P share with MFP < 2.8 A (Guo & Guo) and
< d(Ti-O) = a/2 = 1.948 A (a = 3.8963 A, fixed cell), and the eta distribution of those modes.
Output: data/processed/revision/mfp_vs_eta.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PJ = ROOT / "data/processed/revision"
C_CM = 2.99792458e10  # cm/s
D_TIO = 3.8963 / 2
mt = pd.read_csv(ROOT / "data/processed/sto_production/mode_table_T.csv.gz", comment="#")
km = pd.read_csv(ROOT / "data/processed/sto_production/kappa_mesh_T.csv", comment="#")
rows = []
for T in (300, 900):
    x = mt[(mt.T_K == T) & (mt.mesh == 16) & (mt.omega_cm1 > 1e-3)].copy()
    v = np.sqrt(x.v2_star_m2s2 / x.multiplicity)
    tau = 1 / (2 * 2 * np.pi * C_CM * x.gamma_hwhm_cm1)
    x["mfp_A"] = v * tau * 1e10
    w, k = x.multiplicity.astype(float), x.kappa_mode_W_mK
    kP = float(km[(km.T_K == T) & (km.mesh == 16)].kappa_P.iloc[0])
    for lab, cut in (("2.8 A (Guo & Guo)", 2.8), (f"d(Ti-O) = {D_TIO:.3f} A", D_TIO)):
        s = x.mfp_A < cut
        e = x.eta[s]
        rows.append(dict(T_K=T, mesh=16, mfp_cut=lab, frac_modes=w[s].sum() / w.sum(), kappaP_share=k[s].sum() / kP,
                         kappaP_sum_closure=k.sum() / kP - 1, eta_median=float(np.median(np.repeat(e, w[s]))),
                         eta_max=float(e.max()), frac_eta_gt_0p1=w[s & (x.eta > 0.1)].sum() / w[s].sum(),
                         frac_eta_gt_0p5=w[s & (x.eta > 0.5)].sum() / w[s].sum(),
                         omega_range_cm1=f"{x.omega_cm1[s].min():.0f}-{x.omega_cm1[s].max():.0f}",
                         eta_max_all_modes=float(x.eta.max())))
out = pd.DataFrame(rows)
out.to_csv(PJ / "mfp_vs_eta.csv", index=False, float_format="%.4g")
pd.set_option("display.width", 250)
print(out.to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"mfp_vs_eta.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/processed/sto_production/mode_table_T.csv.gz', 'data/processed/sto_production/kappa_mesh_T.csv'])
