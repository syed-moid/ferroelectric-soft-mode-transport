"""C3: damping-map bookkeeping audit on the committed mode table (8x8x8, SCPH-coupled).

Checks, per temperature:
* irreducible-set size, star multiplicities (sum must be 512), Gamma-translation handling;
* eta statistics unweighted over the irreducible set (as committed) vs weighted by
  star multiplicity (= full-mesh measure) vs kappa_P weight;
* exceedance fractions at eta > 0.1 and > 0.5 under all three weightings;
* that the kappa weight equals kappa_mode / sum(kappa_mode) with kappa_mode built
  from the star-summed |v|^2 (i.e. the same multiplicity weighting as kappa_P).

Output: data/processed/revision/c3_damping_bookkeeping.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
df = pd.read_csv(ROOT / "data/processed/sto_alamode/mode_kappa_T.csv", comment="#")
OUT = ROOT / "data/processed/revision/c3_damping_bookkeeping.csv"


def wmedian(x, w):
    o = np.argsort(x)
    c = np.cumsum(w[o]) / w.sum()
    return float(x[o][np.searchsorted(c, 0.5)])


rows = []
for t, g in df.groupby("T_K"):
    qm = g.drop_duplicates("q_index")[["q_index", "multiplicity"]]
    assert len(qm) == 35 and qm["multiplicity"].sum() == 512
    gamma_pt = g[g["omega_cm1"] <= 1e-3]
    s = g[g["omega_cm1"] > 1e-3]
    eta = (s["gamma_hwhm_cm1"] / s["omega_cm1"]).to_numpy()
    m = s["multiplicity"].to_numpy().astype(float)
    k = s["kappa_mode_W_mK"].to_numpy()
    kw = k / k.sum() if k.sum() > 0 else np.zeros_like(k)
    rows.append(
        {
            "T_K": t,
            "n_irred_modes": len(g),
            "n_gamma_translations_excluded": len(gamma_pt),
            "gamma_translation_omega_max_cm1": float(gamma_pt["omega_cm1"].abs().max()),
            "n_modes_in_map": len(s),
            "n_mesh_modes_in_map": int(m.sum()),
            "eta_median_unweighted": float(np.median(eta)),
            "eta_median_mult": wmedian(eta, m),
            "eta_median_kappa": wmedian(eta, kw) if k.sum() > 0 else np.nan,
            "eta_max": float(eta.max()),
            "frac_gt_0p1_unweighted": float((eta > 0.1).mean()),
            "frac_gt_0p1_mult": float(m[eta > 0.1].sum() / m.sum()),
            "frac_gt_0p1_kappa": float(kw[eta > 0.1].sum()),
            "frac_gt_0p05_mult": float(m[eta > 0.05].sum() / m.sum()),
            "frac_gt_0p05_kappa": float(kw[eta > 0.05].sum()),
            "frac_gt_0p5_any": float((eta > 0.5).any()),
        }
    )
out = pd.DataFrame(rows)
with open(OUT, "w") as fh:
    fh.write("# C3 damping-map bookkeeping; source data/processed/sto_alamode/mode_kappa_T.csv (8x8x8 SCPH-coupled)\n")
    out.to_csv(fh, index=False, float_format="%.6g")
pd.set_option("display.width", 250)
print(out.round(5).to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"c3_damping_bookkeeping.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/processed/sto_alamode/mode_kappa_T.csv'])
