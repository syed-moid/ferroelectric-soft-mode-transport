"""E1: factor decomposition of the coherence pair-sector blocks (8x8x8, submitted dataset).

anphon's coherence element (conductivity.cpp:692-699) for a pair (s, s') at irreducible q:

  kc = H * L * V,
  H  = 2 w1 w2/(w1+w2) * (C1/w1 + C2/w2)                     heat-capacity factor
  L  = 2 (G1+G2) / (4 (w1-w2)^2 + 4 (G1+G2)^2)                Lorentzian overlap (spectral)
  V  = sum over star of v_ss'^x v_s's^x  (from kc / (H L))     off-diagonal velocity factor

We recover V per pair from the solver output (kc) and the recomputed H, L (G from .result
GAMMA_EACH, C from the anphon Cv formula), check the recomposition, and then
  (1) report per block the kappa_C-weighted means of |w1-w2|, G1+G2, H, and the sum of V;
  (2) factor-freezing (labelled numerical experiment): hold ONE factor of each pair at its
      value at T_ref and recompute the blocks at every T (pairs matched by (ik, branch pair)).

Outputs: data/processed/revision/e1_block_factors.csv, e1_factor_freezing.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

from fekappa.alamode_io import read_result
from fekappa.kappa_decomposition import classify_sector, read_kc_elem

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/alamode_sto"
OUT = ROOT / "data/processed/revision"
NK = 512
KB_CM1 = 0.6950348  # k_B in cm^-1/K
TEMPS = list(range(50, 1001, 50))
SHORT = {"acoustic": "A", "soft_optic": "soft", "hard_optic": "hard"}


def cv(w, T):  # per-mode heat capacity in k_B units (anphon Cv / k_B)
    x = w / (KB_CM1 * T)
    return (x / (2 * np.sinh(0.5 * x))) ** 2


def pair_table(t):
    kc = read_kc_elem(RAW / f"STO_RTA_scph_{t}K.kc_elem")
    kc = kc[(kc.i == kc.j) & (kc.ibranch != kc.jbranch)]
    kc = kc.groupby(["ik", "ibranch", "jbranch"], as_index=False).agg(
        w1=("omega1_cm1", "first"), w2=("omega2_cm1", "first"), kc=("kc_real", "sum"))
    res = read_result(RAW / f"STO_RTA_scph_{t}K.result")
    g = res.set_index(["q_index", "branch"])["gamma_hwhm_cm1"]
    kc["G1"] = g.loc[list(zip(kc.ik, kc.ibranch))].to_numpy()
    kc["G2"] = g.loc[list(zip(kc.ik, kc.jbranch))].to_numpy()
    kc["H"] = 2 * kc.w1 * kc.w2 / (kc.w1 + kc.w2) * (cv(kc.w1, t) / kc.w1 + cv(kc.w2, t) / kc.w2)
    Gs = kc.G1 + kc.G2
    kc["L"] = 2 * Gs / (4 * (kc.w1 - kc.w2) ** 2 + 4 * Gs ** 2)
    kc["V"] = kc.kc / (kc.H * kc.L)
    s1 = np.vectorize(SHORT.get)(classify_sector(kc.ibranch, kc.w1))
    s2 = np.vectorize(SHORT.get)(classify_sector(kc.jbranch, kc.w2))
    kc["block"] = ["/".join(sorted(p)) for p in zip(s1, s2)]
    kc["T_K"] = t
    return kc


P = {t: pair_table(t) for t in TEMPS}

rows = []
for t, d in P.items():
    tot = d.kc.sum() / 3 / NK
    for b, g in d.groupby("block"):
        w = g.kc.clip(lower=0)
        rows.append({
            "T_K": t, "block": b, "kappa_block": g.kc.sum() / 3 / NK, "share_of_kC": g.kc.sum() / 3 / NK / tot,
            "n_pairs": len(g),
            "kw_mean_abs_dw_cm1": np.average((g.w1 - g.w2).abs(), weights=w) if w.sum() > 0 else np.nan,
            "kw_mean_G1pG2_cm1": np.average(g.G1 + g.G2, weights=w) if w.sum() > 0 else np.nan,
            "kw_mean_H": np.average(g.H, weights=w) if w.sum() > 0 else np.nan,
            "sum_V": g.V.sum(),
            "frac_pairs_dw_lt_G": float(((g.w1 - g.w2).abs() < (g.G1 + g.G2)).mean()),
            "top_pair_share_of_block": g.kc.max() / g.kc.sum() if g.kc.sum() > 0 else np.nan,
        })
fac = pd.DataFrame(rows)
fac.to_csv(OUT / "e1_block_factors.csv", index=False, float_format="%.5g")

# factor freezing
frz = []
key = ["ik", "ibranch", "jbranch"]
for t_ref in (300, 500):
    ref = P[t_ref].set_index(key)
    for t, d in P.items():
        d = d.set_index(key)
        common = d.index.intersection(ref.index)
        d, r = d.loc[common], ref.loc[common]
        for factor in ("none", "H", "L", "V"):
            H = r.H if factor == "H" else d.H
            L = r.L if factor == "L" else d.L
            V = r.V if factor == "V" else d.V
            k = (H * L * V).groupby(d.block).sum() / 3 / NK
            for b, v in k.items():
                frz.append({"T_ref": t_ref, "frozen": factor, "T_K": t, "block": b, "kappa_block": v})
frz = pd.DataFrame(frz)
frz.to_csv(OUT / "e1_factor_freezing.csv", index=False, float_format="%.5g")

pd.set_option("display.width", 220)
print(fac[fac.T_K.isin([100, 300, 500, 1000])].round(4).to_string(index=False))
for blk in ("hard/hard", "A/soft", "soft/soft", "hard/soft"):
    print(f"\nfactor freezing, block {blk} (T_ref=500):")
    print(frz[(frz.block == blk) & (frz.T_ref == 500)].pivot(index="T_K", columns="frozen", values="kappa_block").round(4).to_string())


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"e1_block_factors.csv", f"e1_factor_freezing.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/raw/alamode_sto (distributed example outputs; IMPORT_NOTE.md)'])
