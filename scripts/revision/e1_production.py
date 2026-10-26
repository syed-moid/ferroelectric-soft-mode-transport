"""E1 on production data (own SCPH q1 = 12^3, PM): factor decomposition of the coherence pair-sector blocks,
kc = H * L * V per pair (anphon conductivity.cpp:692-699; definitions as in e1_factor_decomposition.py):
  H = 2 w1 w2/(w1+w2) (C1/w1 + C2/w2), L = 2 (G1+G2) / (4 (w1-w2)^2 + 4 (G1+G2)^2), V = kc / (H L).
Meshes 12^3 and 16^3, 300-1000 K (main-text range); sector labels under both rules (frequency rule as submitted;
eigenvector-character rule from mode_table_T.csv.gz). Factor freezing (labelled numerical experiment): one factor of
each pair held at its T_ref value (pairs matched by (ik, branch pair) on the same mesh).
Outputs: data/processed/revision/e1_prod_block_factors.csv, e1_prod_factor_freezing.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

from fekappa.alamode_io import read_result
from fekappa.kappa_decomposition import classify_sector, read_kc_elem

ROOT = Path(__file__).resolve().parents[2]
ARCH = ROOT / "data/raw/alamode_sto/revision/pm_q12"
OUT = ROOT / "data/processed/revision"
KB_CM1 = 0.6950348
TEMPS = list(range(300, 1001, 50))
SHORT = {"acoustic": "A", "soft_optic": "soft", "hard_optic": "hard"}
MODES = pd.read_csv(ROOT / "data/processed/sto_production/mode_table_T.csv.gz",
                    usecols=["q_index", "branch", "T_K", "mesh", "sector_char"])


def cv(w, T):
    x = w / (KB_CM1 * T)
    return (x / (2 * np.sinh(0.5 * x))) ** 2


def pair_table(T, mesh):
    job = f"rta_pm_q12_T{T}_m{mesh}"
    kc = read_kc_elem(ARCH / job / f"{job}.kc_elem")
    kc = kc[(kc.i == kc.j) & (kc.ibranch != kc.jbranch)]
    kc = kc.groupby(["ik", "ibranch", "jbranch"], as_index=False).agg(
        w1=("omega1_cm1", "first"), w2=("omega2_cm1", "first"), kc=("kc_real", "sum"))
    res = read_result(ARCH / job / f"{job}.result")
    g = res.set_index(["q_index", "branch"])["gamma_hwhm_cm1"]
    kc["G1"] = g.loc[list(zip(kc.ik, kc.ibranch))].to_numpy()
    kc["G2"] = g.loc[list(zip(kc.ik, kc.jbranch))].to_numpy()
    kc["H"] = 2 * kc.w1 * kc.w2 / (kc.w1 + kc.w2) * (cv(kc.w1, T) / kc.w1 + cv(kc.w2, T) / kc.w2)
    Gs = kc.G1 + kc.G2
    kc["L"] = 2 * Gs / (4 * (kc.w1 - kc.w2) ** 2 + 4 * Gs ** 2)
    kc["V"] = kc.kc / (kc.H * kc.L)
    s1 = np.vectorize(SHORT.get)(classify_sector(kc.ibranch, kc.w1))
    s2 = np.vectorize(SHORT.get)(classify_sector(kc.jbranch, kc.w2))
    kc["block_freq"] = ["/".join(sorted(p)) for p in zip(s1, s2)]
    lab = MODES[(MODES.T_K == T) & (MODES.mesh == mesh)].set_index(["q_index", "branch"]).sector_char
    c1 = lab.loc[list(zip(kc.ik, kc.ibranch))].to_numpy()
    c2 = lab.loc[list(zip(kc.ik, kc.jbranch))].to_numpy()
    kc["block_char"] = ["/".join(sorted(p)) for p in zip(c1, c2)]
    kc["T_K"], kc["mesh"] = T, mesh
    return kc


rows, frz = [], []
for mesh in (12, 16):
    nk = mesh ** 3
    P = {T: pair_table(T, mesh) for T in TEMPS}
    for T, d in P.items():
        tot = d.kc.sum() / 3 / nk
        for rule in ("freq", "char"):
            for b, g in d.groupby(f"block_{rule}"):
                w = g.kc.clip(lower=0)
                ok = w.sum() > 0
                rows.append({"mesh": mesh, "T_K": T, "labels": rule, "block": b,
                             "kappa_block": g.kc.sum() / 3 / nk, "share_of_kC": g.kc.sum() / 3 / nk / tot,
                             "n_pairs": len(g),
                             "kw_mean_abs_dw_cm1": np.average((g.w1 - g.w2).abs(), weights=w) if ok else np.nan,
                             "kw_mean_G1pG2_cm1": np.average(g.G1 + g.G2, weights=w) if ok else np.nan,
                             "kw_mean_H": np.average(g.H, weights=w) if ok else np.nan,
                             "sum_V": g.V.sum(),
                             "frac_pairs_dw_lt_G": float(((g.w1 - g.w2).abs() < (g.G1 + g.G2)).mean())})
    key = ["ik", "ibranch", "jbranch"]
    for t_ref in (300, 600):
        ref = P[t_ref].set_index(key)
        for T, d in P.items():
            d = d.set_index(key)
            common = d.index.intersection(ref.index)
            d, r = d.loc[common], ref.loc[common]
            for factor in ("none", "H", "L", "V"):
                H = r.H if factor == "H" else d.H
                L = r.L if factor == "L" else d.L
                V = r.V if factor == "V" else d.V
                for rule in ("freq", "char"):
                    k = (H * L * V).groupby(d[f"block_{rule}"]).sum() / 3 / nk
                    for b, v in k.items():
                        frz.append({"mesh": mesh, "T_ref": t_ref, "frozen": factor, "labels": rule, "T_K": T,
                                    "block": b, "kappa_block": v})
fac, frz = pd.DataFrame(rows), pd.DataFrame(frz)
fac.to_csv(OUT / "e1_prod_block_factors.csv", index=False, float_format="%.5g")
frz.to_csv(OUT / "e1_prod_factor_freezing.csv", index=False, float_format="%.5g")
pd.set_option("display.width", 220)
s = fac[(fac.mesh == 16) & (fac.T_K.isin([300, 600, 1000]))]
print(s.round(4).to_string(index=False))
for rule in ("freq", "char"):
    x = frz[(frz.mesh == 16) & (frz.T_ref == 300) & (frz.labels == rule) & (frz.block == "hard/hard")]
    print(f"\n16^3 hard/hard, {rule} labels, factors frozen at 300 K:")
    print(x.pivot(index="T_K", columns="frozen", values="kappa_block").round(4).to_string())


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"e1_prod_block_factors.csv", f"e1_prod_factor_freezing.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['zenodo:pm_q12_rta.tar.gz (rta_pm_q12_T*_m16 .kc_elem)', 'data/processed/sto_production/mode_table_T.csv.gz'])
