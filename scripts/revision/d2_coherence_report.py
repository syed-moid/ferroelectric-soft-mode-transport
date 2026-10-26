"""D2: coherence-channel convergence report on production data (q1 = 12^3, PM): 300 K on 12/15/16/17/20^3 (N2 adds
the odd meshes), 500 K on 12/16/20^3 (N5) and 800 K on 12/16/20^3.

For each (T, mesh):
  * kappa_C (pair sum, closes on the solver value) and the six frequency-rule sector blocks;
  * the 20 largest pair contributions (unordered pair, summed over Cartesian diagonal and both orders)
    with q (fractional), branches, frequencies, |dw|, G1+G2, and r = (G1+G2)/dnu_mesh;
  * hypothesis test: kappa_C converges once the pair linewidths exceed the mesh's
    frequency resolution. Per mode, dnu_mesh = |v| / (a n c) (cm^-1): the frequency change across one
    mesh spacing dq = 2 pi/(a n) at group velocity |v| (rms over the q star). Per pair, dnu = mean of the
    two modes. Reported: kappa_C-weighted fraction of pairs with r > 1, and the weighted median r.
Outputs: data/processed/revision/d2_top_pairs.csv, d2_summary.csv, data/processed/sto_production/convergence_report_kC.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

from fekappa.alamode_io import read_result, read_result_velocities
from fekappa.kappa_decomposition import classify_sector, read_kc_elem

ROOT = Path(__file__).resolve().parents[2]
ARCH = ROOT / "data/raw/alamode_sto/revision/pm_q12"
OUT = ROOT / "data/processed/revision"
A_M = 7.363 * 0.529177210903e-10
C_CM = 2.99792458e10
SHORT = {"acoustic": "A", "soft_optic": "soft", "hard_optic": "hard"}


def analyse(T, mesh):
    job = f"rta_pm_q12_T{T}_m{mesh}"
    pre = ARCH / job / job
    nk = mesh ** 3
    kc = read_kc_elem(Path(f"{pre}.kc_elem"))
    kc = kc[(kc.i == kc.j) & (kc.ibranch != kc.jbranch)]
    res = read_result(Path(f"{pre}.result"))
    g = res.set_index(["q_index", "branch"])["gamma_hwhm_cm1"]
    qf = res.drop_duplicates("q_index").set_index("q_index")[["q_frac_x", "q_frac_y", "q_frac_z", "multiplicity"]]
    vel = read_result_velocities(Path(f"{pre}.result"))
    vel["v2"] = vel[["v_x", "v_y", "v_z"]].pow(2).sum(axis=1)
    vrms = np.sqrt(vel.groupby(["q_index", "branch"])["v2"].mean())
    dnu = vrms / (A_M * mesh * C_CM)
    a = kc.assign(lo=np.minimum(kc.ibranch, kc.jbranch), hi=np.maximum(kc.ibranch, kc.jbranch))
    p = a.groupby(["ik", "lo", "hi"]).agg(v=("kc_real", "sum"), w_lo=("omega1_cm1", "first")).reset_index()
    p["v"] = p.v / 3 / nk
    w = res.set_index(["q_index", "branch"])["omega_cm1"]
    p["w1"] = w.loc[list(zip(p.ik, p.lo))].to_numpy()
    p["w2"] = w.loc[list(zip(p.ik, p.hi))].to_numpy()
    p["dw"] = (p.w1 - p.w2).abs()
    p["G1pG2"] = g.loc[list(zip(p.ik, p.lo))].to_numpy() + g.loc[list(zip(p.ik, p.hi))].to_numpy()
    p["dnu_mesh"] = 0.5 * (dnu.loc[list(zip(p.ik, p.lo))].to_numpy() + dnu.loc[list(zip(p.ik, p.hi))].to_numpy())
    p["r"] = p.G1pG2 / p.dnu_mesh
    s1 = np.vectorize(SHORT.get)(classify_sector(p.lo, p.w1))
    s2 = np.vectorize(SHORT.get)(classify_sector(p.hi, p.w2))
    p["block"] = ["/".join(sorted(x)) for x in zip(s1, s2)]
    for c in ("q_frac_x", "q_frac_y", "q_frac_z", "multiplicity"):
        p[c] = qf.loc[p.ik, c].to_numpy()
    p["T_K"], p["mesh"] = T, mesh
    return p


UNITS = [(300, m) for m in (12, 15, 16, 17, 20, 24)] + [(T, m) for T in (400, 500, 650, 800) for m in (12, 16, 20)] + [(T, 16) for T in (900, 950, 1000)] + [(1000, 20)]
rows, tops, pairs = [], [], {}
for T, mesh in UNITS:
        p = analyse(T, mesh)
        pairs[(T, mesh)] = p
        kC = p.v.sum()
        wpos = p.v.clip(lower=0)
        o = np.argsort(p.r.to_numpy())
        cw = np.cumsum(wpos.to_numpy()[o]) / wpos.sum()
        blocks = p.groupby("block").v.sum()
        top = p.sort_values("v", ascending=False).head(20).copy()
        top["share_of_kC"] = top.v / kC
        tops.append(top)
        rows.append({"T_K": T, "mesh": mesh, "kappa_C": kC, "hard_hard_share": blocks.get("hard/hard", 0) / kC,
                     **{f"block_{b.replace('/', '_')}": blocks.get(b, 0.0) for b in ("A/A", "A/hard", "A/soft", "hard/hard", "hard/soft", "soft/soft")},
                     "top1_share": top.share_of_kC.iloc[0], "top20_share": top.share_of_kC.sum(),
                     "kw_median_r": float(p.r.to_numpy()[o][np.searchsorted(cw, 0.5)]),
                     "kw_frac_r_gt_1": float(wpos[p.r > 1].sum() / wpos.sum()),
                     "kw_median_G1pG2": float(p.G1pG2.to_numpy()[np.argsort(p.G1pG2.to_numpy())][np.searchsorted(np.cumsum(wpos.to_numpy()[np.argsort(p.G1pG2.to_numpy())]) / wpos.sum(), 0.5)]),
                     "kw_median_dnu_mesh": float(p.dnu_mesh.to_numpy()[np.argsort(p.dnu_mesh.to_numpy())][np.searchsorted(np.cumsum(wpos.to_numpy()[np.argsort(p.dnu_mesh.to_numpy())]) / wpos.sum(), 0.5)])})


def key(p):  # pair identity across meshes: fractional q (irreducible representative) + branch pair
    return set(zip(p.q_frac_x.round(4), p.q_frac_y.round(4), p.q_frac_z.round(4), p.lo, p.hi))


summ = pd.DataFrame(rows)
for T in (300, 500, 800):
    for m1, m2 in ((12, 16), (16, 20), (12, 20)):
        t1 = key(pairs[(T, m1)].sort_values("v", ascending=False).head(20))
        t2 = key(pairs[(T, m2)].sort_values("v", ascending=False).head(20))
        summ.loc[(summ.T_K == T) & (summ.mesh == m2), f"top20_overlap_with_{m1}"] = len(t1 & t2)
summ.to_csv(OUT / "d2_summary.csv", index=False, float_format="%.5g")
pd.concat(tops).to_csv(OUT / "d2_top_pairs.csv", index=False, float_format="%.5g")
summ.to_csv(ROOT / "data/processed/sto_production/convergence_report_kC.csv", index=False, float_format="%.5g")
pd.set_option("display.width", 250)
print(summ.round(4).to_string(index=False))
for (T, mesh), p in pairs.items():
    t = p.sort_values("v", ascending=False).head(5)
    print(f"\nT={T} mesh={mesh} top-5:")
    print(t[["q_frac_x", "q_frac_y", "q_frac_z", "lo", "hi", "w1", "w2", "dw", "G1pG2", "dnu_mesh", "r", "block", "v"]].round(4).to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"d2_summary.csv", f"d2_top_pairs.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['zenodo:pm_q12_rta.tar.gz (*.kc_elem)', 'data/processed/sto_production/mode_table_T.csv.gz'])
