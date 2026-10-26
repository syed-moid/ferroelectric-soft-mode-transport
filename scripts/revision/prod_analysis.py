"""Production analysis (own SCPH q1 = 12^3, PM; tetrahedron; NONANALYTIC = 3).

Per (T, mesh) with a validated RTA unit and matching eigenvectors (12^3, 16^3 all T; 20^3 at 300/800 K):
  * mode table: omega, Gamma (HWHM), kappa_P per mode (star-summed, closure-checked vs solver),
    eta = Gamma/Omega, multiplicity, two sector labellings:
      - 'freq': branch 1-3 acoustic, else omega < 175 cm^-1 soft, else hard   (as submitted)
      - 'char': eigenvector character relative to the Gamma-point subspaces of the same (T, mesh):
                c_S(q,j) = sum_{i in S(Gamma)} |<e_i(Gamma)|e_j(q)>|^2, S = A (3 translations),
                soft (TO1 triplet = lowest optic triplet at Gamma), hard (rest); label = argmax.
                Exactly degenerate subspaces at q (|dOmega| < 1e-3 cm^-1) are labelled as a whole by argmax of the
                summed c_S, so the label does not depend on the diagonalizer's basis;
                the per-eigenvector label is kept as 'sector_char_eig'
  * validity summary per (T, mesh): max eta; medians (unweighted / mesh-weighted / kappa-weighted);
    fractions eta > 0.05, 0.1, 0.5 (mesh- and kappa-weighted); sector of the max-eta mode
  * particle sector shares and coherence blocks (from .kc_elem pairs) under both labellings
Outputs (data/processed/sto_production/): mode_table_T.csv.gz, eta_validity_T.csv,
  kappa_sectors_T.csv, coherence_blocks_T.csv
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
from fekappa.evec_io import irred_evecs
from fekappa.kappa_decomposition import classify_sector, mode_kappa_table, read_kc_elem  # noqa: E402

ROOT = HERE.parents[1]
ARCH = ROOT / "data/raw/alamode_sto/revision/pm_q12"
OUT = ROOT / "data/processed/sto_production"
SHORT = {"acoustic": "A", "soft_optic": "soft", "hard_optic": "hard"}
DEG_TOL = 1e-3   # cm^-1, exact degeneracy in the mode table
BLOCKS = ["A/A", "A/hard", "A/soft", "hard/hard", "hard/soft", "soft/soft"]
man = json.loads((ROOT / "data/provenance/manifest_final.json").read_text())["jobs"]


def wmedian(x, w):
    o = np.argsort(x)
    c = np.cumsum(w[o]) / w[o].sum()
    return float(x[o][np.searchsorted(c, 0.5)])


def one(T, mesh):
    rj, ej = f"rta_pm_q12_T{T}_m{mesh}", f"evec_pm_q12_T{T}_m{mesh}"
    res = ARCH / rj / f"{rj}.result"
    df = mode_kappa_table(res, nk=mesh ** 3)
    solver = man[rj]["metrics"]["kappa_P"]
    closure = abs(df.kappa_mode_W_mK.sum() - solver) / solver
    assert closure < 5e-3, (rj, closure)
    E, W = irred_evecs(ARCH / ej / f"{ej}.mesh.evec", res)
    # Gamma subspaces: q_index 1 is Gamma; branches sorted by frequency
    w0 = np.sort(df[df.q_index == 1].omega_cm1.to_numpy())
    assert abs(w0[3] - w0[5]) < 1e-3 * w0[3], f"Gamma branches 4-6 not a degenerate triplet at {T} K {mesh}^3"
    G = E[0]
    sub = {"A": G[0:3], "soft": G[3:6], "hard": G[6:15]}
    ch = {S: (np.abs(np.einsum("ic,qjc->qij", v.conj(), E)) ** 2).sum(axis=1) for S, v in sub.items()}  # [q, j]
    qi, bi = df.q_index.to_numpy() - 1, df.branch.to_numpy() - 1
    for S in sub:
        df[f"c_{S}"] = ch[S][qi, bi]
    lab = np.array(["A", "soft", "hard"])
    cols = ["c_A", "c_soft", "c_hard"]
    df["sector_char_eig"] = lab[np.argmax(df[cols].to_numpy(), axis=1)]
    df = df.sort_values(["q_index", "omega_cm1"])
    sub = df.groupby("q_index").omega_cm1.transform(
        lambda w: np.concatenate([[0], np.cumsum(np.diff(w.to_numpy()) > DEG_TOL)]))
    csum = df.groupby([df.q_index, sub])[cols].transform("sum").to_numpy()
    df["sector_char"] = lab[np.argmax(csum, axis=1)]
    df = df.sort_values(["q_index", "branch"])
    df["sector_freq"] = np.vectorize(SHORT.get)(classify_sector(df.branch, df.omega_cm1))
    df["eta"] = np.where(df.omega_cm1 > 1e-3, df.gamma_hwhm_cm1 / df.omega_cm1, np.nan)
    df["mesh"] = mesh
    return df, closure


def blocks(T, mesh, lab):
    rj = f"rta_pm_q12_T{T}_m{mesh}"
    kc = read_kc_elem(ARCH / rj / f"{rj}.kc_elem")
    kc = kc[(kc.i == kc.j) & (kc.ibranch != kc.jbranch)]
    s1 = lab[kc.ik.to_numpy() - 1, kc.ibranch.to_numpy() - 1]
    s2 = lab[kc.ik.to_numpy() - 1, kc.jbranch.to_numpy() - 1]
    b = pd.Series(kc.kc_real.to_numpy() / 3 / mesh ** 3).groupby(["/".join(sorted(p)) for p in zip(s1, s2)]).sum()
    return {k: float(b.get(k, 0.0)) for k in BLOCKS}


units = [(T, m) for T in range(250, 1001, 50) for m in (12, 16)] + [(300, 15), (300, 17)] + [(T, 20) for T in (300, 400, 500, 650, 800, 1000)]
modes, val, sect, coh = [], [], [], []
for T, mesh in units:
    df, closure = one(T, mesh)
    modes.append(df)
    s = df[df.omega_cm1 > 1e-3]
    eta, mult, kap = s.eta.to_numpy(), s.multiplicity.to_numpy().astype(float), s.kappa_mode_W_mK.to_numpy()
    kw = kap / kap.sum()
    imax = int(np.argmax(eta))
    val.append({"T_K": T, "mesh": mesh, "n_irred_modes": len(df), "n_in_map": len(s), "closure_rel": closure,
                "eta_max": float(eta.max()), "eta_max_omega_cm1": float(s.omega_cm1.iloc[imax]),
                "eta_max_q_index": int(s.q_index.iloc[imax]), "eta_max_branch": int(s.branch.iloc[imax]),
                "eta_max_sector_char": s.sector_char.iloc[imax], "eta_max_sector_freq": s.sector_freq.iloc[imax],
                "eta_median_unweighted": float(np.median(eta)), "eta_median_mesh": wmedian(eta, mult),
                "eta_median_kappa": wmedian(eta, kw),
                **{f"frac_eta_gt_{str(x).replace('.', 'p')}_{w}": float(wt[eta > x].sum() / wt.sum())
                   for x in (0.05, 0.1, 0.5) for w, wt in (("mesh", mult), ("kappa", kw))}})
    kP = df.kappa_mode_W_mK.sum()
    for lab in ("char", "freq"):
        g = df.groupby(f"sector_{lab}").kappa_mode_W_mK.sum()
        sect.append({"T_K": T, "mesh": mesh, "labels": lab, "kappa_P": kP,
                     **{f"kappa_P_{S}": float(g.get(S, 0.0)) for S in ("A", "soft", "hard")}})
        L = np.empty((df.q_index.max(), 15), dtype=object)
        L[df.q_index.to_numpy() - 1, df.branch.to_numpy() - 1] = df[f"sector_{lab}"].to_numpy()
        b = blocks(T, mesh, L)
        coh.append({"T_K": T, "mesh": mesh, "labels": lab, "kappa_C_pairsum": sum(b.values()),
                    **{f"kC_{k.replace('/', '_')}": v for k, v in b.items()}})
    print(f"{T} K {mesh}^3: closure {closure:.1e}, eta_max {eta.max():.4f} ({s.sector_char.iloc[imax]}, "
          f"{s.omega_cm1.iloc[imax]:.1f} cm-1)", flush=True)

OUT.mkdir(parents=True, exist_ok=True)
hdr = ("# source: data/raw/alamode_sto/revision/pm_q12 (own SCPH q1=12^3, PM; anphon RTA tetrahedron, NONANALYTIC=3,"
       " KAPPA_COHERENT=2) + MODE=phonons eigenvectors on the same FC2/mesh; script scripts/revision/prod_analysis.py\n")
pd.concat(modes).to_csv(OUT / "mode_table_T.csv.gz", index=False, float_format="%.6g")
for name, rows in (("eta_validity_T.csv", val), ("kappa_sectors_T.csv", sect), ("coherence_blocks_T.csv", coh)):
    with open(OUT / name, "w") as fh:
        fh.write(hdr)
        pd.DataFrame(rows).to_csv(fh, index=False, float_format="%.6g")
print("written")
