#!/usr/bin/env python3
"""A3 step 6: transport quantities for one RTA run (+ optional eigenvector run), as in scripts/v3/prod_analysis.py.

Per run (300 K):
  1. R-mode frequency, half-width (HWHM) and their ratio eta, separately (lowest branch at q = (1/2,1/2,1/2));
  2. kappa_P (mode sum, closure-checked against the solver's .kl value);
  3. share of kappa_P carried by modes with eta > 0.1;
  4. kappa_C (pair sum of .kc_elem; and the solver's .kl_coherent value);
  5. hard/hard share of kappa_C under both classification rules (freq: branches 1-3 acoustic, else soft if
     Omega < 175 cm^-1, else hard; char: projection on the Gamma subspaces of the same run, degenerate subspaces
     labelled as a whole; needs the eigenvector run), with pair counts (ordered branch pairs per irreducible q;
     star-weighted count) and mean contribution per pair (scripts/v3/e1_pair_counts.py convention).
Without an eigenvector run only the freq rule is evaluated.

Option --phase=atom: the char rule projects the eigenvectors in the atom-position phase
convention, e~_j = e_j exp(-2 pi i q.r_j) with q the eigenvector's own wave vector, which makes the labels identical
across each star. Default (--phase=cell) is the published convention (v9.4), unchanged.

Usage: a3_transport_analysis.py [--phase=cell|atom] LABEL RTA_DIR [EVEC_DIR]   (prints one JSON line)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from fekappa.evec_io import irred_evecs, irred_q, kmap, read_evec
from fekappa.kappa_decomposition import classify_sector, mode_kappa_table, read_kc_elem

SHORT = {"acoustic": "A", "soft_optic": "soft", "hard_optic": "hard"}
DEG_TOL = 1e-3
POS = np.array([[0, 0, 0], [.5, .5, .5], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])   # Sr, Ti, O_x, O_y, O_z
ATOM_OF = np.repeat(np.arange(5), 3)


def kl_value(d, suffix, ncomp):
    f = next(Path(d).glob(f"*.{suffix}"))
    for l in f.read_text().splitlines():
        s = l.split()
        if s and not l.startswith("#") and abs(float(s[0]) - 300) < 1e-6:
            v = [float(x) for x in s[1:1 + ncomp]]
            return (v[0] + v[4] + v[8]) / 3 if ncomp == 9 else sum(v) / 3
    raise ValueError(f)


def cell_volume_ang3(rdir):
    """Cell volume from the &cell block of the run's rta.in (cubic: lattice constant in Bohr x unit vectors), so
    that each model's own lattice constant is used (model A's 7.363 Bohr is the library default)."""
    from fekappa.kappa_decomposition import BOHR_TO_ANG
    lines = (Path(rdir) / "rta.in").read_text().splitlines()
    i = next(k for k, l in enumerate(lines) if l.strip().lower() == "&cell")
    a = float(lines[i + 1].split()[0])
    m = np.array([[float(x) for x in lines[i + 2 + r].split()[:3]] for r in range(3)])
    return abs(np.linalg.det(a * m)) * BOHR_TO_ANG ** 3


def main():
    opts = [a for a in sys.argv[1:] if a.startswith("--phase=")]
    argv = [a for a in sys.argv[1:] if not a.startswith("--phase=")]
    phase = opts[-1].split("=", 1)[1] if opts else "cell"
    assert phase in ("cell", "atom"), phase
    label, rdir = argv[0], Path(argv[1])
    edir = Path(argv[2]) if len(argv) > 2 else None
    res = next(rdir.glob("*.result"))
    rin = (rdir / "rta.in").read_text()
    mesh = int([l for l in rin.splitlines() if "KPOINT" not in l and l.split()[:1] and l.split()[0].isdigit()
                and len(l.split()) == 3][0].split()[0])
    df = mode_kappa_table(res, nk=mesh ** 3, V_cell_ang3=cell_volume_ang3(rdir))
    df = df[df.T_K == 300].copy()
    kP_solver = kl_value(rdir, "kl", 9)
    kP = df.kappa_mode_W_mK.sum()
    closure = abs(kP - kP_solver) / kP_solver
    assert closure < 5e-3, closure
    df["eta"] = np.where(df.omega_cm1 > 1e-3, df.gamma_hwhm_cm1 / df.omega_cm1, np.nan)
    df["sector_freq"] = np.vectorize(SHORT.get)(classify_sector(df.branch, df.omega_cm1))
    rules = ["freq"]
    if edir is not None:
        E, W = irred_evecs(next(edir.glob("*.mesh.evec")), res)
        if phase == "atom":
            xk, _, _ = read_evec(next(edir.glob("*.mesh.evec")))
            qv = xk[kmap(xk, irred_q(res))]
            E = E * np.exp(-2j * np.pi * (qv @ POS.T))[:, None, ATOM_OF]
        w0 = np.sort(df[df.q_index == 1].omega_cm1.to_numpy())
        assert abs(w0[3] - w0[5]) < 1e-3 * w0[3], "Gamma branches 4-6 not a degenerate triplet"
        G = E[0]
        sub = {"A": G[0:3], "soft": G[3:6], "hard": G[6:15]}
        ch = {S: (np.abs(np.einsum("ic,qjc->qij", v.conj(), E)) ** 2).sum(axis=1) for S, v in sub.items()}
        qi, bi = df.q_index.to_numpy() - 1, df.branch.to_numpy() - 1
        for S in sub:
            df[f"c_{S}"] = ch[S][qi, bi]
        lab = np.array(["A", "soft", "hard"])
        cols = ["c_A", "c_soft", "c_hard"]
        df = df.sort_values(["q_index", "omega_cm1"])
        grp = df.groupby("q_index").omega_cm1.transform(
            lambda w: np.concatenate([[0], np.cumsum(np.diff(w.to_numpy()) > DEG_TOL)]))
        csum = df.groupby([df.q_index, grp])[cols].transform("sum").to_numpy()
        df["sector_char"] = lab[np.argmax(csum, axis=1)]
        df = df.sort_values(["q_index", "branch"])
        rules.append("char")

    # R mode
    isR = np.all(np.isclose(np.abs(df[["q_frac_x", "q_frac_y", "q_frac_z"]].to_numpy()), 0.5), axis=1)
    r = df[isR].sort_values("omega_cm1").iloc[0]
    s = df[df.omega_cm1 > 1e-3]
    out = {"label": label, "char_phase": phase if edir is not None else None, "mesh": mesh, "kappa_P_solver": kP_solver, "kappa_P_modesum": kP, "closure_rel": closure,
           "R_omega_cm1": float(r.omega_cm1), "R_hwhm_cm1": float(r.gamma_hwhm_cm1), "R_eta": float(r.eta),
           "R_sector_freq": r.sector_freq, "R_sector_char": r.get("sector_char", None),
           "share_kP_eta_gt_0p1": float(s[s.eta > 0.1].kappa_mode_W_mK.sum() / kP),
           "eta_max": float(s.eta.max()),
           "kappa_C_solver": kl_value(rdir, "kl_coherent", 3)}

    kc = read_kc_elem(next(rdir.glob("*.kc_elem")))
    kc = kc[(kc.T_K == 300) & (kc.i == kc.j) & (kc.ibranch != kc.jbranch)]
    pair = kc.groupby(["ik", "ibranch", "jbranch"]).kc_real.sum().reset_index()  # Cartesian diagonal summed
    pair["kc"] = pair.kc_real / 3 / mesh ** 3
    out["kappa_C_pairsum"] = float(pair.kc.sum())
    mult = df.groupby("q_index").multiplicity.first()
    for rule in rules:
        L = np.empty((df.q_index.max(), 15), dtype=object)
        L[df.q_index.to_numpy() - 1, df.branch.to_numpy() - 1] = df[f"sector_{rule}"].to_numpy()
        s1 = L[pair.ik.to_numpy() - 1, pair.ibranch.to_numpy() - 1]
        s2 = L[pair.ik.to_numpy() - 1, pair.jbranch.to_numpy() - 1]
        blk = np.array(["/".join(sorted(p)) for p in zip(s1, s2)])
        hh = blk == "hard/hard"
        kc_hh = float(pair.kc[hh].sum())
        n_hh = int(hh.sum())
        n_hh_w = int(mult.loc[pair.ik[hh]].sum())
        out.update({f"kC_hardhard_{rule}": kc_hh, f"share_kC_hardhard_{rule}": kc_hh / out["kappa_C_pairsum"],
                    f"n_pairs_hardhard_{rule}": n_hh, f"n_pairs_w_hardhard_{rule}": n_hh_w,
                    f"kC_per_pair_hardhard_{rule}": kc_hh / n_hh if n_hh else np.nan,
                    f"kC_per_pair_w_hardhard_{rule}": kc_hh / n_hh_w if n_hh_w else np.nan})
    print(json.dumps(out))


if __name__ == "__main__":
    main()
