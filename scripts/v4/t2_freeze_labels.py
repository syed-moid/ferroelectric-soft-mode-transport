#!/usr/bin/env python3
"""Frozen family labels from a t2_family_census.py table. Same code, no model-specific parameter.

Characters (atom-position phase convention, degenerate groups averaged, as in the census):
  C_pol   plain polar projection on the Gamma TO1 triplet of the same run         (secondary definition)
  C_pol_w = w_pol(q) C_pol, w_pol = cos^2(pi qx) cos^2(pi qy) cos^2(pi qz)         (primary definition)
  C_rot   rotation-field projection with w_a = sin^2(pi q_b) sin^2(pi q_c)          (both definitions)
Label at threshold c: pol if P >= c > C_rot; rot if C_rot >= c > P; mixed if both >= c; else neither.
Translational modes (branches 1-3 at q = 0) are neither.

T3/T4 encoding (one run per model; every definition and threshold recoverable by summing bins):
  level(x) = 0 if x >= 0.6, 1 if 0.5 <= x < 0.6, 2 if 0.4 <= x < 0.5, 3 if x < 0.4
  code = 8 level(C_pol_w) + 2 level(C_rot) + [C_pol >= 0.5]          (0..31; translational -> 30)

Usage: t2_freeze_labels.py CENSUS_CSV MESH_EVEC OUTDIR TAG
Writes OUTDIR/labels_TAG.csv, OUTDIR/t3_code32_TAG.txt, OUTDIR/t3_primary_TAG.txt (0 neither, 1 pol, 2 rot,
3 mixed; primary, c = 0.5) and prints checks.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from fekappa.evec_io import read_evec

CMINS = (0.4, 0.5, 0.6)
CODE = {"neither": 0, "pol": 1, "rot": 2, "mixed": 3}


def label(P, R, c):
    p, r = P >= c, R >= c
    return np.select([p & ~r, r & ~p, p & r], ["pol", "rot", "mixed"], "neither")


def level(x):
    return np.select([x >= 0.6, x >= 0.5, x >= 0.4], [0, 1, 2], 3)


def decode(code):
    """(level of C_pol_w, level of C_rot, plain flag) from the 32-label code."""
    code = np.asarray(code)
    return code // 8, (code % 8) // 2, code % 2


def labels_from_code(code, definition, c):
    """Label recovered from the code alone (used by T4 to sum bins); c in CMINS for primary, 0.5 for secondary."""
    lp, lr, flag = decode(code)
    ok = {0.6: 0, 0.5: 1, 0.4: 2}[c]
    R = lr <= ok
    P = (lp <= ok) if definition == "primary" else (flag == 1)
    assert definition == "primary" or c == 0.5
    return np.select([P & ~R, R & ~P, P & R], ["pol", "rot", "mixed"], "neither")


def main():
    census, evec, out, tag = sys.argv[1], sys.argv[2], Path(sys.argv[3]), sys.argv[4]
    out.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(census)
    xk, _, _ = read_evec(evec)
    assert len(d) == 15 * len(xk) and (d.k.to_numpy() == np.repeat(np.arange(len(xk)), 15)).all()
    d["w_pol"] = np.prod(np.cos(np.pi * xk[d.k.to_numpy()]) ** 2, axis=1)
    d["C_pol_w"] = d.w_pol * d.C_pol
    trans = (d.d_Gamma < 1e-9) & (d.branch <= 3)
    for c in CMINS:
        d[f"primary_c{c}"] = np.where(trans, "neither", label(d.C_pol_w, d.C_rot, c))
        d[f"secondary_c{c}"] = np.where(trans, "neither", label(d.C_pol, d.C_rot, c))
    code = 8 * level(d.C_pol_w) + 2 * level(d.C_rot) + (d.C_pol >= 0.5).astype(int)
    d["code32"] = np.where(trans, 30, code)

    # the code reproduces every label it is meant to carry
    for c in CMINS:
        assert (labels_from_code(d.code32, "primary", c) == d[f"primary_c{c}"]).all(), c
    assert (labels_from_code(d.code32, "secondary", 0.5) == d["secondary_c0.5"]).all()
    # star invariance (members sorted by frequency) and degenerate groups sharing the code
    bad = 0
    for _, g in d.groupby("star"):
        sigs = {tuple(h.sort_values("Omega_cm1").code32) for _, h in g.groupby("k")}
        bad += len(sigs) - 1
    deg = d.groupby(["k", "group"]).code32.nunique().max()
    print(f"{tag}: {len(xk)} q x 15; stars {d.star.nunique()}; star members with a different code sequence: {bad}; "
          f"max distinct codes in a degenerate group: {deg}")
    assert bad == 0 and deg == 1

    cols = ["k", "qx", "qy", "qz", "branch", "group", "gdim", "Omega_cm1", "C_pol", "w_pol", "C_pol_w", "C_rot",
            "code32"] + [f"{x}_c{c}" for x in ("primary", "secondary") for c in CMINS]
    d[cols].to_csv(out / f"labels_{tag}.csv", index=False, float_format="%.10g")
    with open(out / f"t3_code32_{tag}.txt", "w") as f32, open(out / f"t3_primary_{tag}.txt", "w") as fp:
        f32.write("# qx qy qz branch label ; label = 8 level(C_pol_w) + 2 level(C_rot) + [C_pol >= 0.5]\n")
        fp.write("# qx qy qz branch label ; primary, c_min 0.5: 0 neither, 1 pol, 2 rot, 3 mixed\n")
        for (k, b, c32, lp) in zip(d.k, d.branch, d.code32, d["primary_c0.5"]):
            q = xk[k]
            f32.write("%.12f %.12f %.12f %d %d\n" % (*q, b, c32))
            fp.write("%.12f %.12f %.12f %d %d\n" % (*q, b, CODE[lp]))
    for x in ("primary", "secondary"):
        print(x, {c: d[f"{x}_c{c}"].value_counts().to_dict() for c in CMINS})


if __name__ == "__main__":
    main()
