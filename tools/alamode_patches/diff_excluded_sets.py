#!/usr/bin/env python3
"""Compare the excluded diagram-(c) quartets dumped by the full-precision build (alamode_4ph_v3dump) with the
independent count from the MODE = phonons eigenvalues of the same FC2 and mesh (.mesh.evec, 7 significant digits).

Usage: diff_excluded_sets.py <run dir with AUDITDUMP_k*_s*_r*.txt> <mesh.evec> [knum:snum ...]

For every external mode in the dump: the excluded set from the dump ("X" lines, double precision, as the code
decides), the same predicate re-evaluated here on the dumped double-precision frequencies, and the predicate
evaluated on the .mesh.evec frequencies (the earlier independent count). Quartets are matched by the fractional
wave vectors of their legs and their branch indices. Every difference is printed with the frequencies from both
sources.
"""
import glob
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from fekappa.evec_io import read_evec

TOL = 1e-9
run, evec = Path(sys.argv[1]), sys.argv[2]

kx = {}
quart = defaultdict(list)
for f in sorted(glob.glob(str(run / "AUDITDUMP_k*_s*_r*.txt"))):
    tag = Path(f).name.split("_r")[0]  # AUDITDUMP_k<knum>_s<snum>
    for line in open(f):
        if line.startswith("K "):
            _, ik, a, b, c = line.split(); kx[int(ik)] = (float(a), float(b), float(c))
        elif line[:2] in ("X ", "N "):
            t = line.split()
            quart[tag].append((t[0], *map(int, t[1:7]), *map(float, t[7:14]), *map(int, t[14:20])))

xk, ev2, _ = read_evec(evec)
om_py = np.sign(ev2) * np.sqrt(np.abs(ev2))
n = round(len(xk) ** (1 / 3))
key_py = {tuple(np.round((x * n) % n).astype(int) % n): i for i, x in enumerate(xk)}
cpp2py = {ik: key_py[tuple(np.round((np.array(x) * n) % n).astype(int) % n)] for ik, x in kx.items()}


def predicate(knum, w, w1, w2, w3, ik1, ik2, ik3):
    e1 = ik1 == knum and abs(w1 - w) < TOL and abs(w2 - w3) < TOL
    e2 = ik2 == knum and abs(w2 - w) < TOL and abs(w1 - w3) < TOL
    e3 = ik3 == knum and abs(w3 - w) < TOL and abs(w1 - w2) < TOL
    return e1 or e2 or e3


for tag, rows in quart.items():
    knum, snum = int(tag.split("_k")[1].split("_s")[0]), int(tag.split("_s")[1])
    kq = cpp2py[knum]
    wq_py = om_py[kq, snum]
    X_cpp = {r[1:7] for r in rows if r[0] == "X"}
    # Python predicate on the .mesh.evec frequencies, restricted to the dumped candidates (every quartet the
    # code excluded or held within 1e-6 Ry of resonance); the earlier count covered the full loop.
    X_py, disagree = set(), []
    for r in rows:
        ik1, ik2, ik3, is1, is2, is3 = r[1:7]
        w_c, w1_c, w2_c, w3_c = r[7:11]
        p1, p2, p3 = cpp2py[ik1], cpp2py[ik2], cpp2py[ik3]
        w1_p, w2_p, w3_p = om_py[p1, is1], om_py[p2, is2], om_py[p3, is3]
        in_py = predicate(kq, wq_py, w1_p, w2_p, w3_p, p1, p2, p3)
        in_cpp_recomputed = predicate(knum, w_c, w1_c, w2_c, w3_c, ik1, ik2, ik3)
        assert in_cpp_recomputed == (r[0] == "X"), r
        if in_py:
            X_py.add(r[1:7])
        if in_py != (r[0] == "X"):
            disagree.append((r, (wq_py, w1_p, w2_p, w3_p)))
    print(f"mode knum={knum} snum={snum + 1}: excluded C++ (double) = {len(X_cpp)}, "
          f"Python on .mesh.evec (same candidates) = {len(X_py)}, disagreements = {len(disagree)}")
    by_reason = defaultdict(int)
    for r, (wq, w1p, w2p, w3p) in disagree:
        ik1, ik2, ik3, is1, is2, is3 = r[1:7]
        w_c, w1_c, w2_c, w3_c = r[7:11]
        legs = [(ik1, w1_c, w2_c, w3_c, w1p, w2p, w3p), (ik2, w2_c, w1_c, w3_c, w2p, w1p, w3p),
                (ik3, w3_c, w1_c, w2_c, w3p, w1p, w2p)]
        for ik, wi_c, wa_c, wb_c, wi_p, wa_p, wb_p in legs:
            if ik != knum:
                continue
            reason = (f"C++ |w_i-w| = {abs(wi_c - w_c):.3e}, |w_a-w_b| = {abs(wa_c - wb_c):.3e}; "
                      f"evec |w_i-w| = {abs(wi_p - wq):.3e}, |w_a-w_b| = {abs(wa_p - wb_p):.3e}")
            by_reason[("C++ only" if r[0] == "X" else "Python only",
                       "pair split in evec" if abs(wa_p - wb_p) >= TOL else
                       "leg split in evec" if abs(wi_p - wq) >= TOL else "other")] += 1
            if len(by_reason) and sum(by_reason.values()) <= 6:
                print("   ", r[0], r[1:7], reason)
    for k, c in sorted(by_reason.items()):
        print(f"    {k[0]:12s} {k[1]:22s} {c}")
