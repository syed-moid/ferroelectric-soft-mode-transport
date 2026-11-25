#!/usr/bin/env python3
"""Replacement criterion of reconstruction test 3: a T3 label table (qx qy qz branch label, the
full Gamma-centred N^3 mesh in reduced coordinates of the simple cubic cell) is invariant, branch by branch,
under all 48 operations of the cubic point group O_h and under k -> -k. Counts violations (q, branch, op) with
label(S q, branch) != label(q, branch). Labels only; no linewidth is read.

For the simple cubic lattice the 48 point operations act on reduced coordinates as the 48 signed permutation
matrices (the set is closed under inverse and transpose, so the convention of anphon's knum_sym does not matter).

Usage: t3_label_invariance.py TABLE [TABLE ...]
"""
import hashlib
import itertools
import sys

import numpy as np


def ops48():
    out = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            M = np.zeros((3, 3), dtype=int)
            for i, (p, s) in enumerate(zip(perm, signs)):
                M[i, p] = s
            out.append(M)
    return out


def check(path):
    raw = open(path, "rb").read()
    sha = hashlib.sha256(raw).hexdigest()
    d = np.loadtxt(path, comments="#")
    q, br, lab = d[:, :3], d[:, 3].astype(int), d[:, 4].astype(int)
    nq = len(np.unique(np.round(q, 9), axis=0))
    N = round(nq ** (1 / 3))
    assert N ** 3 == nq, f"{path}: {nq} q points is not a full N^3 mesh"
    nb = br.max()
    idx = np.mod(np.rint(q * N).astype(int), N)
    L = np.full((N, N, N, nb + 1), -1, dtype=int)
    L[idx[:, 0], idx[:, 1], idx[:, 2], br] = lab
    assert len(d) == nq * nb and (L[..., 1:] >= 0).all(), f"{path}: mesh x branch not complete / duplicated"
    g = np.stack(np.meshgrid(*[np.arange(N)] * 3, indexing="ij"), -1).reshape(-1, 3)
    base = L[g[:, 0], g[:, 1], g[:, 2], 1:]
    viol, n_ops = 0, 0
    for M in ops48() + [-np.eye(3, dtype=int)]:          # 48 operations of O_h, then k -> -k explicitly
        h = np.mod(g @ M.T, N)
        viol += int((L[h[:, 0], h[:, 1], h[:, 2], 1:] != base).sum())
        n_ops += 1
    return dict(path=path, sha256=sha, N=N, branches=nb, checks=n_ops * N ** 3 * nb, violations=viol)


def main():
    bad = 0
    for p in sys.argv[1:]:
        r = check(p)
        bad += r["violations"]
        print(f"{r['path']}: mesh {r['N']}^3 x {r['branches']} branches; operations 48 + (k -> -k); "
              f"comparisons {r['checks']}; violations {r['violations']}; sha256 {r['sha256']}")
    print("PASS (zero violations)" if bad == 0 else f"FAIL ({bad} violations)")
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
