#!/usr/bin/env python3
"""Neutral label tables for the T3 reconstruction test (nothing related to the mode families).

Writes two T3_LABELS files covering every wave vector of the mesh of MESH_EVEC and every branch:
  labels_branch.txt   : label = branch index - 1 (0..14)
  labels_starhash.txt : label = a hash of the wave vector's star (|q| components folded and sorted) mod 7,
                        combined with branch parity: 2 * (hash mod 7) + (branch mod 2)  (0..13)
Both are invariant under the point-group operations used by anphon's triplet reduction (branch index: eigenvalues
are sorted per q and equal on a star; star hash: by construction), so the patch's assignment of a representative
pair's weight is exact for them (docs/v4/ANALYSIS_STEPS.md, step T3). A hash of the raw wave vector would not be.

Usage: t3_neutral_labels.py MESH_EVEC OUTDIR
"""
import hashlib
import sys
from pathlib import Path

import numpy as np

from fekappa.evec_io import read_evec


def star_key(q):
    f = np.abs(((np.asarray(q) + 0.5) % 1.0) - 0.5)
    return "%.6f_%.6f_%.6f" % tuple(np.sort(f))


def main():
    xk, _, _ = read_evec(sys.argv[1])
    out = Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
    with open(out / "labels_branch.txt", "w") as fb, open(out / "labels_starhash.txt", "w") as fh:
        fb.write("# qx qy qz branch label ; label = branch - 1\n")
        fh.write("# qx qy qz branch label ; label = 2*(sha256(star key) mod 7) + (branch mod 2)\n")
        for q in xk:
            h = int(hashlib.sha256(star_key(q).encode()).hexdigest(), 16) % 7
            for b in range(1, 16):
                fb.write("%.12f %.12f %.12f %d %d\n" % (*q, b, b - 1))
                fh.write("%.12f %.12f %.12f %d %d\n" % (*q, b, 2 * h + b % 2))
    print(f"{len(xk)} wave vectors x 15 branches -> {out}/labels_branch.txt, labels_starhash.txt")


if __name__ == "__main__":
    main()
