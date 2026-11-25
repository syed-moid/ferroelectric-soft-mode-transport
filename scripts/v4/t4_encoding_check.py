#!/usr/bin/env python3
"""T4 encoding check: on the 8^3 mesh, the bins of the 32-label run summed to the primary classes
(c_min 0.5) must equal, per carrier, the bins of a run made directly with the primary labels
(0 neither, 1 pol, 2 rot, 3 mixed), decay and absorption separately, relative to the carrier's total.

Usage: t4_encoding_check.py T3_OUT_CODE32 T3_OUT_PRIMARY
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from t2_freeze_labels import labels_from_code  # noqa: E402
from fekappa.partition import read_t3  # noqa: E402

FAMS = ["neither", "pol", "rot", "mixed"]


def main():
    tot32, b32 = read_t3(sys.argv[1])
    totp, bp = read_t3(sys.argv[2])
    assert set(tot32) == set(totp)
    fam = {c: FAMS.index(f) for c, f in zip(range(32), labels_from_code(np.arange(32), "primary", 0.5))}
    worst_t, worst_b, n = 0.0, 0.0, 0
    for key, T in tot32.items():
        if T == 0:
            continue
        n += 1
        worst_t = max(worst_t, abs(T - totp[key]) / abs(T))
        a = {}
        for l1, l2, dec, ab in b32.get(key, []):
            k = tuple(sorted((fam[l1], fam[l2])))
            a[k] = a.get(k, np.zeros(2)) + (dec, ab)
        d = {(l1, l2): np.array([dec, ab]) for l1, l2, dec, ab in bp.get(key, [])}
        for k in set(a) | set(d):
            worst_b = max(worst_b, np.abs(a.get(k, np.zeros(2)) - d.get(k, np.zeros(2))).max() / abs(T))
    ok = worst_b <= 1e-10 and worst_t <= 1e-10
    print(f"carriers x T (nonzero): {n}; max |total(32) - total(primary)| / total = {worst_t:.3e}; "
          f"max |bin(32 summed) - bin(primary)| / total = {worst_b:.3e} -> {'PASS' if ok else 'FAIL'} (<= 1e-10)")


if __name__ == "__main__":
    main()
