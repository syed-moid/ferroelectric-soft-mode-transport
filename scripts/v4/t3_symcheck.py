#!/usr/bin/env python3
"""T3 weight-assignment check: the same patched RTA run with the triplet symmetry reduction on (TRISYM = 1,
production) and off (TRISYM = 0: every pair explicit except the k1 <-> k2 permutation), same labels.
For every carrier, temperature and bin, the two partitions must agree (relative to the carrier's total).

Usage: t3_symcheck.py T3_OUT_TRISYM1 T3_OUT_TRISYM0
"""
import sys

from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fekappa.partition import read_t3  # noqa: E402


def main():
    tot1, b1 = read_t3(sys.argv[1])
    tot0, b0 = read_t3(sys.argv[2])
    assert set(tot1) == set(tot0), "different carriers"
    worst, worst_tot = 0.0, 0.0
    for key, T in tot1.items():
        if T == 0:
            continue
        worst_tot = max(worst_tot, abs(T - tot0[key]) / abs(T))
        d1 = {(x[0], x[1]): x[2] + 1j * x[3] for x in b1.get(key, [])}
        d0 = {(x[0], x[1]): x[2] + 1j * x[3] for x in b0.get(key, [])}
        for k in set(d1) | set(d0):
            v1, v0 = d1.get(k, 0), d0.get(k, 0)
            worst = max(worst, abs(v1.real - v0.real) / abs(T), abs(v1.imag - v0.imag) / abs(T))
    print(f"carriers x T: {len(tot1)}; max |total(1) - total(0)| / total = {worst_tot:.3e}; "
          f"max |bin(1) - bin(0)| / total (decay and absorption separately) = {worst:.3e}")


if __name__ == "__main__":
    main()
