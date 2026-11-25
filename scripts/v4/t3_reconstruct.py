#!/usr/bin/env python3
"""T3 reconstruction test for one patched run.

Checks, for every carrier (irreducible q, branch) and temperature in T3_OUT:
  1. sum over bins of (decay + absorption) equals the patch's own copy of the base total ("T" line), relative
     difference <= 1e-10 (per carrier, not on average);
  2. that total, converted to cm^-1 (Hz_to_kayser / time_ry, ALAMODE constants), equals the #GAMMA_EACH value of
     the patched run's .result to its printed precision (6 significant digits);
  3. optionally: the patched run's .result equals the base build's .result byte for byte, i.e. the base linewidths
     (and everything else in the file) are unchanged.
Carriers with total 0 (translational modes at Gamma) are reported separately.

Usage: t3_reconstruct.py T3_OUT PATCHED_RESULT [BASE_RESULT]  -> prints a summary, writes T3_OUT.check.csv
"""
import filecmp
import sys

import pandas as pd

from fekappa.partition import HZ_TO_KAYSER, RY_TO_CM, TIME_RY, read_gamma, read_t3  # noqa: E402,F401


def main():
    t3, res = sys.argv[1], sys.argv[2]
    tot, bins = read_t3(t3)
    gam = read_gamma(res)
    rows = []
    for key, T in tot.items():
        b = bins.get(key, [])
        dec = sum(x[2] for x in b); ab = sum(x[3] for x in b)
        s = dec + ab
        rel = abs(s - T) / abs(T) if T != 0 else (abs(s) if s != 0 else 0.0)
        g = gam[key]
        rel_res = abs(T * RY_TO_CM - g) / abs(g) if g != 0 else abs(T * RY_TO_CM)
        rows.append({"iq0": key[0], "branch0": key[1], "T_K": key[2], "total_Ry": T, "sum_bins_Ry": s,
                     "decay_Ry": dec, "absorption_Ry": ab, "rel_bins_vs_total": rel, "result_cm1": g,
                     "rel_total_vs_result": rel_res, "nbins": len(b)})
    df = pd.DataFrame(rows)
    df.to_csv(t3 + ".check.csv", index=False, float_format="%.17g")
    nz = df[df.total_Ry != 0]
    print(f"carriers x T: {len(df)} (zero total: {len(df) - len(nz)}); .result modes: {len(gam)}; "
          f"all .result modes covered: {set(gam) == set(tot)}")
    print(f"1. max per-carrier |sum bins - total| / total = {nz.rel_bins_vs_total.max():.3e} "
          f"-> {'PASS' if nz.rel_bins_vs_total.max() <= 1e-10 else 'FAIL'} (<= 1e-10); "
          f"zero-total carriers with nonzero bins: {int((df[df.total_Ry == 0].sum_bins_Ry != 0).sum())}")
    print(f"2. max |total (cm^-1) - .result| / .result = {nz.rel_total_vs_result.max():.3e} "
          f"-> {'PASS' if nz.rel_total_vs_result.max() <= 5e-6 else 'FAIL'} (printed precision, 6 digits)")
    if len(sys.argv) > 3:
        same = filecmp.cmp(res, sys.argv[3], shallow=False)
        print(f"3. patched .result identical to base .result: {same} -> {'PASS' if same else 'FAIL'}")


if __name__ == "__main__":
    main()
