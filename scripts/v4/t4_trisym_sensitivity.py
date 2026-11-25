#!/usr/bin/env python3
"""T4 sensitivity to TRISYM on 8^3. Compares two t4_analysis summaries of the same model and labels
(TRISYM = 1 and TRISYM = 0): kappa_P-weighted mean share of the carriers' linewidth with 0 / 1 / 2 soft partners,
and by partner family pair, for every definition. Escalation rule (fixed in the analysis plan before the runs): any of these
summary shares differing by more than 2 percentage points -> a 16^3 TRISYM = 0 run for model A after T4.
Carriers: each run's own 50 % selection (the overlap is reported). The rule is applied to the nsoft and family-pair
shares of all four definitions (primary c0.5, secondary c0.5, primary c0.4, c0.6); the absorption share is reported
only. Fixed 21-Nov 14:55Z, before any T4 output was read.

Usage: t4_trisym_sensitivity.py SUMMARY_T1 SUMMARY_T0 MODES_T1 MODES_T0 OUT_JSON
Writes the comparison to OUT_JSON; prints only the largest difference and the verdict.
"""
import json
import sys

import pandas as pd


def main():
    s1, s0 = (json.load(open(p)) for p in sys.argv[1:3])
    m1, m0 = (pd.read_csv(p) for p in sys.argv[3:5])
    c1 = set(map(tuple, m1[m1.carrier & (m1.definition == "primary_c0.5")][["q_index", "branch"]].to_numpy()))
    c0 = set(map(tuple, m0[m0.carrier & (m0.definition == "primary_c0.5")][["q_index", "branch"]].to_numpy()))
    out = {"model": s1["model"], "kappa_P": [s1["kappa_P_modesum"], s0["kappa_P_modesum"]],
           "carriers": [len(c1), len(c0)], "carriers_common": len(c1 & c0), "definitions": {}}
    worst = (0.0, None)
    for name, d1 in s1["definitions"].items():
        d0 = s0["definitions"][name]
        rows = {}
        for k in ("0", "1", "2"):
            rows[f"nsoft{k}"] = [d1["kappa_weighted_nsoft_frac"][k], d0["kappa_weighted_nsoft_frac"][k]]
        for k in d1["kappa_weighted_pair_frac"]:
            rows[f"pair {k}"] = [d1["kappa_weighted_pair_frac"][k], d0["kappa_weighted_pair_frac"][k]]
        rows["absorption"] = [d1["kappa_weighted_absorption_share"], d0["kappa_weighted_absorption_share"]]
        for k, (a, b) in rows.items():
            pp = 100 * (a - b)
            rows[k] = {"trisym1": a, "trisym0": b, "diff_pp": pp}
            if k != "absorption" and abs(pp) > abs(worst[0]):
                worst = (pp, f"{name}: {k}")
        out["definitions"][name] = rows
    out["max_abs_diff_pp"] = abs(worst[0])
    out["max_at"] = worst[1]
    out["escalate"] = abs(worst[0]) > 2.0
    json.dump(out, open(sys.argv[5], "w"), indent=1)
    print(f"{out['model']}: carriers {len(c1)} / {len(c0)} (common {len(c1 & c0)}); "
          f"max |diff| = {abs(worst[0]):.3f} pp -> {'ESCALATE (> 2 pp)' if out['escalate'] else 'no escalation'}")


if __name__ == "__main__":
    main()
