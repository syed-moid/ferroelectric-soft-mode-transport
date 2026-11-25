#!/usr/bin/env python3
"""T5 control labels: label = 4 * family + band, family 0 neither, 1 pol, 2 rot
(frozen primary c_min 0.5), bands Omega < 105, < 145, < 160, >= 160 cm^-1. Degenerate groups share one label (checked).

Usage: t5_control_labels.py LABELS_CSV OUT_TXT [EDGES]   (EDGES comma-separated; default 105,145,160 = T5;
       T5b: 80,95,105, docs/v4/ANALYSIS_STEPS.md, step T5b)
"""
import sys

import numpy as np
import pandas as pd

EDGES = [105.0, 145.0, 160.0]
FAM = {"neither": 0, "pol": 1, "rot": 2}


def main():
    global EDGES
    if len(sys.argv) > 3:
        EDGES = [float(x) for x in sys.argv[3].split(",")]
    L = pd.read_csv(sys.argv[1]).sort_values(["k", "branch"])
    fam = L["primary_c0.5"].map(FAM)
    assert fam.notna().all(), "family outside neither/pol/rot"
    band = np.searchsorted(EDGES, L.Omega_cm1.to_numpy(), side="right")
    L["t5"] = 4 * fam.to_numpy() + band
    assert (L.groupby(["k", "group"]).t5.nunique() == 1).all(), "degenerate group split across labels"
    with open(sys.argv[2], "w") as f:
        f.write("# qx qy qz branch label ; label = 4 family(0 neither,1 pol,2 rot; primary c0.5) + band(" + "/".join(f"{e:g}" for e in EDGES) + " cm^-1)\n")
        for r in L.itertuples():
            f.write(f"{r.qx:.12f} {r.qy:.12f} {r.qz:.12f} {r.branch} {r.t5}\n")
    print(sys.argv[2], dict(sorted(L.t5.value_counts().items())))


if __name__ == "__main__":
    main()
