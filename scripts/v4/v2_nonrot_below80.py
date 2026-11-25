#!/usr/bin/env python3
"""Composition of the non-rotational partner class below 80 cm^-1 (T5b band 0), both models, 300 K, 16^3.

Labels and frequencies only: no linewidth, kappa or scattering quantity is read.

Inputs: the frozen label tables labels_{A,B}_m16.csv (primary c_min 0.5, Omega per mode, full mesh) and the
T5b label tables t5b_labels_{A,B}_m16.txt (label = 4 family + band; band edges 80/95/105 cm^-1). The class
"non-rot" is neither U pol (T5b: label 0 or 4). Counts are on the full 16^3 mesh (star-weighted, 61 440 modes), as N_c
of the T5b analysis output t5b_{A,B}.json.

Distances in reduced units (2 pi / a), folded to the first zone as in scripts/t2_family_census.py:
  d_Gamma = |f|, f = |((q + 1/2) mod 1) - 1/2|;  d_RM = distance to the nearest R-M line (two components = 1/2).
Usage: v2_nonrot_below80.py LABELS_DIR T5B_LABELS_DIR T5B_JSON_DIR
  labels_{A,B}_m16.csv in LABELS_DIR, t5b_labels_{A,B}_m16.txt in T5B_LABELS_DIR, t5b_{A,B}.json in T5B_JSON_DIR.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

EDGES = [80.0, 95.0, 105.0]
FAM = {"neither": 0, "pol": 1, "rot": 2}


def dists(q):
    f = np.abs(((q + 0.5) % 1.0) - 0.5)
    d = 0.5 - f
    dG = np.linalg.norm(f, axis=1)
    dRM = np.min(np.stack([np.hypot(d[:, 1], d[:, 2]), np.hypot(d[:, 0], d[:, 2]), np.hypot(d[:, 0], d[:, 1])]), axis=0)
    return dG, dRM


def stats(x):
    if len(x) == 0:
        return "—"
    return "min %.4f, median %.4f, max %.4f" % (np.min(x), np.median(x), np.max(x))


def main():
    out = {}
    for m in "AB":
        L = pd.read_csv(Path(sys.argv[1]) / f"labels_{m}_m16.csv").sort_values(["k", "branch"]).reset_index(drop=True)
        lab = np.loadtxt(Path(sys.argv[2]) / f"t5b_labels_{m}_m16.txt", comments="#")
        assert len(L) == len(lab) == 61440
        fam = L["primary_c0.5"].map(FAM).to_numpy()
        band = np.searchsorted(EDGES, L.Omega_cm1.to_numpy(), side="right")
        assert (4 * fam + band == lab[:, 4].astype(int)).all(), "labels csv does not reproduce the T5b label table"
        assert np.allclose(L[["qx", "qy", "qz"]].to_numpy(), lab[:, :3], atol=1e-9)
        assert (L.branch.to_numpy() == lab[:, 3].astype(int)).all()
        q = L[["qx", "qy", "qz"]].to_numpy()
        dG, dRM = dists(q)
        om, br = L.Omega_cm1.to_numpy(), L.branch.to_numpy()
        t5b = json.load(open(Path(sys.argv[3]) / f"t5b_{m}.json"))
        res = {}
        for name, sel in (("nonrot", (band == 0) & (fam != 2)), ("rot", (band == 0) & (fam == 2))):
            assert sel.sum() == t5b["bands"]["<80"][name]["N_modes"]
            r = {"N": int(sel.sum()), "pol": int(((fam == 1) & sel).sum())}
            for part, ps in (("branches 1-3", br <= 3), ("branches >= 4", br >= 4)):
                s = sel & ps
                r[part] = {
                    "N": int(s.sum()),
                    "branches": {int(b): int(((br == b) & s).sum()) for b in np.unique(br[s])},
                    "Omega_cm1": stats(om[s]),
                    "d_Gamma": stats(dG[s]),
                    "d_RM": stats(dRM[s]),
                    "at_Gamma (d_Gamma = 0)": int((s & (dG < 1e-9)).sum()),
                    "d_Gamma <= 0.0625": int((s & (dG <= 0.0625 + 1e-9)).sum()),
                    "d_Gamma <= 0.125": int((s & (dG <= 0.125 + 1e-9)).sum()),
                    "d_Gamma <= 0.1875": int((s & (dG <= 0.1875 + 1e-9)).sum()),
                    "d_Gamma <= 0.25": int((s & (dG <= 0.25 + 1e-9)).sum()),
                    "on R-M line (d_RM = 0)": int((s & (dRM < 1e-9)).sum()),
                    "d_RM <= 0.125": int((s & (dRM <= 0.125 + 1e-9)).sum()),
                    "d_RM <= 0.25": int((s & (dRM <= 0.25 + 1e-9)).sum()),
                }
            res[name] = r
        out[m] = res
    print(json.dumps(out, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
