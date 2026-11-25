#!/usr/bin/env python3
"""Tables for the T2 census from t2_family_census.py output. Characters and frequencies only (no kappa,
linewidth or eta).

Proposed label (parameters c_min and the cross-check thresholds are arguments, for the author to fix):
  pol   if C_pol >= c_min and C_rot <  c_min
  rot   if C_rot >= c_min and C_pol <  c_min
  mixed if both  >= c_min
  other otherwise ("neither")
Translational modes (branches 1-3 at q = 0) are excluded.
D1 (harmonically unstable): s_unst > 0.5.   D3 (renormalized): r_ren > r_min.

Polar character: --polar=plain (default; the census's C_pol) or --polar=weighted (the analysis plan primary definition:
C_pol^w = w_pol(q) C_pol, w_pol = cos^2(pi qx) cos^2(pi qy) cos^2(pi qz); the plain value is kept as C_pol_plain).
Same code and no model-specific parameter for every model.

Usage: t2_census_tables.py [--polar=plain|weighted] CENSUS_PREFIX OUT_PREFIX C_MIN R_MIN LABEL [LABEL ...]
Writes OUT_PREFIX_*.csv and prints markdown tables.
"""
import sys

import numpy as np
import pandas as pd

PATH = [("G", (0, 0, 0)), ("X", (.5, 0, 0)), ("M", (.5, .5, 0)), ("G", (0, 0, 0)), ("R", (.5, .5, .5)),
        ("M", (.5, .5, 0))]


def label(d, c_min):
    p, r = d.C_pol >= c_min, d.C_rot >= c_min
    return np.select([p & ~r, r & ~p, p & r], ["pol", "rot", "mixed"], "other")


def path_points(n=16):
    pts = []
    for (a, qa), (b, qb) in zip(PATH[:-1], PATH[1:]):
        qa, qb = np.array(qa), np.array(qb)
        steps = int(round(np.abs(qb - qa).max() * n))
        for s in range(steps + (1 if b == PATH[-1][0] and (a, b) == ("R", "M") else 0)):
            pts.append((f"{a}-{b}", s, qa + (qb - qa) * s / steps))
    return pts


def main():
    opts = [a for a in sys.argv[1:] if a.startswith("--polar=")]
    args = [a for a in sys.argv[1:] if not a.startswith("--polar=")]
    polar = opts[-1].split("=", 1)[1] if opts else "plain"
    assert polar in ("plain", "weighted"), polar
    cpre, opre, c_min, r_min, labs = args[0], args[1], float(args[2]), float(args[3]), args[4:]
    for lab in labs:
        d = pd.read_csv(f"{cpre}_{lab}.csv")
        d["C_pol_plain"] = d.C_pol
        if polar == "weighted":
            d["w_pol"] = np.prod(np.cos(np.pi * d[["qx", "qy", "qz"]].to_numpy()) ** 2, axis=1)
            d["C_pol"] = d.w_pol * d.C_pol_plain
        d = d[~((d.d_Gamma < 1e-9) & (d.branch <= 3))].copy()
        d["label"] = label(d, c_min)
        d["D1"] = d.s_unst > 0.5
        d["D3"] = d.r_ren > r_min
        print(f"\n### model {lab}  (polar = {polar}, c_min = {c_min}, D3: r_ren > {r_min})")

        # symmetry: labels and characters identical across each star (compared after sorting by frequency)
        bad = 0
        for s, g in d.groupby("star"):
            ref = None
            for k, h in g.groupby("k"):
                h = h.sort_values("Omega_cm1")
                sig = (np.round(h.Omega_cm1.to_numpy(), 3), h.label.to_numpy(),
                       h[["C_pol", "C_rot"]].to_numpy())
                if ref is None:
                    ref = sig
                elif not (np.allclose(sig[0], ref[0], atol=2e-3) and np.array_equal(sig[1], ref[1])
                          and np.allclose(sig[2], ref[2], atol=1e-4)):
                    bad += 1
        nstar = d.star.nunique()
        print(f"symmetry check: {nstar} stars, members differing from their star's first member: {bad}")

        # counts: star-weighted (full mesh) and irreducible (one representative per star)
        rep = d[d.k == d.groupby("star").k.transform("min")]
        tab = pd.DataFrame({"irreducible": rep.label.value_counts(), "star_weighted": d.label.value_counts()})
        fr = d.groupby("label").Omega_cm1.agg(["min", "median", "max"]).round(2)
        tab = tab.join(fr).fillna(0)
        print(tab.to_string())
        tab.to_csv(f"{opre}_{lab}_counts.csv")

        # recovery at Gamma and R
        for nm, q in (("Gamma", (0, 0, 0)), ("R", (.5, .5, .5)), ("M", (.5, .5, 0)), ("X", (.5, 0, 0))):
            h = d[np.isclose(d.qx, q[0]) & np.isclose(d.qy, q[1]) & np.isclose(d.qz, q[2])].sort_values("Omega_cm1")
            print(f"{nm}: " + "; ".join(f"{w:.1f}:{l}(p{cp:.2f},r{cr:.2f},rp{crp:.2f})"
                                        for w, l, cp, cr, crp in zip(h.Omega_cm1, h.label, h.C_pol, h.C_rot,
                                                                     h.C_rot_plain)))

        # along the path: per point, labelled branches and their frequencies
        rows = []
        for seg, s, q in path_points():
            h = d[np.isclose(d.qx, q[0]) & np.isclose(d.qy, q[1]) & np.isclose(d.qz, q[2])].sort_values("Omega_cm1")
            if h.empty:  # translational modes removed at Gamma only; points always exist on the 16^3 mesh
                continue
            rows.append({"segment": seg, "step": s, "q": "%.4f %.4f %.4f" % tuple(q),
                         "pol_cm1": " ".join(f"{w:.1f}" for w in h[h.label == "pol"].Omega_cm1),
                         "rot_cm1": " ".join(f"{w:.1f}" for w in h[h.label == "rot"].Omega_cm1),
                         "mixed_cm1": " ".join(f"{w:.1f}" for w in h[h.label == "mixed"].Omega_cm1),
                         "max_C_pol": h.C_pol.max(), "max_C_rot": h.C_rot.max()})
        pdf = pd.DataFrame(rows)
        pdf.to_csv(f"{opre}_{lab}_path.csv", index=False, float_format="%.4g")

        # counts against distance from Gamma and from the nearest R-M line (star-weighted)
        for col, bins in (("d_Gamma", [0, .1, .2, .3, .4, .5, .6, .7, .9]), ("d_RM", [-1e-9, .05, .1, .2, .3, .4, .5, .8])):
            t = pd.crosstab(pd.cut(d[col], bins), d.label)
            t.to_csv(f"{opre}_{lab}_vs_{col}.csv")
            print(f"\nstar-weighted counts against {col}:\n" + t.to_string())

        # mixed / ambiguous: where (near threshold = within 0.1 of c_min for the larger character)
        amb = d[(np.maximum(d.C_pol, d.C_rot) - c_min).abs() < 0.1]
        print(f"\nmixed: {int((d.label == 'mixed').sum())} star-weighted; "
              f"within 0.1 of c_min: {len(amb)} (labels {amb.label.value_counts().to_dict()}); "
              f"mixed at median d_Gamma {d[d.label == 'mixed'].d_Gamma.median():.3f}, "
              f"d_RM {d[d.label == 'mixed'].d_RM.median():.3f}")

        # confusion with D1 and D3 (star-weighted)
        for D in ("D1", "D3"):
            c = pd.crosstab(d.label, d[D])
            c.to_csv(f"{opre}_{lab}_confusion_{D}.csv")
            print(f"\nlabel x {D}:\n" + c.to_string())
        d[["k", "branch", "Omega_cm1", "label", "D1", "D3"]].to_csv(f"{opre}_{lab}_labels_cmin{c_min}.csv",
                                                                    index=False, float_format="%.6g")


if __name__ == "__main__":
    main()
