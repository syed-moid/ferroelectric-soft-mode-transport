#!/usr/bin/env python3
"""Count-only sampling test of the band ratios.

The partition records linewidth shares by class, not by partner mode, so the weight of one symmetry star cannot be
removed exactly. Test: list the stars of each class in a band and every star within 2 cm^-1 of a band edge; move
the nearest rotational star, and separately the nearest non-rotational star, across the edge (counts change,
shares S unchanged) and recompute ratio = (S_rot / N_rot) / (S_nonrot / N_nonrot). If either move takes the ratio
across 2.0 or 1.3, the band's result is "sensitive to sampling".

Unit moved: one star-group = all modes of one degenerate group (same Omega) at the wave vectors of one star of the
16^3 mesh (cubic point group: signed permutations of q, modulo 1); its size is the number of modes it contains.
Family: primary definition c_min 0.5 (frozen labels). Non-rot = neither + pol.

Usage: p3_sampling.py MODEL LABELS_CSV T5_JSON T5B_JSON OUT_JSON
       (bands: '<105' from T5_JSON; '<80', '80-95', '95-105' from T5B_JSON)
"""
import itertools
import json
import sys

import numpy as np
import pandas as pd

THR = (2.0, 1.3)
NEAR = 2.0


def star_key(q):
    q = np.round(np.asarray(q, float) * 16).astype(int) % 16
    best = None
    for perm in itertools.permutations(range(3)):
        for sg in itertools.product((1, -1), repeat=3):
            v = tuple(int((s * q[p]) % 16) for p, s in zip(perm, sg))
            best = v if best is None or v < best else best
    return best


def cls(r):
    return "rot" if r == "rot" else "nonrot"


def band_class(r):
    return None if r is None else ("≥ 2.0" if r >= THR[0] else ("1.3 to < 2.0" if r >= THR[1] else "< 1.3"))


def main():
    model, lab, j5, j5b, out = sys.argv[1:6]
    df = pd.read_csv(lab)
    df["star"] = [star_key(q) for q in df[["qx", "qy", "qz"]].to_numpy()]
    df["fam"] = df["primary_c0.5"].map(cls)
    df["Om"] = df.Omega_cm1.round(4)
    sg = (df.groupby(["star", "Om", "fam"]).size().rename("size").reset_index())   # star-groups
    J = {"<105": (json.load(open(j5))["bands"]["<105"], (None, 105.0)),
         "<80": (json.load(open(j5b))["bands"]["<80"], (None, 80.0)),
         "80-95": (json.load(open(j5b))["bands"]["80-95"], (80.0, 95.0)),
         "95-105": (json.load(open(j5b))["bands"]["95-105"], (95.0, 105.0))}
    R = {"model": model, "unit": "star-group (degenerate group x star)", "bands": {}}
    for b, (e, (lo, hi)) in J.items():
        S = {c: e[c]["S"]["total"] for c in ("rot", "nonrot")}
        N = {c: e[c]["N_modes"] for c in ("rot", "nonrot")}
        inb = sg[((sg.Om >= lo) if lo is not None else True) & (sg.Om < hi)]
        Ncsv = {c: int(inb[inb.fam == c]["size"].sum()) for c in ("rot", "nonrot")}
        ratio = (S["rot"] / N["rot"]) / (S["nonrot"] / N["nonrot"]) if N["rot"] and N["nonrot"] and S["nonrot"] else None
        rec = {"S_total": S, "N": N, "N_from_labels": Ncsv, "ratio": ratio, "class": band_class(ratio),
               "stars_in_band": {c: {"n_star_groups": int((inb.fam == c).sum()),
                                     "sizes": sorted(inb[inb.fam == c]["size"].tolist())} for c in ("rot", "nonrot")},
               "edges": {}, "moves": [], "sensitive_to_sampling": False}
        for edge in [x for x in (lo, hi) if x is not None]:
            near = sg[(sg.Om - edge).abs() <= NEAR].sort_values("Om")
            rec["edges"][str(edge)] = [{"Omega": float(r.Om), "family": r.fam, "size": int(r["size"]),
                                        "side": "inside" if ((lo is None or r.Om >= lo) and r.Om < hi) else "outside"}
                                       for _, r in near.iterrows()]
            for c in ("rot", "nonrot"):
                cand = sg[sg.fam == c].copy()
                cand["d"] = (cand.Om - edge).abs()
                cand = cand.sort_values(["d", "Om"])
                if cand.empty:
                    continue
                r = cand.iloc[0]
                inside = (lo is None or r.Om >= lo) and r.Om < hi
                n2 = dict(N)
                n2[c] += -int(r["size"]) if inside else int(r["size"])
                r2 = ((S["rot"] / n2["rot"]) / (S["nonrot"] / n2["nonrot"])
                      if n2["rot"] > 0 and n2["nonrot"] > 0 and S["nonrot"] else None)
                crosses = (ratio is not None and r2 is not None and
                           any((ratio >= t) != (r2 >= t) for t in THR))
                rec["moves"].append({"edge": edge, "class_moved": c, "star_group_Omega": float(r.Om),
                                     "size": int(r["size"]), "direction": "out of band" if inside else "into band",
                                     "distance_cm1": float(r.d), "N_after": n2, "ratio_after": r2,
                                     "crosses_2.0_or_1.3": bool(crosses)})
                rec["sensitive_to_sampling"] |= bool(crosses)
        R["bands"][b] = rec
    json.dump(R, open(out, "w"), indent=1)
    for b, rec in R["bands"].items():
        print(f"{model} {b:7s} ratio {rec['ratio']:.3f} [{rec['class']}] N {rec['N']} (labels {rec['N_from_labels']}) "
              f"sensitive={rec['sensitive_to_sampling']}")
        for m in rec["moves"]:
            print(f"   move {m['class_moved']:6s} star-group Ω {m['star_group_Omega']:.2f} (size {m['size']}, "
                  f"{m['direction']}, edge {m['edge']:g}, d {m['distance_cm1']:.2f}) -> ratio {m['ratio_after']:.3f}"
                  f"{'  CROSSES' if m['crosses_2.0_or_1.3'] else ''}")


if __name__ == "__main__":
    main()
