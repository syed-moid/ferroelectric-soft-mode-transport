#!/usr/bin/env python3
"""P0: checks of the archived 2015 files against the distributed ones. HOME, no run. Facts only.

Usage: p0_archive_files.py ARCH_XML DIST_FC2_XML DIST_ANHARM_XML[.bz2] ARCH_BORN DIST_BORN TAU_DAT OUR_12x12x12_RTA_DIR OUT_JSON
Prints Markdown sections for docs/v4/ANALYSIS_STEPS.md, step P0 (items 1-4 and 6; item 5 is a node run).
"""
import bz2
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fekappa.modes import TOL_DEG_CM1, cell_volume_ang3, mesh_of  # noqa: E402

from fekappa.kappa_decomposition import mode_kappa_table  # noqa: E402

ORDERS = ["HARMONIC", "ANHARM3", "ANHARM4", "ANHARM5", "ANHARM6"]


def load(path):
    p = str(path)
    if p.endswith(".bz2"):
        with bz2.open(p) as f:
            return ET.parse(f).getroot()
    return ET.parse(p).getroot()


def fc_dict(root, tag):
    node = root.find(f"ForceConstants/{tag}")
    if node is None:
        return None
    d = {}
    for c in node:
        key = tuple(c.attrib[k] for k in sorted(c.attrib) if k.startswith("pair"))
        d[key] = d.get(key, 0.0) + float(c.text)
    return d, len(node)


def structure(root):
    s = root.find("Structure")
    lat = [[float(x) for x in s.find(f"LatticeVector/a{i}").text.split()] for i in (1, 2, 3)]
    pos = [(p.attrib["element"], [float(x) for x in p.text.split()]) for p in s.find("Position")]
    tr = [(t.attrib["tran"], t.attrib["atom"], t.text.strip()) for t in root.find("Symmetry/Translations")]
    return np.array(lat), pos, tr


def cmp_fc(a, b):
    keys = set(a) | set(b)
    only_a, only_b = len(set(a) - set(b)), len(set(b) - set(a))
    va = np.array([a.get(k, 0.0) for k in keys])
    vb = np.array([b.get(k, 0.0) for k in keys])
    d = np.abs(va - vb)
    i = int(np.argmax(d))
    big = np.abs(vb) > 1e-6
    rel = d[big] / np.abs(vb[big])
    return {"n_keys": len(keys), "only_first": only_a, "only_second": only_b, "max_abs_Ry_bohr2": float(d.max()),
            "max_abs_at": list(list(keys)[i]), "value_there": [float(va[i]), float(vb[i])],
            "max_rel_for_abs_gt_1e-6": float(rel.max()) if rel.size else None,
            "largest_abs_fc": float(np.abs(vb).max())}


def read_born(p):
    return np.array([[float(x) for x in l.split()] for l in open(p) if l.strip()])


def blocks_from_tau(tau, n_q, nb=15):
    """Row blocks under the hypothesis: irreducible points in order, branches in order, Γ acoustic (3) omitted."""
    sizes = [nb - 3] + [nb] * (n_q - 1)
    out, i = [], 0
    for s in sizes:
        out.append(tau[i:i + s])
        i += s
    return out, i


def deg_pattern(w):
    pat, start = [], 0
    for j in range(1, len(w) + 1):
        if j == len(w) or abs(w[j] - w[start]) >= max(TOL_DEG_CM1, 1e-6 * abs(w[start])):
            pat.append(j - start)
            start = j
    return pat


def main():
    arch, fc2, anh, born_a, born_d, tau_p, rdir, out = sys.argv[1:9]
    R = {}
    ra, rf, rd = load(arch), load(fc2), load(anh)
    L = []
    w = L.append

    # 1. header
    fit = ra.find("Fitting")
    R["header"] = {"ALM_version": ra.findtext("ALM_version"),
                   "Fitting": {c.tag: c.text for c in fit} if fit is not None else None}
    for nm, r in (("distributed STO_anharm.xml", rd), ("distributed STO222.xml", rf)):
        f = r.find("Fitting")
        R.setdefault("header_others", {})[nm] = {"ALM_version": r.findtext("ALM_version"),
                                                 "Fitting": {c.tag: c.text for c in f} if f is not None else None}

    # structure / symmetry
    la, pa, ta = structure(ra)
    for nm, r in (("STO222.xml", rf), ("STO_anharm.xml", rd)):
        lb, pb, tb = structure(r)
        R.setdefault("structure_vs", {})[nm] = {
            "lattice_max_abs_diff_bohr": float(np.abs(la - lb).max()),
            "n_atoms": [len(pa), len(pb)],
            "elements_equal_in_order": [e for e, _ in pa] == [e for e, _ in pb],
            "positions_max_abs_diff": float(max(np.abs(np.array(x) - np.array(y)).max() for (_, x), (_, y) in zip(pa, pb))),
            "translations_equal": ta == tb, "n_translations": [len(ta), len(tb)]}

    # 2. harmonic
    ha, na = fc_dict(ra, "HARMONIC")
    R["harmonic"] = {}
    for nm, r in (("STO222.xml", rf), ("STO_anharm.xml", rd)):
        hb, nb_ = fc_dict(r, "HARMONIC")
        R["harmonic"][nm] = dict(cmp_fc(ha, hb), entries=[na, nb_])

    # 3. entries per order
    R["entries"] = {}
    for nm, r in (("archived", ra), ("distributed STO_anharm.xml", rd), ("distributed STO222.xml", rf)):
        R["entries"][nm] = {}
        for t in ORDERS:
            x = fc_dict(r, t)
            R["entries"][nm][t] = None if x is None else {"entries": x[1], "distinct_keys": len(x[0])}

    # 4. Born
    ba, bd = read_born(born_a), read_born(born_d)
    R["born"] = {"rows": [ba.tolist(), bd.tolist()], "diff": (ba - bd).tolist(),
                 "max_abs_diff": float(np.abs(ba - bd).max()),
                 "eps_inf_diag": [ba[:3].diagonal().tolist(), bd[:3].diagonal().tolist()]}

    # 6. tau rows against our irreducible 12^3 list
    tau = np.loadtxt(tau_p)
    rdir = Path(rdir)
    res = next(rdir.glob("*.result"))
    df = mode_kappa_table(res, nk=mesh_of(rdir) ** 3, V_cell_ang3=cell_volume_ang3(rdir))
    df = df[df.T_K == 300].sort_values(["q_index", "branch"])
    nq = df.q_index.nunique()
    blocks, used = blocks_from_tau(tau, nq)
    ascending, pat_equal, pat_diff, dw = 0, 0, [], []
    per_q = []
    for (qi, g), b in zip(df.groupby("q_index"), blocks):
        ours = g.omega_cm1.to_numpy()
        if qi == df.q_index.min():
            ours = ours[3:]          # Γ: acoustic omitted in the file
        asc = bool(np.all(np.diff(b[:, 0]) >= -1e-9))
        ascending += asc
        pa_, po = deg_pattern(b[:, 0]), deg_pattern(ours)
        eq = pa_ == po
        pat_equal += eq
        if not eq:
            pat_diff.append({"q_index": int(qi), "q": g.iloc[0][["q_frac_x", "q_frac_y", "q_frac_z"]].tolist(),
                             "file": pa_, "ours": po})
        dw.append(b[:, 0] - ours)
        per_q.append({"q_index": int(qi), "q": [float(x) for x in g.iloc[0][["q_frac_x", "q_frac_y", "q_frac_z"]]],
                      "pattern_file": pa_, "pattern_ours": po, "min_w_file": float(b[0, 0]), "min_w_ours": float(ours[0])})
    # order test: points sharing a degeneracy pattern (ours) cannot be told apart by pattern; within each such set,
    # the optimal one-to-one assignment of file blocks to our points by RMS frequency difference is computed and
    # compared with the positional (identity) assignment. A check of the order hypothesis, not an assignment rule.
    from collections import defaultdict
    from scipy.optimize import linear_sum_assignment
    W = {int(q): g.omega_cm1.to_numpy() for q, g in df.groupby("q_index")}
    q0 = int(df.q_index.min())
    F = {}
    for (qi, _), b in zip(df.groupby("q_index"), blocks):
        F[int(qi)] = np.r_[[np.nan] * 3, b[:, 0]] if int(qi) == q0 else b[:, 0]
    sets = defaultdict(list)
    for p in per_q:
        sets[tuple(p["pattern_ours"])].append(p["q_index"])
    order = []
    for sig, v in sets.items():
        C = np.array([[np.sqrt(np.nanmean((F[i] - W[j]) ** 2)) for j in v] for i in v])
        _, col = linear_sum_assignment(C)
        order.append({"pattern": list(sig), "n_points": len(v), "identity_is_optimal": bool((col == np.arange(len(v))).all()),
                      "moved": int((col != np.arange(len(v))).sum())})
    own = np.array([[np.sqrt(np.nanmean((F[i] - W[j]) ** 2)) for j in sorted(W)] for i in sorted(W)])
    dw = np.concatenate(dw)
    ratio = tau[:, 1] / tau[:, 0]
    k = int(np.argmax(ratio))
    tau_rel = np.abs(tau[:, 2] * 2 * tau[:, 1] * 2.99792458e10 * 2 * np.pi * 1e-12 - 1)  # τ(ps) = 1/(2Γ), Γ in cm^-1 → ω units
    R["tau"] = {"rows": int(len(tau)), "n_irred_q_ours": int(nq), "rows_expected": int(nq * 15 - 3), "rows_used": used,
                "blocks_ascending": int(ascending), "blocks_same_degeneracy_pattern": int(pat_equal),
                "pattern_differences": pat_diff,
                "order_test": {"n_patterns": len(sets), "points_unique_pattern": sum(len(v) == 1 for v in sets.values()),
                               "sets": order, "rms_own_point_median_cm1": float(np.median(np.diag(own))),
                               "rms_own_point_max_cm1": float(np.diag(own).max()),
                               "nearest_point_is_own_point": int((own.argmin(1) == np.arange(len(own))).sum())},
                "freq_diff_file_minus_ours_cm1": {"median_abs": float(np.median(np.abs(dw))),
                                                  "q95_abs": float(np.quantile(np.abs(dw), 0.95)),
                                                  "max_abs": float(np.abs(dw).max())},
                "max_gamma_over_omega": float(ratio[k]), "max_at_row": k + 1,
                "max_at": {"omega_cm1": float(tau[k, 0]), "gamma_cm1": float(tau[k, 1])},
                "tau_eq_1_over_2Gamma_max_rel_dev": float(tau_rel.max()),
                "reference_run": str(rdir), "per_q": per_q}
    json.dump(R, open(out, "w"), indent=1)
    print(json.dumps({k: v for k, v in R.items() if k not in ("tau",)}, indent=1))
    print(json.dumps({k: v for k, v in R["tau"].items() if k != "per_q"}, indent=1))


if __name__ == "__main__":
    main()
