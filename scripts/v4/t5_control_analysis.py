#!/usr/bin/env python3
"""T5 frequency-matched control. Reads the patched 16^3 run made with the
control labels (label = 4 family + band) and the T4 carriers (primary c0.5); writes, per class and per band, the
kappa_P-weighted share of the carriers' group-averaged linewidth from processes with at least one partner in the
class, the star-weighted mode count, the per-mode share, and the rot / non-rot ratio; decay, absorption and total.

Usage: t5_control_analysis.py MODEL T3_OUT RTA_DIR T5_LABEL_TABLE T4_MODES_CSV OUT_JSON [LABELS_CSV]
Band edges and names are read from the label table's header (T5: 105/145/160; T5b: 80/95/105). With LABELS_CSV
(the frozen labels with Omega, same row order as the label table) the mean, min and max frequency of every class and
per-band union are added (unweighted over the full mesh; T5b, the analysis plan). Prints only label-independent checks.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fekappa.modes import cell_volume_ang3, group_average_gamma, mesh_of  # noqa: E402
from fekappa.partition import RY_TO_CM, read_t3  # noqa: E402
from fekappa.modes import deg_groups  # noqa: E402
from t4_analysis import T_K  # noqa: E402

from fekappa.kappa_decomposition import mode_kappa_table  # noqa: E402

FAMS = ["neither", "pol", "rot"]
BANDS = ["<105", "105-145", "145-160", ">=160"]
NONROT = {"nonrot": (0, 1), "neither": (0,)}


def main():
    model, t3, rdir, labtab, t4csv, out = sys.argv[1], sys.argv[2], Path(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6]
    res = next(rdir.glob("*.result"))
    df = mode_kappa_table(res, nk=mesh_of(rdir) ** 3, V_cell_ang3=cell_volume_ang3(rdir))
    df = df[df.T_K == T_K].copy()
    g = group_average_gamma(df)
    with np.errstate(divide="ignore", invalid="ignore"):
        df["kappa_mode_W_mK"] = np.where(df.kappa_mode_W_mK > 0, df.kappa_mode_W_mK * df.gamma_hwhm_cm1 / g, 0.0)
    df["gamma_avg_cm1"] = g
    df = df.sort_values(["q_index", "branch"]).reset_index(drop=True)
    df["grp"] = np.concatenate([deg_groups(h.omega_cm1.to_numpy()) for _, h in df.groupby("q_index")])

    t4 = pd.read_csv(t4csv)
    t4 = t4[(t4.definition == "primary_c0.5") & t4.carrier]
    car = set(zip(t4.q_index, t4.branch))
    df["carrier"] = [(q, b) in car for q, b in zip(df.q_index, df.branch)]
    assert df.carrier.sum() == len(car)
    # kappa weights must match T4 (same widths)
    kk = t4.set_index(["q_index", "branch"]).kappa_mode_W_mK
    kd = df[df.carrier].set_index(["q_index", "branch"]).kappa_mode_W_mK
    assert np.allclose(kd.loc[kk.index].to_numpy(), kk.to_numpy(), rtol=1e-8)

    global BANDS
    edges = [float(x) for x in re.search(r"band\(([\d./]+) cm", open(labtab).readline()).group(1).split("/")]
    BANDS = ([f"<{edges[0]:g}"] + [f"{a:g}-{b:g}" for a, b in zip(edges, edges[1:])] + [f">={edges[-1]:g}"])
    lab = np.loadtxt(labtab, comments="#")[:, 4].astype(int)
    N = np.bincount(lab, minlength=12)                       # star-weighted counts on the full mesh
    omega = None
    if len(sys.argv) > 7:
        L = pd.read_csv(sys.argv[7]).sort_values(["k", "branch"])
        fam = L["primary_c0.5"].map({"neither": 0, "pol": 1, "rot": 2}).to_numpy()
        chk = 4 * fam + np.searchsorted(edges, L.Omega_cm1.to_numpy(), side="right")
        assert len(chk) == len(lab) and (chk == lab).all(), "LABELS_CSV does not reproduce the label table"
        omega = L.Omega_cm1.to_numpy()

    def fstat(sel):
        if omega is None or not sel.any():
            return None
        o = omega[sel]
        return {"mean": float(o.mean()), "min": float(o.min()), "max": float(o.max())}

    tot, bins = read_t3(t3)
    worst, worst_res = 0.0, 0.0
    rows = []
    for (q, gr), h in df[df.carrier].groupby(["q_index", "grp"]):
        members = df[(df.q_index == q) & (df.grp == gr)].branch.to_numpy()
        tl = np.mean([tot[(q - 1, b - 1, T_K)] for b in members])
        acc = {}
        for b in members:
            for l1, l2, dec, ab in bins.get((q - 1, b - 1, T_K), []):
                acc[(l1, l2)] = acc.get((l1, l2), np.zeros(2)) + np.array([dec, ab]) / len(members)
        s = sum(v.sum() for v in acc.values())
        worst = max(worst, abs(s - tl) / abs(tl))
        worst_res = max(worst_res, abs(tl * RY_TO_CM - h.gamma_avg_cm1.iloc[0]) / h.gamma_avg_cm1.iloc[0])
        sd = sum(v[0] for v in acc.values())
        sa = sum(v[1] for v in acc.values())

        def share(pred):
            v = sum((x for (l1, l2), x in acc.items() if pred(l1) or pred(l2)), np.zeros(2))
            return v[0] / s, v[1] / s, v.sum() / s, (v[0] / sd if sd else np.nan), (v[1] / sa if sa else np.nan)

        for _, m in h.iterrows():
            r = {"q_index": q, "branch": int(m.branch), "w": m.kappa_mode_W_mK}
            for c in range(12):
                r[c] = share(lambda l, c=c: l == c)
            for b in range(4):
                for nm, fs in NONROT.items():
                    r[(nm, b)] = share(lambda l, b=b, fs=fs: l % 4 == b and l // 4 in fs)
            rows.append(r)
    w = np.array([r["w"] for r in rows])
    w = w / w.sum()

    def wmean(key):
        a = np.array([r[key] for r in rows])                 # columns: decay/total, abs/total, total, decay/decay, abs/abs
        return {k: float(np.nansum(w * a[:, i])) for i, k in enumerate(
            ["decay_of_total", "absorption_of_total", "total", "decay_of_decay", "absorption_of_absorption"])}

    res_out = {"model": model, "band_edges_cm1": edges, "carriers": int(df.carrier.sum()), "classes": {}, "bands": {},
               "check_bins_vs_total_max_rel": worst, "check_total_vs_result_avg_max_rel": worst_res}
    for c in range(12):
        f, b = divmod(c, 4)
        if f > 2:
            continue
        S = wmean(c)
        res_out["classes"][f"{FAMS[f]}|{BANDS[b]}"] = {"N_modes": int(N[c]), "S": S,
                                                         "per_mode": {k: (v / N[c] if N[c] else None) for k, v in S.items()},
                                                         "omega_cm1": fstat(lab == c)}
    for b in range(4):
        rot = res_out["classes"][f"rot|{BANDS[b]}"]
        e = {"rot": rot}
        for nm, fs in NONROT.items():
            Nn = int(sum(N[4 * f + b] for f in fs))
            S = wmean((nm, b))
            pm = {k: (v / Nn if Nn else None) for k, v in S.items()}
            e[nm] = {"N_modes": Nn, "S": S, "per_mode": pm,
                     "omega_cm1": fstat(np.isin(lab, [4 * f + b for f in fs]))}
            e[f"ratio_rot_over_{nm}"] = {k: (rot["per_mode"][k] / pm[k] if rot["per_mode"][k] is not None and pm[k]
                                             else None) for k in pm}
        res_out["bands"][BANDS[b]] = e
    json.dump(res_out, open(out, "w"), indent=1)
    print(json.dumps({k: res_out[k] for k in ("model", "carriers", "check_bins_vs_total_max_rel",
                                              "check_total_vs_result_avg_max_rel")}))


if __name__ == "__main__":
    main()
