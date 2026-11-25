#!/usr/bin/env python3
"""T4: who scatters the carriers. One patched 16^3 run per model with the frozen 32-label
code; every definition and threshold is obtained by summing bins.

Steps
1. Mode table from the run's .result with the widths averaged over degenerate groups exactly as anphon does before
   forming kappa (kp_eta_distribution.group_average_gamma); kappa per mode rescaled accordingly (sum = solver kappa_P).
2. Labels of the irreducible modes from the frozen table (same q on the full mesh; frequencies checked).
3. Carriers (frozen rule): irreducible modes ranked by their symmetry-weighted kappa_P contribution, taken until
   >= 50 % of kappa_P; the whole degenerate group at the boundary.
4. For each carrier: bins of the partner-label pair (decay, absorption) -> by family pair and by the number of
   soft partners (0 / 1 / 2; soft = pol, rot or mixed); averaged over the carrier's degenerate group as the code
   does. Per-carrier check: sum of averaged bins = averaged total (<= 1e-10), and = group-averaged .result width.
5. kappa_P carried by each family (all modes) and by modes with eta > 0.1 (group-averaged), opened together.

Definitions: primary c_min 0.5 (main), secondary c_min 0.5, primary c_min 0.4 and 0.6 (sensitivity).
Usage: t4_analysis.py MODEL T3_OUT RTA_DIR LABELS_CSV OUTDIR
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fekappa.modes import TOL_DEG_CM1, cell_volume_ang3, deg_groups, group_average_gamma, mesh_of  # noqa: E402,F401
from t2_freeze_labels import labels_from_code  # noqa: E402
from fekappa.partition import RY_TO_CM, read_t3  # noqa: E402

from fekappa.evec_io import irred_q, kmap  # noqa: E402
from fekappa.kappa_decomposition import mode_kappa_table  # noqa: E402

DEFS = [("primary", 0.5), ("secondary", 0.5), ("primary", 0.4), ("primary", 0.6)]
FAMS = ["neither", "pol", "rot", "mixed"]
SOFT = {"pol", "rot", "mixed"}
T_K = float(os.environ.get("T4_T_K", "300"))   # 300 K unless set (T7: 600)


def main():
    model, t3, rdir, labcsv, out = sys.argv[1], sys.argv[2], Path(sys.argv[3]), sys.argv[4], Path(sys.argv[5])
    out.mkdir(parents=True, exist_ok=True)
    res = next(rdir.glob("*.result"))
    mesh = mesh_of(rdir)

    # 1. mode table with group-averaged widths
    df = mode_kappa_table(res, nk=mesh ** 3, V_cell_ang3=cell_volume_ang3(rdir))
    df = df[df.T_K == T_K].copy()
    g = group_average_gamma(df)
    with np.errstate(divide="ignore", invalid="ignore"):
        df["kappa_mode_W_mK"] = np.where(df.kappa_mode_W_mK > 0, df.kappa_mode_W_mK * df.gamma_hwhm_cm1 / g, 0.0)
    df["gamma_raw_cm1"] = df.gamma_hwhm_cm1
    df["gamma_avg_cm1"] = g
    df["eta"] = np.where(df.omega_cm1 > 1e-3, g / df.omega_cm1, np.nan)
    df = df.sort_values(["q_index", "branch"]).reset_index(drop=True)
    df["grp"] = np.concatenate([deg_groups(h.omega_cm1.to_numpy()) for _, h in df.groupby("q_index")])
    kP = df.kappa_mode_W_mK.sum()

    # 2. labels of the irreducible modes
    lab = pd.read_csv(labcsv)
    xk = lab.drop_duplicates("k").sort_values("k")[["qx", "qy", "qz"]].to_numpy()
    km = kmap(xk, irred_q(res))
    L = lab.set_index(["k", "branch"])
    keys = list(zip(km[df.q_index.to_numpy() - 1], df.branch.to_numpy()))
    df["code32"] = L.loc[keys, "code32"].to_numpy()
    dom = np.abs(L.loc[keys, "Omega_cm1"].to_numpy() - df.omega_cm1.to_numpy())
    assert dom[df.omega_cm1 > 1e-3].max() < 2e-3, dom.max()
    fam_of = {}
    for d, c in DEFS:
        name = f"{d}_c{c}"
        fam_of[name] = dict(zip(range(32), labels_from_code(np.arange(32), d, c)))
        df[f"fam_{name}"] = [fam_of[name][x] for x in df.code32]
        if (d, c) != ("secondary", 0.5):
            assert (df[f"fam_{name}"].to_numpy() == L.loc[keys, f"{d}_c{c}"].to_numpy()).all()

    # 3. carriers
    gk = df.groupby(["q_index", "grp"]).kappa_mode_W_mK.sum().sort_values(ascending=False)
    cum = gk.cumsum() / kP
    take = cum.index[: int(np.searchsorted(cum.to_numpy(), 0.5)) + 1]
    sel = set(take)
    df["carrier"] = [(q, gr) in sel for q, gr in zip(df.q_index, df.grp)]
    rank = {k: i + 1 for i, k in enumerate(gk.index)}
    df["group_rank"] = [rank[(q, gr)] for q, gr in zip(df.q_index, df.grp)]
    df["group_cum_share"] = [cum[(q, gr)] for q, gr in zip(df.q_index, df.grp)]

    # 4. partition per carrier (all modes computed; carriers flagged)
    tot, bins = read_t3(t3)
    rows, worst, worst_res = [], 0.0, 0.0
    for (q, gr), h in df.groupby(["q_index", "grp"]):
        members = h.branch.to_numpy()
        tl = np.mean([tot[(q - 1, b - 1, T_K)] for b in members])
        acc = {}
        for b in members:
            for l1, l2, dec, ab in bins.get((q - 1, b - 1, T_K), []):
                acc[(l1, l2)] = acc.get((l1, l2), np.zeros(2)) + np.array([dec, ab]) / len(members)
        s = sum(v.sum() for v in acc.values())
        if tl != 0:
            worst = max(worst, abs(s - tl) / abs(tl))
            worst_res = max(worst_res, abs(tl * RY_TO_CM - h.gamma_avg_cm1.iloc[0]) / h.gamma_avg_cm1.iloc[0])
        for _, m in h.iterrows():
            base = {"model": model, "q_index": q, "branch": int(m.branch), "q": "%.4f %.4f %.4f" % (
                m.q_frac_x, m.q_frac_y, m.q_frac_z), "omega_cm1": m.omega_cm1, "gamma_avg_cm1": m.gamma_avg_cm1,
                "eta": m.eta, "kappa_mode_W_mK": m.kappa_mode_W_mK, "carrier": bool(m.carrier),
                "group_rank": m.group_rank, "group_cum_share": m.group_cum_share, "code32": int(m.code32)}
            for name, fmap in fam_of.items():
                r = dict(base, definition=name, carrier_family=fmap[int(m.code32)])
                parts = {f"{a}/{b}": np.zeros(2) for i, a in enumerate(FAMS) for b in FAMS[i:]}
                ns = np.zeros((3, 2))
                for (l1, l2), v in acc.items():
                    f1, f2 = sorted((fmap[l1], fmap[l2]), key=FAMS.index)
                    parts[f"{f1}/{f2}"] += v
                    ns[(f1 in SOFT) + (f2 in SOFT)] += v
                r["gamma_partition_cm1"] = s * RY_TO_CM
                for k in range(3):
                    r[f"nsoft{k}_decay"], r[f"nsoft{k}_absorption"] = ns[k] * RY_TO_CM
                    r[f"nsoft{k}_frac"] = ns[k].sum() / s if s else np.nan
                for k, v in parts.items():
                    r[f"pair_{k}_decay"], r[f"pair_{k}_absorption"] = v * RY_TO_CM
                    r[f"pair_{k}_frac"] = v.sum() / s if s else np.nan
                rows.append(r)
    P = pd.DataFrame(rows)
    P.to_csv(out / f"t4_{model}_modes.csv", index=False, float_format="%.10g")

    # 5. summaries
    summ = {"model": model, "kappa_P_modesum": kP, "n_irred_modes": len(df),
            "carriers": int(df.carrier.sum()), "carrier_groups": len(sel),
            "carrier_kappa_share": float(df[df.carrier].kappa_mode_W_mK.sum() / kP),
            "check_bins_vs_total_max_rel": worst, "check_total_vs_result_avg_max_rel": worst_res,
            "share_kP_eta_gt_0.1": float(df[df.eta > 0.1].kappa_mode_W_mK.sum() / kP), "definitions": {}}
    for name in fam_of:
        p = P[(P.definition == name) & P.carrier]
        w = p.kappa_mode_W_mK / p.kappa_mode_W_mK.sum()
        fam_k = df.groupby(f"fam_{name}").kappa_mode_W_mK.sum() / kP
        car = p.carrier_family.value_counts().to_dict()
        summ["definitions"][name] = {
            "kappa_P_share_by_family": {f: float(fam_k.get(f, 0.0)) for f in FAMS},
            "carrier_family_counts": {f: int(car.get(f, 0)) for f in FAMS},
            "kappa_weighted_nsoft_frac": {k: float((w * p[f"nsoft{k}_frac"]).sum()) for k in range(3)},
            "median_nsoft_frac": {k: float(p[f"nsoft{k}_frac"].median()) for k in range(3)},
            "kappa_weighted_pair_frac": {c[5:-5]: float((w * p[c]).sum()) for c in p.columns
                                         if c.startswith("pair_") and c.endswith("_frac")},
            "kappa_weighted_absorption_share": float((w * sum(p[f"nsoft{k}_absorption"] for k in range(3))
                                                      / p.gamma_partition_cm1).sum()),
        }
    json.dump(summ, open(out / f"t4_{model}_summary.json", "w"), indent=1)
    # label-independent checks only; family results are opened in the T4 report tables (Markdown)
    print(json.dumps({k: summ[k] for k in ("model", "kappa_P_modesum", "n_irred_modes", "carriers",
                                           "carrier_groups", "carrier_kappa_share", "check_bins_vs_total_max_rel",
                                           "check_total_vs_result_avg_max_rel")}))


if __name__ == "__main__":
    main()
