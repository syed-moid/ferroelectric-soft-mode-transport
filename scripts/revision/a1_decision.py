"""Apply the pre-registered four-phonon mesh rule (docs/recorded_rules.md (R4), written 2026-10-15T17:25Z) to the
spec_A1 results (300 K) or, with A1_T=1000, to the spec_A2 results under R5 (same rule, 1000 K thresholds). Reads data/provenance/manifest_final.json and the production mode table.
Outputs data/processed/revision/a1_per_mesh.csv (counters, Gamma3 calibration, Gamma4/Gamma3 per mesh and eps)
and a1_decision.csv (stability S1/S2, Gamma4_max, outcome per mode). The thresholds are the pre-registered
numbers; they are recomputed from production data only as a consistency check."""
import json
import os
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PJ = ROOT / "data/processed/revision"
man = json.loads(Path(os.environ.get("A1_MANIFEST", ROOT / "data/provenance/manifest_final.json")).read_text())["jobs"]  # override: dry check on a pulled node manifest
modes_tab = pd.read_csv(ROOT / "data/processed/sto_production/mode_table_T.csv.gz", comment="#")

T = int(os.environ.get("A1_T", 300))
MODES = {"0,0,0,4": "Gamma-TO1", "0.5,0.5,0.5,1": "R", "0.5,0.5,0,1": "M", "0,0,0.5,7": "X-b7", "0.5,0.5,0.5,13": "R-b13"}
THR = {300: {"Gamma-TO1": 64.6, "M": 36.5, "R": 8.7},  # pre-registered (R4)
       1000: {"Gamma-TO1": 85.6, "M": 53.2, "R": 39.1, "X-b7": 108.7, "R-b13": 232.1}}[T]  # pre-registered (R5)
SUF = "" if T == 300 else f"_T{T}"
RULE_JOBS = [(6, 2), (6, 5), (6, 10), (8, 5)]
NUM = r"([0-9]+)"
PAT = re.compile(rf"nk = {NUM}, excluded = {NUM} \(external mode itself on the leg: {NUM}\), "
                 rf"denominator check = {NUM}, near-misses kept \(< 1e-6 Ry\) = {NUM}")


def prod_mode(key):
    qx, qy, qz, b = (float(x) for x in key.split(","))
    t = modes_tab[(modes_tab.mesh == 16) & (modes_tab.T_K == T) & (modes_tab.branch == int(b))]
    star = t[["q_frac_x", "q_frac_y", "q_frac_z"]].abs().apply(lambda r: tuple(sorted(r)), axis=1)  # same star
    m = t[star == tuple(sorted((abs(qx), abs(qy), abs(qz))))]
    return float(m.omega_cm1.iloc[0]), float(m.gamma_hwhm_cm1.iloc[0])


rows = []
for jid, j in man.items():
    mm = re.fullmatch(rf"m4ph_audit_q12_T{T}_m(\d+)_e(\d+)(?:_[a-z])?", jid)  # _a/_b: one mesh split over jobs
    if not mm or j["status"] != "done":
        continue
    mesh, eps = int(mm.group(1)), int(mm.group(2))
    audit = j["metrics"].get("audit_4ph_c", [])
    order = [",".join(f"{x:g}" for x in md) for md in j["params"]["modes"]]  # AUDIT lines follow the job's mode order
    for key, v in j["metrics"]["modes"].items():
        i = order.index(key)
        c = PAT.search(audit[i]) if i < len(audit) else None
        om, g3p = prod_mode(key)
        g3, g4 = v["Gamma3_hwhm"][0], v["Gamma4_hwhm"][0]
        rows.append({"mode": MODES[key], "mesh": mesh, "eps": eps, "omega_cm1": v["omega_cm1"], "Gamma3": g3,
                     "Gamma4": g4, "G4_over_G3": g4 / g3, "Gamma3_prod16": g3p,
                     "Gamma3_calib_pct": 100 * (g3 / g3p - 1),
                     "eta3p4_probe": (g3p + g4) / om,
                     "excluded": int(c.group(2)) if c else None, "excluded_self": int(c.group(3)) if c else None,
                     "exact_zero_nonidentity": int(c.group(4)) - int(c.group(2)) if c else None,
                     "near_misses": int(c.group(5)) if c else None, "build": j.get("build", "")})
per = pd.DataFrame(rows).sort_values(["mode", "mesh", "eps"])
per.to_csv(PJ / f"a1_per_mesh{SUF}.csv", index=False, float_format="%.4g")

dec = []
for mode, thr in THR.items():
    d = per[per["mode"] == mode].set_index(["mesh", "eps"]).Gamma4
    om, g3p = prod_mode([k for k, v in MODES.items() if v == mode][0])
    have = all(x in d.index for x in RULE_JOBS)
    rec = {"mode": mode, "Gamma4_thr_prereg": thr, "Gamma4_thr_recomputed": 0.5 * om - g3p, "complete": have}
    if have:
        s1 = 0.5 <= d[(8, 5)] / d[(6, 5)] <= 2
        scan = [d[(6, e)] for e in (2, 5, 10)]
        s2 = max(scan) / min(scan) <= 2
        gmax = max(d[x] for x in RULE_JOBS)
        rec.update({"G4_6_5": d[(6, 5)], "G4_8_5": d[(8, 5)], "S1_ratio_8_6": d[(8, 5)] / d[(6, 5)],
                    "S2_scan_maxmin": max(scan) / min(scan), "S1": s1, "S2": s2, "Gamma4_max": gmax})
        if not (s1 and s2):
            rec["outcome"] = "Gamma4 not determined at feasible meshes (fallback B for this mode)"
        elif gmax <= 0.5 * thr:
            rec["outcome"] = "eta_3+4 < 0.5 with at least a factor-2 margin (diagram (c) included)"
        elif gmax < thr:
            rec["outcome"] = "eta_3+4 < 0.5, within a factor of 2 of the threshold"
        else:
            rec["outcome"] = "eta_3+4 may reach 0.5 for this mode; claim (i) qualified for it"
        if mode == "R":  # the AFD mode only (R-b13 is a hard optic mode)
            rec["outcome"] += " (in this model)"
    else:
        rec["outcome"] = "incomplete: rule jobs missing " + ", ".join(f"{m}^3 e{e}" for m, e in RULE_JOBS
                                                                      if (m, e) not in d.index)
    dec.append(rec)
dec = pd.DataFrame(dec)
dec.to_csv(PJ / f"a1_decision{SUF}.csv", index=False, float_format="%.4g")
pd.set_option("display.width", 250, "display.max_colwidth", 90)
print(per[["mode", "mesh", "eps", "Gamma3", "Gamma3_calib_pct", "Gamma4", "G4_over_G3", "excluded",
           "exact_zero_nonidentity", "near_misses"]].to_string(index=False))
print(dec.to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"a1_per_mesh{SUF}.csv", f"a1_decision{SUF}.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/provenance/manifest_final.json', 'data/processed/sto_production/mode_table_T.csv.gz'])
