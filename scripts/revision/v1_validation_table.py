"""Tabulate the four-phonon validation runs (spec_V1) from the queue manifest.
Per soft mode at 300 K: Gamma4 unmodified / strict-index / enlarged (cm^-1, HWHM), removed amounts, the omitted-sum
check (unmodified - enlarged = sum of omitted terms) and N_k * eps * (enlarged - strict difference).
Output: data/processed/revision/v1_validation.csv (+ printed table)."""
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
m = json.loads((ROOT / "data/provenance/manifest_final.json").read_text())["jobs"]
NAMES = {"0,0,0,4": "Gamma-TO1", "0.5,0.5,0.5,1": "R", "0.5,0.5,0,1": "M"}
NUM = r"([-+0-9.eE]+)"


def audit(jid):
    j = m.get(jid)
    if not j or j["status"] != "done":
        return None
    out = {}
    # the AUDIT lines are printed in the order of params["modes"] (KS_INPUT); the metrics dict may be key-sorted
    keys = [f"{a},{b},{c},{n}" for a, b, c, n in j["params"]["modes"]]
    for key, line in zip(keys, j["metrics"].get("audit_4ph_c", [])):
        d = {"Gamma4_printed": j["metrics"]["modes"][key]["Gamma4_hwhm"][0],
             "omega_cm1": j["metrics"]["modes"][key]["omega_cm1"],
             "Gamma3_same": j["metrics"]["modes"][key]["Gamma3_hwhm"][0]}
        for f in ("excluded", "strict", "strict-not-enlarged"):
            r = re.search(rf"\b{f} = (\d+)", line)
            if r:
                d[f] = int(r.group(1))
        for f in ("ImKept_cm1", "ImOmitted_cm1", "ImEnlargedNotStrict_cm1"):
            r = re.search(rf"{f} = {NUM}", line)
            if r:
                d[f] = float(r.group(1))
        out[NAMES.get(key, key)] = d
    return out


def printed(jid):
    j = m.get(jid)
    if not j or j["status"] != "done":
        return {}
    return {NAMES.get(k, k): v["Gamma4_hwhm"][0] for k, v in j["metrics"]["modes"].items()}


rows = []
for mesh, eps in ((6, 2), (6, 5), (6, 10), (8, 5)):
    nk = mesh ** 3
    strict = audit(f"m4ph_v2strict_q12_T300_m{mesh}_e{eps}")
    enl_audit = printed(f"m4ph_audit_q12_T300_m{mesh}_e{eps}")
    omit = audit(f"m4ph_v1omit_q12_T300_m{mesh}_e{eps}") if (mesh, eps) == (6, 5) else None
    unmod = printed(f"m4ph_pm_q12_T300_m{mesh}_e{eps}") if (mesh, eps) == (6, 5) else {}
    for mode in ("Gamma-TO1", "M", "R"):
        r = {"mesh": mesh, "eps_cm1": eps, "mode": mode, "Gamma4_enlarged_auditbuild": enl_audit.get(mode)}
        if strict and mode in strict:
            s = strict[mode]
            r.update({"omega_cm1": s["omega_cm1"], "Gamma3_same_settings": s["Gamma3_same"],
                      "Gamma4_unmodified": s["ImKept_cm1"] + s["ImOmitted_cm1"],
                      "Gamma4_strict": s["ImKept_cm1"],
                      "Gamma4_enlarged": s["ImKept_cm1"] - s["ImEnlargedNotStrict_cm1"],
                      "removed_strict": s["ImOmitted_cm1"],
                      "removed_enlarged": s["ImOmitted_cm1"] + s["ImEnlargedNotStrict_cm1"],
                      "enlarged_minus_strict": s["ImEnlargedNotStrict_cm1"],
                      "Nk_eps_x_diff": nk * eps * s["ImEnlargedNotStrict_cm1"],
                      "n_strict": s.get("strict"), "n_enlarged": s.get("excluded"),
                      "n_strict_not_enlarged": s.get("strict-not-enlarged")})
        if omit and mode in omit:
            o = omit[mode]
            r.update({"v1_kept": o["ImKept_cm1"], "v1_omitted": o["ImOmitted_cm1"],
                      "Gamma4_unmodified_production_binary": unmod.get(mode)})
            if unmod.get(mode) is not None and enl_audit.get(mode) is not None:
                r["unmod_minus_enlarged_printed"] = unmod[mode] - enl_audit[mode]
                r["omitted_check_rel"] = (unmod[mode] - enl_audit[mode] - o["ImOmitted_cm1"]) / o["ImOmitted_cm1"]
        rows.append(r)
df = pd.DataFrame(rows)
out = ROOT / "data/processed/revision/v1_validation.csv"
with open(out, "w") as f:
    f.write("# validation runs (spec_V1); widths are HWHM in cm^-1 (Gamma = Im Sigma_code); "
            "scripts/v1_validation_table.py\n")
    df.to_csv(f, index=False)
pd.set_option("display.width", 250, "display.max_columns", 40)
print(df.to_string(index=False))

# ------------------------------------------------------------------ removed amount against mesh and epsilon
# Removed amounts (cm^-1, HWHM): unmodified - enlarged and enlarged - strict, per mode, mesh, epsilon.
# 6^3 and 8^3: from the strict-index runs (all three Gamma4 from one loop). 4^3: no strict run; unmodified -
# enlarged from the unmodified production runs and the audit-build runs at the same settings.
rem = []
for mesh, eps in ((4, 2), (4, 5), (4, 10), (6, 2), (6, 5), (6, 10), (8, 5)):
    s = audit(f"m4ph_v2strict_q12_T300_m{mesh}_e{eps}") if mesh > 4 else None
    u, a = printed(f"m4ph_pm_q12_T300_m{mesh}_e{eps}"), printed(f"m4ph_audit_q12_T300_m{mesh}_e{eps}")
    for mode in ("Gamma-TO1", "M", "R"):
        r = {"mesh": mesh, "Nk": mesh ** 3, "eps_cm1": eps, "mode": mode}
        if s and mode in s:
            d = s[mode]
            r.update({"unmod_minus_enlarged": d["ImOmitted_cm1"] + d["ImEnlargedNotStrict_cm1"],
                      "enlarged_minus_strict": d["ImEnlargedNotStrict_cm1"],
                      "unmod_minus_strict": d["ImOmitted_cm1"],
                      "source": f"m4ph_v2strict_q12_T300_m{mesh}_e{eps}"})
        elif mode in u and mode in a:
            r.update({"unmod_minus_enlarged": u[mode] - a[mode],
                      "source": f"m4ph_pm_q12_T300_m{mesh}_e{eps} - m4ph_audit_q12_T300_m{mesh}_e{eps}"})
        else:
            continue
        rem.append(r)
rd = pd.DataFrame(rem)
for col in ("unmod_minus_enlarged", "enlarged_minus_strict"):
    if col in rd:
        rd[f"eps_x_{col}"] = rd["eps_cm1"] * rd[col]
out2 = ROOT / "data/processed/revision/v1_removed_vs_mesh.csv"
with open(out2, "w") as f:
    f.write("# removed diagram-(c) amount (cm^-1, HWHM, 300 K) against mesh and epsilon; measured values, no fitted "
            "law; scripts/v1_validation_table.py\n")
    rd.to_csv(f, index=False)
print(rd.to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"v1_validation.csv", f"v1_removed_vs_mesh.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/provenance/manifest_final.json'])
