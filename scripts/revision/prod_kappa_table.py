"""Production kappa table (q1 = 12^3, PM on): kappa_P, kappa_C, kappa_total per T and transport mesh, from the
queue manifest (validated metrics; solver values, closure-checked) + observed spread across the meshes run.
Output: data/processed/sto_production/kappa_mesh_T.csv (long) and kappa_spread_T.csv (wide, with spreads).
"""
import json
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[2]
m = json.loads((ROOT / "data/provenance/manifest_final.json").read_text())["jobs"]
rows = [{"T_K": j["params"]["T"], "mesh": j["params"]["mesh"], "kappa_P": j["metrics"]["kappa_P"],
         "kappa_C": j["metrics"]["kappa_C"], "kappa_total": j["metrics"]["kappa_P"] + j["metrics"]["kappa_C"],
         "closure_rel": j["metrics"]["closure_rel"], "job": k, "host": j.get("host"), "build": j.get("build")}
        for k, j in m.items() if j["kind"] == "rta" and j["status"] == "done" and j["params"].get("variant") == "pm_q12"]
df = pd.DataFrame(rows).sort_values(["T_K", "mesh"])
out = ROOT / "data/processed/sto_production"
hdr = ("# source: data/provenance/manifest_final.json (validated anphon RTA metrics), raw in data/raw/alamode_sto/revision/pm_q12/\n"
       "# method: own SCPH q1=12^3 with polarization mixing, dfc2 FC2 per T, anphon RTA tetrahedron, NONANALYTIC=3, KAPPA_COHERENT=2\n"
       "# units: W m^-1 K^-1\n")
with open(out / "kappa_mesh_T.csv", "w") as fh:
    fh.write(hdr); df.to_csv(fh, index=False, float_format="%.6g")
w = df.pivot(index="T_K", columns="mesh", values=["kappa_P", "kappa_C"])
w.columns = [f"{a}_m{b}" for a, b in w.columns]
for q in ("kappa_P", "kappa_C"):
    cols = [c for c in w.columns if c.startswith(q) and not c.endswith("_m8")]
    w[f"{q}_spread_12up_pct"] = 100 * (w[cols].max(axis=1) - w[cols].min(axis=1)) / w[f"{q}_m16"]
with open(out / "kappa_spread_T.csv", "w") as fh:
    fh.write(hdr + "# spread = (max - min)/kappa(16^3) over meshes >= 12^3 run at that T (observed spread, not a bound)\n")
    w.reset_index().to_csv(fh, index=False, float_format="%.5g")
pd.set_option("display.width", 250)
print(w.round(3).to_string())
