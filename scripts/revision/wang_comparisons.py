"""Two comparisons with Wang, Zeng & Chen PRB 104, 235205 (2021).
 (a) kappa_P(T) ~ T^-alpha over 300-800 K, least squares in log-log, on every production mesh with >= 3 points;
 (b) share of kappa_P carried by the three lowest branches (frequency ordering, as ALAMODE's branch index) at 300/800 K.
Sources: data/provenance/manifest_final.json (validated solver kappa_P), data/processed/sto_production/mode_table_T.csv.gz.
Output: data/processed/revision/wang_comparisons.csv
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
man = json.loads((ROOT / "data/provenance/manifest_final.json").read_text())["jobs"]
rows = []
for k, j in man.items():
    if k.startswith("rta_pm_q12_T") and j["status"] == "done":
        T, m = k[len("rta_pm_q12_T"):].split("_m")
        rows.append({"T": int(T), "mesh": int(m), "kP": j["metrics"]["kappa_P"]})
k = pd.DataFrame(rows)
out = []
for m, g in k[(k["T"] >= 300) & (k["T"] <= 800)].groupby("mesh"):
    if len(g) < 3:
        continue
    slope, _ = np.polyfit(np.log(g["T"]), np.log(g.kP), 1)
    out.append({"quantity": "alpha_300_800K", "mesh": m, "n_T": len(g), "value": -slope,
                "T_used": " ".join(map(str, sorted(g["T"])))})
modes = pd.read_csv(ROOT / "data/processed/sto_production/mode_table_T.csv.gz")
for (T, m), g in modes[modes.T_K.isin([300, 800])].groupby(["T_K", "mesh"]):
    out.append({"quantity": f"share_lowest3_{T}K", "mesh": m, "n_T": 1,
                "value": g[g.branch <= 3].kappa_mode_W_mK.sum() / g.kappa_mode_W_mK.sum(), "T_used": str(T)})
df = pd.DataFrame(out)
df.to_csv(ROOT / "data/processed/revision/wang_comparisons.csv", index=False, float_format="%.4g")
print(df.round(4).to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"wang_comparisons.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/provenance/manifest_final.json', 'data/processed/sto_production/mode_table_T.csv.gz'])
