"""Per-mode table of the N4 heat-carrier samples (docs/recorded_rules.md (R6)): Γ₃ and Γ₄ (HWHM, cm⁻¹) at the same
settings, from the job metrics; r = Γ₄/Γ₃; the lifetime factor 1/(1 + r); the production Γ₃ of the same mode
(n4_selection.csv, 16³ or 12³ table). Class labels are those of the selection (character rule at 300 K).
Inputs: data/provenance/manifest_final.json (jobs m4ph_n4_*), data/processed/revision/n4_selection.csv.
Output: data/processed/revision/n4_results.csv
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fekappa.provenance import stamp

ROOT = Path(__file__).resolve().parents[2]
PJ = ROOT / "data/processed/revision"
man = json.loads((ROOT / "data/provenance/manifest_final.json").read_text())["jobs"]
sel = pd.read_csv(PJ / "n4_selection.csv", comment="#")
sel["mesh"] = sel["sample"].str.split("^").str[0].astype(int)

rows = []
for jid, j in man.items():  # job order of the manifest
    if not jid.startswith("m4ph_n4_"):
        continue
    p = j["params"]
    mesh, T, eps = p["mesh"], p.get("T_fc2", p.get("tmin")), p["epsilon"]
    block = []
    for key, v in j["metrics"]["modes"].items():
        *q, b = key.split(",")
        q, b = np.array([float(x) for x in q]), int(b)
        s = sel[(sel.mesh == mesh) & (sel.branch == b)]
        s = s[np.all(np.abs(s[["q_frac_x", "q_frac_y", "q_frac_z"]].to_numpy() - q) < 1e-5, axis=1)].iloc[0]
        g3, g4 = v["Gamma3_hwhm"][0], v["Gamma4_hwhm"][0]
        block.append({"sample": mesh, "T": T, "eps": eps, "cls": s["class"],
                      "q": "(" + ",".join(f"{x:g}" for x in (s.q_frac_x, s.q_frac_y, s.q_frac_z)) + ")", "b": b,
                      "om": v["omega_cm1"], "G3": g3, "G4": g4, "r": g4 / g3, "life": 1 / (1 + g4 / g3),
                      "G3prod": s.gamma3_prod_300 if T == 300 else s.gamma3_prod_800})
    rows += block
df = pd.DataFrame(rows).sort_values(["sample", "T", "eps", "cls", "q"], kind="stable")
out = PJ / "n4_results.csv"
df.to_csv(out, index=False, float_format="%.5g")
stamp(out, __file__, ["data/provenance/manifest_final.json (m4ph_n4_*)", "data/processed/revision/n4_selection.csv"])
print(df.to_string(index=False))
