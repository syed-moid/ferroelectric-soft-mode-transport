"""Additive four-phonon margin: eta4 = Gamma4/Omega over every sampled mesh/eps case.
Inputs: a1_per_mesh.csv (soft modes, 300 K; only the 6^3/8^3 meshes of the R4 rule), a1_per_mesh_T1000.csv
(soft modes + X b7 / R b13, 1000 K), the descriptive (1/4,1/2,1/2) b13 jobs m4ph_auditdesc_* (queue manifest),
n4_results.csv (heat carriers, 300/800 K); largest production eta3 per T from eta_validity_T.csv.
Outputs: eta4_margin_cases.csv (one row per mesh/eps case), eta4_margin.csv (per-T margin table)."""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PJ = ROOT / "data/processed/revision"
rd = lambda f: pd.read_csv(f, comment="#")
cases = []
for f, T in (("a1_per_mesh.csv", 300), ("a1_per_mesh_T1000.csv", 1000)):
    a = rd(PJ / f)
    a = a[a.mesh.isin([6, 8])]  # 4^3/5^3 were gate meshes only
    for r in a.itertuples():
        grp = "soft" if r.mode in ("Gamma-TO1", "M", "R") else "hard optic (1000 K check)"
        cases.append(dict(group=grp, mode=r.mode, T_K=T, mesh=r.mesh, eps=r.eps, omega_cm1=r.omega_cm1,
                          Gamma3=r.Gamma3, Gamma4=r.Gamma4))
man = json.loads((ROOT / "data/provenance/manifest_final.json").read_text())["jobs"]
for e in (2, 5, 10):
    (k, v), = man[f"m4ph_auditdesc_q12_T1000_m8_e{e}"]["metrics"]["modes"].items()
    cases.append(dict(group="hard optic (descriptive)", mode="(1/4,1/2,1/2) b13", T_K=1000, mesh=8, eps=e,
                      omega_cm1=v["omega_cm1"], Gamma3=v["Gamma3_hwhm"][0], Gamma4=v["Gamma4_hwhm"][0]))
for r in rd(PJ / "n4_results.csv").itertuples():
    cases.append(dict(group=f"heat carrier ({r.cls})", mode=f"{r.q} b{r.b} [{r.sample}^3 sample]", T_K=r.T,
                      mesh=r.sample, eps=r.eps, omega_cm1=r.om, Gamma3=r.G3, Gamma4=r.G4))
c = pd.DataFrame(cases)
c["eta4"] = c.Gamma4 / c.omega_cm1
c["G4_over_G3"] = c.Gamma4 / c.Gamma3
c.sort_values("eta4", ascending=False).to_csv(PJ / "eta4_margin_cases.csv", index=False, float_format="%.5g")

combos = c.drop_duplicates(["mode", "T_K"])
top = c.sort_values("eta4", ascending=False)
print(f"coverage: {len(combos)} mode-temperature combinations, {len(c)} mesh/eps cases")
print(combos.groupby(["group", "T_K"]).size().to_string())
print("\nlargest eta4:\n" + top.head(8)[["mode", "T_K", "mesh", "eps", "eta4"]].to_string(index=False))
for g, d in c.groupby(c.group.str.startswith("heat")):
    print(("heat carriers" if g else "soft/hard optic"), "max eta4 by T:", d.groupby("T_K").eta4.max().round(4).to_dict())

eta = rd(ROOT / "data/processed/sto_production/eta_validity_T.csv")
emax = c.eta4.max()
rows = []
for T in (300, 800, 950, 1000):
    x = eta[eta.T_K == T].sort_values("eta_max").iloc[-1]
    s = c[c.T_K == T].eta4
    rows.append(dict(T_K=T, eta3_max=x.eta_max, eta3_max_mesh=int(x.mesh), eta3_max_q_index=int(x.eta_max_q_index),
                     eta3_max_branch=int(x.eta_max_branch), eta4_needed=0.5 - x.eta_max,
                     needed_over_eta4_max=(0.5 - x.eta_max) / emax,
                     eta4_max_sampled_at_T=s.max() if len(s) else float("nan")))
m = pd.DataFrame(rows)
with open(PJ / "eta4_margin.csv", "w") as fh:
    fh.write(f"# eta4 needed for eta_3+4 = 0.5 against the largest sampled eta4 = {emax:.4f} "
             f"({top.iloc[0]['mode']}, {top.iloc[0].T_K} K, {top.iloc[0].mesh}^3, eps {top.iloc[0].eps}); "
             "built by scripts/eta4_margin.py\n")
    m.to_csv(fh, index=False, float_format="%.4g")
print("\n" + m.to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"eta4_margin.csv", f"eta4_margin_cases.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/provenance/manifest_final.json', 'data/processed/sto_production/eta_validity_T.csv', 'data/processed/revision/n4_results.csv'])
