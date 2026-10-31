"""Combined damping ratios of SI Tables S23-S25 from unrounded inputs (rounded once, at print).

The tables were printed from s7_widths.csv, s7_margin_inequality.csv and a1_decision*.csv, whose inputs pass through
a1_per_mesh*.csv (4 significant digits) and eta4_margin_cases.csv (5), and whose ratios are stored to 4 decimals and
printed to 3 or 4; X-b7 at 1000 K printed eta_3 + eta_4,max = 0.102 where (Gamma3 + Gamma4,max)/Omega = 0.10146.
Here every ratio is recomputed with the same formulas from the full printed widths and frequencies of the four-phonon
runs (job manifest metrics, as in a1_decision.py) and the unrounded production Gamma3 of the 16^3 mode table (same
star and branch, as in a1_decision.prod_mode; the descriptive (1/4,1/2,1/2) b13 mode as in build_si_tables_v4.py):
  Table S23  eta3 = Gamma3(16^3)/Omega, eta4 = Gamma4/Omega, eta3 + eta4, eta3+4(same) = (Gamma3(same) + Gamma4)/Omega
  Table S24  (i) Gamma4(8^3, eps 5)/Gamma4(6^3, eps 5), (ii) max/min of Gamma4(6^3; eps 2, 5, 10)
  Table S25  1/2(0.5 Omega - Gamma3(16^3)), eta3 + eta4,max = (Gamma3(16^3) + max Gamma4)/Omega
Omega is the frequency printed by the four-phonon run (the column used by the tables). Only stored files are read.
Output: data/processed/revision/s7_ratios_unrounded.csv (one row per printed cell: table, row key, column, unrounded
value, the printed value before and the value rounded once at the same number of decimals).
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SI = ROOT / "data/processed/si_v4"
PJ = ROOT / "data/processed/revision"
man = json.loads((ROOT / "data/provenance/manifest_final.json").read_text())["jobs"]
mt = pd.read_csv(ROOT / "data/processed/sto_production/mode_table_T.csv.gz", comment="#")
rd = lambda f: pd.read_csv(f, comment="#")
MODES = {"0,0,0,4": "Gamma-TO1", "0.5,0.5,0.5,1": "R", "0.5,0.5,0,1": "M", "0,0,0.5,7": "X-b7", "0.5,0.5,0.5,13": "R-b13"}
DESC = "(1/4,1/2,1/2) b13"


def prod_g3(key, T):
    qx, qy, qz, b = (float(x) for x in key.split(","))
    t = mt[(mt.mesh == 16) & (mt.T_K == T) & (mt.branch == int(b))]
    star = t[["q_frac_x", "q_frac_y", "q_frac_z"]].abs().apply(lambda r: tuple(sorted(r)), axis=1)
    return float(t[star == tuple(sorted((abs(qx), abs(qy), abs(qz))))].gamma_hwhm_cm1.iloc[0])


cases = []
for jid, j in man.items():
    mm = re.fullmatch(r"m4ph_audit_q12_T(\d+)_m(\d+)_e(\d+)(?:_[a-z])?", jid)
    if mm and j["status"] == "done" and int(mm.group(2)) in (6, 8):
        T, mesh, eps = (int(x) for x in mm.groups())
        for key, v in j["metrics"]["modes"].items():
            cases.append(dict(T_K=T, mode=MODES[key], rule="decision rule", mesh=mesh, eps=eps, omega=v["omega_cm1"],
                              G3=v["Gamma3_hwhm"][0], G4=v["Gamma4_hwhm"][0], G3p=prod_g3(key, T)))
m16 = mt[(mt.T_K == 1000) & (mt.mesh == 16)]
qs = np.sort(np.abs(m16[["q_frac_x", "q_frac_y", "q_frac_z"]].to_numpy()), axis=1)
off = m16[np.isclose(qs, [0.25, 0.5, 0.5]).all(axis=1) & (m16.branch == 13)].iloc[0]
for e in (2, 5, 10):
    (k, v), = man[f"m4ph_auditdesc_q12_T1000_m8_e{e}"]["metrics"]["modes"].items()
    cases.append(dict(T_K=1000, mode=DESC, rule="descriptive", mesh=8, eps=e, omega=v["omega_cm1"],
                      G3=v["Gamma3_hwhm"][0], G4=v["Gamma4_hwhm"][0], G3p=float(off.gamma_hwhm_cm1)))
c = pd.DataFrame(cases)
c["eta3"] = c.G3p / c.omega
c["eta4"] = c.G4 / c.omega
c["eta34"] = c.eta3 + c.eta4
c["eta34_same"] = (c.G3 + c.G4) / c.omega

out = []


def cell(table, key, col, val, old, nd):
    new = f"{val:.{nd}f}"
    out.append(dict(table=table, row=key, column=col, unrounded=f"{val:.10g}", printed_before=old, rounded_once=new,
                    changed=old != new))


# Table S23: the printed rows (s7_widths.csv, as build_si_v4.t_s7_widths prints them)
w = rd(SI / "s7_widths.csv")
assert len(w) == len(c), (len(w), len(c))
for r in w.itertuples():
    x = c[(c.T_K == r.T_K) & (c["mode"] == r.mode) & (c.mesh == r.mesh) & (c.eps == r.eps)]
    assert len(x) == 1, (r.T_K, r.mode, r.mesh, r.eps)
    x = x.iloc[0]
    key = f"{r.T_K} K, {r.mode}, {r.mesh}^3, eps {r.eps}"
    for col, src, nd in (("η_{3}", "eta3_prod", 4), ("η₄", "eta4", 4), ("η_{3} + η₄", "eta34_prod", 4),
                         ("η_{3+4} same", "eta34_same", 4)):
        cell("S23", key, col, x[{"eta3_prod": "eta3", "eta4": "eta4", "eta34_prod": "eta34",
                                  "eta34_same": "eta34_same"}[src]], f"{getattr(r, src):.{nd}f}", nd)

# Table S24: decision ratios (a1_decision*.csv, printed to 2 decimals)
for f, T in (("a1_decision.csv", 300), ("a1_decision_T1000.csv", 1000)):
    for r in rd(PJ / f).itertuples():
        d = c[(c.T_K == T) & (c["mode"] == r.mode)].set_index(["mesh", "eps"]).G4
        scan = [d[(6, e)] for e in (2, 5, 10)]
        cell("S24", f"{T} K, {r.mode}", "(i) 8³/6³", d[(8, 5)] / d[(6, 5)], f"{r.S1_ratio_8_6:.2f}", 2)
        cell("S24", f"{T} K, {r.mode}", "(ii) 6³ max/min", max(scan) / min(scan), f"{r.S2_scan_maxmin:.2f}", 2)

# Table S25: margin rows (s7_margin_inequality.csv)
for r in rd(SI / "s7_margin_inequality.csv").itertuples():
    x = c[(c.T_K == r.T_K) & (c["mode"] == r.mode) & (c.rule == "decision rule")]
    om = x.omega.iloc[0]
    assert np.allclose(x.omega, om) and np.allclose(x.G3p, x.G3p.iloc[0])
    g3p, g4max = x.G3p.iloc[0], x.G4.max()
    cell("S25", f"{r.T_K} K, {r.mode}", "½(0.5Ω − Γ_{3})", 0.5 * (0.5 * om - g3p), f"{r.allowance_half:.2f}", 2)
    cell("S25", f"{r.T_K} K, {r.mode}", "η_{3} + η_{4,max}", (g3p + g4max) / om, f"{r.eta34_max:.3f}", 3)

df = pd.DataFrame(out)
dst = PJ / "s7_ratios_unrounded.csv"
df.to_csv(dst, index=False)
pd.set_option("display.width", 220, "display.max_colwidth", 60)
print(f"{len(df)} printed ratio cells; changed when rounded once from unrounded inputs: {int(df.changed.sum())}")
print(df[df.changed].to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
stamp(dst, __file__, ["data/provenance/manifest_final.json", "data/processed/sto_production/mode_table_T.csv.gz",
                      "data/processed/si_v4/s7_widths.csv", "data/processed/si_v4/s7_margin_inequality.csv",
                      "data/processed/revision/a1_decision.csv", "data/processed/revision/a1_decision_T1000.csv"])
