"""D3: isolate the SCPH q1-sampling effect and the polarization-mixing (PM) effect (300 K).

Input: data/processed/revision/c1_benchmark_matrix.csv (12 runs; transport 8^3 and 12^3).
(a) q1 = 2^3 -> 4^3 -> 8^3 at fixed PM (on, off): sampling effect isolated.
(b) PM off -> on at fixed q1 (2^3, 4^3, 8^3): mixing effect isolated.
Frequencies in meV at Gamma (TO1), M (lowest), R (lowest) from c1_scph_softmodes.csv.
Outputs: d3_table_q1_sampling.csv, d3_table_pm_mixing.csv
"""
from pathlib import Path
import pandas as pd
D = Path(__file__).resolve().parents[2] / "data/processed/revision"
m = pd.read_csv(D / "c1_benchmark_matrix.csv", comment="#")
w = m.pivot_table(index=["PM", "kmesh_scph"], columns="transport_mesh", values=["kappa_P", "kappa_C"])
w.columns = [f"{a}_{b}3" for a, b in w.columns]
f = m.drop_duplicates(["PM", "kmesh_scph"]).set_index(["PM", "kmesh_scph"])[["Omega_Gamma_TO1_meV", "Omega_M_meV", "Omega_R_meV"]]
t = w.join(f).reset_index()
t["kappa_P_123_vs_8.7_pct"] = 100 * (t["kappa_P_123"] - 8.7) / 8.7
a = t.sort_values(["PM", "kmesh_scph"])
for col in ("kappa_P_123", "Omega_M_meV", "Omega_R_meV", "Omega_Gamma_TO1_meV"):
    a[f"d_{col}_vs_q1_2_pct"] = a.groupby("PM")[col].transform(lambda s: 100 * (s - s.iloc[0]) / s.iloc[0])
a.to_csv(D / "d3_table_q1_sampling.csv", index=False, float_format="%.4g")
off, on = t[~t.PM].set_index("kmesh_scph"), t[t.PM].set_index("kmesh_scph")
b = pd.DataFrame({
    "kappa_P_83_off": off.kappa_P_83, "kappa_P_83_on": on.kappa_P_83,
    "kappa_P_123_off": off.kappa_P_123, "kappa_P_123_on": on.kappa_P_123,
    "dkappa_P_123_pct": 100 * (on.kappa_P_123 - off.kappa_P_123) / off.kappa_P_123,
    "Gamma_off_meV": off.Omega_Gamma_TO1_meV, "Gamma_on_meV": on.Omega_Gamma_TO1_meV,
    "dGamma_pct": 100 * (on.Omega_Gamma_TO1_meV - off.Omega_Gamma_TO1_meV) / off.Omega_Gamma_TO1_meV,
    "dM_pct": 100 * (on.Omega_M_meV - off.Omega_M_meV) / off.Omega_M_meV,
    "dR_pct": 100 * (on.Omega_R_meV - off.Omega_R_meV) / off.Omega_R_meV,
}).reset_index()
b.to_csv(D / "d3_table_pm_mixing.csv", index=False, float_format="%.4g")
pd.set_option("display.width", 250)
print(a.round(3).to_string(index=False)); print(); print(b.round(3).to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"d3_table_q1_sampling.csv", f"d3_table_pm_mixing.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/processed/revision/c1_benchmark_matrix.csv'])
