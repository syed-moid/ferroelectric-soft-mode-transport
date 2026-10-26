"""C1: soft-mode frequencies at Gamma (TO1), R (AFD) and M (M3) from our own MODE=SCPH runs
vs Tadano & Tsuneyuki 2015 Table I (300 K, SCPH q1-grid convergence) and Fig. 5 (Omega_TO1(T)).

Runs: $FEKAPPA_RUNS/scph_{off,on}_k{2,4,8}  (off/on = SELF_OFFDIAG 0/1 = without/with
polarization mixing; k = KMESH_SCPH; KMESH_INTERPOLATE = 2 2 2 fixed by the 2x2x2 harmonic supercell).
Output: data/processed/revision/c1_scph_softmodes.csv
"""

from pathlib import Path

from fekappa.provenance import external_dir
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RUNS = external_dir("FEKAPPA_RUNS", "the pre-queue run directory (Zenodo prequeue_runs)")
CM1_PER_MEV = 8.06554
KX = {"Gamma": 0.0, "M": 0.853346, "R": 2.195772}
TABLE_I_INDEX = {"Gamma": "FE", "R": "AFD", "M": "M3"}

lit = pd.read_csv(ROOT / "data/literature/STO/tadano_2015/scph_convergence_300K.csv")

rows = []
for cfg in ["off_k2", "on_k2", "off_k4", "on_k4", "off_k8", "on_k8"]:
    f = RUNS / f"scph_{cfg}" / f"scph_{cfg}.scph_bands"
    if not f.exists():
        continue
    log = (RUNS / f"scph_{cfg}" / "scph.log").read_text()
    a = np.loadtxt(f)
    for T in sorted(set(a[:, 0])):
        conv = f"Temp = {T:.6e} : convergence achieved" in log
        b = a[a[:, 0] == T]
        for lab, kx in KX.items():
            row = b[np.argmin(np.abs(b[:, 1] - kx))][2:]
            # Gamma: soft TO1 = lowest optic (skip 3 acoustic); R, M: lowest mode
            w = row[3] if lab == "Gamma" else row[0]
            k = int(cfg[-1])
            ref = lit[(lit["branch"] == f"{TABLE_I_INDEX[lab]}_{k}x{k}x{k}")]["omega_meV"]
            ref = float(ref.iloc[0]) if (len(ref) and T == 300) else np.nan
            rows.append({"config": cfg, "PM": cfg.startswith("on"), "kmesh_scph": k, "T_K": T,
                         "scph_converged": conv, "point": lab, "omega_cm1": w,
                         "omega_meV": w / CM1_PER_MEV, "tadano_tableI_meV": ref,
                         "dev_pct": 100 * (w / CM1_PER_MEV - ref) / ref})
df = pd.DataFrame(rows)
df.to_csv(ROOT / "data/processed/revision/c1_scph_softmodes.csv", index=False, float_format="%.5g")
pd.set_option("display.width", 200)
print(df[df["T_K"] == 300].round(3).to_string(index=False))
print(df[(df["point"] == "Gamma") & df["T_K"].isin([0, 100, 200, 300, 400, 500, 1000])]
      .pivot(index="T_K", columns="config", values="omega_meV").round(2))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"c1_scph_softmodes.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['$FEKAPPA_RUNS/scph_{on,off}_k{2,4,8} (zenodo:prequeue_runs.tar.gz)', 'data/literature/STO/tadano_2015/scph_convergence_300K.csv'])
