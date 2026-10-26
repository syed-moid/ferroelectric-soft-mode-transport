"""Lowest SCPH frequencies at Gamma (TO1), X, M, R from every SCPH run on disk,
against Tadano & Tsuneyuki 2015 Table I (300 K; q1 = 2^3..12^3) and Wang, Zeng & Chen 2021 (TDEP, R point:
1.28 THz at 200 K, 2.87 THz at 800 K). Extraction as in d_q1_freq_shift.py (path G-X-M-G-R-M of .scph_bands;
Gamma = 4th branch at k = 0, X/M/R = lowest branch). Only temperatures with 'convergence achieved' are used.
Output: data/processed/revision/r_point_check.csv
"""
from pathlib import Path

from fekappa.provenance import external_dir
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RUNS = external_dir("FEKAPPA_RUNS", "the pre-queue run directory (Zenodo prequeue_runs)")
REV = ROOT / "data/raw/alamode_sto/revision"
SRC = {("off", 2): RUNS / "scph_off_k2", ("off", 4): RUNS / "scph_off_k4", ("off", 8): RUNS / "scph_off_k8",
       ("on", 2): RUNS / "scph_on_k2", ("on", 4): RUNS / "scph_on_k4", ("on", 8): RUNS / "scph_on_k8",
       ("on", 12): REV / "pm_q12/scph_pm_q12", ("on", 16): REV / "pm_q16/scph_pm_q16"}
TADANO_300 = {2: (144, 69, 103), 4: (138, 46, 89), 6: (136, 39, 86), 8: (136, 37, 85), 10: (135, 36, 85),
              12: (135, 35, 85)}  # Gamma15, R25, M3 (cm^-1), Table I p.7 of arXiv:1506.01781v1
THZ = 33.35641  # cm^-1 per THz


def load(d):
    f = next(d.glob("*.scph_bands"))
    txt = f.read_text().splitlines()
    head = [float(x) for x in txt[1].lstrip("#").split()]
    kx = {"Gamma": head[0], "X": head[1], "M": head[2], "R": head[4]}
    a = np.loadtxt(f)
    conv = {float(l.split("=")[1].split(":")[0]) for l in (d / "scph.log").read_text().splitlines()
            if "convergence achieved" in l}
    rows = []
    for T in sorted(set(a[:, 0])):
        b = a[a[:, 0] == T]
        r = {"T_K": T, "converged": T in conv}
        for lab, x in kx.items():
            row = b[np.argmin(np.abs(b[:, 1] - x))][2:]
            r[lab] = row[3] if lab == "Gamma" else row[0]
        rows.append(r)
    return pd.DataFrame(rows)


out = []
for (pm, q1), d in SRC.items():
    if not d.exists():
        print("missing", d); continue
    df = load(d)
    df.insert(0, "q1", q1); df.insert(0, "PM", pm)
    out.append(df)
df = pd.concat(out)
df.to_csv(ROOT / "data/processed/revision/r_point_check.csv", index=False, float_format="%.4g")
c = df[df.converged]
print("lowest converged T per run:", c.groupby(["PM", "q1"]).T_K.min().to_dict())
t300 = c[c.T_K == 300].set_index(["PM", "q1"])
print("\n300 K vs Tadano Table I (ours / his / %):")
for (pm, q1), r in t300.iterrows():
    if q1 in TADANO_300:
        g, R, M = TADANO_300[q1]
        print(f"PM {pm:3s} q1 {q1:2d}: Gamma {r.Gamma:7.2f}/{g} ({100*(r.Gamma/g-1):+.1f}%)  "
              f"R {r.R:6.2f}/{R} ({100*(r.R/R-1):+.1f}%)  M {r.M:6.2f}/{M} ({100*(r.M/M-1):+.1f}%)")
    else:
        print(f"PM {pm:3s} q1 {q1:2d}: Gamma {r.Gamma:7.2f}  R {r.R:6.2f}  M {r.M:6.2f}  (no Table I entry)")
print("\nR point vs Wang TDEP (1.28 THz = %.1f cm^-1 at 200 K; 2.87 THz = %.1f at 800 K):" % (1.28 * THZ, 2.87 * THZ))
for T in (200, 250, 300, 800):
    s = df[(df.T_K == T)][["PM", "q1", "R", "converged"]]
    print(T, s.to_dict("records"))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"r_point_check.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['zenodo:scph_q12.tar.gz (scph_pm_q12)', 'zenodo:scph_q16.tar.gz (scph_pm_q16)', '$FEKAPPA_RUNS/scph_{on,off}_k{2,4,8} (zenodo:prequeue_runs.tar.gz)'])
