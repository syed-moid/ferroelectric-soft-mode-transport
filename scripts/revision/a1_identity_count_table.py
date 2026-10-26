"""Tabulate the independent count of the diagram-(c) excluded set (a1_identity_count.py) for every mesh and mode of
the four-phonon rule cases, so the SI count table is read from a CSV rather than from the audit notes.
Runs a1_identity_count.py (unchanged) on the eigenvalue files of the production archive and parses its output.
Output: data/processed/revision/a1_identity_counts.csv
"""
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RAW = ROOT / "data/raw/alamode_sto/revision/pm_q12"
SOFT = [("Gamma-TO1", "0 0 0 4"), ("R", "0.5 0.5 0.5 1"), ("M", "0.5 0.5 0 1")]
HARD = [("X b7", "0 0 0.5 7"), ("R b13", "0.5 0.5 0.5 13")]
CASES = [(300, 4, SOFT), (300, 6, SOFT), (300, 8, SOFT), (1000, 6, SOFT + HARD), (1000, 8, SOFT + HARD)]
LINE = re.compile(r"branch (\d+), w = \S+ Ry \(([\d.]+) cm-1\), mesh (\d+)\^3: excluded (\d+), exact-zero denominators "
                  r"(\d+), near-misses (\d+): (\{.*\})")

rows = []
for T, m, modes in CASES:
    ev = RAW / f"evec_pm_q12_T{T}_m{m}" / f"evec_pm_q12_T{T}_m{m}.mesh.evec"
    args = [a for _, t in modes for a in t.split()]
    out = subprocess.run([sys.executable, str(HERE / "a1_identity_count.py"), str(ev), *args],
                         capture_output=True, text=True, check=True).stdout
    hits = LINE.findall(out)
    assert len(hits) == len(modes), out
    for (name, _), h in zip(modes, hits):
        cls = eval(h[6])  # dict literal printed by a1_identity_count.py
        rows.append({"T_K": T, "mesh": m, "mode": name, "branch": int(h[0]), "omega_cm1": float(h[1]),
                     "excluded_py": int(h[3]), "exact_zero_py": int(h[4]), "near_misses_py": int(h[5]),
                     "near_identity_like_py": next((v for k, v in cls.items() if k.startswith("leg at q")), 0)})
        print(T, m, name, h[3:6], flush=True)
out = ROOT / "data/processed/revision/a1_identity_counts.csv"
with open(out, "w") as f:
    f.write("# independent count of the diagram-(c) excluded set (scripts/a1_identity_count_table.py -> "
            "a1_identity_count.py on evec_pm_q12_T<T>_m<n>.mesh.evec; tolerance 1e-9 Ry, near-miss window 1e-6 Ry)\n")
    pd.DataFrame(rows).to_csv(f, index=False)
print("wrote", out)


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"a1_identity_counts.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['zenodo:pm_q12_evec.tar.gz (evec_pm_q12_T*_m{4,5,6,8})'])
