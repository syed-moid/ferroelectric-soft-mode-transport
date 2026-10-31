"""SCPH control: our own SCPH run with the ALAMODE example settings against the distributed renormalization file.

Compared files (both written by anphon MODE = SCPH, `<PREFIX>.scph_dfc2`):
  distributed  $FEKAPPA_TOOLCHAIN/sto_example/STO_scph2-2.scph_dfc2 (ALAMODE example/SrTiO3/reference at c5b015c,
               same SHA-256)
  own run      $FEKAPPA_RUNS/scph_off_k2/scph_off_k2.scph_dfc2 (the example scph.in, PREFIX changed; c5b015c build)
Each entry is the SCPH correction to one harmonic force constant, Delta Phi(R; i a, j b), in Ry bohr^-2 (anphon
scph.cpp writes delta_dymat * sqrt(m_i m_j) with masses in Rydberg units; the dfc2 tool adds it to the FC2 of
STO222.xml, which is in Ry bohr^-2). The renormalized force constants are FC2 + Delta Phi with the same FC2, so their
difference equals the difference of these entries.
Per temperature block: number of entries, max |Delta Phi| (distributed), max absolute difference, RMS difference,
and max relative difference over entries with |Delta Phi| >= 1e-6 Ry bohr^-2. Only stored files are read.
Output: data/processed/revision/dfc2_reproduction.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

from fekappa.provenance import external_dir

ROOT = Path(__file__).resolve().parents[2]
DIST = external_dir("FEKAPPA_TOOLCHAIN", "the ALAMODE toolchain directory (sto_example/ holds the force-constant inputs)") \
    / "sto_example/STO_scph2-2.scph_dfc2"
OWN = external_dir("FEKAPPA_RUNS", "the pre-queue run directory (Zenodo prequeue_runs)") / "scph_off_k2/scph_off_k2.scph_dfc2"
REL_FLOOR = 1e-6


def read_dfc2(path):
    blocks, cur = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith("# Temp ="):
            cur = float(line.split("=")[1])
            blocks[cur] = []
        elif cur is not None and line.strip():
            f = line.split()
            blocks[cur].append((tuple(int(x) for x in f[:7]), float(f[7])))
    return blocks


dist, own = read_dfc2(DIST), read_dfc2(OWN)
assert dist.keys() == own.keys(), "temperature blocks differ"
rows = []
for t in dist:
    da, db = dict(dist[t]), dict(own[t])
    assert len(da) == len(dist[t]) and da.keys() == db.keys(), f"entry indices differ at {t} K"
    a = np.array([da[k] for k in da]); b = np.array([db[k] for k in da])  # matched by index (the order may differ)
    d = np.abs(b - a)
    big = np.abs(a) >= REL_FLOOR
    rows.append({"T_K": t, "n_entries": len(a), "max_abs_dPhi_Ry_bohr2": np.abs(a).max(),
                 "max_abs_diff_Ry_bohr2": d.max(), "rms_diff_Ry_bohr2": np.sqrt((d ** 2).mean()),
                 "max_rel_diff_above_1e-6": (d[big] / np.abs(a[big])).max()})
df = pd.DataFrame(rows)
out = ROOT / "data/processed/revision/dfc2_reproduction.csv"
df.to_csv(out, index=False, float_format="%.4g")
pd.set_option("display.width", 200)
print(df.to_string(index=False))
hi = df[df.T_K > 0]
print(f"\n50-1000 K: max abs diff {hi.max_abs_diff_Ry_bohr2.max():.3g} Ry bohr^-2, "
      f"max rel diff {hi['max_rel_diff_above_1e-6'].max():.3g}")


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
stamp(out, __file__, ["$FEKAPPA_TOOLCHAIN/sto_example/STO_scph2-2.scph_dfc2 (ALAMODE example, c5b015c)",
                      "$FEKAPPA_RUNS/scph_off_k2/scph_off_k2.scph_dfc2 (zenodo:prequeue_runs.tar.gz)"])
