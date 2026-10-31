#!/usr/bin/env python3
"""300 K experimental kappa of SrTiO3 read from three published figures, and our gap against it.

All values are "read from figure". Sources (Literature_Data/):
  Tadano_STO.pdf   arXiv v1 p.9 Fig. 7     readings from exp_datasets.py (15-Oct; scale-8 render, auto centroids)
  WangSTO2021.pdf  PRB 104, 235205 p.-4 Fig. 4, scale-8 render, crop (440, 480, 2140, 1920):
                   tick rows 16/12/8/4 W/mK at y = 162/454/745.5/1039 (73.1 px per W/mK); tick columns 200/300/400 K at
                   x = 305.5/516.5/729.5 (2.116 px per K). Hollow markers overlap near 300 K, so centroids were placed by
                   hand on a 4x zoom; uncertainty = marker half-height (~15 px = 0.2 W/mK).
  Gao.pdf          (T. Guo & D. Guo, PLA 467, 128727) p.9 Fig. 10(c), scale-6 render, crop (780, 1220, 1780, 2120):
                   tick rows 20 ... 2 W/mK every 69.3 px from y = 15; tick columns 300 ... 900 K every 117 px from
                   x = 189.5. Marker centroids by colour mask (Ref. 4 purple diamonds, Ref. 5 gold triangles) after a
                   5-px opening; uncertainty = marker half-height (11 px = 0.32 W/mK). In panel (d) the 300 K Ref. 4 marker
                   is partly covered by the "(d)" label; its visible part reads 11.12, the other points agree with (c).
Guo & Guo's Ref. 4 = Muta 2005 and Ref. 5 = Popuri 2014 (their reference list, p.10).
Output: data/processed/revision/exp_300K_crossread.csv; prints the gap of our 300 K kappa against the single-crystal values.
"""
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
PJ = REPO / "data/processed/revision"

wang_y = lambda y: 16 - (y - 162) / 73.1
wang_T = lambda x: 300 + (x - 516.5) / 2.116
guo_y = lambda y: 20 - (y - 15) / 34.65
guo_T = lambda x: 300 + (x - 189.5) / 1.17

rows = [
    # dataset, figure, T, kappa, unc, note
    ("Muta2005", "Tadano 2015 Fig. 7 (arXiv v1)", 296.8, 11.136, 0.13, "exp_datasets.py"),
    ("Muta2005", "Tadano 2015 Fig. 7 (arXiv v1)", 324.3, 10.178, 0.13, "exp_datasets.py"),
    ("Popuri2014", "Tadano 2015 Fig. 7 (arXiv v1)", 321.3, 9.813, 0.10, "exp_datasets.py; lowest T plotted"),
    ("Muta2005", "Wang 2021 Fig. 4", wang_T(517.5), wang_y(517.5), 0.2, "behind the 3ph square; hand-placed centroid"),
    ("Muta2005", "Wang 2021 Fig. 4", wang_T(572.5), wang_y(588.0), 0.2, "hand-placed centroid"),
    ("Popuri2014", "Wang 2021 Fig. 4", wang_T(560.0), wang_y(613.75), 0.2, "lowest T plotted; error bar not read"),
    ("Suemune1965", "Wang 2021 Fig. 4", wang_T(502.5), wang_y(558.75), 0.2, "sample form not established in files on hand"),
    ("Buscaglia2014", "Wang 2021 Fig. 4", wang_T(496.0), wang_y(478.75), 0.2, "La-doped ceramics (Buscaglia et al. 2014)"),
    ("Yu2008", "Wang 2021 Fig. 4", wang_T(492.5), wang_y(517.5), 0.2, "oxygen-deficient (Yu et al. 2008)"),
    ("Yu2008", "Wang 2021 Fig. 4", wang_T(540.0), wang_y(544.5), 0.2, "oxygen-deficient (Yu et al. 2008)"),
    ("Muta2005", "Guo & Guo 2023 Fig. 10(c)", guo_T(189.5), guo_y(321.9), 0.32, "colour-mask centroid"),
    ("Muta2005", "Guo & Guo 2023 Fig. 10(c)", guo_T(217.1), guo_y(352.4), 0.32, "colour-mask centroid"),
    ("Popuri2014", "Guo & Guo 2023 Fig. 10(c)", guo_T(219.0), guo_y(367.8), 0.32, "lowest T plotted"),
]
d = pd.DataFrame(rows, columns=["dataset", "figure", "T_K", "kappa_W_mK", "unc_W_mK", "note"]).round(
    {"T_K": 1, "kappa_W_mK": 2})
d["obtained"] = "read from figure"
d.to_csv(PJ / "exp_300K_crossread.csv", index=False)
print(d.to_string(index=False))

# Muta at 300 K: direct readings at 300 +- 3 K, and the Tadano-figure interpolation (exp_gap_vs_T.csv, 11.025)
m300 = d[(d.dataset == "Muta2005") & d.T_K.between(295, 305)].kappa_W_mK
muta = (min(11.025, m300.min()), max(11.025, m300.max()))
km = pd.read_csv(REPO / "data/processed/sto_production/kappa_mesh_T.csv", comment="#")
k = km[(km.T_K == 300) & (km.mesh >= 15)]
kp, kt = (k.kappa_P.min(), k.kappa_P.max()), (k.kappa_total.min(), k.kappa_total.max())
print(f"\nMuta 300 K: {muta[0]:.2f}-{muta[1]:.2f} (three figures)")
print(f"ours >= 15^3 meshes: kappa_P {kp[0]:.2f}-{kp[1]:.2f} (mean {k.kappa_P.mean():.2f}), "
      f"kappa_total {kt[0]:.2f}-{kt[1]:.2f}")
for lab, (lo, hi) in (("kappa_P", kp), ("kappa_total", kt)):
    g = [100 * (v / e - 1) for v in (lo, hi) for e in muta]
    print(f"{lab} vs Muta: {min(g):.1f} to {max(g):.1f} %")


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"exp_300K_crossread.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/processed/sto_production/kappa_mesh_T.csv', 'figure readings in the script (Muta 2005, Popuri 2014)'])
