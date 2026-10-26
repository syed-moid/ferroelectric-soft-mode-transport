#!/usr/bin/env python3
"""Experimental kappa(T) of SrTiO3 on hand, and our gap against it.

Outputs (data/processed/revision/):
  exp_datasets.csv   every experimental (and the Tadano calculated) kappa point used here, with how it was obtained
  exp_gap_vs_T.csv   our 16^3 kappa_P and kappa_total against the experimental sets (linear interpolation in T, no
                     extrapolation), next to the R-mode deficit (ours q1 = 12^3 PM vs Fauque INS, read from figure)

Figure readings ("read from figure") are made by this script from a render of Tadano & Tsuneyuki arXiv:1506.01781v1
(Literature_Data/Tadano_STO.pdf), page 9, Fig. 7:
  (any PDF renderer: page 9 of the local copy of Tadano_STO.pdf at scale 8, passed with --render)
  .venv/bin/python scripts/revision/exp_datasets.py --render <out>/p9_x8.png
Without --render the readings recorded below (from that render, 15-Oct) are used.
Calibration: frame and tick marks found by pixel scan (x: 200 K at px 139, 1000 K at px 1702; y: 4 W/mK at px 1213,
ticks every 2 W/mK = 220 px). Marker centroids after a morphological opening of the series colour.
Uncertainty: marker half-height (Muta diamonds 14.5 px = 0.13 W/mK; Popuri circles 11 px = 0.10 W/mK; triangles 12 px
= 0.11 W/mK); axis calibration +-1 px = 0.01 W/mK. Muta markers partly hidden behind Popuri markers (visible area
< 300 px) get +-0.25 W/mK. Popuri error bars are read as their visible extent.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PJ = REPO / "data/processed/revision"
CROP = (330 * 8, 55 * 8, 550 * 8, 225 * 8)       # page-9 crop at scale 8 (points * 8)
X0, X1, Y4, PX_PER_W = 139, 1702, 1213, 110.0


def to_T(x):
    return 200 + (x - X0) * 800 / (X1 - X0)


def to_k(y):
    return 4 + (Y4 - y) / PX_PER_W


def digitize(png):
    from PIL import Image
    from scipy import ndimage as nd
    a = np.asarray(Image.open(png).convert("RGB").crop(CROP)).astype(int)
    R, G, B = a[..., 0], a[..., 1], a[..., 2]
    out = []

    def blobs(mask, open_px, min_area, series, unc, legend_T=245):
        mask = mask.copy()
        mask[:700, 850:] = False                      # inset
        lab, n = nd.label(nd.binary_opening(mask, np.ones((open_px, open_px))))
        for i in range(1, n + 1):
            ys, xs = np.nonzero(lab == i)
            if len(ys) < min_area or to_T(xs.mean()) < legend_T:   # legend markers sit left of 245 K
                continue
            out.append(dict(series=series, T_K=round(to_T(xs.mean()), 1), kappa=round(to_k(ys.mean()), 3),
                            area_px=len(ys), unc=unc(len(ys)), ys=(ys.min(), ys.max()), x=int(round(xs.mean()))))

    blobs((R > 100) & (R < 170) & (G < 40) & (B > 100) & (B < 170), 5, 100, "Muta", lambda A: 0.13 if A >= 300 else 0.25)
    blobs((R < 60) & (G > 100) & (G < 170) & (B > 100) & (B < 170), 15, 150, "Popuri", lambda A: 0.10)
    blobs(a.sum(2) < 120, 9, 200, "Tadano_calc", lambda A: 0.11)
    teal = (R < 60) & (G > 100) & (G < 170) & (B > 100) & (B < 170)
    for r in out:
        if r["series"] == "Popuri":
            ys = np.nonzero(teal[:, r["x"] - 1:r["x"] + 2].any(1))[0]
            r["bar_lo"], r["bar_hi"] = round(to_k(ys.max()), 2), round(to_k(ys.min()), 2)
    return out


# readings from the 15-Oct render (identical to digitize(p9_x8.png))
READ = {
    "Muta": [(296.8, 11.136, 362), (324.3, 10.178, 355), (377.4, 9.100, 235), (466.1, 7.653, 138),
             (570.1, 6.778, 351), (678.1, 6.183, 335), (775.7, 5.666, 218)],
    "Popuri": [(321.3, 9.813, 8.80, 10.79), (372.5, 8.923, 8.06, 9.85), (421.9, 8.191, 7.36, 8.97),
               (473.8, 7.741, 6.95, 8.46), (524.3, 7.309, 6.58, 8.04), (573.9, 7.027, 6.28, 7.65),
               (622.8, 6.646, 5.96, 7.27), (673.2, 6.400, 5.75, 6.98), (722.6, 6.200, 5.54, 6.81),
               (772.7, 5.851, 5.24, 6.41), (822.4, 5.591, 5.02, 6.15), (873.3, 5.573, 5.04, 6.16)],
    "Tadano_calc": [(300.3, 8.706), (400.2, 7.582), (500.0, 6.846), (600.2, 6.256), (700.1, 5.786),
                    (800.1, 5.390), (900.1, 5.074), (997.5, 4.794)],
}
SRC = {
    "Muta": ("H. Muta, K. Kurosaki, S. Yamanaka, J. Alloys Compd. 392, 306 (2005)",
             "single crystal (title as cited by Wang 2021: 'reduced and La-doped single-crystalline SrTiO3'); "
             "which sample Tadano plotted is not stated in Tadano", "read from figure: Tadano 2015 arXiv v1 Fig. 7"),
    "Popuri": ("S. R. Popuri et al., RSC Adv. 4, 33720 (2014)",
               "not stated in the files on hand (title: 'Glass-like thermal conductivity in SrTiO3 thermoelectrics "
               "induced by A-site vacancies')", "read from figure: Tadano 2015 arXiv v1 Fig. 7 (with error bars)"),
    "Tadano_calc": ("Tadano & Tsuneyuki, PRB 92, 054301 (2015), arXiv:1506.01781v1",
                    "calculation (SCPH q1 = 8^3, BTE-RTA 12^3 q grid, per p.9 text)",
                    "read from figure: Tadano 2015 arXiv v1 Fig. 7"),
}


def build_datasets():
    rows = []
    for T, k, area in READ["Muta"]:
        rows.append(dict(dataset="Muta2005", T_K=T, kappa_W_mK=k, unc_W_mK=0.13 if area >= 300 else 0.25,
                         bar_lo=np.nan, bar_hi=np.nan, note="" if area >= 300 else "partly hidden behind a Popuri marker"))
    for T, k, lo, hi in READ["Popuri"]:
        rows.append(dict(dataset="Popuri2014", T_K=T, kappa_W_mK=k, unc_W_mK=0.10, bar_lo=lo, bar_hi=hi,
                         note="bar_lo/bar_hi = visible error-bar extent"))
    for T, k in READ["Tadano_calc"]:
        rows.append(dict(dataset="Tadano2015_calc", T_K=T, kappa_W_mK=k, unc_W_mK=0.11, bar_lo=np.nan,
                         bar_hi=np.nan, note="calculation, not experiment"))
    tb = pd.read_csv(REPO / "data/literature/STO/tachibana_2008/kappa_total_T.csv")
    for _, r in tb.iterrows():
        rows.append(dict(dataset="Tachibana2008", T_K=r.T_K, kappa_W_mK=r.kappa_W_per_m_K,
                         unc_W_mK=r.kappa_uncertainty, bar_lo=np.nan, bar_hi=np.nan,
                         note="tracked CSV, method 'exp_single_crystal_approx'; PDF not on hand; origin of value "
                              "and uncertainty not recorded"))
    d = pd.DataFrame(rows)
    src = {"Muta2005": SRC["Muta"], "Popuri2014": SRC["Popuri"], "Tadano2015_calc": SRC["Tadano_calc"],
           "Tachibana2008": ("M. Tachibana, T. Kolodiazhnyi, E. Takayama-Muromachi, APL 93, 092902 (2008)",
                             "single crystal (as recorded in the tracked CSV)", "tracked CSV (data inventory)")}
    d["reference"] = d.dataset.map(lambda s: src[s][0])
    d["sample"] = d.dataset.map(lambda s: src[s][1])
    d["obtained"] = d.dataset.map(lambda s: src[s][2])
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--render", help="page-9 render of Tadano_STO.pdf at scale 8: re-digitize and compare")
    args = ap.parse_args()
    if args.render:
        for r in digitize(args.render):
            print(r)
    d = build_datasets()
    d.to_csv(PJ / "exp_datasets.csv", index=False)

    km = pd.read_csv(REPO / "data/processed/sto_production/kappa_mesh_T.csv", comment="#")
    m16 = km[km.mesh == 16].set_index("T_K").sort_index()
    rp = pd.read_csv(PJ / "r_point_check.csv", comment="#")
    rp = rp[(rp.PM == "on") & (rp.q1 == 12)].set_index("T_K")
    # Fauque INS R-mode, read from figure (read from figure; 1 meV = 8.0655 cm^-1, +-5.0 cm^-1)
    r_exp = {300: 6.61 * 8.0655, 600: 9.85 * 8.0655, 1000: 11.6 * 8.0655}
    out = []
    for T in m16.index:
        if T < 300:
            continue
        row = dict(T_K=T, kP_16=m16.kappa_P[T], ktot_16=m16.kappa_total[T])
        for name in ("Muta2005", "Popuri2014", "Tadano2015_calc"):
            s = d[d.dataset == name].sort_values("T_K")
            inside = s.T_K.min() - 1 <= T <= s.T_K.max() + 1   # within 1 K of the read range (np.interp clamps)
            ke = np.interp(T, s.T_K, s.kappa_W_mK) if inside else np.nan
            row[f"k_{name}"] = round(ke, 3) if inside else np.nan
            if name != "Tadano2015_calc":
                row[f"kP/k_{name}"] = round(m16.kappa_P[T] / ke, 3) if inside else np.nan
                row[f"ktot/k_{name}"] = round(m16.kappa_total[T] / ke, 3) if inside else np.nan
        row["R_ours_cm"] = rp.R.get(T, np.nan)
        row["R_exp_cm_read"] = round(r_exp[T], 1) if T in r_exp else np.nan
        row["R_deficit_pct"] = round(100 * (rp.R[T] / r_exp[T] - 1), 1) if T in r_exp else np.nan
        out.append(row)
    g = pd.DataFrame(out)
    hdr = ("# ours: kappa_mesh_T.csv, 16^3 (own SCPH q1 = 12^3, PM). exp: Muta / Popuri read from Tadano arXiv v1 "
           "Fig. 7, linear interpolation in T between read points, no extrapolation (blank outside the read range).\n"
           "# R: ours q1 = 12^3 PM (r_point_check.csv) vs Fauque Fig. 1(d) INS, read from figure (+-5.0 cm^-1). "
           "Correlation table only.\n")
    (PJ / "exp_gap_vs_T.csv").write_text(hdr + g.to_csv(index=False))
    pd.set_option("display.width", 250)
    print(g.to_string())
    for name in ("Muta2005", "Popuri2014", "Tadano2015_calc"):
        s = d[(d.dataset == name) & d.T_K.between(295, 805)]
        a = -np.polyfit(np.log(s.T_K), np.log(s.kappa_W_mK), 1)[0]
        print(f"{name}: kappa ~ T^-{a:.2f} over {s.T_K.min():.0f}-{s.T_K.max():.0f} K (read points)")
    s = m16.loc[300:800]
    print(f"ours 16^3 kappa_P: T^-{-np.polyfit(np.log(s.index), np.log(s.kappa_P), 1)[0]:.2f} over 300-800 K")


if __name__ == "__main__":
    main()


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"exp_datasets.csv", f"exp_gap_vs_T.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/literature/STO/', 'data/processed/sto_production/kappa_mesh_T.csv', 'data/processed/revision/r_point_check.csv'])
