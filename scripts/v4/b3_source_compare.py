#!/usr/bin/env python3
"""Compare model B's 300 K SCPH soft modes with the source band data (figures/bands, THz) on a matched grid.

The extraction rule is the same for both sides: Gamma = 4th branch at the Gamma point of the path, M and R = the
lowest branch at their path points. The source path is G-X-M-G-R-M with points 0.5, 1.0, 1.7071, 2.5731, 3.0731.

The "disagrees materially" thresholds are |dR| > 5 cm^-1, or |dGamma| or |dM| > 5 %.

Usage: b3_source_compare.py SOFTMODES_CSV RUN=SOURCE_BAND_FILE [...]
"""
import sys

import numpy as np
import pandas as pd

THZ = 33.35641


def src(path):
    blocks, cur = [], []
    for l in open(path):
        if l.startswith("#"):
            continue
        if not l.strip():
            if cur:
                blocks.append(np.array(cur)); cur = []
            continue
        cur.append([float(x) for x in l.split()])
    if cur:
        blocks.append(np.array(cur))
    x, b = blocks[0][:, 0], np.array([bl[:, 1] for bl in blocks]) * THZ
    at = lambda xx: b[:, np.argmin(abs(x - xx))]
    return {"Gamma": at(1.7071)[3], "R": at(2.5731)[0], "M": at(1.0)[0]}


df = pd.read_csv(sys.argv[1])
for spec in sys.argv[2:]:
    run, f = spec.split("=", 1)
    r = df[(df.run == run) & (df.T_K == 300)].iloc[0]
    s = src(f)
    dG, dR, dM = 100 * (r.Gamma / s["Gamma"] - 1), r.R - s["R"], 100 * (r.M / s["M"] - 1)
    mat = abs(dR) > 5 or abs(dG) > 5 or abs(dM) > 5
    print(f"{run} vs {f.split('figures/bands/')[-1]}: Gamma {r.Gamma:.2f}/{s['Gamma']:.2f} ({dG:+.2f} %), "
          f"R {r.R:.2f}/{s['R']:.2f} ({dR:+.2f} cm-1), M {r.M:.2f}/{s['M']:.2f} ({dM:+.2f} %) -> "
          f"{'MATERIAL DISAGREEMENT' if mat else 'no material disagreement'}")
