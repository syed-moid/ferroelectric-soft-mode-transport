"""Readers for the linewidth-partition output of the patched anphon build (tools/alamode_patches/t3_partition.diff)
and for the per-branch linewidths of an anphon .result file.

Partition file format (one line per record, whitespace separated):
  T  iq  ib  _  T_K  _  total_Ry          total three-phonon HWHM of carrier (iq, ib) at T_K, base-code arithmetic
  B  iq  ib  T_K  label1  label2  decay_Ry  absorption_Ry
                                          contribution of partner-label pair (label1, label2)
"""
from __future__ import annotations

import math

import numpy as np

PI = 4.0 * math.atan(1.0)
H_PLANCK, RYD, C_LIGHT = 6.62606896e-34, 4.35974394e-18 / 2.0, 299792458   # ALAMODE's constants
TIME_RY = H_PLANCK / (2.0 * PI) / RYD
HZ_TO_KAYSER = 1.0e-2 / (2.0 * PI * C_LIGHT)
RY_TO_CM = HZ_TO_KAYSER / TIME_RY


def read_t3(path):
    """Return (tot, bins): tot[(iq, ib, T)] = total (Ry); bins[(iq, ib, T)] = [(l1, l2, decay, absorption), ...]."""
    tot, bins = {}, {}
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if p[0] == "T":
                key = (int(p[1]), int(p[2]), float(p[4]))
                assert key not in tot, f"carrier written twice: {key}"
                tot[key] = float(p[6])
            elif p[0] == "B":
                key = (int(p[1]), int(p[2]), float(p[3]))
                bins.setdefault(key, []).append((int(p[4]), int(p[5]), float(p[6]), float(p[7])))
    return tot, bins


def read_gamma(result):
    """Per-branch HWHM (cm^-1) from the #GAMMA_EACH blocks of an anphon .result: {(iq-1, ib-1, T): value}."""
    g, lines = {}, open(result).read().splitlines()
    temps = None
    for i, l in enumerate(lines):
        if l.startswith("#TEMPERATURE"):
            tmin, tmax, dt = map(float, lines[i + 1].split())
            temps = np.arange(tmin, tmax + 0.5 * dt, dt)
        if l == "#GAMMA_EACH":
            iq, ib = map(int, lines[i + 1].split())
            m = int(lines[i + 2])
            vals = [float(lines[i + 3 + m + t]) for t in range(len(temps))]
            for t, v in zip(temps, vals):
                g[(iq - 1, ib - 1, float(t))] = v
    return g
