"""Mode-table helpers shared by the v4 analysis scripts (scripts/v4/).

Degenerate groups and group-averaged linewidths follow anphon's rule: within one irreducible wave vector, in
branch order, modes are joined to the group of the first mode while |Omega - Omega_first| < tol_omega
(1e-7 Ry); the half widths of a group are averaged before kappa is formed.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

TOL_DEG_CM1 = 1.0e-7 * 109737.31568160   # anphon tol_omega (1e-7 Ry) in cm^-1


def deg_groups(om) -> np.ndarray:
    """Group id (index of the first mode of the group) for the branches of one irreducible q, in branch order."""
    om = np.asarray(om, dtype=float)
    gid, start = np.zeros(len(om), int), 0
    for j in range(1, len(om)):
        if abs(om[j] - om[start]) >= TOL_DEG_CM1:
            start = j
        gid[j] = start
    return gid


def group_average_gamma(df):
    """Group-averaged HWHM per mode (df: one temperature, all branches of each irreducible q; columns q_index,
    branch, omega_cm1, gamma_hwhm_cm1). Returns a Series aligned with df."""
    out = df.gamma_hwhm_cm1.copy()
    for _, g in df.sort_values(["q_index", "branch"]).groupby("q_index"):
        w, idx = g.omega_cm1.to_numpy(), g.index.to_numpy()
        start = 0
        for j in range(1, len(w) + 1):
            if j == len(w) or abs(w[j] - w[start]) >= TOL_DEG_CM1:
                out.loc[idx[start:j]] = g.gamma_hwhm_cm1.iloc[start:j].mean()
                start = j
    return out


def mesh_of(rdir) -> int:
    """Transport mesh (first entry of the &kpoint mesh line) of an RTA run directory with rta.in."""
    rin = (Path(rdir) / "rta.in").read_text().splitlines()
    return int([l for l in rin if l.split()[:1] and l.split()[0].isdigit() and len(l.split()) == 3][0].split()[0])


def cell_volume_ang3(rdir) -> float:
    """Cell volume (Angstrom^3) from the &cell block of the run's rta.in (lattice constant in Bohr times the
    three cell vectors), so that each model's own lattice constant is used."""
    from fekappa.kappa_decomposition import BOHR_TO_ANG
    lines = (Path(rdir) / "rta.in").read_text().splitlines()
    i = next(k for k, l in enumerate(lines) if l.strip().lower() == "&cell")
    a = float(lines[i + 1].split()[0])
    m = np.array([[float(x) for x in lines[i + 2 + r].split()[:3]] for r in range(3)])
    return abs(np.linalg.det(a * m)) * BOHR_TO_ANG ** 3
