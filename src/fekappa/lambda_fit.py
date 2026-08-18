"""Small-q dispersion stiffness Lambda^2(T) of the renormalized soft branch.

Fits ``Omega^2(q) = Omega_0^2 + Lambda^2 q^2`` to the SCPH-renormalized TO1
branch along Gamma -> X, computed by ``anphon`` in phonons mode with the
dfc2-renormalized fc2 at each temperature (``STO_band_scph_<T>K.bands``).
This replaces the effective Verdi Gamma->X curvature proxy used previously
(old caveat C2): here both Omega_0(T) and Lambda^2(T) come from the same
renormalized dynamical matrix at each T.

Conventions:

* ``.bands`` k-axis is in bohr^-1 (2*pi included); converted to Angstrom^-1
  here, once.
* Frequencies in the ``.bands`` file are in cm^-1; Lambda^2 is reported in
  meV^2 Angstrom^2 to match the kernel's units (old Eq. 4.4 machinery).
* The TO1 branch is the transverse doublet starting at the renormalized
  zone-center soft frequency; along the [100] path it is branches 4 and 5
  (1-indexed, ascending order) and stays below every crossing candidate in
  the fitted q-range. Degeneracy of the doublet is asserted, not assumed.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

# cm^-1 -> meV (1 meV = 8.06554 cm^-1).
CM1_TO_MEV: float = 1.0 / 8.06554
# bohr^-1 -> Angstrom^-1.
BOHR1_TO_ANG1: float = 1.0 / 0.529177210903

# 1-indexed columns of the TO1 transverse doublet in the [100] band file.
TO1_BRANCHES: tuple[int, int] = (4, 5)


def read_bands(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read an ``anphon`` ``.bands`` file.

    Returns:
        ``(q_ang1, omega_cm1)``: the k-axis converted to Angstrom^-1 and the
        branch frequencies (n_q, n_branch) in cm^-1, ascending per row as
        written by anphon.
    """
    data = np.loadtxt(path, comments="#")
    q_ang1 = data[:, 0] * BOHR1_TO_ANG1
    return q_ang1, data[:, 1:]


def fit_lambda_sq(
    q_ang1: np.ndarray,
    omega_cm1: np.ndarray,
    *,
    q_max_ang1: float,
    branches: tuple[int, int] = TO1_BRANCHES,
    degeneracy_rtol: float = 1e-3,
) -> dict[str, float]:
    """Fit ``Omega^2 = Omega_0^2 + Lambda^2 q^2`` to the TO1 doublet.

    Args:
        q_ang1: k-axis in Angstrom^-1.
        omega_cm1: branch frequencies in cm^-1, shape (n_q, n_branch).
        q_max_ang1: upper end of the fit window (lower end is Gamma).
        branches: 1-indexed columns of the transverse doublet.
        degeneracy_rtol: allowed relative split of the doublet in the window.

    Returns:
        Dict with ``Omega0_cm1``, ``Omega0_meV``, ``lambda_sq_meV2_ang2``,
        ``qmax_fit_ang1``, ``npoints`` and ``R2`` (coefficient of
        determination of the linear fit in the (q^2, Omega^2) plane).

    Raises:
        ValueError: If the doublet is split beyond ``degeneracy_rtol`` or the
            fit window holds fewer than 5 points.
    """
    mask = q_ang1 <= q_max_ang1
    if mask.sum() < 5:
        raise ValueError(f"only {mask.sum()} points at q <= {q_max_ang1} A^-1")
    b1, b2 = (b - 1 for b in branches)
    w1, w2 = omega_cm1[mask, b1], omega_cm1[mask, b2]
    if np.any(np.abs(w1 - w2) > degeneracy_rtol * np.abs(w1)):
        raise ValueError("TO1 doublet is split beyond tolerance in the fit window")

    q = q_ang1[mask]
    omega_meV = 0.5 * (w1 + w2) * CM1_TO_MEV
    omega_sq = omega_meV**2

    slope, intercept = np.polyfit(q**2, omega_sq, 1)
    pred = intercept + slope * q**2
    ss_res = float(np.sum((omega_sq - pred) ** 2))
    ss_tot = float(np.sum((omega_sq - omega_sq.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    omega0_meV = float(np.sqrt(intercept))
    return {
        "Omega0_cm1": omega0_meV / CM1_TO_MEV,
        "Omega0_meV": omega0_meV,
        "lambda_sq_meV2_ang2": float(slope),
        "qmax_fit_ang1": float(q_max_ang1),
        "npoints": int(mask.sum()),
        "R2": r2,
    }


_HEADER_LAMBDA = (
    "# source: STO_band_scph_<T>K.bands (anphon MODE=phonons, "
    "FCSXML = dfc2-renormalized fc2 per T, NONANALYTIC=3), "
    "Gamma->(0.15,0,0) r.l.u., 61 points\n"
    "# method: linear fit of Omega^2 vs q^2 for the TO1 transverse doublet "
    "(branches 4-5) over q <= qmax_fit_ang1; Omega_0 from the same fit\n"
    "# units: Omega0 in cm^-1 / meV; lambda_sq in meV^2 Angstrom^2; "
    "q in Angstrom^-1\n"
    "# T-grid: 0-400 K step 50. Above ~450 K the renormalized TO1 anticrosses "
    "TO2 at the zone center (188.98 vs 189.65 cm^-1 at 450 K), so a "
    "single-branch quadratic small-q fit is no longer defined; the soft-mode "
    "kernel consumes this table only in the 0-400 K regime\n"
    "# provenance: ext-FP-alamode (ALAMODE tutorial force constants -- see "
    "data/raw/alamode_sto/IMPORT_NOTE.md)\n"
)

# TO1 is a well-separated soft branch only below the ~450 K TO1-TO2
# anticrossing; the quadratic fit is emitted for this range alone.
LAMBDA_TEMPS_K: tuple[int, ...] = tuple(range(0, 401, 50))


def emit_lambda_sq_csv(
    band_dir: Path,
    out_path: Path,
    temps_K: tuple[int, ...] = LAMBDA_TEMPS_K,
    *,
    q_max_ang1: float = 0.15,
) -> pd.DataFrame:
    """Fit every temperature's band file and write ``lambda_sq_T.csv``.

    Returns the emitted DataFrame (one row per temperature).
    """
    band_dir, out_path = Path(band_dir), Path(out_path)
    rows = []
    for t in temps_K:
        q, omega = read_bands(band_dir / f"STO_band_scph_{t}K.bands")
        rows.append({"T_K": float(t), **fit_lambda_sq(q, omega, q_max_ang1=q_max_ang1)})
    df = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(_HEADER_LAMBDA)
        df.to_csv(fh, index=False)
    return df
