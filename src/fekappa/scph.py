"""Self-consistent phonon (SCPH) fits and soft-mode renormalization.

The renormalized soft-mode frequency in a displacive ferroelectric follows a
Cochran-like temperature dependence on either side of the transition. SCPH
(Tadano & Tsuneyuki, Phys. Rev. B 92, 054301 (2015)) is the primary scheme for
producing the temperature-dependent force constants; this module provides the
phenomenological two-branch fit used to summarize the soft-mode data and the
full self-consistency loop used to generate it.

A scalar-SCPH iteration with the Ehsan 2021 two-anchor calibration was drafted
during a first Turn-4 attempt and aborted (the calibration disagreed with the
hybrid BTO dataset by >2800% at T_c). That snapshot is preserved in
``archive/scph_calibration/scph_calibration_snapshot.py`` and is not part of
the active package; v1 instead works from literature-extracted
benchmark curves loaded via :mod:`fekappa.literature_io` (added in v1).
"""

from __future__ import annotations

import numpy as np
from lmfit import Model
from numpy.typing import ArrayLike


def _two_branch_omega_sq(
    T: np.ndarray, T_c: float, A_below: float, A_above: float
) -> np.ndarray:
    """Two-branch Cochran model for Omega^2(T) (continuous, zero at T_c)."""
    T = np.asarray(T, dtype=float)
    below = A_below * (T_c - T)
    above = A_above * (T - T_c)
    return np.where(T < T_c, below, above)


def fit_soft_mode_two_branch(
    T: ArrayLike,
    omega_meV: ArrayLike,
    T_c_init: float = 400.0,
    A_init: float = 0.3,
) -> dict:
    """Joint two-branch Cochran-like fit of the soft mode with a shared T_c.

    Models the squared soft-mode frequency with a V-shaped form that vanishes
    at the transition temperature ``T_c``::

        Omega^2(T) = A_below * (T_c - T)    for T < T_c
        Omega^2(T) = A_above * (T - T_c)    for T > T_c

    Both branches give 0 at ``T = T_c`` by construction, so the model is
    continuous there. The three parameters ``T_c``, ``A_below`` and
    ``A_above`` are optimized jointly (lmfit, least squares) against the data
    ``Omega^2 = omega_meV**2``. Parameter bounds: ``T_c`` in (100, 700),
    ``A_below`` and ``A_above`` in (0, 10).

    Args:
        T: Temperatures (K).
        omega_meV: Soft-mode frequencies (meV) at each temperature.
        T_c_init: Initial guess / split point for the transition temperature.
        A_init: Initial guess for both branch slopes (meV^2 / K).

    Returns:
        A dict with keys ``T_c``, ``T_c_stderr``, ``A_below``,
        ``A_below_stderr``, ``A_above``, ``A_above_stderr``, ``chi_squared``,
        ``reduced_chi_squared``, ``n_data``, ``n_params`` and ``fit_report``.

    Raises:
        ValueError: If fewer than two data points fall on either side of
            ``T_c_init`` (a branch slope cannot then be constrained).
    """
    T = np.asarray(T, dtype=float)
    omega_meV = np.asarray(omega_meV, dtype=float)
    omega_sq = omega_meV**2

    n_below = int(np.sum(T < T_c_init))
    n_above = int(np.sum(T > T_c_init))
    if n_below < 2 or n_above < 2:
        raise ValueError(
            "Need at least two data points on each side of T_c_init="
            f"{T_c_init} K to fit both branches; got {n_below} below and "
            f"{n_above} above. Adjust T_c_init or provide more data."
        )

    model = Model(_two_branch_omega_sq)
    params = model.make_params(T_c=T_c_init, A_below=A_init, A_above=A_init)
    params["T_c"].set(min=100.0, max=700.0)
    params["A_below"].set(min=0.0, max=10.0)
    params["A_above"].set(min=0.0, max=10.0)

    result = model.fit(omega_sq, params, T=T)

    def _stderr(name: str) -> float:
        err = result.params[name].stderr
        return float(err) if err is not None else float("nan")

    return {
        "T_c": float(result.params["T_c"].value),
        "T_c_stderr": _stderr("T_c"),
        "A_below": float(result.params["A_below"].value),
        "A_below_stderr": _stderr("A_below"),
        "A_above": float(result.params["A_above"].value),
        "A_above_stderr": _stderr("A_above"),
        "chi_squared": float(result.chisqr),
        "reduced_chi_squared": float(result.redchi),
        "n_data": int(result.ndata),
        "n_params": int(result.nvarys),
        "fit_report": result.fit_report(),
    }


def scph_iterate(*args, **kwargs):
    """Run the SCPH self-consistency loop for the renormalized frequencies.

    Iterates the self-consistent phonon equations (Tadano-Tsuneyuki) until the
    renormalized frequencies converge, given harmonic frequencies and quartic
    force constants at a target temperature.

    Note:
        A scalar-SCPH iteration calibrated from Ehsan 2021 anchors was drafted
        and aborted in v1 (snapshot kept under
        ``archive/scph_calibration/``); v1 takes a
        literature-benchmark approach instead. A full SCPH self-consistency
        loop will be implemented in a later turn against new first-principles
        force constants.
    """
    raise NotImplementedError(
        "SCPH self-consistency loop: deferred. Turn-4 scalar-SCPH attempt "
        "aborted (see archive/scph_calibration/)."
    )
