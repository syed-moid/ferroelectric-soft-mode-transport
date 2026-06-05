"""Self-consistent phonon (SCPH) fits and soft-mode renormalization.

The renormalized soft-mode frequency in a displacive ferroelectric follows a
Cochran-like temperature dependence on either side of the transition. SCPH
(Tadano & Tsuneyuki, Phys. Rev. B 92, 054301 (2015)) is the primary scheme for
producing the temperature-dependent force constants; this module provides:

* The phenomenological two-branch Cochran fit (:func:`fit_soft_mode_two_branch`).
* An externally calibrated **scalar SCPH** iteration (Turn 4): the soft-mode
  self-consistency loop reduced to a single mode whose quartic constant ``A``
  is fixed from two literature anchors (Ehsan *et al.* PRB 103, 094108 (2021)
  two-anchor scheme, option (b)). Section 3.3-3.7 of the manuscript.
"""

from __future__ import annotations

import numpy as np
from lmfit import Model
from numpy.typing import ArrayLike
from scipy.constants import physical_constants

# Boltzmann constant in meV/K. T (K) -> kT (meV) is just KB_MEV_K * T;
# the inverse 1 / KB_MEV_K = 11.604518 K/meV is the conversion the brief quotes.
KB_MEV_K: float = 1e3 * physical_constants["Boltzmann constant in eV/K"][0]


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


def _bose_factor(E_meV: float, T_K: float) -> float:
    """Bose-Einstein occupation ``n_B(E, T)`` with ``E`` in meV, ``T`` in K.

    When ``E`` is negative (transient during SCPH iteration), the magnitude
    ``|E|`` is used, per the brief. Returns 0 at ``T <= 0`` or for arguments
    that would overflow ``exp``.
    """
    if T_K <= 0.0:
        return 0.0
    x = abs(E_meV) / (KB_MEV_K * T_K)
    if x > 700.0:  # exp(700) overflows; n_B effectively 0
        return 0.0
    return 1.0 / (np.exp(x) - 1.0)


def calibrate_A_from_ehsan(
    E0_sq_meV2: float = -424.0,
    E_anchor_meV: float = 21.08,
    T_anchor_K: float = 500.0,
) -> float:
    """Solve the scalar-SCPH quartic constant ``A`` from a single high-T anchor.

    Equation (3.5) of the manuscript, option (b) of the Ehsan two-anchor
    scheme. Given the bare squared frequency ``E_0^2`` (at 0 K, Gamma point;
    negative for an unstable soft mode) and the renormalized soft-mode energy
    ``E_anchor`` at temperature ``T_anchor``, the scalar SCPH self-consistency

        E_anchor^2 = E_0^2 + A * (2 * n_B(E_anchor, T_anchor) + 1)

    yields

        A = (E_anchor^2 - E_0^2) / (2 * n_B(E_anchor, T_anchor) + 1).

    Defaults reproduce the Ehsan *et al.* PRB 103, 094108 (2021) BaTiO3 anchors
    (Lock 1 + Lock 2 of Turn 4): ``E_0^2 = -424`` meV^2 (166i cm^-1) and
    ``hbar Omega(500 K) = 21.08`` meV.

    Args:
        E0_sq_meV2: Bare squared frequency in meV^2 (negative for soft mode).
        E_anchor_meV: Renormalized soft-mode energy at the high-T anchor (meV).
        T_anchor_K: Anchor temperature (K).

    Returns:
        Quartic renormalization constant ``A`` in meV^2.
    """
    n_B = _bose_factor(E_anchor_meV, T_anchor_K)
    return float((E_anchor_meV**2 - E0_sq_meV2) / (2.0 * n_B + 1.0))


def scph_iterate(
    E0_sq_meV2: float,
    A_meV2: float,
    T_K: float,
    E_init_meV: float = 5.0,
    alpha: float = 0.2,
    tol_meV2: float = 1e-4,
    max_iter: int = 500,
) -> tuple[float, int, bool]:
    """Iterate the scalar SCPH self-consistency at temperature ``T_K``.

    Equation (3.1) of the manuscript:

        E^2(T) = E_0^2 + A * (2 * n_B(E, T) + 1).

    Solved with linear mixing

        E_{m+1}^2 = (1 - alpha) * E_m^2 + alpha * [E_0^2 + A * (2 n_B(E_m, T) + 1)],

    starting from ``E_init_meV``. The Bose factor uses ``sqrt(|E^2|)`` when the
    proposal becomes momentarily negative (per the Turn 4 brief). If the
    iteration converges to a negative ``E^2``, that value is returned as-is;
    the caller interprets the sign.

    Args:
        E0_sq_meV2: Bare squared frequency (meV^2).
        A_meV2: Quartic renormalization constant (meV^2).
        T_K: Temperature (K).
        E_init_meV: Initial guess for ``E`` (meV); ``E_init^2`` seeds the loop.
        alpha: Linear-mixing factor in [0.1, 0.3]; smaller is more conservative.
        tol_meV2: Convergence tolerance on ``|Delta E^2|`` (meV^2).
        max_iter: Maximum iterations.

    Returns:
        ``(E_sq_meV2, n_iter, converged)``.
    """
    if not 0.0 < alpha <= 1.0:
        raise ValueError(f"alpha must be in (0, 1]; got {alpha}")

    E_sq = float(E_init_meV) ** 2
    for k in range(1, max_iter + 1):
        E_for_bose = float(np.sqrt(abs(E_sq)))
        n_B = _bose_factor(E_for_bose, T_K)
        E_sq_rhs = E0_sq_meV2 + A_meV2 * (2.0 * n_B + 1.0)
        E_sq_new = (1.0 - alpha) * E_sq + alpha * E_sq_rhs
        if abs(E_sq_new - E_sq) < tol_meV2:
            return float(E_sq_new), k, True
        E_sq = E_sq_new
    return float(E_sq), max_iter, False


def scph_curve(
    T_grid_K: ArrayLike,
    E0_sq_meV2: float,
    A_meV2: float,
    **kwargs,
) -> dict:
    """Run :func:`scph_iterate` across a temperature grid.

    Args:
        T_grid_K: Temperatures (K).
        E0_sq_meV2: Bare squared frequency (meV^2).
        A_meV2: Quartic renormalization constant (meV^2).
        **kwargs: Forwarded to :func:`scph_iterate` (``E_init_meV``, ``alpha``,
            ``tol_meV2``, ``max_iter``).

    Returns:
        Dict with keys ``T_K``, ``E_sq_meV2``, ``n_iter``, ``converged``,
        ``A_meV2``, ``E0_sq_meV2``.
    """
    T = np.asarray(T_grid_K, dtype=float)
    E_sq = np.zeros_like(T)
    n_iter = np.zeros(T.size, dtype=int)
    converged = np.zeros(T.size, dtype=bool)
    for i, t in enumerate(T):
        E_sq[i], n_iter[i], converged[i] = scph_iterate(
            E0_sq_meV2, A_meV2, float(t), **kwargs
        )
    return {
        "T_K": T,
        "E_sq_meV2": E_sq,
        "n_iter": n_iter,
        "converged": converged,
        "A_meV2": float(A_meV2),
        "E0_sq_meV2": float(E0_sq_meV2),
    }


def scph_residual(
    T_K: ArrayLike,
    E_sq_scph: ArrayLike,
    E_sq_dfpt: ArrayLike,
) -> np.ndarray:
    """Dimensionless residual ``delta(T)`` between scalar SCPH and DFPT.

    Equation (3.7) of the manuscript::

        delta(T) = (Omega^2_SCPH(T) - Omega^2_DFPT(T)) / Omega^2_DFPT(T).

    Args:
        T_K: Temperatures (K). (Carried for interface symmetry; not used in the
            residual itself.)
        E_sq_scph: SCPH-predicted squared frequencies (meV^2).
        E_sq_dfpt: DFPT-measured squared frequencies (meV^2).

    Returns:
        Dimensionless residual array (same shape as inputs).
    """
    _ = np.asarray(T_K)  # interface contract
    s = np.asarray(E_sq_scph, dtype=float)
    d = np.asarray(E_sq_dfpt, dtype=float)
    return (s - d) / d


def scph_tc_estimate(
    T_grid_K: ArrayLike,
    E_sq_meV2: ArrayLike,
) -> float | None:
    """Estimate the operational SCPH transition temperature.

    Finds the temperature at which the scalar-SCPH ``Omega^2(T)`` reaches its
    minimum value on the grid, using a parabolic-vertex interpolation through
    the three points bracketing the discrete minimum. For a scalar SCPH with
    ``E_0^2 < 0`` and a stabilizing quartic the minimum is the operational
    transition temperature. Returns ``None`` if the minimum sits at a grid
    endpoint (i.e. ``Omega^2`` is monotonic on the grid, so no interior
    extremum exists and the "Tc" is not resolved by the grid).

    Args:
        T_grid_K: Temperatures (K), strictly increasing.
        E_sq_meV2: SCPH-predicted squared frequencies (meV^2).

    Returns:
        Estimated transition temperature in K, or ``None``.
    """
    T = np.asarray(T_grid_K, dtype=float)
    E2 = np.asarray(E_sq_meV2, dtype=float)
    if T.size < 3:
        return None
    i_min = int(np.argmin(E2))
    if i_min == 0 or i_min == T.size - 1:
        return None  # boundary minimum -> no resolved interior extremum
    x1, x2, x3 = T[i_min - 1], T[i_min], T[i_min + 1]
    y1, y2, y3 = E2[i_min - 1], E2[i_min], E2[i_min + 1]
    denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
    if abs(denom) < 1e-12:
        return float(x2)
    a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
    b = (x3**2 * (y1 - y2) + x2**2 * (y3 - y1) + x1**2 * (y2 - y3)) / denom
    if abs(a) < 1e-12:
        return float(x2)
    return float(-b / (2.0 * a))
