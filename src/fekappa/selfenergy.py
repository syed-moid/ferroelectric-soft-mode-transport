"""Phonon self-energy: three- and four-phonon damping decomposition.

The imaginary part of the phonon self-energy gives the linewidth (damping)
Gamma(T). Near the soft-mode transition the four-phonon channel becomes
important, so the total damping is decomposed into three-phonon and
four-phonon contributions. This module provides the empirical Gamma(T)
polynomial fit plus stubs for the microscopic channel calculations.
"""

from __future__ import annotations

import numpy as np
from lmfit import Model
from lmfit.models import LorentzianModel
from numpy.typing import ArrayLike
from scipy.signal import find_peaks


def _seed_centers(
    freq_cm1: np.ndarray, intensity: np.ndarray, n_peaks: int
) -> list[float]:
    """Seed ``n_peaks`` Lorentzian centers from an intensity spectrum.

    Uses :func:`scipy.signal.find_peaks` with a prominence threshold scaled to
    the signal range and a minimum inter-peak distance. If fewer than
    ``n_peaks`` peaks are found, the list is padded with evenly spaced
    positions across the frequency range.
    """
    span = float(intensity.max() - intensity.min())
    prominence = 0.02 * span if span > 0 else None
    distance = max(1, len(freq_cm1) // (2 * n_peaks))
    idx, props = find_peaks(intensity, prominence=prominence, distance=distance)

    if len(idx) >= n_peaks:
        order = np.argsort(props["prominences"])[::-1][:n_peaks]
        chosen = np.sort(idx[order])
        return freq_cm1[chosen].tolist()

    centers = list(freq_cm1[idx])
    for extra in np.linspace(freq_cm1.min(), freq_cm1.max(), n_peaks + 2)[1:-1]:
        if len(centers) >= n_peaks:
            break
        centers.append(float(extra))
    return sorted(centers[:n_peaks])


def fit_ir_lorentzians(
    freq_cm1: np.ndarray,
    intensity: np.ndarray,
    n_peaks: int,
    init_centers_cm1: list[float] | None = None,
) -> list[dict]:
    """Fit a sum of ``n_peaks`` Lorentzians to an IR spectrum.

    Each component is an :class:`lmfit.models.LorentzianModel`, whose ``sigma``
    parameter is the half-width at half-maximum (HWHM). Centers are constrained
    to be positive and HWHM/amplitudes non-negative.

    Args:
        freq_cm1: Frequency grid (cm^-1).
        intensity: IR intensity sampled on ``freq_cm1``.
        n_peaks: Number of Lorentzians to fit.
        init_centers_cm1: Optional initial center guesses (cm^-1). When
            ``None``, centers are seeded with :func:`scipy.signal.find_peaks`.

    Returns:
        A list of ``n_peaks`` dicts (sorted by ``center_cm1``), each with keys
        ``center_cm1``, ``hwhm_cm1``, ``amplitude`` (lmfit integrated area),
        ``stderr_center`` and ``stderr_hwhm``.

    Raises:
        ValueError: If ``init_centers_cm1`` is given but its length does not
            match ``n_peaks``.
    """
    freq_cm1 = np.asarray(freq_cm1, dtype=float)
    intensity = np.asarray(intensity, dtype=float)

    if init_centers_cm1 is None:
        centers = _seed_centers(freq_cm1, intensity, n_peaks)
    else:
        if len(init_centers_cm1) != n_peaks:
            raise ValueError(
                f"init_centers_cm1 has {len(init_centers_cm1)} entries but "
                f"n_peaks={n_peaks}"
            )
        centers = sorted(float(c) for c in init_centers_cm1)

    span = float(freq_cm1.max() - freq_cm1.min())
    sigma0 = max(span / 100.0, 1e-3)

    model = None
    params = None
    for k, center in enumerate(centers):
        prefix = f"p{k}_"
        component = LorentzianModel(prefix=prefix)
        height = float(intensity[np.argmin(np.abs(freq_cm1 - center))])
        comp_params = component.make_params()
        comp_params[f"{prefix}center"].set(value=center, min=0.0)
        comp_params[f"{prefix}sigma"].set(value=sigma0, min=1e-6)
        comp_params[f"{prefix}amplitude"].set(
            value=max(height, 1e-9) * np.pi * sigma0, min=0.0
        )
        if model is None:
            model, params = component, comp_params
        else:
            model = model + component
            params.update(comp_params)

    result = model.fit(intensity, params, x=freq_cm1)

    peaks: list[dict] = []
    for k in range(n_peaks):
        prefix = f"p{k}_"
        center = result.params[f"{prefix}center"]
        sigma = result.params[f"{prefix}sigma"]
        amplitude = result.params[f"{prefix}amplitude"]
        peaks.append(
            {
                "center_cm1": float(center.value),
                "hwhm_cm1": float(sigma.value),
                "amplitude": float(amplitude.value),
                "stderr_center": (
                    float(center.stderr) if center.stderr is not None else float("nan")
                ),
                "stderr_hwhm": (
                    float(sigma.stderr) if sigma.stderr is not None else float("nan")
                ),
            }
        )
    peaks.sort(key=lambda d: d["center_cm1"])
    return peaks


def _gamma_quadratic(T: np.ndarray, a: float, b: float, c: float) -> np.ndarray:
    """Quadratic damping model Gamma(T) = a + b*T + c*T^2."""
    T = np.asarray(T, dtype=float)
    return a + b * T + c * T**2


def fit_gamma_polynomial(
    T: ArrayLike,
    gamma: ArrayLike,
    gamma_err: ArrayLike | None = None,
    use_lmfit: bool = True,
) -> dict:
    """Fit ``Gamma(T) = a + b*T + c*T^2`` to phonon HWHM data.

    Physical interpretation of the three coefficients:

    * ``a`` -- boundary / impurity / zero-T contribution (small, positive).
    * ``b`` -- three-phonon scattering-rate coefficient (cubic bubble, linear
      in T at high T).
    * ``c`` -- four-phonon scattering-rate coefficient (quartic loop, quadratic
      in T).

    All three coefficients are constrained non-negative (the damping and both
    physical anharmonic channels are positive). The 3-phonon and 4-phonon terms
    are equal at the crossover temperature ``T = b / c``.

    Args:
        T: Temperatures (K).
        gamma: Linewidth / damping (HWHM) at each temperature.
        gamma_err: Optional per-point uncertainties on ``gamma``; when ``None``,
            uniform unit weights are used.
        use_lmfit: If True (default), fit with a bounded lmfit model; if False,
            fall back to an unconstrained ``numpy`` least-squares fit with
            errors from the covariance matrix.

    Returns:
        A dict with keys ``a``, ``a_stderr``, ``b``, ``b_stderr``, ``c``,
        ``c_stderr``, ``chi_squared``, ``reduced_chi_squared``, ``n_data``,
        ``n_params``, ``gamma_3ph_500K`` (``b * 500``), ``gamma_4ph_500K``
        (``c * 500**2``), ``crossover_temperature_K`` (``b / c`` if ``c > 0``
        else ``inf``) and ``fit_report`` (str).
    """
    T = np.asarray(T, dtype=float)
    gamma = np.asarray(gamma, dtype=float)
    err = None if gamma_err is None else np.asarray(gamma_err, dtype=float)

    if use_lmfit:
        model = Model(_gamma_quadratic)
        params = model.make_params(a=max(gamma.min(), 0.0), b=0.0, c=0.0)
        for name in ("a", "b", "c"):
            params[name].set(min=0.0)
        weights = None if err is None else 1.0 / err
        result = model.fit(gamma, params, T=T, weights=weights)

        def _stderr(name: str) -> float:
            e = result.params[name].stderr
            return float(e) if e is not None else float("nan")

        a = float(result.params["a"].value)
        b = float(result.params["b"].value)
        c = float(result.params["c"].value)
        a_err, b_err, c_err = _stderr("a"), _stderr("b"), _stderr("c")
        chi_squared = float(result.chisqr)
        reduced_chi_squared = float(result.redchi)
        n_data = int(result.ndata)
        n_params = int(result.nvarys)
        fit_report = result.fit_report()
    else:
        # Unconstrained weighted polynomial least squares (degree 2).
        w = None if err is None else 1.0 / err
        coeffs, cov = np.polyfit(T, gamma, deg=2, w=w, cov=True)
        c, b, a = (float(v) for v in coeffs)  # polyfit returns highest power first
        c_err, b_err, a_err = (float(np.sqrt(cov[i, i])) for i in range(3))
        resid = gamma - _gamma_quadratic(T, a, b, c)
        weights_sq = np.ones_like(gamma) if err is None else 1.0 / err**2
        chi_squared = float(np.sum(weights_sq * resid**2))
        n_data = int(T.size)
        n_params = 3
        reduced_chi_squared = (
            chi_squared / (n_data - n_params) if n_data > n_params else float("nan")
        )
        fit_report = (
            f"numpy polyfit (unconstrained): a={a:.6g}, b={b:.6g}, c={c:.6g}"
        )

    crossover = float(b / c) if c > 0 else float("inf")
    return {
        "a": a,
        "a_stderr": a_err,
        "b": b,
        "b_stderr": b_err,
        "c": c,
        "c_stderr": c_err,
        "chi_squared": chi_squared,
        "reduced_chi_squared": reduced_chi_squared,
        "n_data": n_data,
        "n_params": n_params,
        "gamma_3ph_500K": float(b * 500.0),
        "gamma_4ph_500K": float(c * 500.0**2),
        "crossover_temperature_K": crossover,
        "fit_report": fit_report,
    }


def gamma_three_phonon(*args, **kwargs):
    """Three-phonon contribution to the soft-mode damping.

    Raises:
        NotImplementedError: Always, for now.
    """
    raise NotImplementedError("Turn N: three-phonon damping channel")


def gamma_four_phonon(*args, **kwargs):
    """Four-phonon contribution to the soft-mode damping.

    Raises:
        NotImplementedError: Always, for now.
    """
    raise NotImplementedError("Turn N: four-phonon damping channel")
