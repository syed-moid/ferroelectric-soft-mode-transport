"""Smoke tests for the fekappa.selfenergy damping fits."""

from __future__ import annotations

import numpy as np

from fekappa.selfenergy import fit_gamma_polynomial


def test_fit_gamma_polynomial_synthetic() -> None:
    """Recover (a, b, c) within ~10% and crossover b/c = 600 K within 50 K."""
    rng = np.random.default_rng(0)
    a_true, b_true, c_true = 0.5, 0.003, 5e-6
    T = np.linspace(50.0, 700.0, 40)
    gamma = a_true + b_true * T + c_true * T**2
    gamma = gamma + rng.normal(0.0, 0.03, size=T.size)

    res = fit_gamma_polynomial(T, gamma)
    assert abs(res["a"] - a_true) / a_true < 0.10
    assert abs(res["b"] - b_true) / b_true < 0.10
    assert abs(res["c"] - c_true) / c_true < 0.10
    # Crossover where b*T = c*T^2  ->  T = b/c = 600 K.
    assert abs(res["crossover_temperature_K"] - 600.0) < 50.0
