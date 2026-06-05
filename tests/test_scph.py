"""Smoke tests for the fekappa.scph soft-mode fits."""

from __future__ import annotations

import numpy as np
import pytest

from fekappa.scph import fit_soft_mode_two_branch


def test_fit_soft_mode_two_branch_synthetic() -> None:
    """Recover T_c (within 5 K) and both slopes (within 5%) from clean data."""
    rng = np.random.default_rng(0)
    T = np.array(
        [100, 150, 200, 250, 300, 350, 450, 500, 550, 600, 650, 700], dtype=float
    )
    T_c_true, a_below_true, a_above_true = 400.0, 0.4, 0.2
    omega_sq = np.where(
        T < T_c_true,
        a_below_true * (T_c_true - T),
        a_above_true * (T - T_c_true),
    )
    omega_sq = omega_sq + rng.normal(0.0, 0.2, size=T.size)
    omega = np.sqrt(np.clip(omega_sq, 1e-6, None))

    res = fit_soft_mode_two_branch(T, omega, T_c_init=400.0)
    assert abs(res["T_c"] - T_c_true) < 5.0
    assert abs(res["A_below"] - a_below_true) / a_below_true < 0.05
    assert abs(res["A_above"] - a_above_true) / a_above_true < 0.05


def test_fit_soft_mode_two_branch_requires_data_on_both_sides() -> None:
    """A side with fewer than two points raises a clear ValueError."""
    T = np.array([100.0, 200.0, 300.0, 350.0, 450.0])  # only one point above 400
    omega = np.array([10.0, 9.0, 7.0, 5.0, 3.0])
    with pytest.raises(ValueError, match="both branches"):
        fit_soft_mode_two_branch(T, omega, T_c_init=400.0)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-v"]))
