"""Tests for the small-q Lambda^2(T) fit and the Phase-2 convergence tables."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fekappa.lambda_fit import CM1_TO_MEV, fit_lambda_sq

REPO = Path(__file__).resolve().parents[1]
TIDY_DIR = REPO / "data" / "processed" / "sto_alamode"


def _synthetic_bands(omega0_meV: float, lambda_sq: float, n_q: int = 40):
    """Exact quadratic TO1 doublet plus dummy acoustic/upper branches."""
    q = np.linspace(0.0, 0.2, n_q)
    to1 = np.sqrt(omega0_meV**2 + lambda_sq * q**2) / CM1_TO_MEV
    omega = np.zeros((n_q, 6))
    omega[:, 0] = 30000.0 * q  # acoustic stand-ins, far below TO1 here
    omega[:, 1] = 30000.0 * q
    omega[:, 2] = 40000.0 * q
    omega[:, 3] = to1
    omega[:, 4] = to1
    omega[:, 5] = to1 * 3.0
    return q, omega


def test_fit_recovers_synthetic_parameters() -> None:
    """The (q^2, Omega^2) fit recovers exact synthetic Omega_0 and Lambda^2."""
    q, omega = _synthetic_bands(omega0_meV=20.0, lambda_sq=850.0)
    fit = fit_lambda_sq(q, omega, q_max_ang1=0.15)
    assert fit["Omega0_meV"] == pytest.approx(20.0, rel=1e-10)
    assert fit["lambda_sq_meV2_ang2"] == pytest.approx(850.0, rel=1e-10)
    assert fit["R2"] == pytest.approx(1.0, abs=1e-12)


def test_fit_rejects_split_doublet() -> None:
    """A split TO1 doublet must raise, not silently average."""
    q, omega = _synthetic_bands(omega0_meV=20.0, lambda_sq=850.0)
    omega[:, 4] *= 1.05
    with pytest.raises(ValueError, match="doublet"):
        fit_lambda_sq(q, omega, q_max_ang1=0.15)


def test_fit_requires_enough_points() -> None:
    q, omega = _synthetic_bands(omega0_meV=20.0, lambda_sq=850.0, n_q=40)
    with pytest.raises(ValueError, match="points"):
        fit_lambda_sq(q, omega, q_max_ang1=0.01)


# --- committed Phase-2 tables ----------------------------------------------
def test_lambda_sq_csv_full_grid_and_quality() -> None:
    """lambda_sq_T.csv: 0-400 K (below the TO1-TO2 anticrossing), R^2 > 0.99."""
    df = pd.read_csv(TIDY_DIR / "lambda_sq_T.csv", comment="#")
    assert len(df) == 9
    assert df["T_K"].min() == 0.0 and df["T_K"].max() == 400.0
    assert (df["R2"] > 0.99).all()
    row300 = df[df["T_K"] == 300.0].iloc[0]
    # Zone-center SCPH TO1 at 300 K (fit intercept must match the .result value).
    assert row300["Omega0_cm1"] == pytest.approx(167.13, abs=0.1)
    assert row300["lambda_sq_meV2_ang2"] == pytest.approx(842.0, abs=5.0)


def test_lambda_sq_monotonic_hardening() -> None:
    """The renormalized soft mode hardens monotonically with T (SCPH signature)."""
    df = pd.read_csv(TIDY_DIR / "lambda_sq_T.csv", comment="#").sort_values("T_K")
    omega0 = df["Omega0_cm1"].to_numpy()
    assert (np.diff(omega0) > 0).all()


def test_mesh_convergence_table() -> None:
    """mesh_convergence.csv holds both meshes with the locked 8x8x8 kappa."""
    df = pd.read_csv(TIDY_DIR / "mesh_convergence.csv", comment="#")
    assert sorted(df["mesh"]) == ["12x12x12", "8x8x8"]
    k8 = df.loc[df["mesh"] == "8x8x8", "kappa_W_mK"].item()
    assert k8 == pytest.approx(9.5283, abs=1e-4)
    # Coherent channel is present, positive, and small relative to kappa_P.
    assert (df["kappa_coherent_W_mK"] > 0).all()
    assert (df["kappa_coherent_W_mK"] < 0.2 * df["kappa_W_mK"]).all()
