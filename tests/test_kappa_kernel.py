"""Tests for fekappa.kappa_kernel -- the K_soft particle-like soft-mode kernel."""

import numpy as np
import pytest

from fekappa import kappa_kernel as kk


# --- building blocks -------------------------------------------------------
def test_n_bose_limits():
    # High-T classical limit: n_B -> k_B T / Omega.
    assert kk.n_bose(1.0, 1e4) == pytest.approx(8.617333262e-2 * 1e4 / 1.0, rel=1e-3)
    # Frozen-out (Omega >> k_B T) -> ~0; T <= 0 -> exactly 0.
    assert kk.n_bose(100.0, 1.0) == pytest.approx(0.0, abs=1e-12)
    assert kk.n_bose(10.0, 0.0) == 0.0
    assert kk.n_bose(10.0, -5.0) == 0.0


def test_heat_capacity_limits():
    # Dulong-Petit per mode: C -> k_B as Omega/k_B T -> 0.
    assert kk.heat_capacity_mode(1e-4, 300.0) == pytest.approx(kk.KB_J_K, rel=1e-4)
    # Frozen-out at low T (no overflow), and exactly 0 at T <= 0.
    c_low = kk.heat_capacity_mode(50.0, 1.0)
    assert 0.0 <= c_low < 1e-30
    assert kk.heat_capacity_mode(10.0, 0.0) == 0.0


def test_debye_qmax():
    assert kk.debye_qmax(59.55) == pytest.approx((6 * np.pi**2 / 59.55) ** (1 / 3))
    assert kk.debye_qmax(59.55) == pytest.approx(0.9981, abs=1e-3)


# --- K_soft core -----------------------------------------------------------
def _sto():
    # Omega_0(300 K) ~ 10.5 meV (Verdi exp); Lambda^2 = 889; eta = 0.109.
    return dict(omega0_sq_meV2=10.5**2, lambda_sq_meV2_ang2=889.0,
                eta=0.109, T_K=300.0, V_cell_ang3=59.55)


def test_K_soft_scalar_positive_finite_and_physical():
    k = kk.K_soft(**_sto())
    assert np.isscalar(k) or np.ndim(k) == 0
    assert np.isfinite(k) and k > 0
    # Soft-mode-only contribution: small but nonzero, far below total STO kappa(~10).
    assert 1e-4 < k < 1e1


def test_K_soft_vectorized_over_T():
    p = _sto()
    p.pop("T_K")
    p.pop("omega0_sq_meV2")
    T = np.array([100.0, 200.0, 300.0])
    k = kk.K_soft(omega0_sq_meV2=10.5**2, T_K=T, **p)
    assert k.shape == T.shape
    assert np.all(np.isfinite(k)) and np.all(k > 0)


def test_kappa_scales_inverse_eta():
    # tau ~ 1/Gamma ~ 1/eta  =>  kappa ~ 1/eta at fixed Omega_0, Lambda^2.
    base = _sto()
    k1 = kk.K_soft(**base)
    k2 = kk.K_soft(**{**base, "eta": 2 * base["eta"]})
    assert k2 == pytest.approx(0.5 * k1, rel=1e-6)


def test_zero_lambda_no_transport():
    # Dispersionless soft mode (Lambda^2 = 0) -> v = 0 -> kappa = 0.
    assert kk.K_soft(**{**_sto(), "lambda_sq_meV2_ang2": 0.0}) == pytest.approx(0.0, abs=1e-30)


def test_n_q_converged():
    base = _sto()
    k512 = kk.K_soft(**{**base, "n_q": 512})
    k2048 = kk.K_soft(**{**base, "n_q": 2048})
    assert k512 == pytest.approx(k2048, rel=1e-2)


# --- guards: unstable mode and overdamping ---------------------------------
def test_unstable_mode_skipped_with_warning():
    # Omega_0^2 < 0 (below T_c / Schmidt ST-18O / Ehsan 0 K) -> NaN + warning.
    with pytest.warns(RuntimeWarning, match="Omega_0\\^2 < 0"):
        k = kk.K_soft(**{**_sto(), "omega0_sq_meV2": -424.0})
    assert np.isnan(k)


def test_unstable_row_in_array_skipped():
    with pytest.warns(RuntimeWarning):
        k = kk.K_soft(omega0_sq_meV2=np.array([-100.0, 110.25, 200.0]),
                      lambda_sq_meV2_ang2=889.0, eta=0.109,
                      T_K=np.array([50.0, 300.0, 500.0]), V_cell_ang3=59.55)
    assert np.isnan(k[0]) and np.all(np.isfinite(k[1:]))


def test_overdamped_warns_but_returns_value():
    with pytest.warns(kk.OverdampedModeWarning, match="(?i)order-of-magnitude"):
        k = kk.K_soft(omega0_sq_meV2=15.0**2, lambda_sq_meV2_ang2=972.0,
                      eta=2.3, T_K=400.0, V_cell_ang3=64.63)
    assert np.isfinite(k) and k > 0


def test_underdamped_no_overdamped_warning():
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")          # any warning -> test failure
        kk.K_soft(**_sto())


# --- K_soft_sweep: scalar vs directional dict ------------------------------
def test_sweep_scalar_matches_K_soft():
    base = _sto()
    L2 = base.pop("lambda_sq_meV2_ang2")
    assert kk.K_soft_sweep(base["omega0_sq_meV2"], L2, base["eta"],
                           base["T_K"], base["V_cell_ang3"]) == pytest.approx(kk.K_soft(**_sto()))


def test_sweep_dict_runs_soft_stiff_bracket():
    out = kk.K_soft_sweep(15.0**2, {"soft_100": 972.0, "stiff_110_111": 4750.0},
                          0.109, 400.0, 64.63)   # eta forced underdamped to isolate Lambda^2
    assert set(out) == {"soft_100", "stiff_110_111"}
    # Stiffer branch -> larger group velocity -> larger kappa (the upper bracket).
    assert out["stiff_110_111"] > out["soft_100"] > 0


# --- material registry -----------------------------------------------------
def test_material_registry_consistency():
    for m, d in kk.MATERIALS.items():
        assert d["V_cell_ang3"] == pytest.approx(d["a_ang"] ** 3, rel=1e-3), m


def test_K_soft_material_sto_kto_bto():
    ksto = kk.K_soft_material("STO", 10.5**2, 300.0)
    assert np.isfinite(ksto) and ksto > 0
    kkto = kk.K_soft_material("KTO", 10.3**2, 300.0)
    assert np.isfinite(kkto) and kkto > 0
    # BTO -> directional dict + overdamped warning.
    with pytest.warns(kk.OverdampedModeWarning):
        kbto = kk.K_soft_material("BTO", 15.0**2, 400.0)
    assert set(kbto) == {"soft_100", "stiff_110_111"}
    assert kbto["stiff_110_111"] > kbto["soft_100"] > 0


def test_K_soft_material_unknown_raises():
    with pytest.raises(ValueError, match="unknown material"):
        kk.K_soft_material("PbTiO3", 100.0, 300.0)
