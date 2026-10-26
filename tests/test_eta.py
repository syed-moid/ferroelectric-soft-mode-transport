"""Tests for fekappa.eta -- the eta = Gamma_HWHM/Omega damping ratio module.

Locks the v1 results:
  * eta_STO ~ 0.109 and eta_KTaO3 ~ 0.126 from the dense Vogt 1995 sources,
  * the two quantum-paraelectrics cluster (|delta| small, both << 0.5),
  * the Vogt FWHM->HWHM /2 is applied exactly once (regression vs double/none),
  * Yamada-Shirane gamma is consumed directly as HWHM (no /2).
"""

import numpy as np
import pytest

from fekappa import eta


# --- convention conversion -------------------------------------------------
def test_gamma_to_hwhm_factors():
    assert eta.gamma_to_hwhm(3.04, "FWHM") == pytest.approx(1.52)   # Vogt: halved
    assert eta.gamma_to_hwhm(0.84, "HWHM") == pytest.approx(0.84)   # YS: direct
    np.testing.assert_allclose(
        eta.gamma_to_hwhm([2.0, 4.0], "FWHM"), [1.0, 2.0]
    )


def test_gamma_convention_rejects_unknown():
    # A missing/unknown convention must fail loudly, never silently use a default.
    with pytest.raises(ValueError):
        eta.gamma_to_hwhm(1.0, "HALF")
    with pytest.raises(ValueError):
        eta.gamma_to_hwhm(1.0, "")


# --- primary eta_fit on the committed Vogt sources -------------------------
def test_eta_sto_vogt_reproduces_0p109():
    _T, e = eta.eta_sto_vogt()
    # ~300 K endpoint is the locked value (dense marker-center, ±0.02).
    assert e[-1] == pytest.approx(0.109, abs=0.02)
    assert eta.eta_underdamped(e)


def test_eta_kto_vogt_reproduces_0p126():
    _T, e = eta.eta_kto_vogt()
    assert e[-1] == pytest.approx(0.126, abs=0.02)
    assert eta.eta_underdamped(e)


def test_two_qpe_cluster():
    _, es = eta.eta_sto_vogt()
    _, ek = eta.eta_kto_vogt()
    # Two underdamped quantum paraelectrics cluster; both far below threshold.
    assert abs(es[-1] - ek[-1]) < 0.05
    assert es[-1] < 0.5 and ek[-1] < 0.5


# --- regression: the /2 is applied exactly once ----------------------------
def test_fwhm_half_applied_exactly_once():
    o = __import__("fekappa.literature_io", fromlist=["x"])
    so = o.load_vogt1995_sto_omega_T()
    so = so[so["phase_mode"] == "cubic_F1u"]
    sg = o.load_vogt1995_sto_gamma_T()
    sg = sg[sg["phase_mode"] == "cubic_F1u"]

    _, e_fwhm = eta.eta_fit(so["T_K"], so["Omega_meV"], sg["T_K"], sg["gamma_meV"],
                            gamma_convention="FWHM")
    _, e_hwhm = eta.eta_fit(so["T_K"], so["Omega_meV"], sg["T_K"], sg["gamma_meV"],
                            gamma_convention="HWHM")
    # Treating the same FWHM data as HWHM yields exactly 2x (no /2). The FWHM path
    # must be exactly half of that -- i.e. the factor is applied once, not zero or
    # two times.
    np.testing.assert_allclose(e_fwhm, 0.5 * e_hwhm, rtol=1e-12)
    assert e_fwhm[-1] == pytest.approx(0.109, abs=0.02)         # single application
    assert e_fwhm[-1] != pytest.approx(0.109 / 2, abs=0.005)    # not double-applied
    assert e_fwhm[-1] != pytest.approx(0.109 * 2, abs=0.005)    # not un-applied


# --- Yamada-Shirane cross-check: gamma used directly as HWHM ---------------
def test_ys_gamma_is_hwhm_direct():
    from fekappa import literature_io as lit

    yo = lit.load_yamada_shirane_omega_T()
    yg = lit.load_yamada_shirane_gamma_T()
    _T, e = eta.eta_fit(yo["T_K"], yo["omega_meV"], yg["T_K"], yg["gamma_meV"],
                       gamma_convention="HWHM")
    # YS is underdamped and the value must be the raw gamma/Omega (no /2 baked in):
    # halving it (wrong convention) would shift every point by exactly 2x.
    assert eta.eta_underdamped(e)
    _, e_wrong = eta.eta_fit(yo["T_K"], yo["omega_meV"], yg["T_K"], yg["gamma_meV"],
                             gamma_convention="FWHM")
    np.testing.assert_allclose(e_wrong, 0.5 * e, rtol=1e-12)


# --- bracket fallback ------------------------------------------------------
def test_eta_bracket_fallback():
    lo, hi = eta.eta_bracket(10.0, 1.0, 3.0, gamma_convention="FWHM")
    assert (lo, hi) == pytest.approx((0.05, 0.15))  # (0.5/10, 1.5/10)


def test_eta_underdamped_threshold():
    assert eta.eta_underdamped([0.1, 0.2, 0.49])
    assert not eta.eta_underdamped([0.1, 0.6])
    assert eta.eta_underdamped(2.3, threshold=3.0)   # BTO overdamped vs custom thr
    assert not eta.eta_underdamped(2.3)              # BTO overdamped at default 0.5
