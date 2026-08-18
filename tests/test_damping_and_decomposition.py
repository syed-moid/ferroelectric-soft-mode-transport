"""Tests for the damping map, kappa decomposition, and recomputed soft kernel.

All assertions run against the committed CSVs under
``data/processed/sto_alamode/`` -- the data of record -- so they hold without
the local-only raw ALAMODE files.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fekappa.kappa_kernel import K_soft

REPO = Path(__file__).resolve().parents[1]
TIDY = REPO / "data" / "processed" / "sto_alamode"


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(TIDY / name, comment="#")


# --- units gate: mode-resolved kappa_P closes on the anphon value ------------
def test_mode_sum_reproduces_anphon_kappa() -> None:
    """Sum of committed per-mode kappa contributions matches the .kl kappa.

    This is the dimensional gate: C [J/K] * tau [s] * v^2 [m^2/s^2] / V [m^3]
    summed over the mesh must land on the anphon W/mK value (<= 0.5%; the
    observed closure is +0.07-0.27%, the floor set by the 6-significant-figure
    printed .result inputs).
    """
    modes = _read("mode_kappa_T.csv")
    kappa = _read("kappa_rta_T.csv")
    for t, group in modes.groupby("T_K"):
        ref = kappa.loc[kappa["T_K"] == t, "kappa_W_mK"].item()
        total = group["kappa_mode_W_mK"].sum()
        if ref == 0.0:
            assert total == 0.0
        else:
            assert total == pytest.approx(ref, rel=5e-3), f"T={t}"


def test_five_way_partition_sums_and_closure() -> None:
    """Sectors are exhaustive: P sectors sum to kappa_P; total = P + C."""
    df = _read("kappa_five_way_T.csv")
    p_sum = df["kappa_acoustic_P"] + df["kappa_soft_optic_P"] + df["kappa_hard_optic_P"]
    assert np.allclose(p_sum, df["kappa_P_modesum"], rtol=1e-12)
    assert np.allclose(
        df["kappa_total"], df["kappa_P_modesum"] + df["kappa_C"], rtol=1e-12
    )
    # Mode sum closes on the anphon .kl value at every T and boundary.
    nonzero = df[df["kappa_P_anphon"] > 0]
    assert np.allclose(nonzero["kappa_P_modesum"], nonzero["kappa_P_anphon"], rtol=5e-3)


def test_five_way_boundary_sensitivity_rows() -> None:
    """All three sector boundaries present; the acoustic sector ignores them."""
    df = _read("kappa_five_way_T.csv")
    assert sorted(df["boundary_cm1"].unique()) == [150.0, 175.0, 200.0]
    for _, group in df.groupby("T_K"):
        assert group["kappa_acoustic_P"].nunique() == 1
        assert group["kappa_C"].nunique() == 1


def test_mode_partition_disjoint_and_complete() -> None:
    """Every mode carries exactly one sector; branches 1-3 are the acoustic tag."""
    modes = _read("mode_kappa_T.csv")
    assert set(modes["sector"].unique()) == {"acoustic", "soft_optic", "hard_optic"}
    assert (modes.loc[modes["branch"] <= 3, "sector"] == "acoustic").all()
    assert (modes.loc[modes["branch"] > 3, "sector"] != "acoustic").all()
    # 525 modes per temperature, each counted once.
    assert modes.groupby("T_K").size().eq(525).all()


# --- eta map and validity ----------------------------------------------------
def test_eta_map_matches_gamma_over_omega() -> None:
    eta_map = _read("eta_map_T.csv")
    recomputed = eta_map["gamma_hwhm_cm1"] / eta_map["omega_cm1"]
    assert np.allclose(eta_map["eta"], recomputed, rtol=1e-12)
    assert (eta_map["eta"] >= 0).all()


def test_narrow_width_holds_zone_wide() -> None:
    """Headline result: no mode beyond eta = 0.5 at any T, 0-1000 K."""
    val = _read("eta_validity_T.csv")
    assert len(val) == 21
    assert (val["frac_eta_gt_0p5_count"] == 0).all()
    assert (val["frac_eta_gt_1_count"] == 0).all()
    assert (val["frac_eta_gt_0p5_kappa"] == 0).all()
    assert (val["eta_max"] < 0.15).all()


def test_zone_center_eta_dual_reporting() -> None:
    """eta_0(300 K) both ways, and the ratio to Vogt's literature 0.109."""
    val = _read("eta_validity_T.csv")
    r = val[val["T_K"] == 300.0].iloc[0]
    assert r["eta0_computed"] == pytest.approx(0.0225, abs=5e-4)
    assert r["eta0_vogt_anchored"] == pytest.approx(0.0426, abs=5e-4)
    # Computed eta_0 sits ~5x below the literature value -- reported, not tuned.
    assert r["eta0_computed"] / 0.109 == pytest.approx(0.207, abs=0.01)
    # The fixed-cell SCPH Omega overestimate at 300 K (C8): Vogt/computed ~ 0.53.
    assert r["vogt_over_computed"] == pytest.approx(0.529, abs=0.005)


def test_eta_sensitivity_conclusions_robust() -> None:
    """Vogt-rescaled soft-branch eta moves no mode across 0.5 (150-300 K).

    Covers all three rescalings: Omega-side (eta/r_omega), Gamma-side
    (r_gamma*eta with the ~2.6x 300 K HWHM ratio), and both combined.
    """
    s = _read("eta_validity_sensitivity.csv")
    assert list(s["T_K"]) == [150.0, 200.0, 250.0, 300.0]
    assert (s["rescale_ratio"] < 1).all()  # Omega rescaling can only RAISE soft eta
    assert (s["frac_eta_gt_0p5_count_rescaled"] == 0).all()
    assert (s["frac_eta_gt_1_kappa_rescaled"] == 0).all()
    assert (s["eta_max_rescaled"] < 0.15).all()
    # Gamma-side: disclosed Gamma_0 underestimate is ~2.6x vs Vogt at 300 K.
    assert s["gamma_ratio_300K"].iloc[0] == pytest.approx(2.57, abs=0.02)
    assert (s["frac_eta_gt_0p5_count_gamma"] == 0).all()
    assert (s["frac_eta_gt_0p5_kappa_gamma"] == 0).all()
    # Combined (strictest): still nowhere near the narrow-width boundary.
    assert (s["frac_eta_gt_0p5_count_both"] == 0).all()
    assert (s["frac_eta_gt_0p5_kappa_both"] == 0).all()
    assert (s["eta_max_both_rescaled"] < 0.2).all()


# --- soft kernel -------------------------------------------------------------
def test_kernel_reduces_to_old_literature_path() -> None:
    """K_soft with the old literature inputs reproduces the committed triple."""
    triples = pd.read_csv(REPO / "data" / "processed" / "section4_eta_kappa_triples.csv")
    sto = triples[triples["material"] == "STO"].iloc[0]
    kappa = K_soft(
        sto["omega0_meV"] ** 2, sto["lambda_sq_meV2ang2"], sto["eta"],
        sto["T_K"], 59.55, warn=False,
    )
    assert kappa == pytest.approx(sto["kappa_soft_W_m_K"], abs=5e-4)


def test_soft_kernel_variants_at_300K() -> None:
    """Computed and Vogt-anchored kernels at 300 K, vs the old 0.137 W/mK."""
    df = _read("soft_kernel_T.csv")
    at300 = df[(df["T_K"] == 300.0) & (df["q_max"] == "debye")]
    computed = at300[at300["variant"] == "computed"]["kappa_soft_W_mK"].item()
    anchored = at300[at300["variant"] == "vogt_anchored"]["kappa_soft_W_mK"].item()
    assert computed == pytest.approx(0.3217, abs=2e-3)
    assert anchored == pytest.approx(0.3382, abs=2e-3)
    # Gate-3 sanity: within ~3x of the old literature-input kernel...
    assert computed / 0.13712 < 3.0
    # ...and still a small fraction of kappa_total (order 1%, here ~3%).
    kappa_total = _read("kappa_rta_T.csv")
    k300 = kappa_total.loc[kappa_total["T_K"] == 300.0, "kappa_W_mK"].item()
    assert computed / k300 < 0.05


def test_soft_kernel_qc_envelope_ordering() -> None:
    """The q_c-truncated kernel grows with q_c and stays below the Debye value."""
    df = _read("soft_kernel_T.csv")
    for (_t, _v), group in df.groupby(["T_K", "variant"]):
        by_q = group.set_index("q_max")["kappa_soft_W_mK"]
        assert by_q["0.10"] < by_q["0.15"] < by_q["0.20"] <= by_q["debye"]


# --- bare-vs-SCPH negative control ------------------------------------------
def test_bare_vs_scph_control_table() -> None:
    df = _read("bare_vs_scph_control.csv")
    kap = df[df["record"] == "kappa_T"]
    assert len(kap) == 21
    r300 = kap[kap["T_K"] == 300.0].iloc[0]
    assert r300["bare"] == pytest.approx(3.3439, abs=1e-4)
    assert r300["scph"] == pytest.approx(9.5283, abs=1e-4)
    assert r300["scph_over_bare"] == pytest.approx(2.85, abs=0.01)
    tau = df[df["record"] == "tau_band_300K"]
    assert not tau.empty
    # SCPH lengthens every band-median lifetime; low-omega bands by the most.
    assert (tau["scph_over_bare"] >= 1.0).all()
    low = tau[tau["band_hi_cm1"] <= 200.0]["scph_over_bare"]
    assert (low > 3.0).all()