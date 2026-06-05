"""Tests for ``fekappa.literature_io`` loaders (Turn 4B §2.1).

One test per loader, verifying shape, expected columns (implicit -- the loader
raises if columns disagree), and key sentinel values matching the headline
assertions in ``notebooks/05b_data_review.ipynb`` (Turn 4B brief §1).
"""

from __future__ import annotations

import math
from types import SimpleNamespace

import pandas as pd
import pytest

from fekappa import literature_io as lit


# ---------------------------------------------------------------------------
# BTO Ehsan 2021
# ---------------------------------------------------------------------------
def test_load_ehsan_anchors() -> None:
    ns = lit.load_ehsan_anchors()
    assert isinstance(ns, SimpleNamespace)
    assert math.isclose(ns.T0_Omega_sq_meV2, -424.0, abs_tol=0.5)
    assert math.isclose(ns.T500_Omega_meV, 21.08, abs_tol=0.02)
    assert math.isclose(ns.T500_Omega_sq_meV2, 444.4, abs_tol=0.5)
    assert "Ehsan2021" in ns.provenance


def test_load_ehsan_scp_omega_T() -> None:
    df = lit.load_ehsan_scp_omega_T()
    assert isinstance(df, pd.DataFrame)
    assert df.empty  # by design: no Omega(T) curve in Ehsan paper
    assert "Ehsan2021" in df.attrs["provenance"]


def test_load_ehsan_free_energy_T() -> None:
    df = lit.load_ehsan_free_energy_T()
    assert isinstance(df, pd.DataFrame)
    assert df.shape == (20, 4)
    # Ehsan sets F = 0 at the cubic-tetra crossing (T ≈ 455 K), not at 0 K.
    crossing = df[df["T_K"] == 455]
    assert not crossing.empty
    assert math.isclose(float(crossing["F_cubic_meV"].iloc[0]), 0.0, abs_tol=1.0)
    assert math.isclose(float(crossing["F_tetra_meV"].iloc[0]), 0.0, abs_tol=1.0)
    assert "Ehsan2021" in df.attrs["provenance"]


def test_load_ehsan_structural() -> None:
    df = lit.load_ehsan_structural()
    cubic = df[df["phase"] == "cubic"].iloc[0]
    assert math.isclose(float(cubic["a_angstrom"]), 4.013, abs_tol=0.001)
    assert math.isclose(float(cubic["V_cell_A3"]), 64.63, abs_tol=0.02)
    tetra = df[df["phase"] == "tetragonal"].iloc[0]
    assert math.isclose(float(tetra["c_over_a"]), 1.007, abs_tol=0.001)


def test_load_ehsan_cubic_dispersion() -> None:
    df = lit.load_ehsan_cubic_dispersion()
    # Turn 4B unit fix: q_inv_angstrom is now Å^-1, X point should be ~0.78 Å^-1.
    x_rows = df[df["q_label"] == "X"]
    assert not x_rows.empty
    q_X = float(x_rows["q_inv_angstrom"].iloc[0])
    assert math.isclose(q_X, 0.7829, abs_tol=0.001), (
        "Ehsan cubic_dispersion q(X) should be 0.7829 Å^-1 after Turn 4B "
        f"unit fix; got {q_X}"
    )
    # Soft mode at Γ matches the anchors.
    soft_g_0K = df[
        (df["branch_label"] == "Ferroelectric_Soft_Mode")
        & (df["q_label"] == "Gamma")
        & (df["condition"] == "nonthermalized")
    ]
    assert math.isclose(float(soft_g_0K["omega_cm_inv"].iloc[0]), -166.0, abs_tol=1.0)


def test_load_ehsan_born_charges() -> None:
    df = lit.load_ehsan_born_charges()
    cubic_Ti = df[(df["phase"] == "cubic") & (df["atom"] == "Ti")].iloc[0]
    # Z*(Ti) in cubic BaTiO3 is the large polar-active charge.
    assert math.isclose(float(cubic_Ti["Z_xx"]), 7.26, abs_tol=0.01)
    assert math.isclose(float(cubic_Ti["epsilon_inf_xx"]), 6.00, abs_tol=0.05)


# ---------------------------------------------------------------------------
# BTO Tadano 2022
# ---------------------------------------------------------------------------
def test_load_tadano2022_free_energy_T() -> None:
    df = lit.load_tadano2022_free_energy_T()
    assert df.shape[1] == 5  # T_K + 4 phases
    assert df.shape[0] >= 5
    assert "Masuki-Tadano2022" in df.attrs["provenance"]


def test_load_tadano2022_scp_dispersion_plus1pct() -> None:
    df = lit.load_tadano2022_scp_dispersion_plus1pct()
    # +1% lattice flag should appear in provenance.
    assert "+1% lattice" in df.attrs["provenance"]
    assert df.shape[0] > 0


def test_load_tadano2022_transition_temperatures() -> None:
    df = lit.load_tadano2022_transition_temperatures()
    ctt = df[df["transition"] == "cubic-tetra"].iloc[0]
    assert math.isclose(float(ctt["T_c_calc_K"]), 606.0, abs_tol=1.0)
    assert math.isclose(float(ctt["T_c_exp_K"]), 390.0, abs_tol=1.0)


def test_load_tadano2022_kappa_total_T() -> None:
    df = lit.load_tadano2022_kappa_total_T()
    calc = df[df["method"] == "Tadano_calc"]
    k_200 = float(calc[calc["T_K"] == 200]["kappa_W_per_m_K"].iloc[0])
    k_800 = float(calc[calc["T_K"] == 800]["kappa_W_per_m_K"].iloc[0])
    assert math.isclose(k_200, 10.3, abs_tol=0.1)
    assert math.isclose(k_800, 5.4, abs_tol=0.1)


# ---------------------------------------------------------------------------
# BTO Shirane 1970
# ---------------------------------------------------------------------------
def test_load_shirane_Lambda_sq() -> None:
    ns = lit.load_shirane_Lambda_sq()
    assert isinstance(ns, SimpleNamespace)
    assert math.isclose(ns.central_meV2_A2, 5100.0, abs_tol=1.0)
    assert ns.T_K == 295
    assert "Shirane1970" in ns.provenance


def test_load_shirane_cubic_dispersion() -> None:
    df = lit.load_shirane_cubic_dispersion()
    # 295 K tetragonal, both Σ (Fig. 3) and Δ (Fig. 6) branches present.
    assert all(df["T_K"] == 295)
    assert all(df["phase"] == "tetragonal")
    assert set(df["source_figure"].unique()) == {"Fig_3", "Fig_6"}


# ---------------------------------------------------------------------------
# STO Tadano 2015
# ---------------------------------------------------------------------------
def test_load_tadano2015_scph_omega_T() -> None:
    df = lit.load_tadano2015_scph_omega_T()
    methods = set(df["method"].unique())
    expected = {"SCPH_with_PM", "SCPH_without_PM", "TDEP",
                "exp_Servoin", "exp_YamadaShirane"}
    assert methods == expected
    modes = set(df["mode"].unique())
    assert modes == {"TO1", "TO2"}


def test_load_tadano2015_scph_dispersion() -> None:
    df = lit.load_tadano2015_scph_dispersion()
    # SCPH soft optical at Γ_mid ≈ 12.4 meV at 300 K.
    soft = df[
        (df["branch"] == "Soft_Optical")
        & (df["method"] == "SCPH")
    ].iloc[0]
    assert math.isclose(float(soft["omega_meV"]), 12.4, abs_tol=0.2)


def test_load_tadano2015_scph_convergence_300K() -> None:
    df = lit.load_tadano2015_scph_convergence_300K()
    # FE soft converges from 17.85 (2x2x2) toward 16.74 meV at fine grids.
    fe_8x8x8 = df[df["branch"] == "FE_8x8x8"]
    assert not fe_8x8x8.empty
    val = float(fe_8x8x8["omega_meV"].iloc[0])
    assert math.isclose(val, 16.86, abs_tol=0.2)


def test_load_tadano2015_kappa_total_T() -> None:
    df = lit.load_tadano2015_kappa_total_T()
    assert df.empty  # template only; not digitized
    assert "Tadano2015" in df.attrs["provenance"]


# ---------------------------------------------------------------------------
# STO Verdi 2023
# ---------------------------------------------------------------------------
def test_load_verdi_sscha_omega_T() -> None:
    df = lit.load_verdi_sscha_omega_T()
    mlff = df[df["method"] == "MLFF_SSCHA"]
    # Shift_meV is constant 4.93 for the MLFF curve.
    assert set(mlff["shift_meV"].unique()) == {4.93}
    ins = df[df["method"] == "INS_Ref39"]
    assert set(ins["shift_meV"].unique()) == {0.0}


def test_load_verdi_dispersion() -> None:
    df = lit.load_verdi_dispersion()
    # Soft mode at Γ, 0 K, DFT_harmonic should be imaginary (-160 cm^-1).
    soft_g_harm = df[
        (df["q_label"] == "Gamma")
        & (df["T_K"] == 0)
        & (df["method"] == "DFT_harmonic")
        & (df["branch"] == "Soft_Optical")
    ]
    assert not soft_g_harm.empty
    assert math.isclose(float(soft_g_harm["omega_cm_inv"].iloc[0]), -160.0,
                        abs_tol=1.0)


def test_load_verdi_inverse_dielectric_T() -> None:
    df = lit.load_verdi_inverse_dielectric_T()
    mlff = df[df["method"] == "MLFF_SSCHA"]
    assert math.isclose(
        float(mlff[mlff["T_K"] == 10]["inverse_epsilon"].iloc[0]), 0.0006,
        abs_tol=5e-5,
    )
    assert math.isclose(
        float(mlff[mlff["T_K"] == 300]["inverse_epsilon"].iloc[0]), 0.0243,
        abs_tol=5e-4,
    )


# ---------------------------------------------------------------------------
# STO Yamada-Shirane 1969 (Table II q=0 Omega(T) pre-filled Turn 4B;
# gamma_T still an empty placeholder pending Fig. 6 digitization)
# ---------------------------------------------------------------------------
def test_load_yamada_shirane_omega_T() -> None:
    df = lit.load_yamada_shirane_omega_T()
    assert not df.empty
    assert "YamadaShirane1969" in df.attrs["provenance"]
    # Table II q=0 soft-mode T-series sentinels (SrTiO3, q_0_TO_soft).
    q0 = df[
        (df["material"] == "SrTiO3")
        & (df["condition_direction"] == "q_0_TO_soft")
    ]
    assert len(q0) == 19
    assert math.isclose(
        float(q0[q0["T_K"] == 293]["omega_meV"].iloc[0]), 10.8, abs_tol=1e-6
    )
    assert math.isclose(
        float(q0[q0["T_K"] == 4.5]["omega_meV"].iloc[0]), 1.9, abs_tol=1e-6
    )


def test_load_yamada_shirane_gamma_T() -> None:
    """STO Yamada-Shirane 1969 Fig. 6 intrinsic Γ(T) at q=0 (Wave-2 verified).

    Endpoints sentinel-checked against Fig. 6 (rendered image, 2026-05-31):
    Γ rises monotonically from 0.24 meV at 20 K to 0.84 meV at 250 K. This is
    the canonical low-T experimental STO soft-mode linewidth -- strongly
    *underdamped* (η = Γ/Ω ≈ 0.09 at 250 K with Ω ≈ 9.6 meV), in contrast to
    the overdamped cubic-BTO soft mode (HAS 1971, η ≈ 2.3).
    """
    df = lit.load_yamada_shirane_gamma_T()
    assert not df.empty
    assert math.isclose(
        float(df[df["T_K"] == 20]["gamma_meV"].iloc[0]), 0.24, abs_tol=1e-6
    )
    assert math.isclose(
        float(df[df["T_K"] == 250]["gamma_meV"].iloc[0]), 0.84, abs_tol=1e-6
    )
    g = df.sort_values("T_K")["gamma_meV"].to_numpy()
    assert (g[1:] >= g[:-1]).all()  # monotonically non-decreasing
    assert "YamadaShirane1969" in df.attrs["provenance"]


# ---------------------------------------------------------------------------
# BTO Harada-Axe-Shirane 1971 (cubic phase; Wave 2)
# ---------------------------------------------------------------------------
def test_load_has1971_cubic_dispersion_small_q_soft_TO() -> None:
    """Soft-TO sampled at small-to-moderate q (usable for the Λ² fit), cubic."""
    df = lit.load_harada_axe_shirane_1971_cubic_dispersion()
    assert (df["phase"] == "cubic").all()
    to = df[df["branch"] == "TO"].sort_values("q_inv_angstrom")
    # q = 0, 0.313, 0.478 Å⁻¹ (ζ = 0, 0.2, 0.305) -- NOT zone-boundary only.
    assert list(to["q_inv_angstrom"]) == [0.0, 0.313, 0.478]
    assert to["q_inv_angstrom"].max() < 0.5  # below the [100] ZB (0.784 Å⁻¹)
    assert "HaradaAxeShirane1971" in df.attrs["provenance"]


def test_load_has1971_mode_parameters_overdamped() -> None:
    """Verified overdamping: Γ_TO ≳ 2 Ω_TO, η = Γ/Ω ≈ 2.3 (cubic BTO)."""
    df = lit.load_harada_axe_shirane_1971_mode_parameters()
    assert (df["phase"] == "cubic").all()
    assert set(df["T_K"].unique()) == {423, 523, 723}  # 150/250/450 °C
    to = df[df["mode"] == "TO"]
    # Every TO row is overdamped: Γ ≥ 2Ω is the DHO overdamping threshold.
    assert (to["gamma_meV"] >= 2.0 * to["omega_meV"]).all()
    hi = to[to["T_K"] == 723].iloc[0]
    assert math.isclose(float(hi["omega_meV"]), 15.0, abs_tol=1e-6)
    assert math.isclose(float(hi["gamma_meV"]), 34.0, abs_tol=1e-6)
    assert math.isclose(float(hi["gamma_meV"]) / float(hi["omega_meV"]),
                        2.3, abs_tol=0.05)
    assert "HaradaAxeShirane1971" in df.attrs["provenance"]


def test_load_has1971_anisotropy_signed_angles() -> None:
    """Fig. 7 abscissa: signed angle from [100]=0°, within (-90, 90)."""
    df = lit.load_harada_axe_shirane_1971_anisotropy()
    assert df["angle_deg"].between(-90, 90).all()
    assert (df["angle_deg"] == 0.0).any()       # [100] reference present
    assert df["angle_deg"].min() < 0            # [111] side negative
    assert df["angle_deg"].max() > 0            # [110] side positive
    assert "HaradaAxeShirane1971" in df.attrs["provenance"]


def test_load_has1971_Lambda_sq_cubic() -> None:
    """Cubic-phase Λ² = A(100) = 972 ± 50 meV²Å² (HAS p.160 text)."""
    ns = lit.load_harada_axe_shirane_1971_Lambda_sq()
    assert isinstance(ns, SimpleNamespace)
    assert math.isclose(ns.central_meV2_A2, 972.0, abs_tol=1.0)
    assert math.isclose(ns.uncertainty_meV2_A2, 50.0, abs_tol=1.0)
    assert ns.T_K == 423  # 150 °C, cubic
    assert "HaradaAxeShirane1971" in ns.provenance


# ---------------------------------------------------------------------------
# STO Fauqué 2022 (Wave 1)
# ---------------------------------------------------------------------------
def test_load_fauque2022_omega_T_no_gamma() -> None:
    """Regression lock: Fauqué Fig. 1(d) reports NO DHO linewidth Γ.

    Both gamma columns must be entirely NaN. If a future re-extraction
    populates them, this test fails on purpose -- forcing a re-confirmation
    that the source really does report Γ before η switches from a bracket to
    a fit.
    """
    df = lit.load_fauque2022_omega_T()
    assert df["gamma_meV"].isna().all()
    assert df["gamma_uncertainty_meV"].isna().all()
    assert "Fauque2022" in df.attrs["provenance"]


def test_load_fauque2022_omega_T_fe_soft_mode_endpoints() -> None:
    """The FE soft mode is the Γ-point TO branch; sentinel low/high-T Ω."""
    df = lit.load_fauque2022_omega_T()
    to = df[(df["mode"] == "TO") & (df["method"] == "experimental")]
    # Γ-point (2,0,0) FE soft mode, 4-1400 K endpoints.
    assert math.isclose(
        float(to[to["T_K"] == 4]["omega_meV"].iloc[0]), 2.0, abs_tol=1e-6
    )
    assert math.isclose(
        float(to[to["T_K"] == 1400]["omega_meV"].iloc[0]), 17.3, abs_tol=1e-6
    )
    assert (to["q_point_rlu"] == "(2,0,0)").all()


def test_load_fauque2022_omega_T_nonfe_modes_tagged() -> None:
    """R-AFD and M modes are present but tagged distinctly from the FE TO."""
    df = lit.load_fauque2022_omega_T()
    assert set(df["mode"].unique()) == {"TO", "AFD", "M_mode"}
    # The non-FE branches sit at the R and M points, never (2,0,0).
    nonfe = df[df["mode"].isin(["AFD", "M_mode"])]
    assert (nonfe["q_label"].isin(["R", "M"])).all()


def test_load_fauque2022_tdep_dispersion_300K() -> None:
    df = lit.load_fauque2022_tdep_dispersion_300K()
    assert (df["T_K"] == 300).all()
    soft = df[(df["q_label"] == "Gamma") & (df["branch"] == "TO_Soft_Mode")]
    assert math.isclose(float(soft["omega_meV"].iloc[0]), 9.5, abs_tol=1e-6)
    assert "Fauque2022" in df.attrs["provenance"]


# ---------------------------------------------------------------------------
# STO Schmidt-Spaldin 2025 (Wave 1)
# ---------------------------------------------------------------------------
def test_load_schmidt2025_omega_T_isotopes_and_imaginary() -> None:
    df = lit.load_schmidt2025_omega_T()
    assert {"ST16O", "ST18O"} <= set(df["isotope"].unique())
    # ST18O SSCHA is imaginary (Ω² < 0) at low T -> negative Omega convention.
    st18 = df[(df["isotope"] == "ST18O") & (df["method"] == "SSCHA")]
    low_T = st18[st18["T_K"] == 0]
    assert float(low_T["Omega_meV"].iloc[0]) < 0
    assert float(low_T["Omega_sq_meV2"].iloc[0]) < 0
    # Sign of Omega and Omega_sq agree row-by-row (imaginary <-> negative).
    same_sign = (df["Omega_meV"] < 0) == (df["Omega_sq_meV2"] < 0)
    assert same_sign.all()
    assert "Schmidt-Spaldin2025" in df.attrs["provenance"]


def test_load_schmidt2025_energy_polarization_unit_and_isotopes() -> None:
    df = lit.load_schmidt2025_energy_polarization()
    # Provenance-only; column carries the µC/cm² unit suffix.
    assert "polarization_uC_per_cm2" in df.columns
    assert {"SrTi16O", "SrTi18O"} <= set(df["isotope"].unique())
    assert {"RPA", "SSCHA"} <= set(df["method"].unique())
    # µC/cm²-plausible magnitudes (well minimum ~0.22 µC/cm²), not C/m².
    assert df["polarization_uC_per_cm2"].abs().max() < 1.0
    assert "Schmidt-Spaldin2025" in df.attrs["provenance"]


# ---------------------------------------------------------------------------
# KTO Meier 2022 (Wave 1)
# ---------------------------------------------------------------------------
def test_load_meier2022_omega_T_plateau_and_endpoints() -> None:
    df = lit.load_meier2022_omega_T()
    assert {"Calculated", "Shirane_exp"} <= set(df["method"].unique())
    calc = df[df["method"] == "Calculated"]
    # Quantum-paraelectric low-T plateau: Ω(0 K) == Ω(10 K).
    assert math.isclose(
        float(calc[calc["T_K"] == 0]["Omega_meV"].iloc[0]),
        float(calc[calc["T_K"] == 10]["Omega_meV"].iloc[0]),
        abs_tol=1e-6,
    )
    assert math.isclose(
        float(calc[calc["T_K"] == 300]["Omega_meV"].iloc[0]), 10.86, abs_tol=1e-6
    )
    assert "Meier2022" in df.attrs["provenance"]


def test_load_meier2022_inverse_dielectric_T() -> None:
    df = lit.load_meier2022_inverse_dielectric_T()
    assert {"Calculated", "Fujishita_exp", "Barrett_fit"} <= set(df["method"].unique())
    assert "Meier2022" in df.attrs["provenance"]


def test_load_meier2022_dispersion_T_gamma_only_soft_mode() -> None:
    """Meier Fig. 2 soft mode is sampled at Γ only -> no small-q Λ² fit."""
    df = lit.load_meier2022_dispersion_T()
    soft = df[df["branch"] == "Soft_Optical"]
    # Soft branch lives at Γ exclusively; Λ²_KTO is not fittable from here.
    assert (soft["q_label"] == "Gamma").all()
    # Harmonic small-displacement reference is imaginary at Γ.
    harm = soft[soft["method"] == "harmonic_small_disp"]
    assert float(harm["omega_meV"].iloc[0]) < 0
    sscha = soft[soft["method"] == "Calculated_SSCHA"]
    assert math.isclose(
        float(sscha[sscha["T_K"] == 0]["omega_meV"].iloc[0]), 3.31, abs_tol=1e-6
    )
    assert "Meier2022" in df.attrs["provenance"]


# ---------------------------------------------------------------------------
# KTO Ranalli 2023 (2nd KTaO3 source; KTaO3 wave)
# ---------------------------------------------------------------------------
def test_load_ranalli2023_omega_T_shift_and_plateau() -> None:
    """Verdi-family down-shift on the computed curve; QPE low-T plateau."""
    df = lit.load_ranalli2023_omega_T()
    assert {"MLFF+SSCHA", "expt", "fit_paraelectric"} <= set(df["method"].unique())
    ml = df[df["method"] == "MLFF+SSCHA"]
    # Computed curve shifted DOWN by 34.71 cm⁻¹ ≈ 4.30 meV (negative = down).
    assert (ml["shift_meV"] == -4.30).all()
    # Experimental + √-fit rows are unshifted.
    assert (df[df["method"] != "MLFF+SSCHA"]["shift_meV"] == 0.0).all()
    # Quantum-paraelectric low-T plateau: Ω(5 K) == Ω(20 K).
    assert math.isclose(
        float(ml[ml["T_K"] == 5]["Omega_meV"].iloc[0]),
        float(ml[ml["T_K"] == 20]["Omega_meV"].iloc[0]), abs_tol=1e-6,
    )
    # Endpoints (shifted MLFF and exp at 300 K).
    assert math.isclose(float(ml[ml["T_K"] == 300]["Omega_meV"].iloc[0]),
                        9.05, abs_tol=0.05)
    ex = df[df["method"] == "expt"]
    assert math.isclose(float(ex[ex["T_K"] == 300]["Omega_meV"].iloc[0]),
                        10.29, abs_tol=0.05)
    assert "Ranalli2023" in df.attrs["provenance"]


def test_ranalli_vs_meier_two_source_consistency() -> None:
    """Two-source KTaO3 cross-check: experiments agree; computed curves spread.

    The INS experimental anchors agree across sources (Ranalli expt vs Meier
    Shirane_exp at 300 K within ~1 meV). The computed curves carry a larger
    ~1-2 meV spread from differing shift conventions / methods -- so we assert
    the experimental agreement but only sanity-bound (do NOT pin) the computed
    gap, flagging it as a real inter-source spread rather than a digitization
    error.
    """
    ran = lit.load_ranalli2023_omega_T()
    mei = lit.load_meier2022_omega_T()
    ran_exp_300 = float(ran[(ran["method"] == "expt") & (ran["T_K"] == 300)]["Omega_meV"].iloc[0])
    mei_exp_300 = float(mei[(mei["method"] == "Shirane_exp") & (mei["T_K"] == 300)]["Omega_meV"].iloc[0])
    assert abs(ran_exp_300 - mei_exp_300) < 1.0  # experiments agree (~0.4 meV)
    # Low-T plateaus both lie in the QPE ~2-4 meV band (same ballpark).
    ran_plat = float(ran[(ran["method"] == "MLFF+SSCHA") & (ran["T_K"] == 5)]["Omega_meV"].iloc[0])
    mei_plat = float(mei[(mei["method"] == "Calculated") & (mei["T_K"] == 0)]["Omega_meV"].iloc[0])
    assert 2.0 <= ran_plat <= 4.0 and 2.0 <= mei_plat <= 4.0
    # Computed 300 K spread is real but bounded (< 3 meV); not asserted equal.
    ran_calc_300 = float(ran[(ran["method"] == "MLFF+SSCHA") & (ran["T_K"] == 300)]["Omega_meV"].iloc[0])
    mei_calc_300 = float(mei[(mei["method"] == "Calculated") & (mei["T_K"] == 300)]["Omega_meV"].iloc[0])
    assert abs(ran_calc_300 - mei_calc_300) < 3.0  # ~1.8 meV shift/model spread


def test_load_ranalli2023_inverse_dielectric_T() -> None:
    df = lit.load_ranalli2023_inverse_dielectric_T()
    assert {"MLFF+SSCHA", "expt", "Curie-Weiss_fit"} <= set(df["method"].unique())
    assert "Ranalli2023" in df.attrs["provenance"]


# ---------------------------------------------------------------------------
# Vogt 1995 (experimental hyper-Raman Ω(T) + γ(T) for both QPEs; Vogt wave)
# ---------------------------------------------------------------------------
# γ-convention adopted for tests: plotted Vogt γ = FWHM -> Γ_HWHM = γ/2 (see
# the literature_io module note). η = (γ/2)/Ω.
ETA_GAMMA_TO_HWHM = 0.5  # γ (FWHM) -> Γ_HWHM factor


def test_load_vogt1995_kto_experimental_qzero() -> None:
    o = lit.load_vogt1995_kto_omega_T()
    g = lit.load_vogt1995_kto_gamma_T()
    assert {"exp_HyperRaman_Vogt1995", "exp_HyperRaman_Vogt1988"} <= set(o["method"].unique())
    assert (o["shift_meV"] == 0.0).all()                 # experimental, no shift
    assert list(g.columns) == ["T_K", "gamma_meV", "gamma_uncertainty_meV"]  # single-phase 3-col
    assert "Vogt1995" in o.attrs["provenance"] and "Vogt1995" in g.attrs["provenance"]


def test_load_vogt1995_sto_phase_split_and_cubic_filter() -> None:
    o = lit.load_vogt1995_sto_omega_T()
    g = lit.load_vogt1995_sto_gamma_T()
    assert (o["shift_meV"] == 0.0).all()
    # Two-variant gamma: STO carries phase_mode (4-col).
    assert "phase_mode" in g.columns
    assert {"tetragonal_A2u", "tetragonal_Eu", "cubic_F1u"} <= set(o["phase_mode"].unique())
    # Kernel-relevant cubic single mode is filterable and >= 110 K only.
    cub = o[o["phase_mode"] == "cubic_F1u"]
    assert cub["T_K"].min() >= 110   # sub-105 K merged points re-tagged transition_merged
    # dense extraction (non-round T): F1u top ~299 K -> ~10.96 meV (η denominator)
    assert 10.0 < float(cub.sort_values("T_K")["Omega_meV"].iloc[-1]) < 11.5


def test_vogt_eta_cluster_underdamped() -> None:
    """Both QPEs cluster as moderately underdamped (η<0.5) on the FWHM basis."""
    so = lit.load_vogt1995_sto_omega_T()
    sg = lit.load_vogt1995_sto_gamma_T()
    ko = lit.load_vogt1995_kto_omega_T()
    kg = lit.load_vogt1995_kto_gamma_T()
    # Both are dense WebPlotDigitizer extractions (non-round T) -> use high-T
    # endpoints, not hardcoded T. STO from cubic_F1u; KTaO3 Ω from the 1995
    # curve (excludes the 1988 Table III point).
    so_f1u = so[so["phase_mode"] == "cubic_F1u"].sort_values("T_K")
    sg_f1u = sg[sg["phase_mode"] == "cubic_F1u"].sort_values("T_K")
    sO = float(so_f1u["Omega_meV"].iloc[-1])           # ~10.96 meV @ 299 K
    sG = float(sg_f1u["gamma_meV"].iloc[-1])           # ~2.40 meV @ 298 K (dense)
    ko95 = ko[ko["method"] == "exp_HyperRaman_Vogt1995"].sort_values("T_K")
    kO = float(ko95["Omega_meV"].iloc[-1])             # ~10.15 meV @ 294 K
    kG = float(kg.sort_values("T_K")["gamma_meV"].iloc[-1])  # ~2.57 meV @ 296 K
    eta_sto = (ETA_GAMMA_TO_HWHM * sG) / sO            # ~0.11 (dense STO gamma)
    eta_kto = (ETA_GAMMA_TO_HWHM * kG) / kO            # ~0.127
    assert eta_sto < 0.5 and eta_kto < 0.5            # both underdamped
    assert abs(eta_sto - eta_kto) < 0.05              # cluster (≈0.11–0.13)
    # Sanity vs the overdamped end (BTO HAS η≈2.3): QPEs are far below it.
    assert eta_sto < 0.5 < 2.3


def test_ys_vs_vogt_sto_damping_documented_spread() -> None:
    """YS and Vogt STO damping DISAGREE by design -- documented content, NOT to
    be reconciled to equality (Amendment-2 §2/§15). With the dense STO γ the raw
    spread is ~2.4× at 250 K (was ~3.3× on the round γ); the qualitative point
    (Vogt > YS, q-res + central-peak physics) is unchanged."""
    ys = lit.load_yamada_shirane_gamma_T()           # HWHM, use directly
    sg = lit.load_vogt1995_sto_gamma_T()
    ys_250 = float(ys[ys["T_K"] == 250]["gamma_meV"].iloc[0])
    # dense Vogt has no exactly-250 K row -> nearest cubic_F1u point (~249 K).
    sgf = sg[sg["phase_mode"] == "cubic_F1u"].copy()
    vg_250 = float(sgf.iloc[(sgf["T_K"] - 250).abs().argmin()]["gamma_meV"])  # ~2.05 meV
    raw_ratio = vg_250 / ys_250                       # ~2.4 (dense; vs YS HWHM)
    assert raw_ratio > 2.0                            # genuinely differ (documented)
    # After FWHM->HWHM /2, Vogt still exceeds YS (residual > 1).
    residual = (ETA_GAMMA_TO_HWHM * vg_250) / ys_250  # ~1.2
    assert residual > 1.0


# ---------------------------------------------------------------------------
# Umbrella
# ---------------------------------------------------------------------------
def test_load_all_keys_and_count() -> None:
    data = lit.load_all()
    assert len(data) == 38
    # Sample keys covering the {material}.{source}.{quantity} convention.
    assert "BTO.ehsan2021.anchors" in data
    assert "BTO.shirane1970.Lambda_sq" in data
    assert "STO.tadano2015.scph_omega_T" in data
    assert "STO.verdi2023.sscha_omega_T" in data
    assert "STO.yamada_shirane1969.gamma_T" in data
    # Wave 1 arXiv sources.
    assert "STO.fauque2022.omega_T" in data
    assert "STO.schmidt2025.omega_T" in data
    assert "KTO.meier2022.omega_T" in data
    assert "KTO.ranalli2023.omega_T" in data
    assert "KTO.vogt1995.gamma_T" in data
    assert "STO.vogt1995.gamma_T" in data
    # Wave 2 cubic-BTO library source.
    assert "BTO.harada_axe_shirane1971.cubic_dispersion" in data
    assert "BTO.harada_axe_shirane1971.mode_parameters" in data
    assert isinstance(data["BTO.harada_axe_shirane1971.Lambda_sq"], SimpleNamespace)
    # Anchors and Λ² return SimpleNamespace; everything else is DataFrame.
    assert isinstance(data["BTO.ehsan2021.anchors"], SimpleNamespace)
    assert isinstance(data["BTO.shirane1970.Lambda_sq"], SimpleNamespace)
    assert isinstance(data["BTO.tadano2022.kappa_total_T"], pd.DataFrame)
    assert isinstance(data["STO.verdi2023.dispersion"], pd.DataFrame)


@pytest.mark.parametrize("rel_path", sorted(lit._SCHEMAS))
def test_schema_dict_matches_test_inputs(rel_path: str) -> None:
    """Cross-check that ``literature_io._SCHEMAS`` matches the canonical
    schemas in ``tests.test_literature_inputs.EXPECTED_SCHEMAS``."""
    from tests.test_literature_inputs import EXPECTED_SCHEMAS

    assert lit._SCHEMAS[rel_path] == EXPECTED_SCHEMAS[rel_path]
