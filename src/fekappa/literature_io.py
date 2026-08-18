"""Literature CSV loaders for the Turn 4B reduced-kernel pipeline.

One loader per CSV under ``data/literature/{BTO,STO}/``. Each loader:

* Reads the file with :func:`pandas.read_csv`.
* Asserts the expected column list, which encodes the unit-suffix conventions
  (e.g. ``omega_meV``, ``omega_cm_inv``, ``Lambda_sq_meV2_A2``, ``shift_meV``).
* Returns a typed object: :class:`pandas.DataFrame` for tabular data,
  :class:`types.SimpleNamespace` for the two scalar/anchor files
  (Ehsan anchors, Shirane Λ²).
* Attaches a provenance string ``"<CitationKey> :: <figure-or-table>"`` -- as
  ``DataFrame.attrs["provenance"]`` for tables, or as the ``provenance``
  attribute for namespaces.

Loaders follow ``load_<source>_<quantity>()`` naming; the umbrella
:func:`load_all` returns everything keyed by ``"{material}.{source}.{quantity}"``.

Schema reference: ``data/literature/digitization_plan.md`` (kept in sync with
``tests/test_literature_inputs.py::EXPECTED_SCHEMAS``).
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd

LIT_ROOT: Path = Path(__file__).resolve().parents[2] / "data" / "literature"

# ---------------------------------------------------------------------------
# Schemas (mirror tests/test_literature_inputs.EXPECTED_SCHEMAS)
# ---------------------------------------------------------------------------
_SCHEMAS: dict[str, list[str]] = {
    # BTO Ehsan 2021
    "BTO/ehsan_2021/anchors.csv":
        ["T_K", "Omega_meV", "Omega_sq_meV2", "note"],
    "BTO/ehsan_2021/scp_omega_T.csv":
        ["T_K", "Omega_meV", "Omega_sq_meV2"],
    "BTO/ehsan_2021/free_energy_T.csv":
        ["T_K", "F_cubic_meV", "F_tetra_meV", "source_figure"],
    "BTO/ehsan_2021/structural_properties.csv":
        ["phase", "a_angstrom", "c_angstrom", "c_over_a", "V_cell_A3",
         "atom", "Wyckoff_x", "Wyckoff_y", "Wyckoff_z", "source"],
    "BTO/ehsan_2021/cubic_dispersion.csv":
        ["T_K", "condition", "branch_label", "q_label", "q_inv_angstrom",
         "omega_cm_inv"],
    "BTO/ehsan_2021/born_charges_epsilon_infty.csv":
        ["phase", "atom", "Z_xx", "Z_yy", "Z_zz",
         "epsilon_inf_xx", "epsilon_inf_yy", "epsilon_inf_zz", "source"],
    # BTO Tadano 2022
    "BTO/tadano_2022/free_energy_T.csv":
        ["T_K", "F_cubic_meV", "F_tetra_meV", "F_ortho_meV", "F_rhombo_meV"],
    "BTO/tadano_2022/scp_dispersion_plus1pct.csv":
        ["T_K", "branch_label", "q_label", "omega_cm_inv"],
    "BTO/tadano_2022/transition_temperatures.csv":
        ["transition", "T_c_calc_K", "T_c_exp_K"],
    "BTO/tadano_2022/kappa_total_T.csv":
        ["T_K", "kappa_W_per_m_K", "kappa_uncertainty", "method", "source"],
    # BTO Shirane 1970
    "BTO/shirane_1970/Lambda_sq.csv":
        ["Lambda_sq_meV2_A2", "Lambda_sq_uncertainty_meV2_A2", "T_K",
         "reference_q_range"],
    "BTO/shirane_1970/cubic_dispersion.csv":
        ["q_label", "q_inv_angstrom", "omega_meV", "omega_uncertainty_meV",
         "branch", "polarization_e", "wavevector_q", "phase", "T_K",
         "source_figure", "note"],
    # BTO Harada-Axe-Shirane 1971 (cubic-phase; Wave 2)
    "BTO/harada_axe_shirane_1971/cubic_dispersion.csv":
        ["q_label", "q_inv_angstrom", "omega_meV", "omega_uncertainty_meV",
         "branch", "polarization_e", "wavevector_q", "phase", "T_K",
         "source_figure", "note"],
    "BTO/harada_axe_shirane_1971/mode_parameters.csv":
        ["q_label", "q_inv_angstrom", "T_K", "phase", "mode",
         "omega_meV", "omega_uncertainty_meV",
         "gamma_meV", "gamma_uncertainty_meV",
         "structure_factor", "coupling_param", "source"],
    "BTO/harada_axe_shirane_1971/anisotropy.csv":
        ["angle_deg", "q_inv_angstrom", "omega_meV", "omega_uncertainty_meV",
         "T_K", "phase", "source_figure", "note"],
    "BTO/harada_axe_shirane_1971/Lambda_sq.csv":
        ["Lambda_sq_meV2_A2", "Lambda_sq_uncertainty_meV2_A2", "T_K",
         "reference_q_range"],
    # STO Tadano 2015
    "STO/tadano_2015/scph_omega_T.csv":
        ["T_K", "mode", "method", "Omega_meV", "Omega_sq_meV2", "source_figure"],
    "STO/tadano_2015/scph_dispersion.csv":
        ["T_K", "q_label", "q_inv_angstrom", "omega_meV", "branch",
         "method", "source_figure"],
    "STO/tadano_2015/scph_convergence_300K.csv":
        ["T_K", "q_label", "q_inv_angstrom", "omega_meV", "branch",
         "method", "source_figure"],
    "STO/tadano_2015/kappa_total_T.csv":
        ["T_K", "kappa_W_per_m_K", "kappa_uncertainty", "method", "source"],
    "STO/tachibana_2008/kappa_total_T.csv":
        ["T_K", "kappa_W_per_m_K", "kappa_uncertainty", "method", "source"],
    # STO Verdi 2023
    "STO/verdi_2023/sscha_omega_T.csv":
        ["T_K", "Omega_meV", "Omega_sq_meV2", "method", "shift_meV",
         "source_figure"],
    "STO/verdi_2023/dispersion.csv":
        ["T_K", "q_label", "q_inv_angstrom", "omega_cm_inv", "branch",
         "method", "supercell", "source_figure"],
    "STO/verdi_2023/inverse_dielectric_T.csv":
        ["T_K", "inverse_epsilon", "method", "source_figure"],
    # STO Yamada-Shirane 1969 (and Shirane 1970 Fig. 1 STO-panel fallback)
    "STO/yamada_shirane_1969/omega_T.csv":
        ["material", "condition_direction", "q_inv_angstrom", "T_K",
         "omega_meV", "omega_uncertainty_meV"],
    "STO/yamada_shirane_1969/gamma_T.csv":
        ["T_K", "gamma_meV", "gamma_uncertainty_meV"],
    # STO Fauqué 2022 (Wave 1)
    "STO/fauque_2022/dho_fits_T.csv":
        ["T_K", "q_label", "q_point_rlu", "mode",
         "omega_meV", "omega_uncertainty_meV",
         "gamma_meV", "gamma_uncertainty_meV",
         "method", "source_figure"],
    "STO/fauque_2022/tdep_dispersion_300K.csv":
        ["T_K", "q_label", "q_inv_angstrom", "omega_meV", "branch",
         "method", "source_figure"],
    # STO Schmidt-Spaldin 2025 (Wave 1)
    "STO/schmidt_2025/sscha_omega_T.csv":
        ["T_K", "isotope", "Omega_meV", "Omega_sq_meV2", "method",
         "source_figure"],
    "STO/schmidt_2025/energy_polarization.csv":
        ["T_K", "polarization_uC_per_cm2", "E_meV_per_fu", "isotope",
         "method", "source_figure"],
    # KTO Meier 2022 (Wave 1; third material)
    "KTO/meier_2022/omega_T.csv":
        ["T_K", "Omega_meV", "Omega_sq_meV2", "method", "source_figure"],
    "KTO/meier_2022/inverse_dielectric_T.csv":
        ["T_K", "inverse_epsilon", "method", "source_figure"],
    "KTO/meier_2022/dispersion_T.csv":
        ["T_K", "q_label", "q_inv_angstrom", "omega_meV", "branch",
         "method", "source_figure"],
    # KTO Ranalli 2023 (2nd KTaO3 source; schema identical to Verdi sscha_omega_T)
    "KTO/ranalli_2023/omega_T.csv":
        ["T_K", "Omega_meV", "Omega_sq_meV2", "method", "shift_meV",
         "source_figure"],
    "KTO/ranalli_2023/inverse_dielectric_T.csv":
        ["T_K", "inverse_epsilon", "method", "source_figure"],
    # KTO + STO Vogt 1995 (experimental Ω(T) + γ(T) hyper-Raman; Vogt wave).
    # omega uses Omega_uncertainty_meV (not Omega_sq); shift_meV=0 (exp).
    # gamma schema is intentionally TWO-VARIANT: STO-Vogt carries phase_mode
    # (Fig. 6 straddles the 105 K AFD transition), KTaO3-Vogt does not.
    "KTO/vogt_1995/omega_T.csv":
        ["T_K", "Omega_meV", "Omega_uncertainty_meV", "method", "shift_meV",
         "source_figure"],
    "KTO/vogt_1995/gamma_T.csv":
        ["T_K", "gamma_meV", "gamma_uncertainty_meV"],
    "STO/vogt_1995/omega_T.csv":
        ["T_K", "Omega_meV", "Omega_uncertainty_meV", "method", "shift_meV",
         "source_figure", "phase_mode"],
    "STO/vogt_1995/gamma_T.csv":
        ["T_K", "gamma_meV", "gamma_uncertainty_meV", "phase_mode"],
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _provenance(citation: str, source: str) -> str:
    """Build a ``'CitationKey :: source-locator'`` provenance string."""
    return f"{citation} :: {source}"


def _load_csv(rel_path: str, provenance: str) -> pd.DataFrame:
    """Read a CSV under :data:`LIT_ROOT` and validate its column schema.

    Attaches ``df.attrs["provenance"]`` and ``df.attrs["source_path"]``.
    """
    path = LIT_ROOT / rel_path
    expected = _SCHEMAS[rel_path]
    df = pd.read_csv(path)
    cols = list(df.columns)
    if cols != expected:
        raise ValueError(
            f"{rel_path}: column schema mismatch. "
            f"Expected {expected}, got {cols}."
        )
    df.attrs["provenance"] = provenance
    df.attrs["source_path"] = str(path)
    return df


# ===========================================================================
# BTO -- Ehsan 2021
# ===========================================================================
def load_ehsan_anchors() -> SimpleNamespace:
    """BTO Ehsan 2021 SCP anchors at 0 K (bare imaginary) and 500 K (thermalized).

    Returns:
        SimpleNamespace with attributes ``T0_Omega_sq_meV2``,
        ``T500_Omega_meV``, ``T500_Omega_sq_meV2``, ``provenance``.
    """
    df = _load_csv("BTO/ehsan_2021/anchors.csv",
                   _provenance("Ehsan2021", "Sec. IV C anchors"))
    by_T = {int(float(r["T_K"])): r for _, r in df.iterrows()}
    a0, a500 = by_T[0], by_T[500]
    return SimpleNamespace(
        T0_Omega_sq_meV2=float(a0["Omega_sq_meV2"]),
        T500_Omega_meV=float(a500["Omega_meV"]),
        T500_Omega_sq_meV2=float(a500["Omega_sq_meV2"]),
        provenance=df.attrs["provenance"],
    )


def load_ehsan_scp_omega_T() -> pd.DataFrame:
    """BTO Ehsan 2021 SCP Omega(T) curve placeholder.

    Empty by design: the Ehsan paper contains no T-dependent Omega(T) curve
    (Figs. 4 and 5 are U/TS and free energies respectively). Loader exists for
    schema completeness; callers should check ``df.empty``.
    """
    return _load_csv("BTO/ehsan_2021/scp_omega_T.csv",
                     _provenance("Ehsan2021", "n/a (no Omega(T) curve in paper)"))


def load_ehsan_free_energy_T() -> pd.DataFrame:
    """BTO Ehsan 2021 Fig. 5 free energies F_cubic(T) and F_tetra(T)."""
    return _load_csv("BTO/ehsan_2021/free_energy_T.csv",
                     _provenance("Ehsan2021", "Fig. 5"))


def load_ehsan_structural() -> pd.DataFrame:
    """BTO Ehsan 2021 Table I cubic + tetragonal structural properties."""
    return _load_csv("BTO/ehsan_2021/structural_properties.csv",
                     _provenance("Ehsan2021", "Table I"))


def load_ehsan_cubic_dispersion() -> pd.DataFrame:
    """BTO Ehsan 2021 Fig. 2 cubic phonon band dispersion (0 K + 500 K).

    Notes:
        ``q_inv_angstrom`` carries Å⁻¹ values after the Turn 4B pre-flight
        unit fix (converted from r.l.u. via 2π/4.013, Ehsan Table I cubic a).
        Zone-boundary points (X, M, R) are *not* used for Λ² refits in
        NB 06 -- they lie outside the small-q regime where
        Ω²(q) = Ω₀² + Λ² q² is valid.
    """
    return _load_csv("BTO/ehsan_2021/cubic_dispersion.csv",
                     _provenance("Ehsan2021", "Fig. 2"))


def load_ehsan_born_charges() -> pd.DataFrame:
    """BTO Ehsan 2021 Table II Born effective charges Z* and ε_∞."""
    return _load_csv("BTO/ehsan_2021/born_charges_epsilon_infty.csv",
                     _provenance("Ehsan2021", "Table II"))


# ===========================================================================
# BTO -- Masuki/Tadano 2022
# ===========================================================================
def load_tadano2022_free_energy_T() -> pd.DataFrame:
    """BTO Tadano 2022 Fig. 4 free energies (cubic, tetra, ortho, rhombo)."""
    return _load_csv("BTO/tadano_2022/free_energy_T.csv",
                     _provenance("Masuki-Tadano2022", "Fig. 4"))


def load_tadano2022_scp_dispersion_plus1pct() -> pd.DataFrame:
    """BTO Tadano 2022 Fig. 15 SCP dispersion. **At +1% lattice constant.**"""
    return _load_csv("BTO/tadano_2022/scp_dispersion_plus1pct.csv",
                     _provenance("Masuki-Tadano2022", "Fig. 15 (+1% lattice)"))


def load_tadano2022_transition_temperatures() -> pd.DataFrame:
    """BTO Tadano 2022 Table I cubic-tetra / tetra-ortho / ortho-rhombo T_c."""
    return _load_csv("BTO/tadano_2022/transition_temperatures.csv",
                     _provenance("Masuki-Tadano2022", "Table I"))


def load_tadano2022_kappa_total_T() -> pd.DataFrame:
    """BTO Tadano 2022 Fig. 7 total κ(T) + Muta/Popuri experimental overlays.

    **Scale-overlay only.** Turn 4 explicitly does not validate against total κ.
    """
    return _load_csv("BTO/tadano_2022/kappa_total_T.csv",
                     _provenance("Masuki-Tadano2022", "Fig. 7"))


# ===========================================================================
# BTO -- Shirane 1970
# ===========================================================================
def load_shirane_Lambda_sq() -> SimpleNamespace:
    """BTO Shirane 1970 Λ² central lock (prior-chat fit).

    Returns:
        SimpleNamespace with ``central_meV2_A2``, ``uncertainty_meV2_A2``
        (may be ``nan``), ``T_K``, ``reference_q_range``, ``provenance``.
    """
    df = _load_csv("BTO/shirane_1970/Lambda_sq.csv",
                   _provenance("Shirane1970", "Δ-branch prior fit"))
    row = df.iloc[0]
    return SimpleNamespace(
        central_meV2_A2=float(row["Lambda_sq_meV2_A2"]),
        uncertainty_meV2_A2=(
            float(row["Lambda_sq_uncertainty_meV2_A2"])
            if pd.notna(row["Lambda_sq_uncertainty_meV2_A2"])
            else float("nan")
        ),
        T_K=int(row["T_K"]),
        reference_q_range=str(row["reference_q_range"]),
        provenance=df.attrs["provenance"],
    )


def load_shirane_cubic_dispersion() -> pd.DataFrame:
    """BTO Shirane 1970 Figs. 3 + 6 soft-TO dispersion at 295 K.

    **Tetragonal-phase** measurements (295 K, below T_c ≈ 403 K). The cubic-
    phase kernel inherits this approximation; the ±20% Λ² envelope absorbs
    the phase-shift uncertainty.
    """
    return _load_csv("BTO/shirane_1970/cubic_dispersion.csv",
                     _provenance("Shirane1970", "Figs. 3 + 6 (Σ, Δ branches)"))


# ===========================================================================
# BTO -- Harada-Axe-Shirane 1971 (cubic phase; Wave 2)
# ===========================================================================
def load_harada_axe_shirane_1971_cubic_dispersion() -> pd.DataFrame:
    """BTO Harada-Axe-Shirane 1971 Figs. 5 + 8 soft-TO + TA dispersion (cubic).

    **The first genuine cubic-phase (T > T_c ≈ 130 °C) BTO dispersion in the
    project** -- it removes the tetragonal-phase caveat that Shirane 1970
    (295 K) carries. Measured at 150 °C (423 K), just above T_c.

    The soft-TO branch is sampled at q = 0, 0.313, 0.478 Å⁻¹ (ζ = 0, 0.2,
    0.305 along [ζ00]); the paper fits these to the small-q form
    Ω²(q) = Ω₀² + A q² (Fig. 5 dotted line), giving the cubic Λ² in
    :func:`load_harada_axe_shirane_1971_Lambda_sq`. The TO mode is heavily
    overdamped (see :func:`load_harada_axe_shirane_1971_mode_parameters`), so
    these Ω values are the *quasiharmonic* (damping-projected) frequencies,
    not resolvable peak positions.
    """
    return _load_csv("BTO/harada_axe_shirane_1971/cubic_dispersion.csv",
                     _provenance("HaradaAxeShirane1971", "Figs. 5 + 8"))


def load_harada_axe_shirane_1971_mode_parameters() -> pd.DataFrame:
    """BTO Harada-Axe-Shirane 1971 Table I quasiharmonic Ω, damping Γ, coupling.

    Long format, ``mode`` ∈ {``TO``, ``TA``}, cubic phase at 423/523/723 K
    (150/250/450 °C; the paper quotes temperatures in °C). ``omega_meV`` and
    ``gamma_meV`` are the renormalized (quasiharmonic) Ω_j and the DHO
    mode-damping Γ_j of the coupled-mode response, Eq. (2).

    **Overdamping (paper finding, verified against Table I + p.159 text):**
    the soft TO satisfies Γ_TO ≳ 2 Ω_TO at every (q, T); the paper reports the
    empirical η = Γ_TO/Ω_TO ≈ 2.3 (neutron, finite q) and cites Γ/Ω = 4-6 at
    q = 0 (optical) and 9 at q = 0.08 a*. Γ_TO ≈ 34 meV vs Ω_TO ≈ 15 meV at
    q = 0.2 a* is therefore a *real* cubic-BTO linewidth (the DHO damping
    constant, not a FWHM/HWHM artifact), and direct evidence that cubic BTO
    is overdamped -- the underdamped soft-mode framing is STO-specific. This
    is surfaced in the NB 06 BTO methodological note.

    Auxiliary Table-I columns (context, not kernel inputs):
    ``structure_factor`` = the F_TO/F_TA inelastic structure-factor ratio
    (e.g. -2.34); ``coupling_param`` = ℏ²λ² mode-coupling parameter (meV²,
    e.g. 61.0).
    """
    return _load_csv("BTO/harada_axe_shirane_1971/mode_parameters.csv",
                     _provenance("HaradaAxeShirane1971", "Table I"))


def load_harada_axe_shirane_1971_anisotropy() -> pd.DataFrame:
    """BTO Harada-Axe-Shirane 1971 Fig. 7 soft-TO anisotropy near q = 0.

    The abscissa ``angle_deg`` is the *signed* angular orientation of q along
    arcs connecting the principal symmetry axes ([100] = 0°, [111] ≈ -54.7°,
    [110] = +45°). The soft-TO curvature A(q̂) is strongly anisotropic
    (A(100) = 972 vs A(110)₁ = A(111) = 4750 meV²Å²; Eq. 5, p.160), which
    drives the (100)-plane diffuse-scattering sheets discussed by Comès/Hüller.
    """
    return _load_csv("BTO/harada_axe_shirane_1971/anisotropy.csv",
                     _provenance("HaradaAxeShirane1971", "Fig. 7"))


def load_harada_axe_shirane_1971_Lambda_sq() -> SimpleNamespace:
    """BTO Harada-Axe-Shirane 1971 cubic-phase Λ² (soft-TO curvature, [100]).

    Returns:
        SimpleNamespace with ``central_meV2_A2``, ``uncertainty_meV2_A2``,
        ``T_K``, ``reference_q_range``, ``provenance``.

    The value is the paper's self-consistent A(100) = (972 ± 50) meV²Å²
    (p.160 text), the [ζ00] soft-TO curvature in Ω²(q) = Ω₀² + A q² with
    Ω₀ = (3 ± 1) meV at 150 °C. **This is the cubic-phase counterpart to the
    tetragonal Shirane 1970 Λ² = 5100.** The two are NOT in conflict: the
    soft TO is ~5× stiffer along [110]/[111] (A = 4750 ≈ Shirane's steep
    direction) than along [100] (A = 972), so the apparent gap is the
    soft-mode anisotropy, not a discrepancy. (A 3-point Fig. 5 digitization
    over ζ ≤ 0.3 gave 859; the text value 972 ± 50 is authoritative.)
    """
    df = _load_csv("BTO/harada_axe_shirane_1971/Lambda_sq.csv",
                   _provenance("HaradaAxeShirane1971", "Fig. 5 / p.160 A(100)"))
    row = df.iloc[0]
    return SimpleNamespace(
        central_meV2_A2=float(row["Lambda_sq_meV2_A2"]),
        uncertainty_meV2_A2=(
            float(row["Lambda_sq_uncertainty_meV2_A2"])
            if pd.notna(row["Lambda_sq_uncertainty_meV2_A2"])
            else float("nan")
        ),
        T_K=int(row["T_K"]),
        reference_q_range=str(row["reference_q_range"]),
        provenance=df.attrs["provenance"],
    )


# ===========================================================================
# STO -- Tadano 2015
# ===========================================================================
def load_tadano2015_scph_omega_T() -> pd.DataFrame:
    """STO Tadano 2015 Fig. 5 Ω(T): TO1 + TO2 × {SCPH±PM, TDEP, exp_Servoin, exp_YamadaShirane}."""
    return _load_csv("STO/tadano_2015/scph_omega_T.csv",
                     _provenance("Tadano2015", "Fig. 5"))


def load_tadano2015_scph_dispersion() -> pd.DataFrame:
    """STO Tadano 2015 Fig. 4 cubic phonon dispersion at 300 K (SCPH + harm + exp)."""
    return _load_csv("STO/tadano_2015/scph_dispersion.csv",
                     _provenance("Tadano2015", "Fig. 4"))


def load_tadano2015_scph_convergence_300K() -> pd.DataFrame:
    """STO Tadano 2015 Table I SCPH convergence at 300 K (FE, AFD, M₃ vs q₁ grid)."""
    return _load_csv("STO/tadano_2015/scph_convergence_300K.csv",
                     _provenance("Tadano2015", "Table I"))


def load_tadano2015_kappa_total_T() -> pd.DataFrame:
    """STO Tadano 2015 κ(T): the digitized 300 K value (8.7 W/mK, 12×12×12 mesh).

    Scale-overlay / benchmark-marker use only (single digitized point; the
    full curve remains undigitized).
    """
    return _load_csv("STO/tadano_2015/kappa_total_T.csv",
                     _provenance("Tadano2015", "κ(T) figure, 300 K point"))


def load_tachibana2008_kappa_total_T() -> pd.DataFrame:
    """STO Tachibana 2008 (APL 93, 092902) experimental κ(300 K) ≈ 11 W/mK.

    Approximate anchor as recorded in the data inventory; benchmark-marker
    use only.
    """
    return _load_csv("STO/tachibana_2008/kappa_total_T.csv",
                     _provenance("Tachibana2008", "κ(T), 300 K anchor"))


# ===========================================================================
# STO -- Verdi 2023
# ===========================================================================
def load_verdi_sscha_omega_T() -> pd.DataFrame:
    """STO Verdi 2023 Fig. 4(a) Ω_TO1(T): MLFF_SSCHA + INS_Ref39 + sqrt_fit.

    **Critical caveat:** the MLFF_SSCHA curve as plotted in Fig. 4(a) is
    shifted DOWN by 39.8 cm⁻¹ ≈ 4.93 meV to match the experimental low-T
    value (per Verdi caption). Recorded in the ``shift_meV`` column:
    ``4.93`` for MLFF_SSCHA and sqrt_fit, ``0`` for INS_Ref39. Subtract the
    shift back out (i.e. add 4.93 to Omega_meV) for the unshifted MLFF+SSCHA
    Cochran fit. The unshifted reference is Fig. S1 of the SI (not yet
    captured).
    """
    return _load_csv("STO/verdi_2023/sscha_omega_T.csv",
                     _provenance("Verdi2023", "Fig. 4(a)"))


def load_verdi_dispersion() -> pd.DataFrame:
    """STO Verdi 2023 Figs. 3 + 5 phonon dispersion at 0 K and 300 K."""
    return _load_csv("STO/verdi_2023/dispersion.csv",
                     _provenance("Verdi2023", "Figs. 3 + 5"))


def load_verdi_inverse_dielectric_T() -> pd.DataFrame:
    """STO Verdi 2023 Fig. 4(b) 1/ε(T) MLFF_SSCHA + Curie-Weiss fit.

    Used in NB 07 for the LST cross-check against the SSCHA Ω²(T).
    """
    return _load_csv("STO/verdi_2023/inverse_dielectric_T.csv",
                     _provenance("Verdi2023", "Fig. 4(b)"))


# ===========================================================================
# STO -- Yamada-Shirane 1969 (placeholders; PDF outstanding)
# ===========================================================================
def load_yamada_shirane_omega_T() -> pd.DataFrame:
    """STO Yamada-Shirane 1969 Ω(q, T).

    Holds the Table II q=0 soft-mode Ω(T) series (19 points, 4.5-293 K;
    ``condition_direction == "q_0_TO_soft"``) plus the 4.5 K [100]/[110]
    dispersion (Figs. 2/3). **Full Wave-2 verification of this source is
    still pending.**
    """
    return _load_csv("STO/yamada_shirane_1969/omega_T.csv",
                     _provenance("YamadaShirane1969", "Table II + Figs. 2/3"))


def load_yamada_shirane_gamma_T() -> pd.DataFrame:
    """STO Yamada-Shirane 1969 Fig. 6 intrinsic soft-mode Γ(T) at q=0.

    The experimental STO Γ(T) curve (≈ 0.24 → 0.84 meV over 20-250 K) that
    NB 08 will use to replace ``eta.eta_bracket`` with ``eta.eta_fit``.
    **Full Wave-2 verification of this source is still pending.**
    """
    return _load_csv("STO/yamada_shirane_1969/gamma_T.csv",
                     _provenance("YamadaShirane1969", "Fig. 6"))


# ===========================================================================
# STO -- Fauqué 2022 (Wave 1; arXiv:2203.15495)
# ===========================================================================
def load_fauque2022_omega_T() -> pd.DataFrame:
    """STO Fauqué 2022 Fig. 1(d) DHO mode energies Ω(T), 4-1400 K.

    The **FE soft mode** is ``mode == "TO"`` at the Γ point
    (``q_point_rlu == "(2,0,0)"``); this is the curve that extends the STO
    Ω(T) benchmark to 1400 K. The file also carries the **non-FE** R-point
    (``mode == "AFD"``) and M-point (``mode == "M_mode"``) branches for
    context -- those must NOT be used for η or the κ kernel.

    **No Γ(T) in this source.** Fauqué Fig. 1 plots only mode energies E (=ω);
    neither the panels nor the SI tabulate the DHO linewidth Γ
    (user-verified 2026-05-31). ``gamma_meV`` and ``gamma_uncertainty_meV``
    are therefore entirely NaN, by design -- a fact locked by a regression
    test. STO Γ(T) must come from another source (Yamada-Shirane 1969 Fig. 6,
    Wave 2).

    Note: Fig. 1(d)'s inset reports the Γ-vs-R splitting E₀ ≈ 4.1 meV; do not
    confuse the Γ-TO soft mode with the R-AFD mode when reading the figure.
    """
    return _load_csv("STO/fauque_2022/dho_fits_T.csv",
                     _provenance("Fauque2022", "Fig. 1(d)"))


def load_fauque2022_tdep_dispersion_300K() -> pd.DataFrame:
    """STO Fauqué 2022 Fig. 1(c) TDEP dispersion at 300 K along M-Γ-R.

    The Γ-point ``TO_Soft_Mode`` (≈ 9.5 meV at 300 K) is the FE soft mode;
    ``AFD_Soft_Mode`` at R is the non-FE octahedral-rotation branch.
    """
    return _load_csv("STO/fauque_2022/tdep_dispersion_300K.csv",
                     _provenance("Fauque2022", "Fig. 1(c)"))


# ===========================================================================
# STO -- Schmidt-Spaldin 2025 (Wave 1; arXiv:2508.10735)
# ===========================================================================
def load_schmidt2025_omega_T() -> pd.DataFrame:
    """STO Schmidt-Spaldin 2025 Fig. 3 (top left) Ω_Eu(T) for ST¹⁶O + ST¹⁸O.

    Third independent SSCHA Ω(T) for STO, with the ¹⁶O/¹⁸O isotope pair.

    **Imaginary-frequency convention:** the unstable harmonic reference
    gives Ω² < 0 at low T; these are tabulated as **negative** ``Omega_meV``
    (and negative ``Omega_sq_meV2``), matching the Ehsan/Verdi convention.
    ``kappa_kernel`` already skips Ω² < 0 rows, so the negative entries are
    expected and must not be clipped here. (e.g. ST¹⁸O SSCHA at 0 K:
    Ω = -2.11 meV.)
    """
    return _load_csv("STO/schmidt_2025/sscha_omega_T.csv",
                     _provenance("Schmidt-Spaldin2025", "Fig. 3 (top left)"))


def load_schmidt2025_energy_polarization() -> pd.DataFrame:
    """STO Schmidt-Spaldin 2025 Fig. 2 SSCHA/RPA energy vs polarization (0 K).

    **Provenance-only -- never a kernel input.** Illustrates the
    double-well anharmonic potential (RPA) vs the SSCHA-renormalized single
    well for ¹⁶O and ¹⁸O.

    ``polarization_uC_per_cm2`` is in **µC/cm²** (Schmidt Fig. 2 axis,
    user-confirmed 2026-05-31). The column was renamed from the original
    ``polarization_C_per_m2`` to keep the unit-suffix-in-name invariant.
    """
    return _load_csv("STO/schmidt_2025/energy_polarization.csv",
                     _provenance("Schmidt-Spaldin2025", "Fig. 2"))


# ===========================================================================
# KTO -- Meier 2022 (Wave 1; arXiv:2206.08296; third material)
# ===========================================================================
def load_meier2022_omega_T() -> pd.DataFrame:
    """KTaO₃ Meier 2022 Fig. 3(a) TO soft-mode Ω(T), 0-300 K.

    Primary KTaO₃ benchmark (``method`` ∈ {``Calculated`` = SSCHA+MLIP,
    ``Shirane_exp`` = Shirane INS [Ref 10]). The low-T quantum-paraelectric
    plateau holds Ω ≈ 3.3 meV up to ~25 K (Ω(0 K) = Ω(10 K) = 3.33 meV).

    **Lattice lock for the κ kernel:** a = 4.007 Å (Meier SSCHA+MLIP, 0 K) →
    V_cell = 64.34 Å³. The experimental value is a = 3.990 Å (Vousden, via
    Meier Ref. 39) → V_cell = 63.52 Å³; the <1.5% V difference sits inside
    the kernel's sensitivity envelope. Use the Meier calc value for internal
    consistency with the SSCHA Ω(T).
    """
    return _load_csv("KTO/meier_2022/omega_T.csv",
                     _provenance("Meier2022", "Fig. 3(a)"))


def load_meier2022_inverse_dielectric_T() -> pd.DataFrame:
    """KTaO₃ Meier 2022 Fig. 3(b) 1/ε(T): Calculated + Fujishita exp + Barrett fit.

    LST cross-check on the Fig. 3(a) soft-mode Ω²(T) (NB 07, Wave 2).
    """
    return _load_csv("KTO/meier_2022/inverse_dielectric_T.csv",
                     _provenance("Meier2022", "Fig. 3(b)"))


def load_meier2022_dispersion_T() -> pd.DataFrame:
    """KTaO₃ Meier 2022 Fig. 2 phonon spectrum at 0 K and 300 K.

    **Γ-point soft mode only -- NOT a small-q dispersion.** The soft optical
    branch is sampled at Γ alone (harmonic small-displacement Ω = -12.41 meV
    imaginary; SSCHA Ω = 3.31 meV at 0 K, 10.75 meV at 300 K); the X/M/R
    points carry only acoustic and highest-optical branches. There is no
    near-Γ soft-TO curvature here, so **Λ²_KTO is NOT fittable from Meier** --
    it would need another source (e.g. a dedicated KTaO₃ INS dispersion).
    """
    return _load_csv("KTO/meier_2022/dispersion_T.csv",
                     _provenance("Meier2022", "Fig. 2"))


# ===========================================================================
# KTO -- Ranalli 2023 (second KTaO3 source; MLFF+SSCHA, Verdi-lineage)
# ===========================================================================
def load_ranalli2023_omega_T() -> pd.DataFrame:
    """KTaO₃ Ranalli 2023 Fig. 4(a) TO1 Γ-point Ω(T), 0-300 K.

    Second KTaO₃ benchmark (after Meier 2022); same MLFF+SSCHA lineage as
    Verdi STO (Verdi co-authors both). ``method`` ∈ {``MLFF+SSCHA`` (computed,
    incl. bubble term — the curve NB 07b triangulates against Meier),
    ``expt`` (INS experimental anchor), ``fit_paraelectric`` (the
    ω²_TO ∝ C(T−Tc) √-fit line, showing the <30 K quantum-paraelectric
    plateau breakdown)}.

    **Native unit / shift convention (same family as Verdi −4.93 meV and
    Schmidt negative-Ω):** Fig. 4(a)'s y-axis is **cm⁻¹**; ``Omega_meV`` stores
    the converted value (×0.1239842). The ``MLFF+SSCHA`` rows are plotted
    **shifted DOWN by 34.71 cm⁻¹ ≈ 4.30 meV** to match the experimental low-T
    soft mode, recorded as ``shift_meV = -4.30`` (negative = down-shift);
    ``expt`` and ``fit_paraelectric`` rows are unshifted (``shift_meV = 0``).
    To recover the unshifted computed Cochran frequency, subtract the shift:
    Ω_unshifted = Omega_meV − shift_meV (i.e. + 4.30 meV). Do not pre-correct.

    **phase = cubic throughout** (KTaO₃ is cubic at all T — its appeal as a
    clean quantum-paraelectric control; not carried as a column, matching the
    Verdi/Meier omega_T schema).

    Cross-check vs Meier (Ω(0)=3.33, Ω(300)=10.86 meV): the *experimental*
    endpoints agree well (Ranalli expt 300 K ≈ 10.3 meV ≈ Meier exp 10.7); the
    *computed* curves carry a ~1-2 meV model/shift-convention spread (Ranalli
    shifted MLFF+SSCHA 300 K ≈ 9.05 meV; low-T plateau ≈ 2.48 meV vs Meier
    3.33) — a genuine inter-source spread, not a digitization error.
    """
    return _load_csv("KTO/ranalli_2023/omega_T.csv",
                     _provenance("Ranalli2023", "Fig. 4(a) [native cm⁻¹]"))


def load_ranalli2023_inverse_dielectric_T() -> pd.DataFrame:
    """KTaO₃ Ranalli 2023 Fig. 4(b) 1/ε(T) + Curie-Weiss line.

    ``method`` ∈ {``MLFF+SSCHA`` (computed from the shifted Fig. 4(a)
    frequencies, per the caption), ``expt``, ``Curie-Weiss_fit`` (the blue
    classical reference line)}. Used for the KTaO₃ LST cross-check (1/ε vs Ω²)
    parallel to the STO one, and shows the Barrett low-T saturation.
    """
    return _load_csv("KTO/ranalli_2023/inverse_dielectric_T.csv",
                     _provenance("Ranalli2023", "Fig. 4(b)"))


# ===========================================================================
# KTO + STO -- Vogt 1995 (PRB 51, 9426; experimental hyper-Raman, Vogt wave)
# ===========================================================================
# γ-CONVENTION (gating verdict, recorded once for all four loaders):
# Vogt's Fig. 5/6 fit is "a damped harmonic oscillator response and a δ function
# at the center" (Sec. III) -- the central peak is SEPARATED, so the plotted γ is
# the intrinsic zone-center DHO damping (kernel-appropriate). The plotted symbol
# is lowercase γ (the standard DHO damping constant), distinct from the
# theoretical Eq. (13) self-energy Γ (which is HWHM). Reading γ as the standard
# DHO constant => γ = FWHM => Γ_HWHM = γ/2. This matches the brief's §1 numbers
# (STO 300 K: Ω=10.54, γ=3.04 → η=0.15 on the FWHM basis) and turns the raw
# YS-vs-Vogt 3.3× into a ~1.7× residual physical difference.
# **Factor-of-2 systematic FLAG:** if Vogt's fit actually used his Eq.-13 Γ
# (HWHM), use γ directly and η doubles (0.15→0.29). Refs. 29-30 (PRB 34,410
# (1986); 41,1184 (1990)) -- not in repo -- would settle it. BOTH readings keep
# STO and KTaO3 underdamped (η<0.5) and clustered; only the absolute η
# calibration carries the factor. eta.py applies Γ_HWHM = γ/2 by default (later
# wave). These loaders return the plotted γ AS-IS (FWHM); conversion is at the
# η-consumption layer, never silent.
def load_vogt1995_kto_omega_T() -> pd.DataFrame:
    """KTaO₃ Vogt 1995 Fig. 5 experimental TO soft-mode Ω(T), 5-300 K (q=0).

    Hyper-Raman, true q=0, central peak separated. ``method`` ∈
    {``exp_HyperRaman_Vogt1995`` (Fig. 5 curve), ``exp_HyperRaman_Vogt1988``
    (the single 300 K = 10.54 meV point from Table III)}. ``shift_meV = 0``
    (experimental). The experimental tie-breaker for the ~1.8 meV Meier-Ranalli
    computed-curve spread (supports the higher high-T frequency).
    """
    return _load_csv("KTO/vogt_1995/omega_T.csv",
                     _provenance("Vogt1995", "Fig. 5"))


def load_vogt1995_kto_gamma_T() -> pd.DataFrame:
    """KTaO₃ Vogt 1995 Fig. 5 experimental soft-mode damping γ(T), 5-300 K.

    Single-phase (KTaO₃ cubic at all T) -> 3-col schema, no ``phase_mode``.
    Plotted γ is the DHO damping constant (= FWHM; see module γ-convention note);
    η = (γ/2)/Ω. **Primary η source for KTaO₃** (Amendment-2 §4) -> eta_fit.
    """
    return _load_csv("KTO/vogt_1995/gamma_T.csv",
                     _provenance("Vogt1995", "Fig. 5"))


def load_vogt1995_sto_omega_T() -> pd.DataFrame:
    """STO Vogt 1995 Fig. 6 experimental TO soft-mode Ω(T), 5-300 K (q=0).

    Hyper-Raman, true q=0. Fig. 6 straddles the 105 K AFD transition, so rows
    carry ``phase_mode`` ∈ {``tetragonal_A2u``, ``tetragonal_Eu`` (< 105 K),
    ``cubic_F1u`` (≥ 110 K)}. **The kernel/η must filter to ``cubic_F1u``**
    (the sub-105 K A₂ᵤ/Eᵤ split is low-T QPE discussion only). ``shift_meV=0``.
    Adds the 6th source to the STO Ω(T) triangulation (Fig. C / NB 07).
    """
    return _load_csv("STO/vogt_1995/omega_T.csv",
                     _provenance("Vogt1995", "Fig. 6"))


def load_vogt1995_sto_gamma_T() -> pd.DataFrame:
    """STO Vogt 1995 Fig. 6 experimental soft-mode damping γ(T), 5-300 K.

    Phase-split -> 4-col schema with ``phase_mode``. Below 105 K Vogt fit the
    A₂ᵤ/Eᵤ split as two DHOs of different Ω but the SAME γ, so both tetragonal
    rows share one damping; ``cubic_F1u`` (≥ 110 K) is the kernel-relevant
    single mode. **Primary η source for STO** (Amendment-2 §4) -> eta_fit, using
    ``cubic_F1u``. Plotted γ = FWHM (see module note); η = (γ/2)/Ω. The ~3.3×
    spread vs YS 1969 (HWHM, finite-q, central-peak-contaminated) is documented
    content (Amendment-2 §2), NOT a data error to reconcile.
    """
    return _load_csv("STO/vogt_1995/gamma_T.csv",
                     _provenance("Vogt1995", "Fig. 6"))


# ===========================================================================
# Umbrella
# ===========================================================================
def load_all() -> dict[str, Any]:
    """Load every literature CSV; key by ``"{material}.{source}.{quantity}"``.

    Returns:
        A dict of 38 entries (21 Turn-4A + 7 Wave-1 arXiv + 4 Wave-2
        Harada-Axe-Shirane 1971 cubic-BTO + 2 Ranalli 2023 KTaO3 + 4 Vogt 1995
        experimental KTaO3/STO Ω(T)+γ(T); the Yamada-Shirane 1969 loaders were
        already present). Values are :class:`pandas.DataFrame`
        except the three scalar/anchor entries (``BTO.ehsan2021.anchors``,
        ``BTO.shirane1970.Lambda_sq``, ``BTO.harada_axe_shirane1971.Lambda_sq``),
        which are :class:`types.SimpleNamespace`.
    """
    return {
        "BTO.ehsan2021.anchors":                load_ehsan_anchors(),
        "BTO.ehsan2021.scp_omega_T":            load_ehsan_scp_omega_T(),
        "BTO.ehsan2021.free_energy_T":          load_ehsan_free_energy_T(),
        "BTO.ehsan2021.structural":             load_ehsan_structural(),
        "BTO.ehsan2021.cubic_dispersion":       load_ehsan_cubic_dispersion(),
        "BTO.ehsan2021.born_charges":           load_ehsan_born_charges(),
        "BTO.tadano2022.free_energy_T":         load_tadano2022_free_energy_T(),
        "BTO.tadano2022.scp_dispersion_plus1pct": load_tadano2022_scp_dispersion_plus1pct(),
        "BTO.tadano2022.transition_temperatures": load_tadano2022_transition_temperatures(),
        "BTO.tadano2022.kappa_total_T":         load_tadano2022_kappa_total_T(),
        "BTO.shirane1970.Lambda_sq":            load_shirane_Lambda_sq(),
        "BTO.shirane1970.cubic_dispersion":     load_shirane_cubic_dispersion(),
        "BTO.harada_axe_shirane1971.cubic_dispersion": load_harada_axe_shirane_1971_cubic_dispersion(),
        "BTO.harada_axe_shirane1971.mode_parameters":  load_harada_axe_shirane_1971_mode_parameters(),
        "BTO.harada_axe_shirane1971.anisotropy":       load_harada_axe_shirane_1971_anisotropy(),
        "BTO.harada_axe_shirane1971.Lambda_sq":        load_harada_axe_shirane_1971_Lambda_sq(),
        "STO.tadano2015.scph_omega_T":          load_tadano2015_scph_omega_T(),
        "STO.tadano2015.scph_dispersion":       load_tadano2015_scph_dispersion(),
        "STO.tadano2015.scph_convergence_300K": load_tadano2015_scph_convergence_300K(),
        "STO.tadano2015.kappa_total_T":         load_tadano2015_kappa_total_T(),
        "STO.verdi2023.sscha_omega_T":          load_verdi_sscha_omega_T(),
        "STO.verdi2023.dispersion":             load_verdi_dispersion(),
        "STO.verdi2023.inverse_dielectric_T":   load_verdi_inverse_dielectric_T(),
        "STO.yamada_shirane1969.omega_T":       load_yamada_shirane_omega_T(),
        "STO.yamada_shirane1969.gamma_T":       load_yamada_shirane_gamma_T(),
        # Wave 1 -- new arXiv sources
        "STO.fauque2022.omega_T":               load_fauque2022_omega_T(),
        "STO.fauque2022.tdep_dispersion_300K":  load_fauque2022_tdep_dispersion_300K(),
        "STO.schmidt2025.omega_T":              load_schmidt2025_omega_T(),
        "STO.schmidt2025.energy_polarization":  load_schmidt2025_energy_polarization(),
        "KTO.meier2022.omega_T":                load_meier2022_omega_T(),
        "KTO.meier2022.inverse_dielectric_T":   load_meier2022_inverse_dielectric_T(),
        "KTO.meier2022.dispersion_T":           load_meier2022_dispersion_T(),
        "KTO.ranalli2023.omega_T":              load_ranalli2023_omega_T(),
        "KTO.ranalli2023.inverse_dielectric_T": load_ranalli2023_inverse_dielectric_T(),
        # Vogt 1995 wave -- experimental Ω(T) + γ(T), both QPEs
        "KTO.vogt1995.omega_T":                 load_vogt1995_kto_omega_T(),
        "KTO.vogt1995.gamma_T":                 load_vogt1995_kto_gamma_T(),
        "STO.vogt1995.omega_T":                 load_vogt1995_sto_omega_T(),
        "STO.vogt1995.gamma_T":                 load_vogt1995_sto_gamma_T(),
    }
