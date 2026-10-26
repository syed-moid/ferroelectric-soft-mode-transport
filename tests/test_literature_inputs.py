"""Validation tests for the v1 literature CSV templates.

These tests are designed to be run *after* the data are digitized into the
empty templates under ``data/literature/``. Several of them are expected to
fail until the corresponding CSVs are filled; that is by design -- they
serve as a checklist for the digitization pass.

Pre-filled CSVs (Ehsan 2021 anchors, Shirane Lambda_sq central value) and
the existence/schema tests pass immediately. The remaining tests pass as the
required CSVs are populated.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
LIT = REPO_ROOT / "data" / "literature"

# Header schemas declared in data/literature/digitization_plan.md.
EXPECTED_SCHEMAS: dict[str, list[str]] = {
    # --- BTO Ehsan 2021 ---
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
    # --- BTO Tadano 2022 ---
    "BTO/tadano_2022/free_energy_T.csv":
        ["T_K", "F_cubic_meV", "F_tetra_meV", "F_ortho_meV", "F_rhombo_meV"],
    "BTO/tadano_2022/scp_dispersion_plus1pct.csv":
        ["T_K", "branch_label", "q_label", "omega_cm_inv"],
    "BTO/tadano_2022/transition_temperatures.csv":
        ["transition", "T_c_calc_K", "T_c_exp_K"],
    "BTO/tadano_2022/kappa_total_T.csv":
        ["T_K", "kappa_W_per_m_K", "kappa_uncertainty", "method", "source"],
    # --- BTO Shirane 1970 ---
    "BTO/shirane_1970/Lambda_sq.csv":
        ["Lambda_sq_meV2_A2", "Lambda_sq_uncertainty_meV2_A2", "T_K",
         "reference_q_range"],
    "BTO/shirane_1970/cubic_dispersion.csv":
        ["q_label", "q_inv_angstrom", "omega_meV", "omega_uncertainty_meV",
         "branch", "polarization_e", "wavevector_q", "phase", "T_K",
         "source_figure", "note"],
    # --- BTO Harada-Axe-Shirane 1971 (cubic-phase counterpart to the
    #     tetragonal Shirane 1970; removes the tetragonal-phase caveat on the
    #     BTO soft mode via Fig. 5/8 dispersion + Table I damping) ---
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
    # --- STO Tadano 2015 (schemas extended v1: mode + method on
    #     scph_omega_T; method + source_figure on scph_dispersion) ---
    "STO/tadano_2015/scph_omega_T.csv":
        ["T_K", "mode", "method", "Omega_meV", "Omega_sq_meV2",
         "source_figure"],
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
    # --- STO Verdi 2023 (sscha_omega_T extended with shift_meV column;
    #     dispersion and inverse_dielectric_T are new v1 additions) ---
    "STO/verdi_2023/sscha_omega_T.csv":
        ["T_K", "Omega_meV", "Omega_sq_meV2", "method", "shift_meV",
         "source_figure"],
    "STO/verdi_2023/dispersion.csv":
        ["T_K", "q_label", "q_inv_angstrom", "omega_cm_inv", "branch",
         "method", "supercell", "source_figure"],
    "STO/verdi_2023/inverse_dielectric_T.csv":
        ["T_K", "inverse_epsilon", "method", "source_figure"],
    # --- STO Yamada-Shirane 1969 (new 6-col omega schema includes
    #     material + q so the Shirane 1970 Fig. 1 STO panel fits cleanly) ---
    "STO/yamada_shirane_1969/omega_T.csv":
        ["material", "condition_direction", "q_inv_angstrom", "T_K",
         "omega_meV", "omega_uncertainty_meV"],
    "STO/yamada_shirane_1969/gamma_T.csv":
        ["T_K", "gamma_meV", "gamma_uncertainty_meV"],
    # --- STO Fauque 2022 (v1 HIGHEST-impact STO source: multi-T DHO
    #     fits at Γ, R, M over 4-1400 K) ---
    "STO/fauque_2022/dho_fits_T.csv":
        ["T_K", "q_label", "q_point_rlu", "mode",
         "omega_meV", "omega_uncertainty_meV",
         "gamma_meV", "gamma_uncertainty_meV",
         "method", "source_figure"],
    "STO/fauque_2022/tdep_dispersion_300K.csv":
        ["T_K", "q_label", "q_inv_angstrom", "omega_meV", "branch",
         "method", "source_figure"],
    # --- STO Schmidt 2025 (third independent SSCHA Ω(T) for STO; isotope
    #     16O/18O cross-check) ---
    "STO/schmidt_2025/sscha_omega_T.csv":
        ["T_K", "isotope", "Omega_meV", "Omega_sq_meV2", "method",
         "source_figure"],
    "STO/schmidt_2025/energy_polarization.csv":
        ["T_K", "polarization_uC_per_cm2", "E_meV_per_fu", "isotope",
         "method", "source_figure"],
    # --- KTO Meier 2022 (third material: KTaO3 SSCHA+MLIP, parallel to
    #     STO Verdi 2023; feeds Section 4 parametric κ_soft/κ_0 fit) ---
    "KTO/meier_2022/omega_T.csv":
        ["T_K", "Omega_meV", "Omega_sq_meV2", "method", "source_figure"],
    "KTO/meier_2022/inverse_dielectric_T.csv":
        ["T_K", "inverse_epsilon", "method", "source_figure"],
    "KTO/meier_2022/dispersion_T.csv":
        ["T_K", "q_label", "q_inv_angstrom", "omega_meV", "branch",
         "method", "source_figure"],
    # --- KTO Ranalli 2023 (second KTaO3 source; MLFF+SSCHA, same lineage as
    #     Verdi STO -> omega_T schema is identical to Verdi's, incl. shift_meV.
    #     Native Fig. 4 axis is cm^-1; CSV stores meV (x0.1239842). KTaO3 is
    #     cubic at all T, so phase is constant and lives in provenance, not a
    #     column -- matching the Verdi/Meier omega_T pattern.) ---
    "KTO/ranalli_2023/omega_T.csv":
        ["T_K", "Omega_meV", "Omega_sq_meV2", "method", "shift_meV",
         "source_figure"],
    "KTO/ranalli_2023/inverse_dielectric_T.csv":
        ["T_K", "inverse_epsilon", "method", "source_figure"],
    # --- Vogt 1995 (experimental hyper-Raman Ω(T) + γ(T) for BOTH QPEs) ---
    # omega: Omega_uncertainty_meV (not Omega_sq); shift_meV = 0 (experimental).
    # gamma is INTENTIONALLY two-variant: STO-Vogt carries phase_mode (Fig. 6
    # straddles the 105 K AFD transition: tetragonal_A2u/_Eu < 105 K, cubic_F1u
    # ≥ 110 K); KTaO3-Vogt is single-phase (3-col). Do NOT reconcile to uniform.
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

# CSVs that MUST have at least one data row before v1 can start.
# Pre-filled and to-be-digitized CSVs are both in here. CSVs not listed here
# are optional and may stay empty as long as a note is added to the
# extraction log.
REQUIRED_NONEMPTY: set[str] = {
    # Pre-filled
    "BTO/ehsan_2021/anchors.csv",
    "BTO/shirane_1970/Lambda_sq.csv",
    # v1-extracted (Ehsan high/medium priority)
    "BTO/ehsan_2021/free_energy_T.csv",
    "BTO/ehsan_2021/structural_properties.csv",
    "BTO/ehsan_2021/cubic_dispersion.csv",
    "BTO/ehsan_2021/born_charges_epsilon_infty.csv",
    # v1-extracted (Tadano 2022)
    "BTO/tadano_2022/free_energy_T.csv",
    "BTO/tadano_2022/scp_dispersion_plus1pct.csv",
    "BTO/tadano_2022/transition_temperatures.csv",
    "BTO/tadano_2022/kappa_total_T.csv",
    # v1-extracted (Shirane 1970 Fig. 3 + Fig. 6)
    "BTO/shirane_1970/cubic_dispersion.csv",
    # v1-extracted (Yamada-Shirane 1969 Table II pre-filled: q=0 Ω(T))
    "STO/yamada_shirane_1969/omega_T.csv",
    # v1-extracted (STO Tadano 2015 high + medium priority)
    "STO/tadano_2015/scph_omega_T.csv",         # Fig. 5
    "STO/tadano_2015/scph_dispersion.csv",      # Fig. 4
    "STO/tadano_2015/scph_convergence_300K.csv", # Table I
    # v1-extracted (STO Verdi 2023 high + medium priority)
    "STO/verdi_2023/sscha_omega_T.csv",         # Fig. 4(a)
    "STO/verdi_2023/dispersion.csv",            # Figs. 3 + 5
    "STO/verdi_2023/inverse_dielectric_T.csv",  # Fig. 4(b)
    # v1 HIGH-priority TODOs (templates created; awaiting digitization)
    "STO/fauque_2022/dho_fits_T.csv",           # Fig. 1(d) Γ/R/M + Figs. 2-3 TA
    "STO/schmidt_2025/sscha_omega_T.csv",       # Fig. 3 top left (16O + 18O)
    "KTO/meier_2022/omega_T.csv",               # Fig. 3(a) KTaO3 Ω²(T)
    # v1 HIGH-priority TODOs (new library acquisitions)
    "STO/yamada_shirane_1969/gamma_T.csv",      # Fig. 6 experimental STO Γ(T)
    "BTO/harada_axe_shirane_1971/cubic_dispersion.csv",  # Fig. 5 + Fig. 8 cubic Λ²
    "BTO/harada_axe_shirane_1971/mode_parameters.csv",   # Table I ω + Γ (cubic)
    # v1 HIGH TODO: second KTaO3 source (template; awaiting digitization)
    "KTO/ranalli_2023/omega_T.csv",             # Fig. 4(a) TO1 Ω(T), 0-300 K
    # Vogt wave: experimental Ω(T)+γ(T), primary η source for both QPEs
    "KTO/vogt_1995/omega_T.csv",                # Fig. 5 KTaO3 Ω(T)
    "KTO/vogt_1995/gamma_T.csv",                # Fig. 5 KTaO3 γ(T)
    "STO/vogt_1995/omega_T.csv",                # Fig. 6 STO Ω(T) (phase-split)
    "STO/vogt_1995/gamma_T.csv",                # Fig. 6 STO γ(T) (phase-split)
}

OMEGA_T_CSVS: dict[str, dict] = {
    # path -> {"required_monotonic": bool, "T_col": str, "Omega_col": str}.
    # All present Omega(T) curves mix multiple methods or modes in one file
    # (TO1/TO2, SCPH_with_PM/without_PM/TDEP/exp, MLFF_SSCHA/INS/sqrt_fit),
    # so a naive T-sort + monotonicity check is ill-defined. Strict
    # monotonicity is deferred to per-method notebook checks.
    "BTO/ehsan_2021/scp_omega_T.csv": {
        "required_monotonic": False, "T_col": "T_K", "Omega_col": "Omega_meV",
    },
    "STO/tadano_2015/scph_omega_T.csv": {
        "required_monotonic": False, "T_col": "T_K", "Omega_col": "Omega_meV",
    },
    "STO/verdi_2023/sscha_omega_T.csv": {
        "required_monotonic": False, "T_col": "T_K", "Omega_col": "Omega_meV",
    },
    "STO/yamada_shirane_1969/omega_T.csv": {
        "required_monotonic": False, "T_col": "T_K", "Omega_col": "omega_meV",
    },
}

# Column-wise sanity ranges. Frequencies (Omega_meV, omega_meV, omega_cm_inv)
# allow negative values: imaginary frequencies of unstable soft modes are
# conventionally tabulated as -|omega|.
UNIT_RANGES: dict[str, tuple[float, float]] = {
    "T_K": (0.0, 1500.0),
    "Omega_meV": (-200.0, 200.0),
    "omega_meV": (-200.0, 200.0),
    "Omega_sq_meV2": (-1000.0, 10000.0),
    "omega_sq_meV2": (-1000.0, 10000.0),
    "q_inv_angstrom": (0.0, 3.0),
    "Lambda_sq_meV2_A2": (100.0, 50000.0),
    "Lambda_sq_uncertainty_meV2_A2": (0.0, 50000.0),
    # Bound set deliberately from HAS 1971 (Wave-2 figure review, not a
    # placeholder): cubic-BTO soft TO is OVERDAMPED. Table I max Γ_TO = 34±6
    # meV (q=0.2 a*, 450 °C; Γ_TO/Ω_TO ≈ 2.3, i.e. Γ ≳ 2Ω); the cited optical
    # Γ/Ω = 9 at q = 0.08 a* with Ω₀ ≈ 3 meV implies Γ up to ~44 meV at small
    # q. 50 meV brackets that overdamped regime with margin. STO modes are far
    # narrower (Yamada-Shirane Fig. 6: Γ ≤ 0.84 meV over 20-250 K).
    "gamma_meV": (0.0, 50.0),
    "omega_cm_inv": (-500.0, 1500.0),
    "T_c_calc_K": (0.0, 1500.0),
    "T_c_exp_K": (0.0, 1500.0),
    "F_cubic_meV": (-1000.0, 1000.0),
    "F_tetra_meV": (-1000.0, 1000.0),
    "F_ortho_meV": (-1000.0, 1000.0),
    "F_rhombo_meV": (-1000.0, 1000.0),
    "omega_uncertainty_meV": (0.0, 100.0),
    "Omega_uncertainty_meV": (0.0, 100.0),  # Vogt omega files (capital Ω)
    "gamma_uncertainty_meV": (0.0, 20.0),
    # Structural / electronic properties (v1 extensions)
    "a_angstrom": (3.0, 6.0),
    "c_angstrom": (3.0, 6.0),
    "c_over_a": (0.9, 1.2),
    "V_cell_A3": (20.0, 200.0),
    "Wyckoff_x": (0.0, 1.0),
    "Wyckoff_y": (0.0, 1.0),
    "Wyckoff_z": (0.0, 1.0),
    "Z_xx": (-10.0, 15.0),
    "Z_yy": (-10.0, 15.0),
    "Z_zz": (-10.0, 15.0),
    "epsilon_inf_xx": (1.0, 20.0),
    "epsilon_inf_yy": (1.0, 20.0),
    "epsilon_inf_zz": (1.0, 20.0),
    # Verdi 2023 / Tadano 2015 STO extensions
    "shift_meV": (-100.0, 100.0),
    "inverse_epsilon": (0.0, 0.1),
    "kappa_W_per_m_K": (0.0, 1000.0),
    "kappa_uncertainty": (0.0, 1000.0),
    # v1 Schmidt 2025 extensions (polarization axis is µC/cm² per
    # Schmidt Fig. 2, user-confirmed 2026-05-31; column renamed from the
    # original _C_per_m2 to keep the unit-suffix-in-name invariant correct)
    "polarization_uC_per_cm2": (-2.0, 2.0),
    "E_meV_per_fu": (-200.0, 200.0),
    # v1 Harada-Axe-Shirane 1971 extension (Fig. 7 soft-TO anisotropy
    # arc): the abscissa is a *signed* angular orientation of q measured from
    # a reference symmetry axis ([100]=0°, [111]≈-54.7°, [110]=+45°).
    "angle_deg": (-90.0, 90.0),
}


def _read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Return ``(fieldnames, rows)`` for a CSV file. Empty CSV -> ``([], [])``."""
    if not path.exists():
        return [], []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames or [])
        rows = [dict(r) for r in reader]
    return fieldnames, rows


# ---------------------------------------------------------------------------
# 1. Existence
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("rel_path", sorted(EXPECTED_SCHEMAS))
def test_all_required_csvs_exist(rel_path: str) -> None:
    """Every CSV path declared in digitization_plan.md exists (even empty)."""
    assert (LIT / rel_path).exists(), (
        f"Missing template CSV: data/literature/{rel_path}. "
        "Re-run v1 scaffolding."
    )


# ---------------------------------------------------------------------------
# 2. Schema match
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("rel_path", sorted(EXPECTED_SCHEMAS))
def test_csv_schemas_match_plan(rel_path: str) -> None:
    """Header row matches the schema declared in the plan."""
    expected = EXPECTED_SCHEMAS[rel_path]
    fieldnames, _ = _read_rows(LIT / rel_path)
    assert fieldnames == expected, (
        f"Schema mismatch in {rel_path}: expected {expected}, got {fieldnames}"
    )


# ---------------------------------------------------------------------------
# 3. Required CSVs must have at least one data row (initial-failure check)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("rel_path", sorted(REQUIRED_NONEMPTY))
def test_required_csv_has_data(rel_path: str) -> None:
    """Required CSVs must have >= 1 data row beyond the header.

    Pre-filled CSVs pass immediately. To-be-digitized CSVs fail until digitized.
    """
    _, rows = _read_rows(LIT / rel_path)
    assert rows, (
        f"{rel_path} is empty (header only). Digitize per "
        "data/literature/digitization_plan.md, or move it to OPTIONAL with a "
        "note in extraction_log.md."
    )


# ---------------------------------------------------------------------------
# 4. Monotonicity of Omega(T) where it is required (STO SCPH/SSCHA)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("rel_path", sorted(OMEGA_T_CSVS))
def test_omega_T_curves_monotonicity(rel_path: str) -> None:
    """For STO SCPH/SSCHA curves, Omega(T) should be monotonically increasing.

    BTO curves are exempt (the soft mode can be V-shaped through T_c).
    Empty CSVs are skipped here -- they're caught by ``test_required_csv_has_data``
    when required, and tolerated otherwise.
    """
    cfg = OMEGA_T_CSVS[rel_path]
    _, rows = _read_rows(LIT / rel_path)
    if not rows:
        pytest.skip(f"{rel_path} is empty; monotonicity check deferred")
    if not cfg["required_monotonic"]:
        return  # BTO: V-shape allowed

    pairs = []
    for r in rows:
        try:
            T_val = float(r[cfg["T_col"]])
            Omega_val = float(r[cfg["Omega_col"]])
        except (KeyError, ValueError):
            continue
        pairs.append((T_val, Omega_val))
    pairs.sort()
    omegas = [o for _, o in pairs]
    # Allow tiny numerical wobble (1%); flag larger drops.
    for i in range(1, len(omegas)):
        prev, cur = omegas[i - 1], omegas[i]
        if cur < prev - 0.01 * max(prev, 1.0):
            raise AssertionError(
                f"{rel_path}: Omega(T) decreasing at T={pairs[i][0]} "
                f"({prev:.3f} -> {cur:.3f} meV); expected monotone increase"
            )


# ---------------------------------------------------------------------------
# 5. Unit ranges on populated rows
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("rel_path", sorted(EXPECTED_SCHEMAS))
def test_units_in_range(rel_path: str) -> None:
    """Every populated numeric cell falls inside a sanity range for its unit."""
    fieldnames, rows = _read_rows(LIT / rel_path)
    if not rows:
        pytest.skip(f"{rel_path} is empty; unit checks deferred")
    for i, row in enumerate(rows, start=2):  # row 1 = header
        for col in fieldnames:
            if col not in UNIT_RANGES:
                continue
            raw = (row.get(col) or "").strip()
            if not raw or raw.lower() in ("nan", "none", "null", "n/a"):
                continue  # explicit missing markers are allowed
            try:
                v = float(raw)
            except ValueError:
                raise AssertionError(
                    f"{rel_path} row {i} col {col!r}: cannot parse {raw!r} as float"
                ) from None
            if math.isnan(v):
                continue
            lo, hi = UNIT_RANGES[col]
            assert lo <= v <= hi, (
                f"{rel_path} row {i} col {col!r}: value {v} outside "
                f"sanity range [{lo}, {hi}]"
            )


# ---------------------------------------------------------------------------
# 6. Ehsan anchors present and unmodified
# ---------------------------------------------------------------------------
def test_ehsan_anchors_present() -> None:
    """The pre-filled Ehsan two anchors must be present and unmodified."""
    _, rows = _read_rows(LIT / "BTO/ehsan_2021/anchors.csv")
    assert len(rows) == 2, "Ehsan anchors must have exactly two rows"
    by_T = {int(float(r["T_K"])): r for r in rows}
    assert 0 in by_T and 500 in by_T, "Anchors must be at T=0 K and T=500 K"
    assert float(by_T[0]["Omega_sq_meV2"]) == pytest.approx(-424.0, abs=0.1)
    assert float(by_T[500]["Omega_meV"]) == pytest.approx(21.08, abs=0.01)
    assert float(by_T[500]["Omega_sq_meV2"]) == pytest.approx(444.4, abs=0.5)


# ---------------------------------------------------------------------------
# 7. Shirane Lambda^2 central value
# ---------------------------------------------------------------------------
def test_shirane_lambda_sq_central_value() -> None:
    """Lambda^2 central value must be 5100 meV^2 * A^2 (Shirane 1970)."""
    _, rows = _read_rows(LIT / "BTO/shirane_1970/Lambda_sq.csv")
    assert rows, "Lambda_sq.csv must be pre-filled with the central value"
    val = float(rows[0]["Lambda_sq_meV2_A2"])
    assert val == pytest.approx(5100.0, rel=0.001), (
        f"Lambda^2 central value must be 5100 meV^2 * A^2, got {val}"
    )


# ---------------------------------------------------------------------------
# 8. Extraction log exists
# ---------------------------------------------------------------------------
def test_extraction_log_exists() -> None:
    """Extraction log must exist with notes on each source."""
    log = LIT / "extraction_log.md"
    assert log.exists(), "data/literature/extraction_log.md missing"
    assert log.stat().st_size > 0, "extraction_log.md is empty"
