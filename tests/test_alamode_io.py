"""Tests for the ALAMODE RTA parsers and the committed SrTiO3 tidy CSVs.

The committed CSVs under ``data/processed/sto_alamode/`` are tested
unconditionally (they are the data of record). Tests that re-parse the raw
``.result`` files skip when ``data/raw/alamode_sto/`` is absent -- the raw
directory is local-only (gitignored).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from fekappa.alamode_io import (
    KL_COLUMNS,
    RESULT_COLUMNS,
    SCPH_TEMPS_K,
    emit_tidy_csvs,
    kappa_isotropic,
    read_kl,
    read_result,
    read_result_velocities,
    tau_ps_from_gamma_hwhm_cm1,
)

REPO = Path(__file__).resolve().parents[1]
RAW_DIR = REPO / "data" / "raw" / "alamode_sto"
TIDY_DIR = REPO / "data" / "processed" / "sto_alamode"

needs_raw = pytest.mark.skipif(
    not (RAW_DIR / "STO_RTA_production.result").is_file(),  # the directory itself holds the tracked IMPORT_NOTE.md
    reason="raw ALAMODE dataset not present (local-only)"
)


def _read_tidy(name: str) -> pd.DataFrame:
    return pd.read_csv(TIDY_DIR / name, comment="#")


# --- committed tidy CSVs ----------------------------------------------------
def test_scph_csv_row_counts() -> None:
    """SCPH table: 35 irreducible q x 15 branches = 525 modes at each of 21 T."""
    df = _read_tidy("omega_gamma_scph_T.csv")
    assert sorted(df["T_K"].unique()) == list(SCPH_TEMPS_K)
    assert df.groupby("T_K").size().eq(525).all()
    assert len(df) == 525 * 21


def test_bare_csv_row_counts() -> None:
    """Bare table: 525 modes at each of 21 temperatures (0-1000 K)."""
    df = _read_tidy("omega_gamma_bare_T.csv")
    assert df["T_K"].nunique() == 21
    assert df["T_K"].min() == 0.0 and df["T_K"].max() == 1000.0
    assert df.groupby("T_K").size().eq(525).all()
    assert len(df) == 525 * 21


def test_kappa_300K_locked_value() -> None:
    """kappa(300 K) from the SCPH .kl files is the locked 9.5283 W/mK."""
    kappa = _read_tidy("kappa_rta_T.csv")
    assert len(kappa) == 21
    row300 = kappa[kappa["T_K"] == 300.0].iloc[0]
    assert row300["kappa_W_mK"] == pytest.approx(9.5283, abs=1e-4)
    # Direct Wigner coherent channel from KAPPA_COHERENT (not a bound).
    assert row300["kappa_coherent_W_mK"] == pytest.approx(0.3610, abs=1e-4)


def test_unstable_mode_counts_scph_vs_bare() -> None:
    """SCPH stabilizes every mode at 300 K; the bare surface has 11/525 unstable."""
    scph = _read_tidy("omega_gamma_scph_T.csv")
    scph300 = scph[scph["T_K"] == 300.0]
    assert (scph300["omega_cm1"] < 0).sum() == 0
    bare = _read_tidy("omega_gamma_bare_T.csv")
    bare300 = bare[bare["T_K"] == 300.0]
    assert (bare300["omega_cm1"] < 0).sum() == 11


def test_acoustic_54p9_mode_lifetimes() -> None:
    """The nearest-to-Gamma acoustic modes (54.9 cm^-1) give tau 16.19/17.54 ps."""
    scph = _read_tidy("omega_gamma_scph_T.csv")
    mode = scph[(scph["T_K"] == 300.0) & ((scph["omega_cm1"] - 54.9011).abs() < 1e-3)]
    assert len(mode) == 2  # the two degenerate TA branches at q = (0, 0, 1/8)
    tau = sorted(tau_ps_from_gamma_hwhm_cm1(mode["gamma_hwhm_cm1"].to_numpy()))
    assert tau[0] == pytest.approx(16.19, abs=0.01)
    assert tau[1] == pytest.approx(17.54, abs=0.01)


def test_provenance_headers_present() -> None:
    """Every tidy CSV declares its provenance; the bare one carries the warning."""
    for name in ("omega_gamma_scph_T.csv", "omega_gamma_bare_T.csv", "kappa_rta_T.csv"):
        header = (TIDY_DIR / name).read_text()[:2000]
        assert "provenance: ext-FP-alamode" in header, name
    bare_header = (TIDY_DIR / "omega_gamma_bare_T.csv").read_text()[:2000]
    assert "SUPERSEDED" in bare_header
    scph_header = (TIDY_DIR / "omega_gamma_scph_T.csv").read_text()[:2000]
    assert "HWHM" in scph_header


# --- HWHM-once regression ---------------------------------------------------
def test_tau_hwhm_applied_exactly_once() -> None:
    """tau = 1/(2*Gamma_HWHM): the factor 2 lives in the helper, exactly once.

    Feeding a FWHM (= 2*Gamma_HWHM) must give tau/2 -- if a consumer pre-halved
    the linewidth AND used this helper, the spot-checked 16.19/17.54 ps values
    above would come out doubled, so the pair of tests guards both directions.
    """
    tau_hwhm = float(tau_ps_from_gamma_hwhm_cm1(1.0))
    tau_fwhm_fed = float(tau_ps_from_gamma_hwhm_cm1(2.0))
    # 1/(2 * 2*pi*c * 1 cm^-1) = 2.6544 ps: the absolute anchor.
    assert tau_hwhm == pytest.approx(2.65442, abs=1e-4)
    assert tau_fwhm_fed == pytest.approx(tau_hwhm / 2.0, rel=1e-12)


# --- raw-file parsers (local-only data) -------------------------------------
@needs_raw
def test_read_result_scph_shape_and_columns() -> None:
    """A per-T SCPH .result parses to 525 rows with the documented schema."""
    df = read_result(RAW_DIR / "STO_RTA_scph_300K.result")
    assert list(df.columns) == list(RESULT_COLUMNS)
    assert len(df) == 525
    assert df["T_K"].unique().tolist() == [300.0]
    # Star sizes sum to the full 8x8x8 mesh for each branch.
    assert df[df["branch"] == 1]["multiplicity"].sum() == 512


@needs_raw
def test_read_result_bare_full_grid() -> None:
    """The bare production .result holds all 21 temperatures x 525 modes."""
    df = read_result(RAW_DIR / "STO_RTA_production.result")
    assert len(df) == 525 * 21
    assert df["T_K"].nunique() == 21


@needs_raw
def test_read_result_velocities_full_star() -> None:
    """Velocities cover every star member: 512 k-points x 15 branches."""
    vel = read_result_velocities(RAW_DIR / "STO_RTA_scph_300K.result")
    assert len(vel) == 512 * 15


@needs_raw
def test_read_kl_tensor_and_isotropy() -> None:
    """.kl parses to the 9-component tensor and collapses isotropically."""
    kl = read_kl(RAW_DIR / "STO_RTA_scph_300K.kl")
    assert list(kl.columns) == list(KL_COLUMNS)
    kappa = kappa_isotropic(kl)
    assert kappa["kappa_W_mK"].item() == pytest.approx(9.5283, abs=1e-4)


@needs_raw
def test_emit_round_trip_matches_committed(tmp_path: Path) -> None:
    """Re-emitting the tidy CSVs from raw reproduces the committed files."""
    paths = emit_tidy_csvs(RAW_DIR, tmp_path)
    for name, path in paths.items():
        committed = (TIDY_DIR / path.name).read_text()
        assert path.read_text() == committed, name
