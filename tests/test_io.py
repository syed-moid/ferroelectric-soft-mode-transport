"""Smoke tests for the fekappa.io data parsers and IR Lorentzian fitter."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fekappa.io import (
    BTO_PHONON_COLUMNS,
    load_bto_phonon_modes,
    load_ir_spectra,
    load_shirane,
)
from fekappa.selfenergy import fit_ir_lorentzians

DATA_RAW = Path(__file__).resolve().parents[1] / "data" / "raw"
BTO_DIR = DATA_RAW / "BTO"


def test_load_bto_phonon_modes_temps_and_columns() -> None:
    """The BTO loader returns all six temperatures and the expected columns."""
    df = load_bto_phonon_modes(BTO_DIR)
    assert isinstance(df, pd.DataFrame)
    assert sorted(df["temperature_K"].unique()) == [100, 200, 300, 400, 500, 600]
    assert set(BTO_PHONON_COLUMNS).issubset(df.columns)


def test_load_shirane_fig4a_shape() -> None:
    """Shirane Fig 4a parses to a nonempty two-column DataFrame."""
    df = load_shirane(DATA_RAW / "ShiraneFig4a.csv")
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert list(df.columns) == ["energy_loss_meV", "counts"]


def test_soft_mode_rows_sorted_by_temperature() -> None:
    """polar_soft rows are temperature-ordered when sorted on temperature_K."""
    df = load_bto_phonon_modes(BTO_DIR)
    soft = df[df["mode_type"] == "polar_soft"].sort_values("temperature_K")
    temps = soft["temperature_K"].tolist()
    assert temps == sorted(temps)
    # One soft mode per temperature: 100..600 K.
    assert temps == [100, 200, 300, 400, 500, 600]


def test_load_ir_spectra_keys_and_ordering() -> None:
    """load_ir_spectra returns the expected keys, sorted spectra, and grid."""
    ir = load_ir_spectra(DATA_RAW / "STO_ir_data.json")
    assert set(ir.keys()) == {"material", "freq_THz", "freq_cm1", "peaks", "spectra"}
    # Frequency grid monotonically increasing.
    assert np.all(np.diff(ir["freq_THz"]) > 0)
    # Five spectra, temperatures parse to int K and ascend.
    temps = [int(k.removeprefix("T_").removesuffix("K")) for k in ir["spectra"]]
    assert len(temps) == 5
    assert temps == sorted(temps)
    # Peak count matches the raw JSON.
    raw = json.loads((DATA_RAW / "STO_ir_data.json").read_text())
    assert len(ir["peaks"]) == len(raw["peaks"])


def test_fit_ir_lorentzians_synthetic() -> None:
    """fit_ir_lorentzians recovers two synthetic peak centers within 1%."""
    x = np.linspace(0.0, 200.0, 2000)

    def lorentzian(x: np.ndarray, center: float, hwhm: float, amp: float) -> np.ndarray:
        return amp * (hwhm / np.pi) / ((x - center) ** 2 + hwhm**2)

    y = lorentzian(x, 50.0, 3.0, 100.0) + lorentzian(x, 120.0, 5.0, 80.0)
    peaks = fit_ir_lorentzians(x, y, n_peaks=2)
    centers = sorted(p["center_cm1"] for p in peaks)
    assert abs(centers[0] - 50.0) / 50.0 < 0.01
    assert abs(centers[1] - 120.0) / 120.0 < 0.01


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-v"]))
