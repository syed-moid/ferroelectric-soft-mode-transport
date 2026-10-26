"""Parsers for the raw BaTiO3 / SrTiO3 data sets under ``data/raw/``.

The raw inputs are:

* ``BTO/phonon_modes_<T>K.csv`` -- per-temperature Quantum-ESPRESSO phonon
  mode tables for BaTiO3 (six temperatures, 100-600 K).
* ``BTO_ir_data.json`` / ``STO_ir_data.json`` -- infrared spectra, peak
  parameterizations, and a shared frequency grid.
* ``ShiraneFig4a.csv`` / ``ShiraneFig4b.csv`` -- digitized inelastic-neutron
  energy-loss spectra of single-domain **BaTiO3** at 22 C (~295 K) from two
  different (q, e) configurations in Fig. 4 of Shirane, Axe, Harada & Linz,
  Phys. Rev. B 2, 3651 (1970). (Note: corrected from SrTiO3 in v1.)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

# THz -> cm^-1 conversion factor (1 THz = 33.3564 cm^-1).
THZ_TO_CM1: float = 33.3564

# Columns expected in every BTO/phonon_modes_<T>K.csv file.
BTO_PHONON_COLUMNS: tuple[str, ...] = (
    "mode_name",
    "label",
    "omega_Q_meV",
    "Gamma_Q_HWHM_meV",
    "Q_x",
    "Q_y",
    "Q_z",
    "expected_peak_barns",
    "mode_type",
    "temperature_K",
)


def load_bto_phonon_modes(data_dir: Path) -> pd.DataFrame:
    """Load and concatenate all BTO per-temperature phonon-mode tables.

    Reads every ``phonon_modes_*K.csv`` file in ``data_dir`` and stacks them
    into a single tidy :class:`pandas.DataFrame`. Each row is one phonon mode
    at one temperature. A ``temperature_K`` column is guaranteed; when a file
    omits it, the value is recovered from the file name (e.g. ``300`` from
    ``phonon_modes_300K.csv``).

    Args:
        data_dir: Directory containing the ``phonon_modes_*K.csv`` files
            (typically ``data/raw/BTO``).

    Returns:
        A DataFrame with the columns listed in :data:`BTO_PHONON_COLUMNS`,
        sorted by ``temperature_K`` then ``mode_name`` and re-indexed.

    Raises:
        FileNotFoundError: If no matching CSV files are found in ``data_dir``.
    """
    data_dir = Path(data_dir)
    paths = sorted(data_dir.glob("phonon_modes_*K.csv"))
    if not paths:
        raise FileNotFoundError(
            f"No 'phonon_modes_*K.csv' files found in {data_dir!s}"
        )

    frames: list[pd.DataFrame] = []
    for path in paths:
        frame = pd.read_csv(path)
        if "temperature_K" not in frame.columns:
            # Recover temperature from the file name, e.g. phonon_modes_300K.csv
            stem = path.stem  # phonon_modes_300K
            temp_token = stem.split("_")[-1].rstrip("Kk")
            frame["temperature_K"] = int(temp_token)
        frames.append(frame)

    df = pd.concat(frames, ignore_index=True)
    df["temperature_K"] = df["temperature_K"].astype(int)
    df = df.sort_values(["temperature_K", "mode_name"]).reset_index(drop=True)
    return df


def load_ir_spectra(json_path: Path) -> dict:
    """Parse a ``*_ir_data.json`` infrared-spectra file.

    The frequency grid in cm^-1 is recomputed from the THz grid using
    :data:`THZ_TO_CM1` (it agrees with the ``frequency_grid_cm1`` stored in the
    file). The ``spectra`` block is returned sorted by ascending temperature,
    with each ``T_<T>K`` key validated as an integer temperature in K.

    Args:
        json_path: Path to ``BTO_ir_data.json`` or ``STO_ir_data.json``.

    Returns:
        A dict with keys:

        * ``material`` (str): chemical formula as stored in the file.
        * ``freq_THz`` (numpy.ndarray): frequency grid in THz.
        * ``freq_cm1`` (numpy.ndarray): frequency grid in cm^-1
          (``freq_THz * THZ_TO_CM1``).
        * ``peaks`` (list[dict]): per-mode peak parameterizations, verbatim
          from the file (keys include ``frequency_cm1``, ``A_k``, ``B_k``,
          ``C_k``).
        * ``spectra`` (dict[str, numpy.ndarray]): ``T_<T>K`` -> intensity grid,
          ordered by ascending temperature.

    Raises:
        ValueError: If a key in the ``spectra`` block cannot be parsed as an
            integer temperature in K.
    """
    json_path = Path(json_path)
    with open(json_path) as fh:
        data = json.load(fh)

    freq_THz = np.asarray(data["frequency_grid_THz"], dtype=float)
    freq_cm1 = freq_THz * THZ_TO_CM1

    # Validate and sort the spectra block by ascending temperature.
    temp_keyed: list[tuple[int, str]] = []
    for key in data["spectra"]:
        try:
            temp_k = int(key.removeprefix("T_").removesuffix("K"))
        except ValueError as exc:
            raise ValueError(
                f"Cannot parse temperature from spectra key {key!r} in "
                f"{json_path!s}"
            ) from exc
        temp_keyed.append((temp_k, key))

    spectra = {
        key: np.asarray(data["spectra"][key], dtype=float)
        for _, key in sorted(temp_keyed)
    }

    return {
        "material": data["material"],
        "freq_THz": freq_THz,
        "freq_cm1": freq_cm1,
        "peaks": data["peaks"],
        "spectra": spectra,
    }


def load_shirane(csv_path: Path) -> pd.DataFrame:
    """Parse a digitized Shirane inelastic-neutron energy-loss spectrum.

    The source CSVs (``ShiraneFig4a.csv``, ``ShiraneFig4b.csv``) hold two
    columns -- neutron energy loss in meV and a count rate -- under verbose
    digitizer headers. These are renamed to the canonical ``energy_loss_meV``
    and ``counts``.

    Args:
        csv_path: Path to ``ShiraneFig4a.csv`` or ``ShiraneFig4b.csv``.

    Returns:
        A two-column DataFrame ``[energy_loss_meV, counts]`` in file order.
    """
    df = pd.read_csv(csv_path)
    if df.shape[1] != 2:
        raise ValueError(
            f"Expected 2 columns in {csv_path!s}, got {df.shape[1]}: "
            f"{list(df.columns)}"
        )
    df.columns = ["energy_loss_meV", "counts"]
    df = df.astype(float)
    return df
