"""Parsers for ALAMODE ``anphon`` RTA output under ``data/raw/alamode_sto/``.

The raw inputs are the ``.result`` / ``.kl`` files produced by ``anphon`` in
``MODE = RTA`` for cubic SrTiO3 on the 8x8x8 mesh (35 irreducible q-points x
15 branches = 525 modes):

* ``STO_RTA_scph_<T>K.result`` -- SCPH-coupled runs (one temperature each;
  ``FC2XML`` = dfc2-renormalized harmonic fc2 at that temperature).
* ``STO_RTA_production.result`` -- bare-RTA run (21 temperatures, 0-1000 K).
  SUPERSEDED at low frequency by the SCPH-coupled files; kept as the
  negative-control dataset only.
* ``*.kl`` -- thermal-conductivity tensor vs temperature.

Conventions (pinned by the Tadano cross-check recorded in
``data/raw/alamode_sto/PROVENANCE.md``):

* ``#GAMMA_EACH`` values are anharmonic linewidths **Gamma = HWHM in cm^-1**.
  This module stores them verbatim (``gamma_hwhm_cm1``); any unit conversion
  happens exactly once, at the consumption layer.
* ``tau = 1/(2*Gamma)``: the lifetime helper here restores the angular-
  frequency factor (``2*pi*c``) and is the ONLY place that conversion lives.
* Negative ``omega_cm1`` marks an unstable (imaginary-frequency) mode, as
  written by ``anphon`` (bare file: 11/525 at low T; SCPH files: none).

``.result`` layout parsed here: a ``#KPOINT`` block (irreducible q list with
fractional coordinates and weights), a ``#TEMPERATURE`` block (``TMIN TMAX
DT``), a ``##Phonon Frequency`` section (``iq branch omega_cm1``), and one
``#GAMMA_EACH`` block per mode: ``iq branch`` / ``multiplicity`` /
``multiplicity`` lines of group velocity (m/s, one per member of the q star) /
one Gamma value per temperature in the grid.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

# Speed of light in cm/s: converts a wavenumber in cm^-1 to an angular
# frequency via omega = 2*pi*c*nu.
_C_CM_PER_S: float = 2.99792458e10

RESULT_COLUMNS: tuple[str, ...] = (
    "q_index",
    "q_frac_x",
    "q_frac_y",
    "q_frac_z",
    "branch",
    "omega_cm1",
    "gamma_hwhm_cm1",
    "T_K",
    "multiplicity",
)

KL_COLUMNS: tuple[str, ...] = (
    "T_K",
    "k_xx", "k_xy", "k_xz",
    "k_yx", "k_yy", "k_yz",
    "k_zx", "k_zy", "k_zz",
)


def tau_ps_from_gamma_hwhm_cm1(gamma_hwhm_cm1):
    """Lifetime ``tau = 1/(2*Gamma)`` in ps from a HWHM linewidth in cm^-1.

    The wavenumber -> angular-frequency factor (``2*pi*c``) is restored HERE,
    exactly once. Do not pre-convert the input (that would double-apply it).
    Zero linewidths (e.g. the T = 0 rows of the bare file) map to ``inf``.
    """
    gamma = np.asarray(gamma_hwhm_cm1, dtype=float)
    gamma_rad_per_s = 2.0 * np.pi * _C_CM_PER_S * gamma
    with np.errstate(divide="ignore"):
        return 1e12 / (2.0 * gamma_rad_per_s)


def _parse_blocks(lines: list[str]) -> dict[str, list[str]]:
    """Split ``#TAG ... #END TAG`` blocks of the General-information header."""
    blocks: dict[str, list[str]] = {}
    tag: str | None = None
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("##"):
            continue
        if stripped.startswith("#END"):
            tag = None
        elif stripped.startswith("#"):
            tag = stripped.lstrip("#").strip()
            blocks[tag] = []
        elif tag is not None:
            blocks[tag].append(stripped)
    return blocks


def _temperature_grid(tmin: float, tmax: float, dt: float) -> np.ndarray:
    """Temperature grid as anphon builds it from ``TMIN TMAX DT`` (inclusive)."""
    n = round((tmax - tmin) / dt) + 1 if dt > 0 else 1
    return tmin + dt * np.arange(n)


def read_result(path: Path) -> pd.DataFrame:
    """Parse an ``anphon`` RTA ``.result`` file into a tidy mode table.

    Args:
        path: Path to a ``.result`` file (SCPH-coupled per-T file or the bare
            21-temperature production file).

    Returns:
        A DataFrame with the columns in :data:`RESULT_COLUMNS`; one row per
        (irreducible q-point, branch, temperature). ``omega_cm1`` is the
        (SCPH-renormalized or bare) frequency, ``gamma_hwhm_cm1`` the
        anharmonic HWHM linewidth, both verbatim in cm^-1. ``multiplicity``
        is the size of the q star (weights in the file satisfy
        ``weight = multiplicity / nk_total``).

    Raises:
        ValueError: If the mode count, star weights, or per-block Gamma count
            are inconsistent with the file's own header.
    """
    path = Path(path)
    lines = path.read_text().splitlines()

    # --- General information header ---------------------------------------
    end_general = lines.index("##END General information")
    blocks = _parse_blocks(lines[:end_general])

    kpoint = blocks["KPOINT"]
    nk_grid = tuple(int(v) for v in kpoint[0].split())
    nk_total = int(np.prod(nk_grid))
    n_irred = int(kpoint[1])
    q_frac = np.empty((n_irred, 3))
    q_weight = np.empty(n_irred)
    for row, line in enumerate(kpoint[2 : 2 + n_irred]):
        _, coords = line.split(":")
        vals = [float(v) for v in coords.split()]
        q_frac[row] = vals[:3]
        q_weight[row] = vals[3]

    tmin, tmax, dt = (float(v) for v in blocks["TEMPERATURE"][0].split())
    temps = _temperature_grid(tmin, tmax, dt)

    # --- ##Phonon Frequency section ----------------------------------------
    i_freq = lines.index("##Phonon Frequency")
    i_freq_end = lines.index("##END Phonon Frequency")
    omega = {}
    for line in lines[i_freq + 1 : i_freq_end]:
        if line.strip().startswith("#"):
            continue
        iq_s, branch_s, omega_s = line.split()
        omega[(int(iq_s), int(branch_s))] = float(omega_s)
    n_branch = max(branch for _, branch in omega)
    if len(omega) != n_irred * n_branch:
        raise ValueError(
            f"{path.name}: expected {n_irred} x {n_branch} frequency rows, got {len(omega)}"
        )

    # --- #GAMMA_EACH blocks -------------------------------------------------
    records: list[tuple] = []
    i = i_freq_end + 1
    n_lines = len(lines)
    while i < n_lines:
        if lines[i].strip() != "#GAMMA_EACH":
            i += 1
            continue
        iq, branch = (int(v) for v in lines[i + 1].split())
        multiplicity = int(lines[i + 2])
        expected_mult = round(q_weight[iq - 1] * nk_total)
        if multiplicity != expected_mult:
            raise ValueError(
                f"{path.name}: q {iq} star size {multiplicity} != header weight "
                f"{q_weight[iq - 1]} x {nk_total} k-points"
            )
        i_gamma = i + 3 + multiplicity  # skip the per-star velocity lines
        gammas = [float(lines[i_gamma + t]) for t in range(len(temps))]
        if lines[i_gamma + len(temps)].strip() != "#END GAMMA_EACH":
            raise ValueError(
                f"{path.name}: mode ({iq},{branch}) does not hold exactly "
                f"{len(temps)} Gamma values"
            )
        for t_k, gamma in zip(temps, gammas):
            records.append(
                (iq, *q_frac[iq - 1], branch, omega[(iq, branch)], gamma, t_k, multiplicity)
            )
        i = i_gamma + len(temps) + 1

    df = pd.DataFrame.from_records(records, columns=RESULT_COLUMNS)
    expected_rows = n_irred * n_branch * len(temps)
    if len(df) != expected_rows:
        raise ValueError(
            f"{path.name}: parsed {len(df)} rows, expected {expected_rows} "
            f"({n_irred} q x {n_branch} branches x {len(temps)} temperatures)"
        )
    df = df.sort_values(["T_K", "q_index", "branch"]).reset_index(drop=True)
    df["q_index"] = df["q_index"].astype(int)
    df["branch"] = df["branch"].astype(int)
    df["multiplicity"] = df["multiplicity"].astype(int)
    return df


def read_result_velocities(path: Path) -> pd.DataFrame:
    """Extract per-mode group velocities from a ``.result`` file.

    Each ``#GAMMA_EACH`` block carries one velocity vector (m/s) per member of
    the irreducible q-point's star; all members are returned so tensor sums
    (e.g. ``sum v (x) v``) can run over the full star.

    Returns:
        A DataFrame ``[q_index, branch, star_index, v_x, v_y, v_z]`` with
        ``star_index`` counting the star members of each q-point from 1.
    """
    path = Path(path)
    lines = path.read_text().splitlines()
    records: list[tuple] = []
    i = 0
    n_lines = len(lines)
    while i < n_lines:
        if lines[i].strip() != "#GAMMA_EACH":
            i += 1
            continue
        iq, branch = (int(v) for v in lines[i + 1].split())
        multiplicity = int(lines[i + 2])
        for member in range(multiplicity):
            vx, vy, vz = (float(v) for v in lines[i + 3 + member].split())
            records.append((iq, branch, member + 1, vx, vy, vz))
        i += 3 + multiplicity
    return pd.DataFrame.from_records(
        records, columns=["q_index", "branch", "star_index", "v_x", "v_y", "v_z"]
    )


def read_kl(path: Path) -> pd.DataFrame:
    """Parse an ``anphon`` ``.kl`` thermal-conductivity file.

    Returns:
        A DataFrame with the columns in :data:`KL_COLUMNS`: temperature and
        the nine Cartesian components of the kappa tensor in W/mK, one row
        per temperature.
    """
    df = pd.read_csv(
        path, sep=r"\s+", comment="#", header=None, names=KL_COLUMNS
    )
    return df.astype(float)


def kappa_isotropic(kl: pd.DataFrame, rtol: float = 1e-6) -> pd.DataFrame:
    """Collapse a cubic kappa tensor table to ``[T_K, kappa_W_mK]``.

    Verifies the tensor is diagonal and isotropic (as cubic symmetry demands)
    before collapsing; raises ``ValueError`` otherwise rather than averaging
    over a symmetry violation.
    """
    diag = kl[["k_xx", "k_yy", "k_zz"]].to_numpy()
    off = kl[["k_xy", "k_xz", "k_yx", "k_yz", "k_zx", "k_zy"]].to_numpy()
    scale = np.maximum(np.abs(diag).max(axis=1), 1.0)
    if np.any(np.abs(off).max(axis=1) > rtol * scale):
        raise ValueError("kappa tensor has non-zero off-diagonal components")
    if np.any((diag.max(axis=1) - diag.min(axis=1)) > rtol * scale):
        raise ValueError("kappa tensor diagonal is not isotropic")
    return pd.DataFrame({"T_K": kl["T_K"], "kappa_W_mK": diag.mean(axis=1)})


# ---------------------------------------------------------------------------
# Tidy-CSV emission (data/raw/alamode_sto -> data/processed/sto_alamode)
# ---------------------------------------------------------------------------
SCPH_TEMPS_K: tuple[int, ...] = tuple(range(0, 1001, 50))

_HEADER_COMMON = (
    "# method: ALAMODE 1.5.0 anphon MODE=RTA, 8x8x8 mesh "
    "(35 irreducible q x 15 branches = 525 modes per temperature)\n"
    "# units: omega_cm1 and gamma_hwhm_cm1 in cm^-1; Gamma is HWHM "
    "(tau = 1/(2*Gamma), angular factor restored at consumption); "
    "q_frac fractional; negative omega_cm1 = unstable mode\n"
    "# provenance: ext-FP-alamode (ALAMODE tutorial force constants, "
    "cell 7.363 bohr -- see data/raw/alamode_sto/IMPORT_NOTE.md)\n"
)

_HEADER_SCPH = (
    "# source: data/raw/alamode_sto/STO_RTA_scph_{0..1000}K.result "
    "(SCPH-coupled: FC2XML = dfc2-renormalized fc2 per T, "
    "FCSXML = STO_anharm.xml)\n" + _HEADER_COMMON +
    "# T-grid: 0-1000 K step 50 (21 temperatures, one run each)\n"
)

_HEADER_BARE = (
    "# source: data/raw/alamode_sto/STO_RTA_production.result "
    "(bare-RTA: harmonic fc2 from STO_anharm.xml, no SCPH renormalization)\n"
    "# WARNING -- SUPERSEDED at low omega: the bare harmonic surface is "
    "quantitatively wrong in the low-frequency manifold (kappa(300 K) 3.34 vs "
    "9.53 W/mK; low-omega lifetimes 2-85x short; 11/525 unstable modes). "
    "Use omega_gamma_scph_T.csv for paper-facing quantities; this table is "
    "the negative-control dataset only.\n" + _HEADER_COMMON +
    "# T-grid: 0-1000 K step 50 (21 temperatures, single sweep)\n"
)

_HEADER_KAPPA = (
    "# source: data/raw/alamode_sto/STO_RTA_scph_{0..1000}K.kl and "
    ".kl_coherent (SCPH-coupled anphon RTA, KAPPA_COHERENT = 2)\n"
    "# method: ALAMODE 1.5.0 anphon MODE=RTA, 8x8x8 mesh; kappa tensor "
    "verified diagonal and isotropic before collapsing; kappa_coherent is "
    "the Wigner coherent (interbranch) part\n"
    "# units: kappa_W_mK and kappa_coherent_W_mK in W m^-1 K^-1\n"
    "# provenance: ext-FP-alamode (ALAMODE tutorial force constants -- see "
    "data/raw/alamode_sto/IMPORT_NOTE.md)\n"
    "# T-grid: 0-1000 K step 50 (21 temperatures)\n"
)


def _write_csv(df: pd.DataFrame, path: Path, header: str) -> None:
    with open(path, "w") as fh:
        fh.write(header)
        df.to_csv(fh, index=False)


def emit_tidy_csvs(raw_dir: Path, out_dir: Path) -> dict[str, Path]:
    """Emit the committed tidy CSVs from the raw ALAMODE files.

    Writes ``omega_gamma_scph_T.csv`` (7 temperatures), ``omega_gamma_bare_T.csv``
    (21 temperatures, superseded-at-low-omega warning in the header) and
    ``kappa_rta_T.csv`` (SCPH-coupled kappa) into ``out_dir``.

    Returns:
        Mapping of logical name -> written path.
    """
    raw_dir, out_dir = Path(raw_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    scph = pd.concat(
        [read_result(raw_dir / f"STO_RTA_scph_{t}K.result") for t in SCPH_TEMPS_K],
        ignore_index=True,
    )
    bare = read_result(raw_dir / "STO_RTA_production.result")
    kappa_parts = []
    for t in SCPH_TEMPS_K:
        part = kappa_isotropic(read_kl(raw_dir / f"STO_RTA_scph_{t}K.kl"))
        kl_c = read_kl_coherent(raw_dir / f"STO_RTA_scph_{t}K.kl_coherent")
        part["kappa_coherent_W_mK"] = kl_c[["k_xx", "k_yy", "k_zz"]].mean(axis=1)
        kappa_parts.append(part)
    kappa = pd.concat(kappa_parts, ignore_index=True)

    paths = {
        "scph": out_dir / "omega_gamma_scph_T.csv",
        "bare": out_dir / "omega_gamma_bare_T.csv",
        "kappa": out_dir / "kappa_rta_T.csv",
    }
    _write_csv(scph, paths["scph"], _HEADER_SCPH)
    _write_csv(bare, paths["bare"], _HEADER_BARE)
    _write_csv(kappa, paths["kappa"], _HEADER_KAPPA)
    return paths


_HEADER_MESH = (
    "# source: data/raw/alamode_sto/STO_RTA_scph_300K{,_mesh12,_mesh12_iso}."
    "{result,kl,kl_coherent} (SCPH-coupled anphon RTA at 300 K, 8x8x8 vs "
    "12x12x12 q-mesh; _iso adds ISOTOPE=1 natural-abundance scattering)\n"
    "# method: kappa from .kl (isotropic collapse); kappa_C from .kl_coherent; "
    "eta = gamma_hwhm_cm1 / omega_cm1 per mode, fractions unweighted by kappa\n"
    "# units: kappa in W m^-1 K^-1; eta dimensionless\n"
    "# provenance: ext-FP-alamode (ALAMODE tutorial force constants -- see "
    "data/raw/alamode_sto/IMPORT_NOTE.md)\n"
    "# T-grid: 300 K only (convergence check)\n"
)


def read_kl_coherent(path: Path) -> pd.DataFrame:
    """Parse an ``anphon`` ``.kl_coherent`` file to ``[T_K, k_xx, k_yy, k_zz]``."""
    df = pd.read_csv(
        path, sep=r"\s+", comment="#", header=None, names=["T_K", "k_xx", "k_yy", "k_zz"]
    )
    return df.astype(float)


def _eta_summary(df: pd.DataFrame) -> dict[str, float]:
    """Damping-ratio summary of one .result table (positive-frequency modes)."""
    stable = df[df["omega_cm1"] > 0]
    eta = stable["gamma_hwhm_cm1"] / stable["omega_cm1"]
    return {
        "n_modes": len(df),
        "frac_eta_gt_0p5": float((eta > 0.5).mean()),
        "frac_eta_gt_1": float((eta > 1.0).mean()),
        "median_eta": float(eta.median()),
        "max_eta": float(eta.max()),
    }


def emit_mesh_convergence_csv(raw_dir: Path, out_path: Path) -> pd.DataFrame:
    """Write the 8x8x8-vs-12x12x12 SCPH-coupled 300 K convergence table."""
    raw_dir, out_path = Path(raw_dir), Path(out_path)
    rows = []
    stems = [("8x8x8", "STO_RTA_scph_300K"), ("12x12x12", "STO_RTA_scph_300K_mesh12")]
    if (raw_dir / "STO_RTA_scph_300K_mesh12_iso.kl").exists():
        # optional isotope-scattering check (ISOTOPE = 1, natural abundances)
        stems.append(("12x12x12_iso", "STO_RTA_scph_300K_mesh12_iso"))
    for mesh, stem in stems:
        kappa = kappa_isotropic(read_kl(raw_dir / f"{stem}.kl"))["kappa_W_mK"].item()
        kl_c = read_kl_coherent(raw_dir / f"{stem}.kl_coherent")
        kappa_c = float(kl_c[["k_xx", "k_yy", "k_zz"]].mean(axis=1).item())
        rows.append(
            {
                "mesh": mesh,
                "T_K": 300.0,
                "kappa_W_mK": kappa,
                "kappa_coherent_W_mK": kappa_c,
                **_eta_summary(read_result(raw_dir / f"{stem}.result")),
            }
        )
    df = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(_HEADER_MESH)
        df.to_csv(fh, index=False)
    return df


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[2]
    raw = repo_root / "data" / "raw" / "alamode_sto"
    processed = repo_root / "data" / "processed" / "sto_alamode"
    written = emit_tidy_csvs(raw, processed)
    for name, out_path in written.items():
        print(f"{name}: {out_path.relative_to(repo_root)}")
    emit_mesh_convergence_csv(raw, processed / "mesh_convergence.csv")
    print(f"mesh: {(processed / 'mesh_convergence.csv').relative_to(repo_root)}")
