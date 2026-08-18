"""Mode-resolved RTA conductivity and the five-way sector decomposition.

Rebuilds the particle-like (Peierls/populations) conductivity mode-by-mode from
the SCPH-coupled ALAMODE dataset,

    kappa_P = (1 / (N_k V_cell)) sum_{q in star, s} C_qs |v_qs|^2/3 tau_qs ,

and checks it against the ``anphon`` ``.kl`` value at every temperature -- this
closure is the units gate for everything downstream (C from
:func:`fekappa.kappa_kernel.heat_capacity_mode`, tau = 1/(2*Gamma_HWHM) via
:func:`fekappa.alamode_io.tau_ps_from_gamma_hwhm_cm1`, star velocities from the
``.result`` file).

Sector split (stated per the brief): **omega-band partition with a branch-index
acoustic tag** -- ``acoustic`` = the three lowest branches at each q,
``soft-optic`` = branches >= 4 with Omega below the boundary (default 175 cm^-1,
the partition validated in the companion lattice-dynamics project),
``hard-optic`` = the rest. Eigenvector character is not available in the
``.result`` files; at zone-boundary q where renormalized soft-optic modes dip
below the acoustic branches the index tag can swap labels *within* the low-omega
manifold, so the acoustic+soft-optic **sum** is the partition-stable quantity and
the boundary sensitivity is emitted alongside (150/175/200 cm^-1).

The coherence channel kappa_C comes directly from the ``KAPPA_COHERENT`` runs
(``.kl_coherent``); the five-way table is the three particle sectors + kappa_C
+ their total (Eq. 2.67 of the manuscript).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from fekappa.alamode_io import (
    SCPH_TEMPS_K,
    kappa_isotropic,
    read_kl,
    read_kl_coherent,
    read_result,
    read_result_velocities,
    tau_ps_from_gamma_hwhm_cm1,
)
from fekappa.kappa_kernel import heat_capacity_mode

# Run cell of the ALAMODE dataset: cubic, celldm(1) = 7.363 bohr.
BOHR_TO_ANG: float = 0.529177210903
A_CELL_ANG: float = 7.363 * BOHR_TO_ANG
V_CELL_ANG3: float = A_CELL_ANG**3
NK_MESH: int = 512  # 8x8x8

CM1_TO_MEV: float = 1.0 / 8.06554

# Sector boundary (cm^-1) between soft-optic and hard-optic; validated at ~175
# in the companion project, sensitivity emitted for 150/200.
SECTOR_BOUNDARY_CM1: float = 175.0
SECTOR_BOUNDARIES_CM1: tuple[float, ...] = (150.0, 175.0, 200.0)


def classify_sector(branch, omega_cm1, boundary_cm1: float = SECTOR_BOUNDARY_CM1):
    """Sector label per mode: ``acoustic`` / ``soft_optic`` / ``hard_optic``."""
    branch = np.asarray(branch)
    omega = np.asarray(omega_cm1, dtype=float)
    return np.where(
        branch <= 3,
        "acoustic",
        np.where(omega < boundary_cm1, "soft_optic", "hard_optic"),
    )


def mode_kappa_table(
    result_path: Path,
    *,
    V_cell_ang3: float = V_CELL_ANG3,
    nk: int = NK_MESH,
) -> pd.DataFrame:
    """Per-mode isotropic kappa_P contributions from one ``.result`` file.

    Returns the tidy mode table plus ``v2_star_m2s2`` (star-summed |v|^2) and
    ``kappa_mode_W_mK`` such that ``kappa_P = kappa_mode.sum()`` per temperature.
    """
    df = read_result(result_path)
    vel = read_result_velocities(result_path)
    vel["v2"] = vel[["v_x", "v_y", "v_z"]].pow(2).sum(axis=1)
    v2_star = vel.groupby(["q_index", "branch"])["v2"].sum().rename("v2_star_m2s2")
    df = df.merge(v2_star, on=["q_index", "branch"], validate="many_to_one")

    omega_meV = df["omega_cm1"].to_numpy() * CM1_TO_MEV
    C = heat_capacity_mode(omega_meV, df["T_K"].to_numpy())  # J/K
    tau_s = tau_ps_from_gamma_hwhm_cm1(df["gamma_hwhm_cm1"].to_numpy()) * 1e-12
    v2 = df["v2_star_m2s2"].to_numpy()

    V_m3 = V_cell_ang3 * 1e-30
    with np.errstate(invalid="ignore"):
        contrib = C * tau_s * (v2 / 3.0) / (nk * V_m3)
    # Excluded from the sum: the three Gamma-point translations (omega
    # numerically ~1e-10 cm^-1, Gamma = 0 -> tau infinite), unstable modes
    # (omega < 0, bare file only), and modes with C = 0 or zero velocity.
    excluded = (
        (C == 0.0)
        | (v2 == 0.0)
        | (df["omega_cm1"].to_numpy() < 1e-3)
    )
    contrib = np.where(excluded, 0.0, contrib)
    if not np.all(np.isfinite(contrib)):
        bad = df.loc[~np.isfinite(contrib), ["q_index", "branch", "T_K"]]
        raise ValueError(f"non-finite kappa contributions:\n{bad.head()}")
    df["kappa_mode_W_mK"] = contrib
    df["sector"] = classify_sector(df["branch"], df["omega_cm1"])
    return df


def five_way_row(
    mode_df: pd.DataFrame,
    kappa_coherent: float,
    *,
    boundary_cm1: float = SECTOR_BOUNDARY_CM1,
) -> dict[str, float]:
    """Collapse one temperature's mode table to the five-way decomposition."""
    sector = classify_sector(mode_df["branch"], mode_df["omega_cm1"], boundary_cm1)
    kappa = mode_df["kappa_mode_W_mK"]
    parts = {
        s: float(kappa[sector == s].sum())
        for s in ("acoustic", "soft_optic", "hard_optic")
    }
    kappa_p = float(kappa.sum())
    return {
        "kappa_acoustic_P": parts["acoustic"],
        "kappa_soft_optic_P": parts["soft_optic"],
        "kappa_hard_optic_P": parts["hard_optic"],
        "kappa_C": kappa_coherent,
        "kappa_total": kappa_p + kappa_coherent,
        "kappa_P_modesum": kappa_p,
    }


_HEADER_MODE_KAPPA = (
    "# source: data/raw/alamode_sto/STO_RTA_scph_{0..1000}K.result "
    "(SCPH-coupled anphon RTA, 8x8x8)\n"
    "# method: kappa_mode = C(Omega,T) * tau * (star-summed |v|^2)/3 / (N_k V_cell); "
    "C per-mode heat capacity, tau = 1/(2*Gamma_HWHM); "
    "V_cell = (7.363 bohr)^3, N_k = 512; sum over modes reproduces the anphon .kl "
    "kappa (verified in tests)\n"
    "# units: omega_cm1, gamma_hwhm_cm1 in cm^-1; v2_star_m2s2 in m^2/s^2; "
    "kappa_mode_W_mK in W m^-1 K^-1\n"
    "# sector: acoustic = branches 1-3; soft_optic = branch >= 4 and omega < 175 cm^-1; "
    "hard_optic = rest (boundary sensitivity in kappa_five_way_T.csv)\n"
    "# provenance: ext-FP-alamode (see data/raw/alamode_sto/IMPORT_NOTE.md)\n"
)

_HEADER_FIVE_WAY = (
    "# source: mode_kappa_T.csv (particle sectors) + "
    "data/raw/alamode_sto/STO_RTA_scph_*K.kl_coherent (kappa_C, KAPPA_COHERENT=1 "
    "direct evaluation, not a bound)\n"
    "# method: five-way decomposition kappa_total = kappa_acoustic_P + "
    "kappa_soft_optic_P + kappa_hard_optic_P + kappa_C; kappa_P_modesum and "
    "kappa_P_anphon record the closure of the mode sum against the .kl value; "
    "rows with boundary_cm1 != 175 give the sector-boundary sensitivity\n"
    "# units: kappa in W m^-1 K^-1; boundary_cm1 in cm^-1\n"
    "# provenance: ext-FP-alamode (see data/raw/alamode_sto/IMPORT_NOTE.md)\n"
)


def emit_decomposition_csvs(raw_dir: Path, out_dir: Path) -> dict[str, Path]:
    """Emit ``mode_kappa_T.csv`` and ``kappa_five_way_T.csv`` (all temperatures).

    The five-way table carries one row per (T, boundary) with the central
    boundary 175 cm^-1 plus the 150/200 sensitivity rows.
    """
    raw_dir, out_dir = Path(raw_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    mode_frames = []
    rows = []
    for t in SCPH_TEMPS_K:
        mode_df = mode_kappa_table(raw_dir / f"STO_RTA_scph_{t}K.result")
        mode_frames.append(mode_df)
        kl_c = read_kl_coherent(raw_dir / f"STO_RTA_scph_{t}K.kl_coherent")
        kappa_c = float(kl_c[["k_xx", "k_yy", "k_zz"]].mean(axis=1).item())
        kappa_anphon = kappa_isotropic(read_kl(raw_dir / f"STO_RTA_scph_{t}K.kl"))[
            "kappa_W_mK"
        ].item()
        for boundary in SECTOR_BOUNDARIES_CM1:
            rows.append(
                {
                    "T_K": float(t),
                    "boundary_cm1": boundary,
                    **five_way_row(mode_df, kappa_c, boundary_cm1=boundary),
                    "kappa_P_anphon": kappa_anphon,
                }
            )

    paths = {
        "mode_kappa": out_dir / "mode_kappa_T.csv",
        "five_way": out_dir / "kappa_five_way_T.csv",
    }
    with open(paths["mode_kappa"], "w") as fh:
        fh.write(_HEADER_MODE_KAPPA)
        pd.concat(mode_frames, ignore_index=True).to_csv(fh, index=False)
    with open(paths["five_way"], "w") as fh:
        fh.write(_HEADER_FIVE_WAY)
        pd.DataFrame(rows).to_csv(fh, index=False)
    return paths


_HEADER_CONTROL = (
    "# source: data/raw/alamode_sto/STO_RTA_production.kl (bare-RTA) vs "
    "STO_RTA_scph_*K.kl (SCPH-coupled); tau bands from "
    "data/processed/alamode_crosscheck/alamode_sto_tau300K_bands_{bare,scph}.csv\n"
    "# method: controlled negative result -- identical anharmonic (scattering) "
    "force constants, identical mesh; the ONLY difference is the harmonic part "
    "(bare STO_anharm.xml fc2 vs dfc2-renormalized fc2 per T). record='kappa_T' "
    "rows compare kappa(T); record='tau_band_300K' rows compare band-median "
    "tau(300 K) (the 2-85x low-omega failure of the bare surface)\n"
    "# units: kappa in W m^-1 K^-1; tau in ps; bands in cm^-1\n"
    "# provenance: ext-FP-alamode (see data/raw/alamode_sto/IMPORT_NOTE.md)\n"
)


def emit_bare_vs_scph_csv(
    raw_dir: Path, crosscheck_dir: Path, out_path: Path
) -> pd.DataFrame:
    """Emit ``bare_vs_scph_control.csv``: kappa(T) both ways + tau band ratios."""
    raw_dir, crosscheck_dir = Path(raw_dir), Path(crosscheck_dir)
    bare_kl = kappa_isotropic(read_kl(raw_dir / "STO_RTA_production.kl"))
    rows: list[dict] = []
    for t in SCPH_TEMPS_K:
        scph = kappa_isotropic(read_kl(raw_dir / f"STO_RTA_scph_{t}K.kl"))[
            "kappa_W_mK"
        ].item()
        bare = bare_kl.loc[bare_kl["T_K"] == float(t), "kappa_W_mK"].item()
        rows.append(
            {
                "record": "kappa_T",
                "T_K": float(t),
                "bare": bare,
                "scph": scph,
                "scph_over_bare": scph / bare if bare > 0 else np.nan,
            }
        )
    bands_bare = pd.read_csv(
        crosscheck_dir / "alamode_sto_tau300K_bands_bare.csv", comment="#"
    )
    bands_scph = pd.read_csv(
        crosscheck_dir / "alamode_sto_tau300K_bands_scph.csv", comment="#"
    )
    merged = bands_bare.merge(
        bands_scph, on=["band_lo_cm1", "band_hi_cm1"], suffixes=("_bare", "_scph")
    )
    for _, r in merged.iterrows():
        if not (np.isfinite(r["tau_median_ps_bare"]) and np.isfinite(r["tau_median_ps_scph"])):
            continue
        rows.append(
            {
                "record": "tau_band_300K",
                "band_lo_cm1": r["band_lo_cm1"],
                "band_hi_cm1": r["band_hi_cm1"],
                "bare": r["tau_median_ps_bare"],
                "scph": r["tau_median_ps_scph"],
                "scph_over_bare": r["tau_median_ps_scph"] / r["tau_median_ps_bare"],
            }
        )
    df = pd.DataFrame(rows)
    out_path = Path(out_path)
    with open(out_path, "w") as fh:
        fh.write(_HEADER_CONTROL)
        df.to_csv(fh, index=False)
    return df


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[2]
    written = emit_decomposition_csvs(
        repo_root / "data" / "raw" / "alamode_sto",
        repo_root / "data" / "processed" / "sto_alamode",
    )
    for name, out_path in written.items():
        print(f"{name}: {out_path.relative_to(repo_root)}")
    control_path = repo_root / "data" / "processed" / "sto_alamode" / "bare_vs_scph_control.csv"
    emit_bare_vs_scph_csv(
        repo_root / "data" / "raw" / "alamode_sto",
        repo_root / "data" / "processed" / "alamode_crosscheck",
        control_path,
    )
    print(f"control: {control_path.relative_to(repo_root)}")