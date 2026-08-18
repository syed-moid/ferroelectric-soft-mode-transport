"""Zone-resolved damping map eta_qs(T) = Gamma/Omega and its validity summary.

Consumes the committed ``mode_kappa_T.csv`` (which carries Omega, Gamma_HWHM,
sector labels, and the per-mode kappa_P contributions) and emits:

* ``eta_map_T.csv`` -- eta per mode per temperature, with each mode's
  kappa_P weight;
* ``eta_validity_T.csv`` -- per temperature: fractions of modes beyond the
  narrow-width criterion (eta > 0.5) and overdamped (eta > 1), both
  count-weighted and kappa_P-weight-weighted; min/median/max eta; the
  zone-center soft-mode value reported BOTH ways (computed-Omega and
  Vogt-anchored) where Vogt's cubic data overlaps;
* ``eta_validity_sensitivity.csv`` -- the same validity fractions after
  rescaling the soft-branch Omega by the Vogt/computed zone-center ratio at
  each overlapping temperature (eta -> eta / r). kappa_P weights are kept
  fixed (first-order sensitivity of the classification metric only).

Conventions: eta = Gamma_HWHM / Omega (both in cm^-1, ratio dimensionless);
eta < 0.5 narrow-width, eta > 1 overdamped, in between "beyond the
narrow-width regime". The zone-wide map stays on computed inputs; the
fixed-cell SCPH Omega overestimate (METHODS_AND_CAVEATS.md C8) biases the
soft-sector eta DOWNWARD, which is exactly what the sensitivity table
quantifies.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

CM1_TO_MEV: float = 1.0 / 8.06554

# Literature zone-center benchmark (Vogt 1995, cubic F1u): eta(300 K) = 0.109.
ETA_STO_VOGT_300K: float = 0.109


def load_mode_table(mode_kappa_csv: Path) -> pd.DataFrame:
    """Load the committed mode table and attach eta and kappa weights."""
    df = pd.read_csv(mode_kappa_csv, comment="#")
    stable = df["omega_cm1"] > 1e-3
    df = df[stable].copy()
    df["eta"] = df["gamma_hwhm_cm1"] / df["omega_cm1"]
    kappa_tot = df.groupby("T_K")["kappa_mode_W_mK"].transform("sum")
    with np.errstate(invalid="ignore"):
        df["kappa_weight"] = np.where(
            kappa_tot > 0, df["kappa_mode_W_mK"] / kappa_tot, 0.0
        )
    return df


def zone_center_soft(df_t: pd.DataFrame) -> tuple[float, float]:
    """(Omega_0, Gamma_0) of the zone-center TO1 triplet (q 1, branches 4-6).

    The ``.result`` zone-center rows carry no non-analytic correction, so the
    soft F1u mode is the lowest optic triplet; Omega is triply degenerate and
    Gamma is averaged over the triplet.
    """
    g = df_t[(df_t["q_index"] == 1) & df_t["branch"].between(4, 6)]
    if len(g) != 3:
        raise ValueError(f"expected the TO1 triplet at q=1, got {len(g)} rows")
    omega = g["omega_cm1"].to_numpy()
    if not np.allclose(omega, omega[0], rtol=1e-6):
        raise ValueError("zone-center branches 4-6 are not degenerate (not TO1?)")
    return float(omega[0]), float(g["gamma_hwhm_cm1"].mean())


def _validity_fractions(eta: np.ndarray, weight: np.ndarray) -> dict[str, float]:
    return {
        "frac_eta_gt_0p5_count": float((eta > 0.5).mean()),
        "frac_eta_gt_1_count": float((eta > 1.0).mean()),
        "frac_eta_gt_0p5_kappa": float(weight[eta > 0.5].sum()),
        "frac_eta_gt_1_kappa": float(weight[eta > 1.0].sum()),
    }


def vogt_omega_interp(T_K, *, edge_tol_K: float = 5.0):
    """Vogt 1995 cubic-F1u Omega(T) [cm^-1] interpolated onto ``T_K``.

    Returns NaN outside Vogt's cubic temperature range (>= 110 K, up to
    ~299 K); ``edge_tol_K`` admits grid points a few kelvin past the last
    datum (e.g. 300 K vs Vogt's 299.3 K endpoint) as endpoint values.
    """
    from fekappa import literature_io as lit

    v = lit.load_vogt1995_sto_omega_T()
    v = v[v["phase_mode"] == "cubic_F1u"].sort_values("T_K")
    T = np.asarray(T_K, dtype=float)
    omega = np.interp(T, v["T_K"], v["Omega_meV"] / CM1_TO_MEV)
    in_range = (T >= v["T_K"].min() - edge_tol_K) & (T <= v["T_K"].max() + edge_tol_K)
    out = np.where(in_range, omega, np.nan)
    return out if out.ndim else float(out)


def validity_table(df: pd.DataFrame) -> pd.DataFrame:
    """Per-temperature validity summary (computed inputs; dual zone-center eta)."""
    rows = []
    for t, df_t in df.groupby("T_K"):
        eta = df_t["eta"].to_numpy()
        weight = df_t["kappa_weight"].to_numpy()
        omega0, gamma0 = zone_center_soft(df_t)
        omega_vogt = vogt_omega_interp(t)
        eta0 = gamma0 / omega0 if t > 0 else np.nan
        rows.append(
            {
                "T_K": t,
                "n_modes": len(df_t),
                **_validity_fractions(eta, weight),
                "eta_min": float(eta.min()),
                "eta_median": float(np.median(eta)),
                "eta_max": float(eta.max()),
                "Omega0_cm1": omega0,
                "Gamma0_hwhm_cm1": gamma0,
                "eta0_computed": eta0,
                "Omega0_vogt_cm1": omega_vogt,
                "eta0_vogt_anchored": gamma0 / omega_vogt,
                "vogt_over_computed": omega_vogt / omega0,
            }
        )
    return pd.DataFrame(rows).sort_values("T_K").reset_index(drop=True)


def sensitivity_table(df: pd.DataFrame) -> pd.DataFrame:
    """Validity fractions with soft-branch Omega rescaled to the Vogt anchor.

    At each temperature where Vogt's cubic data exists, soft-optic-sector modes
    get Omega -> r * Omega with r = Omega_Vogt / Omega_0_computed (zone-center
    ratio), i.e. eta -> eta / r; acoustic and hard-optic modes are untouched.
    """
    rows = []
    for t, df_t in df.groupby("T_K"):
        omega0, _ = zone_center_soft(df_t)
        omega_vogt = vogt_omega_interp(t)
        if not np.isfinite(omega_vogt):
            continue
        r = omega_vogt / omega0
        soft = (df_t["sector"] == "soft_optic").to_numpy()
        eta = df_t["eta"].to_numpy().copy()
        eta_rescaled = np.where(soft, eta / r, eta)
        weight = df_t["kappa_weight"].to_numpy()
        base = _validity_fractions(eta, weight)
        resc = _validity_fractions(eta_rescaled, weight)
        rows.append(
            {
                "T_K": t,
                "rescale_ratio": r,
                "n_soft_modes": int(soft.sum()),
                **{f"{k}_computed": v for k, v in base.items()},
                **{f"{k}_rescaled": v for k, v in resc.items()},
                "eta_max_rescaled": float(eta_rescaled.max()),
            }
        )
    return pd.DataFrame(rows).sort_values("T_K").reset_index(drop=True)


_HEADER_MAP = (
    "# source: data/processed/sto_alamode/mode_kappa_T.csv (SCPH-coupled "
    "ALAMODE dataset, 8x8x8)\n"
    "# method: eta = gamma_hwhm_cm1 / omega_cm1 per mode; kappa_weight = mode "
    "share of kappa_P at its temperature; Gamma-point translations excluded\n"
    "# units: omega_cm1 in cm^-1; eta and kappa_weight dimensionless\n"
    "# provenance: ext-FP-alamode (see data/raw/alamode_sto/IMPORT_NOTE.md)\n"
)

_HEADER_VALIDITY = (
    "# source: eta_map_T.csv inputs (computed Omega, Gamma); Vogt 1995 "
    "cubic-F1u Omega(T) for the anchored zone-center variant\n"
    "# method: narrow-width criterion eta < 0.5, overdamped eta > 1; fractions "
    "count-weighted and kappa_P-weighted; eta0 reported both ways per "
    "METHODS_AND_CAVEATS.md C8 (computed Omega_0 carries the fixed-cell SCPH "
    "overestimate; Vogt columns NaN outside Vogt's cubic range 110-300 K); "
    "above the ~450 K TO1-TO2 anticrossing the zone-center columns describe "
    "the lowest optic triplet, no longer a pure TO1 soft mode\n"
    "# units: Omega, Gamma in cm^-1; eta dimensionless\n"
    "# provenance: ext-FP-alamode + literature benchmark (derived)\n"
)

_HEADER_SENSITIVITY = (
    "# source: eta_map_T.csv inputs; Vogt 1995 cubic-F1u Omega(T)\n"
    "# method: soft-optic-sector Omega rescaled by r = Omega_Vogt/Omega_0_computed "
    "(zone-center ratio) -> eta/r for soft modes only; kappa_P weights held fixed; "
    "quantifies the downward eta bias of the fixed-cell SCPH Omega overestimate "
    "(METHODS_AND_CAVEATS.md C8) on the validity fractions\n"
    "# units: dimensionless\n"
    "# provenance: derived (ext-FP-alamode + literature anchor)\n"
)


def emit_damping_csvs(mode_kappa_csv: Path, out_dir: Path) -> dict[str, Path]:
    """Emit ``eta_map_T.csv``, ``eta_validity_T.csv``, ``eta_validity_sensitivity.csv``."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = load_mode_table(mode_kappa_csv)
    paths = {
        "map": out_dir / "eta_map_T.csv",
        "validity": out_dir / "eta_validity_T.csv",
        "sensitivity": out_dir / "eta_validity_sensitivity.csv",
    }
    with open(paths["map"], "w") as fh:
        fh.write(_HEADER_MAP)
        df[
            ["T_K", "q_index", "branch", "omega_cm1", "gamma_hwhm_cm1",
             "eta", "sector", "kappa_weight"]
        ].to_csv(fh, index=False)
    with open(paths["validity"], "w") as fh:
        fh.write(_HEADER_VALIDITY)
        validity_table(df).to_csv(fh, index=False)
    with open(paths["sensitivity"], "w") as fh:
        fh.write(_HEADER_SENSITIVITY)
        sensitivity_table(df).to_csv(fh, index=False)
    return paths


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[2]
    processed = repo_root / "data" / "processed" / "sto_alamode"
    written = emit_damping_csvs(processed / "mode_kappa_T.csv", processed)
    for name, out_path in written.items():
        print(f"{name}: {out_path.relative_to(repo_root)}")