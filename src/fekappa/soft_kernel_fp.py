"""Reduced soft-mode kernel kappa_soft^P(T) from first-principles inputs.

Re-evaluates the manuscript's reduced kernel (old Eq. 4.4, machinery in
:func:`fekappa.kappa_kernel.K_soft`) with the computed SCPH quantities --
Omega_0(T) and Lambda^2(T) from the small-q fit (``lambda_sq_T.csv``),
Gamma_0(T) from the zone-center TO1 triplet of the ``.result`` files
(``eta_validity_T.csv``) -- replacing the previous literature-only inputs
(Vogt Omega_0/eta + Verdi effective curvature; committed record in
``section4_eta_kappa_triples.csv``: 0.137 [0.102, 0.167] W/mK at ~300 K).

Two variants are emitted per METHODS_AND_CAVEATS.md C8 (no rescaling anywhere;
the second evaluation is labeled, not a correction):

* ``computed``           -- all inputs computed (fixed-cell SCPH Omega_0);
                            cell volume = the run cell (7.363 bohr)^3.
* ``vogt_anchored``      -- Vogt 1995 cubic-F1u Omega_0(T) with computed
                            Gamma_0(T) and Lambda^2(T); exp cell volume
                            (59.55 Ang^3) as in the old kernel. Only where
                            Vogt's cubic data overlaps (150-300 K grid).

Each variant is evaluated with the Debye-sphere cutoff (the old kernel's
convention, bracketed by the Lambda^2 +/-30% envelope) and with the small-q
cutoffs q_c in {0.10, 0.15, 0.20} 1/Ang.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from fekappa.damping_map import vogt_omega_interp
from fekappa.kappa_decomposition import V_CELL_ANG3
from fekappa.kappa_kernel import K_soft

CM1_TO_MEV: float = 1.0 / 8.06554

# Old kernel's cell volume (experimental a = 3.905 Ang) and Lambda^2 envelope.
V_CELL_EXP_ANG3: float = 59.55
LAMBDA_SQ_FRAC_UNC: float = 0.30

Q_C_ENVELOPE_ANG1: tuple[float, ...] = (0.10, 0.15, 0.20)

_HEADER_KERNEL = (
    "# source: data/processed/sto_alamode/lambda_sq_T.csv (Omega_0, Lambda^2), "
    "eta_validity_T.csv (Gamma_0, Vogt anchor); kernel machinery "
    "fekappa.kappa_kernel.K_soft (old Eq. 4.4)\n"
    "# method: variant 'computed' = fully computed fixed-cell SCPH inputs, "
    "V_cell = (7.363 bohr)^3; variant 'vogt_anchored' = Vogt 1995 cubic-F1u "
    "Omega_0(T) + computed Gamma_0(T), Lambda^2(T), V_cell = 59.55 Ang^3 "
    "(labeled second evaluation per METHODS_AND_CAVEATS.md C8, NOT a "
    "correction); q_max = 'debye' (old convention, with Lambda^2 +/-30% "
    "lo/hi bracket) or the q_c cutoff in 1/Ang\n"
    "# units: kappa in W m^-1 K^-1; Omega_0, Gamma_0 in meV; Lambda^2 in "
    "meV^2 Ang^2\n"
    "# literature-input reference for the before/after comparison: "
    "0.13712 [0.10217, 0.16662] W/mK at 299.3 K "
    "(data/processed/section4_eta_kappa_triples.csv, STO row)\n"
    "# provenance: derived (ext-FP-alamode inputs; vogt_anchored variant "
    "mixes a literature Omega_0)\n"
)


def kernel_rows(
    T_K: float,
    omega0_meV: float,
    gamma0_meV: float,
    lambda_sq: float,
    variant: str,
    V_cell_ang3: float,
) -> list[dict]:
    """Kernel evaluations (Debye + q_c envelope) for one (T, variant)."""
    eta0 = gamma0_meV / omega0_meV
    common = {
        "T_K": T_K,
        "variant": variant,
        "Omega0_meV": omega0_meV,
        "Gamma0_meV": gamma0_meV,
        "eta0": eta0,
        "lambda_sq_meV2_ang2": lambda_sq,
    }
    rows = []
    kw = {"warn": False}
    kappa = K_soft(omega0_meV**2, lambda_sq, eta0, T_K, V_cell_ang3, **kw)
    lo = K_soft(
        omega0_meV**2, lambda_sq * (1 - LAMBDA_SQ_FRAC_UNC), eta0, T_K, V_cell_ang3, **kw
    )
    hi = K_soft(
        omega0_meV**2, lambda_sq * (1 + LAMBDA_SQ_FRAC_UNC), eta0, T_K, V_cell_ang3, **kw
    )
    rows.append(
        {**common, "q_max": "debye", "kappa_soft_W_mK": kappa,
         "kappa_soft_lo": min(lo, hi), "kappa_soft_hi": max(lo, hi)}
    )
    for q_c in Q_C_ENVELOPE_ANG1:
        rows.append(
            {
                **common,
                "q_max": f"{q_c:.2f}",
                "kappa_soft_W_mK": K_soft(
                    omega0_meV**2, lambda_sq, eta0, T_K, V_cell_ang3,
                    q_max_ang1=q_c, **kw
                ),
                "kappa_soft_lo": np.nan,
                "kappa_soft_hi": np.nan,
            }
        )
    return rows


def emit_soft_kernel_csv(
    lambda_csv: Path, validity_csv: Path, out_path: Path
) -> pd.DataFrame:
    """Emit ``soft_kernel_T.csv`` with both variants over their valid T-grids."""
    lam = pd.read_csv(lambda_csv, comment="#")
    val = pd.read_csv(validity_csv, comment="#")
    rows: list[dict] = []
    for _, r in lam.iterrows():
        t = r["T_K"]
        if t <= 0:
            continue  # kernel vanishes at T = 0 (C = 0); skip the trivial row
        gamma0_meV = float(
            val.loc[val["T_K"] == t, "Gamma0_hwhm_cm1"].item()
        ) * CM1_TO_MEV
        rows += kernel_rows(
            t, r["Omega0_meV"], gamma0_meV, r["lambda_sq_meV2_ang2"],
            "computed", V_CELL_ANG3,
        )
        omega_vogt_cm1 = vogt_omega_interp(t)
        if np.isfinite(omega_vogt_cm1):
            rows += kernel_rows(
                t, float(omega_vogt_cm1) * CM1_TO_MEV, gamma0_meV,
                r["lambda_sq_meV2_ang2"], "vogt_anchored", V_CELL_EXP_ANG3,
            )
    df = pd.DataFrame(rows)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(_HEADER_KERNEL)
        df.to_csv(fh, index=False)
    return df


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[2]
    processed = repo_root / "data" / "processed" / "sto_alamode"
    out = emit_soft_kernel_csv(
        processed / "lambda_sq_T.csv",
        processed / "eta_validity_T.csv",
        processed / "soft_kernel_T.csv",
    )
    print(f"soft_kernel: {(processed / 'soft_kernel_T.csv').relative_to(repo_root)}")
    at300 = out[(out["T_K"] == 300.0) & (out["q_max"] == "debye")]
    print(at300.to_string(index=False))