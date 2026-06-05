"""Soft-mode particle-like thermal-conductivity kernel ``K_soft`` (Turn 4B §9).

The diagonal (particle-like / populations) contribution of the single soft
optical branch, in the relaxation-time approximation -- the ``K_P`` term of the
Wigner decomposition (Simoncelli, Marzari & Mauri). For an isotropic soft branch

    Omega(q)^2 = Omega_0^2 + Lambda^2 q^2                         (meV^2)

the mode-resolved RTA conductivity is the Hardy energy-current expression

    kappa = (1 / V_crystal) sum_q C_q v_q^2 tau_q
          = (1 / (2 pi)^3) integral_BZ d^3q  C(Omega_q) <v^2>_dir tau_q ,

with the angular average ``<v_x^2> = v^2 / 3`` (isotropic branch) and the BZ
replaced by a Debye sphere holding exactly one mode per primitive cell,

    (4/3) pi q_max^3 = (2 pi)^3 / V_cell   =>   q_max = (6 pi^2 / V_cell)^(1/3).

Collapsing the angular integral leaves a 1-D radial quadrature

    kappa_soft = (1 / 6 pi^2) integral_0^{q_max} q^2 C(Omega_q) v_q^2 tau_q dq .

Per-factor units (the unit-check the brief asks for; SI unless noted):

    q           [1/m]    radial wavevector (quadrature done in 1/Ang, * 1e30 -> 1/m^3)
    C(Omega_q)  [J/K]    mode heat capacity  k_B x^2 e^x/(e^x-1)^2,  x = Omega/k_B T
    v_q         [m/s]    group velocity  (1/hbar) dOmega/dq,  dOmega/dq = Lambda^2 q / Omega_q
    tau_q       [s]      lifetime  hbar / (2 Gamma_q),  Gamma_q = eta * Omega_q  (HWHM)
    ------------------------------------------------------------------------------
    q^2 dq C v^2 tau                          [1/m^3][J/K][m^2/s^2][s] = W/(m K)  OK

``Gamma_q = eta * Omega_q`` uses the dimensionless damping ``eta`` from
:mod:`fekappa.eta` (HWHM/Omega). ``tau = hbar/(2 Gamma_HWHM)`` is the RTA lifetime
(spectral-function FWHM = 2 Gamma_HWHM, tau = 1/FWHM_angular).

Validity:
  * STO, KTaO3 -- underdamped (eta ~ 0.11-0.13 << 0.5): K_P is meaningful.
  * BTO -- overdamped (eta ~ 2.3): the quasiparticle (well-defined v_q, tau_q)
    assumption FAILS. ``K_soft`` still returns a number but emits an
    ``OverdampedModeWarning`` -- treat it as order-of-magnitude only, pending the
    Wigner-coherent (``K_C``) treatment.
  * Omega_0^2 < 0 (unstable soft mode -- below T_c, Schmidt ST-18O, Ehsan 0 K):
    that temperature is skipped (kappa = NaN) with a warning; no imaginary v_q.
"""

from __future__ import annotations

import warnings

import numpy as np

# --- physical constants (CODATA) -------------------------------------------
HBAR_J_S = 1.054571817e-34      # reduced Planck constant [J s]
KB_J_K = 1.380649e-23           # Boltzmann constant [J/K]
KB_MEV_K = 8.617333262e-2       # Boltzmann constant [meV/K]
MEV_TO_J = 1.602176634e-22      # 1 meV in J
ANG_TO_M = 1.0e-10              # 1 Angstrom in m

# Group-velocity conversion: (Lambda^2 q / Omega) is in meV*Ang; divide by hbar
# after going to SI to get m/s.  v[m/s] = (dOmega/dq)[meV*Ang] * _V_CONV.
_V_CONV = MEV_TO_J * ANG_TO_M / HBAR_J_S          # ~151.93 (m/s)/(meV*Ang)
# Radial measure q^2 dq is evaluated in Ang^-3; convert to m^-3.
_Q3_ANG_TO_SI = 1.0 / ANG_TO_M**3                 # 1e30

# np.trapz was renamed np.trapezoid in NumPy 2.0; support both.
_trapz = getattr(np, "trapezoid", None) or np.trapz


class OverdampedModeWarning(UserWarning):
    """The soft mode is overdamped (eta >= threshold); K_P is order-of-magnitude only."""


def n_bose(omega_meV, T_K):
    """Bose-Einstein occupation ``n_B = 1/(exp(Omega/k_B T) - 1)``, analytic & vectorized.

    ``omega_meV`` and ``T_K`` broadcast. Returns 0 where ``T <= 0`` or the argument
    overflows ``exp`` (n_B is then negligible). No pre-tabulation.
    """
    omega = np.asarray(omega_meV, dtype=float)
    T = np.asarray(T_K, dtype=float)
    x = np.where(T > 0.0, np.abs(omega) / (KB_MEV_K * np.where(T > 0.0, T, 1.0)), np.inf)
    x = np.clip(x, 1e-300, 700.0)              # 700 -> exp overflow guard
    out = 1.0 / np.expm1(x)
    return np.where(np.asarray(T) > 0.0, out, 0.0)


def heat_capacity_mode(omega_meV, T_K):
    """Single-mode heat capacity ``C = k_B x^2 e^x/(e^x-1)^2`` [J/K], ``x = Omega/k_B T``.

    Uses the numerically stable form ``x^2 e^{-x}/(1 - e^{-x})^2`` so that
    ``C -> k_B`` as ``x -> 0`` (Dulong-Petit per mode) and ``C -> 0`` as ``x -> inf``
    (frozen-out) without overflow. Vectorized; ``C = 0`` at ``T <= 0``.
    """
    omega = np.asarray(omega_meV, dtype=float)
    T = np.asarray(T_K, dtype=float)
    safeT = np.where(T > 0.0, T, 1.0)
    x = np.abs(omega) / (KB_MEV_K * safeT)
    x = np.clip(x, 1e-300, 700.0)
    ex = np.exp(-x)
    C = KB_J_K * x**2 * ex / np.expm1(-x) ** 2   # expm1(-x) = e^{-x}-1 = -(1-e^{-x})
    return np.where(T > 0.0, C, 0.0)


def debye_qmax(V_cell_ang3: float) -> float:
    """Debye-sphere radius ``q_max = (6 pi^2 / V_cell)^(1/3)`` [1/Ang] (one mode/cell)."""
    return float((6.0 * np.pi**2 / V_cell_ang3) ** (1.0 / 3.0))


def K_soft(
    omega0_sq_meV2,
    lambda_sq_meV2_ang2: float,
    eta,
    T_K,
    V_cell_ang3: float,
    *,
    n_q: int = 512,
    overdamped_threshold: float = 0.5,
    warn: bool = True,
):
    """Particle-like soft-mode thermal conductivity ``kappa_soft(T)`` [W/(m K)].

    Args:
        omega0_sq_meV2: Soft-mode squared energy ``Omega_0^2(T)`` [meV^2], scalar or
            aligned with ``T_K``. Rows with ``Omega_0^2 < 0`` (unstable mode) yield
            ``NaN`` with a warning.
        lambda_sq_meV2_ang2: Soft-branch curvature ``Lambda^2`` [meV^2 Ang^2] (float;
            a single isotropic branch). For BTO's directional dict use
            :func:`K_soft_sweep`.
        eta: Dimensionless damping ``Gamma_HWHM/Omega`` (scalar or per-``T``). From
            :func:`fekappa.eta.eta_fit`. ``eta >= overdamped_threshold`` warns.
        T_K: Temperature(s) [K]. Vectorized.
        V_cell_ang3: Primitive-cell volume [Ang^3]; sets the Debye cutoff ``q_max``.
        n_q: Radial quadrature points on ``[0, q_max]``.
        overdamped_threshold: ``eta`` above which the K_P picture is flagged invalid.
        warn: Emit overdamped / skipped-row warnings (set ``False`` to silence).

    Returns:
        ``kappa_soft`` [W/(m K)] -- scalar if ``T_K`` is scalar, else an array over ``T``.
    """
    T = np.atleast_1d(np.asarray(T_K, dtype=float))
    O2 = np.broadcast_to(np.asarray(omega0_sq_meV2, dtype=float), T.shape).astype(float)
    eta_arr = np.broadcast_to(np.asarray(eta, dtype=float), T.shape).astype(float)
    scalar_in = np.ndim(T_K) == 0

    if warn and np.any(eta_arr >= overdamped_threshold):
        emax = float(np.max(eta_arr))
        warnings.warn(
            f"soft mode overdamped (eta_max={emax:.2f} >= {overdamped_threshold}): "
            "the K_P quasiparticle assumption (well-defined v_q, tau_q) fails -- "
            "kappa_soft is ORDER-OF-MAGNITUDE ONLY (use Wigner K_C for the real value).",
            OverdampedModeWarning,
            stacklevel=2,
        )

    unstable = O2 < 0.0
    if warn and np.any(unstable):
        warnings.warn(
            f"skipping {int(unstable.sum())} of {T.size} temperature(s) with "
            "Omega_0^2 < 0 (unstable soft mode -> imaginary frequency): kappa = NaN.",
            RuntimeWarning,
            stacklevel=2,
        )

    qmax = debye_qmax(V_cell_ang3)
    q = np.linspace(0.0, qmax, n_q)                       # [1/Ang]
    O2_safe = np.where(unstable, np.nan, O2)

    # 2-D grids (T, q):
    Omega = np.sqrt(O2_safe[:, None] + lambda_sq_meV2_ang2 * q[None, :] ** 2)   # [meV]
    dOmega_dq = lambda_sq_meV2_ang2 * q[None, :] / Omega                        # [meV*Ang]
    v = dOmega_dq * _V_CONV                                                     # [m/s]
    Gamma = eta_arr[:, None] * Omega                                           # [meV] HWHM
    tau = HBAR_J_S / (2.0 * Gamma * MEV_TO_J)                                  # [s]
    C = heat_capacity_mode(Omega, T[:, None])                                 # [J/K]

    integrand = q[None, :] ** 2 * C * v**2 * tau                              # [1/Ang^2 ...]
    kappa = (1.0 / (6.0 * np.pi**2)) * _Q3_ANG_TO_SI * _trapz(integrand, q, axis=1)
    kappa = np.where(unstable, np.nan, kappa)
    return float(kappa[0]) if scalar_in else kappa


def K_soft_sweep(
    omega0_sq_meV2,
    lambda_sq,
    eta,
    T_K,
    V_cell_ang3: float,
    **kwargs,
):
    """``K_soft`` over a Lambda^2 that may be a scalar or a directional dict.

    * ``lambda_sq`` float -> single isotropic branch; returns ``kappa_soft(T)``
      (delegates to :func:`K_soft`).
    * ``lambda_sq`` dict ``{direction: Lambda^2}`` -> adds a direction axis and runs
      the soft/stiff bracket (BTO ``{"soft_100": 972, "stiff_110_111": 4750}``);
      returns ``{direction: kappa_soft(T)}``. Higher Lambda^2 (stiff) gives the
      larger group velocity and hence the upper kappa bound; soft is the lower.

    Extra ``kwargs`` pass through to :func:`K_soft`.
    """
    if isinstance(lambda_sq, dict):
        return {
            direction: K_soft(omega0_sq_meV2, float(L2), eta, T_K, V_cell_ang3, **kwargs)
            for direction, L2 in lambda_sq.items()
        }
    return K_soft(omega0_sq_meV2, float(lambda_sq), eta, T_K, V_cell_ang3, **kwargs)


# ---------------------------------------------------------------------------
# Per-material registry (Turn 4B §9 locked inputs)
# ---------------------------------------------------------------------------
#   lambda_sq        [meV^2 Ang^2]  soft-branch curvature
#   lambda_sq_frac_unc             fractional 1-sigma envelope on Lambda^2
#   eta                            dimensionless damping (fekappa.eta cluster)
#   a_ang, V_cell_ang3             cubic lattice constant [Ang] and cell volume [Ang^3]
MATERIALS: dict[str, dict] = {
    "STO": {
        "a_ang": 3.905,
        "V_cell_ang3": 59.55,
        # Verdi 2023 soft branch, 3x3x3 MLFF-SSCHA Gamma(22 cm^-1) -> X-TO1(190 cm^-1):
        # Lambda^2 = (Omega_X^2 - Omega_Gamma^2)/q_X^2 = 888.5 meV^2 Ang^2. This is an
        # EFFECTIVE BRANCH CURVATURE, *not* a rigorous small-q Omega^2 = Omega_0^2 + Lambda^2 q^2
        # fit: the Gamma->X span (q_X = 0.785 1/Ang) is outside the small-q regime (the same
        # reason Ehsan zone-boundary points were rejected for the BTO refit in Wave 1).
        # Verdi resolves only Gamma/X/M/R -- no points in (0, 0.15 1/Ang) -- so no rigorous
        # small-q value exists; STO's soft branch is smooth enough that 889 is a reasonable
        # scale. The curvature is ~T-independent; Omega_0(T) carries the T-dependence
        # (Omega_0(300 K) ~ 85 cm^-1 = 10.5 meV exp). +/-30% covers the small-q-vs-full-branch
        # ambiguity; the 861 (2x2x2) <-> 889 (3x3x3) grid spread is smaller and sits inside it.
        "lambda_sq": 889.0,
        "lambda_sq_frac_unc": 0.30,
        "eta": 0.109,
    },
    "KTO": {
        "a_ang": 4.007,
        "V_cell_ang3": 64.34,
        # No KTaO3 finite-q soft dispersion available -> STO-proxy curvature, wider
        # +/-50% envelope (on top of STO's own +/-30%). eta from Vogt 1995 (clusters w/ STO).
        "lambda_sq": 889.0,
        "lambda_sq_frac_unc": 0.50,
        "eta": 0.126,
    },
    "BTO": {
        "a_ang": 4.013,
        "V_cell_ang3": 64.63,
        # Harada-Axe-Shirane 1971 cubic phase (NB 06): strongly anisotropic.
        # soft [100] = 972 (paper-text A(100)); stiff [110]/[111] = 4750. Each +/-20%.
        # Overdamped (eta ~ 2.3) -> K_soft is order-of-magnitude only (warns).
        "lambda_sq": {"soft_100": 972.0, "stiff_110_111": 4750.0},
        "lambda_sq_frac_unc": 0.20,
        "eta": 2.3,
    },
}


def K_soft_material(material: str, omega0_sq_meV2, T_K, **kwargs):
    """Convenience wrapper: run :func:`K_soft_sweep` with the :data:`MATERIALS` entry.

    Pulls ``lambda_sq``, ``eta`` and ``V_cell_ang3`` for ``material`` (``"STO"`` /
    ``"KTO"`` / ``"BTO"``); the caller supplies ``Omega_0^2(T)`` and ``T``. Explicit
    ``kwargs`` (e.g. ``eta=...``) override the registry default.
    """
    try:
        m = MATERIALS[material]
    except KeyError:
        raise ValueError(f"unknown material {material!r}; choose from {sorted(MATERIALS)}") from None
    eta = kwargs.pop("eta", m["eta"])
    lambda_sq = kwargs.pop("lambda_sq", m["lambda_sq"])
    return K_soft_sweep(omega0_sq_meV2, lambda_sq, eta, T_K, m["V_cell_ang3"], **kwargs)
