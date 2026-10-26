"""Soft-mode damping ratio eta = Gamma_HWHM / Omega.

Canonical convention: ``eta = Gamma_HWHM / Omega``, underdamped threshold
``eta < 0.5`` (HWHM < Omega/2  <=>  FWHM < Omega).

**Per-source Gamma basis -- the FWHM->HWHM conversion happens HERE, exactly once,
and is never silent:**

- Yamada-Shirane 1969: the plotted/quoted Gamma is already HWHM (Eq. 7 Lorentzian
  half-width) -> ``gamma_convention="HWHM"`` (used directly).
- Vogt 1995: the plotted "damping constant gamma" is the DHO FWHM -> ``gamma_convention="FWHM"``
  so ``Gamma_HWHM = gamma / 2``.

The literature CSVs store gamma **as plotted** (FWHM for Vogt, HWHM for YS). This
module applies the ``/2`` for FWHM sources once, at consumption. Do NOT pre-halve the
CSV, and do NOT pass an already-HWHM value with ``convention="FWHM"`` (that would
double-apply the factor) -- a regression test guards against double-application.

API
---
``eta_fit``        primary: multi-point eta(T) from dense gamma(T) + Omega(T)
                   (STO and KTaO3 via Vogt 1995).
``eta_bracket``    fallback/cross-check: (eta_lo, eta_hi) from a Gamma range at one Omega.
``eta_underdamped`` bool, threshold 0.5.
``eta_sto_vogt`` / ``eta_kto_vogt`` convenience wrappers wiring ``literature_io``.
"""

from __future__ import annotations

import numpy as np

# Factor that converts a source's plotted gamma to Gamma_HWHM. The ONLY place the
# FWHM->HWHM /2 lives, so it is applied exactly once per consumption.
_HWHM_FACTOR: dict[str, float] = {"HWHM": 1.0, "FWHM": 0.5}


def gamma_to_hwhm(gamma, gamma_convention: str):
    """Convert a plotted ``gamma`` to ``Gamma_HWHM`` per the source's convention.

    ``gamma_convention`` is ``"HWHM"`` (use directly, e.g. Yamada-Shirane) or
    ``"FWHM"`` (halve, e.g. Vogt). Raises ``ValueError`` on any other value so a
    missing/unknown convention can never be silently mis-applied.
    """
    try:
        factor = _HWHM_FACTOR[gamma_convention]
    except KeyError:
        raise ValueError(
            f"gamma_convention must be one of {sorted(_HWHM_FACTOR)}, got {gamma_convention!r}"
        ) from None
    return factor * np.asarray(gamma, dtype=float)


def eta_fit(T_omega, omega_meV, T_gamma, gamma_meV, *, gamma_convention: str):
    """Multi-point ``eta(T) = Gamma_HWHM(T) / Omega(T)`` -- the primary path.

    ``gamma(T)`` is converted to HWHM once (per ``gamma_convention``) and linearly
    interpolated onto ``Omega``'s temperature grid over the T-overlap.

    Returns:
        ``(T, eta)`` numpy arrays over the overlapping temperature range.
    """
    To = np.asarray(T_omega, dtype=float)
    Om = np.asarray(omega_meV, dtype=float)
    Tg = np.asarray(T_gamma, dtype=float)
    Gh = gamma_to_hwhm(gamma_meV, gamma_convention)  # HWHM, conversion applied ONCE

    oi = np.argsort(To)
    To, Om = To[oi], Om[oi]
    gi = np.argsort(Tg)
    Tg, Gh = Tg[gi], Gh[gi]

    mask = (To >= Tg.min()) & (To <= Tg.max())
    To, Om = To[mask], Om[mask]
    if To.size == 0:
        raise ValueError("no temperature overlap between Omega(T) and gamma(T)")
    Gh_on_To = np.interp(To, Tg, Gh)
    return To, Gh_on_To / Om


def eta_bracket(omega_meV: float, gamma_lo: float, gamma_hi: float, *, gamma_convention: str):
    """Fallback ``(eta_lo, eta_hi)`` from a Gamma range at a single ``Omega``.

    Retained for sources without dense gamma(T); the same one-time convention
    conversion is applied to both endpoints.
    """
    glo = float(gamma_to_hwhm(gamma_lo, gamma_convention))
    ghi = float(gamma_to_hwhm(gamma_hi, gamma_convention))
    return glo / float(omega_meV), ghi / float(omega_meV)


def eta_underdamped(eta, threshold: float = 0.5) -> bool:
    """True iff every eta value is below the underdamped threshold (default 0.5)."""
    return bool(np.all(np.asarray(eta, dtype=float) < threshold))


# ---------------------------------------------------------------------------
# Convenience wrappers wiring literature_io (the committed dense Vogt sources)
# ---------------------------------------------------------------------------
def eta_sto_vogt():
    """eta_STO(T) from Vogt 1995 cubic_F1u (gamma = FWHM -> /2). Primary STO eta."""
    from fekappa import literature_io as lit

    o = lit.load_vogt1995_sto_omega_T()
    o = o[o["phase_mode"] == "cubic_F1u"]          # cubic phase only (>=110 K)
    g = lit.load_vogt1995_sto_gamma_T()
    g = g[g["phase_mode"] == "cubic_F1u"]
    return eta_fit(o["T_K"], o["Omega_meV"], g["T_K"], g["gamma_meV"],
                   gamma_convention="FWHM")


def eta_kto_vogt():
    """eta_KTaO3(T) from Vogt 1995 (gamma = FWHM -> /2). Primary KTaO3 eta."""
    from fekappa import literature_io as lit

    o = lit.load_vogt1995_kto_omega_T()
    o = o[o["method"] == "exp_HyperRaman_Vogt1995"]  # exclude the 1988 cross-ref pt
    g = lit.load_vogt1995_kto_gamma_T()
    return eta_fit(o["T_K"], o["Omega_meV"], g["T_K"], g["gamma_meV"],
                   gamma_convention="FWHM")
