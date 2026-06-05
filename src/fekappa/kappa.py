"""Lattice thermal-conductivity integrators.

Implements the three contributions to the lattice thermal conductivity:

* ``kappa_bte_rta`` -- the particle-like Boltzmann transport contribution in
  the relaxation-time approximation (Eq. (32) of the manuscript), split into
  acoustic (kA) and optical (kO) channels.
* ``kappa_wigner`` -- the wave-like coherent contribution of
  Simoncelli, Marzari & Mauri.
* ``kappa_total`` -- their sum.
"""

from __future__ import annotations


def kappa_bte_rta(*args, **kwargs):
    """Boltzmann transport thermal conductivity in the RTA.

    Evaluates the particle-like (populations) contribution to the lattice
    thermal conductivity following Eq. (32) of the manuscript, accumulating the
    acoustic (kA) and optical (kO) channels.

    Raises:
        NotImplementedError: Always, for now.
    """
    raise NotImplementedError("Turn N: Boltzmann RTA kappa, Eq. (32)")


def kappa_wigner(*args, **kwargs):
    """Wigner (coherent) thermal-conductivity contribution.

    Evaluates the wave-like coherent term of Simoncelli, Marzari & Mauri,
    important when interband spacings are comparable to phonon linewidths --
    as for the overdamped soft mode near the transition.

    Raises:
        NotImplementedError: Always, for now.
    """
    raise NotImplementedError("Turn N: Wigner coherent kappa term")


def kappa_total(*args, **kwargs):
    """Total lattice thermal conductivity (RTA + Wigner).

    Raises:
        NotImplementedError: Always, for now.
    """
    raise NotImplementedError("Turn N: total kappa = kappa_bte_rta + kappa_wigner")
