"""fekappa: lattice thermal conductivity of anharmonic ferroelectric crystals.

Support code for the paper:

    "Thermal conductivity of anharmonic ferroelectric crystals across the
     soft-mode transition: SCPH validation and Wigner-corrected predictions
     for BaTiO3 and SrTiO3"

The package modernizes a 1987 Green's-function derivation of lattice thermal
conductivity in displacive ferroelectrics, validates it against DFPT /
Quantum-ESPRESSO data and Shirane neutron data, and develops a low-parameter
predictive fit. The primary theoretical framework is self-consistent phonon
theory (SCPH; Tadano & Tsuneyuki, Phys. Rev. B 92, 054301 (2015)), with the
stochastic self-consistent harmonic approximation (SSCHA) used as a
variational benchmark.

Subpackages / modules
---------------------
io
    Parsers for the raw BTO/STO data sets.
scph
    SCPH fits and soft-mode renormalization.
selfenergy
    Three- and four-phonon damping (self-energy) decomposition.
kappa
    Boltzmann RTA, Wigner coherent, and total thermal-conductivity integrators.
plots
    Publication matplotlib style.
"""

__version__ = "3.0.0"
