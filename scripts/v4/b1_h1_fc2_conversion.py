#!/usr/bin/env python3
"""H1 conversion: periodic 2x2x2 supercell FC2 from QE DFPT dynamical matrices at the 8 commensurate q.

Steps
1. Read C(q) (Ry/bohr^2, not mass-weighted) for every q of every star in DYN_PREFIX1..N (ph.x output).
   At Gamma the file holds the analytic part only, which is what a periodic supercell calculation gives.
2. Phi(kappa L; kappa' L') = (1/8) sum_q C_{kappa kappa'}(q) exp(i 2pi q.(R_L' - R_L)), q in 2pi/alat,
   R in alat. For the commensurate q of a 2x2x2 cell the phases are +-1. The result is the full harmonic
   force constant matrix of the periodic supercell (all periodic images summed), as finite displacements
   would give. The q2r.x real-space file is NOT used (it has the dipole part subtracted).
3. Checks written to stdout: |Im| of the inverse transform, symmetry Phi_ij = Phi_ji, translational sum
   rule sum_j Phi_ij per row before any correction; frequencies at the 8 q from the supercell matrix
   (back-transform) against the frequencies printed by ph.x.
4. Writes a synthetic ALM DFSET: NCONF configurations of random displacements (all atoms, components
   uniform in [-AMP, AMP] bohr, fixed seed) with exact harmonic forces F = -Phi u, in the atom order of
   the ALM &position list REF_ALM_IN. ALM then projects Phi onto its symmetry/ASR basis (NORDER = 1,
   least squares), which yields an FC2XML in ALM's own format. This is a format conversion with
   synthetic forces, not a fit to the dataset.

Usage: b1_h1_fc2_conversion.py DYN_PREFIX NSTAR REF_ALM_IN OUT_DFSET [NCONF AMP SEED]
"""
import re
import sys
from pathlib import Path

import numpy as np

RY_TO_CM = 109737.31568160  # 1 Ry in cm^-1
AMU_RY = 911.44424310865645  # 1 amu in Ry atomic mass units (m_e/2)


def read_dyn(path):
    t = Path(path).read_text()
    head = t.split("Dynamical  Matrix in cartesian axes")[0].splitlines()
    ntyp, nat = int(head[2].split()[0]), int(head[2].split()[1])
    alat = float(head[2].split()[3])
    i0 = [k for k, l in enumerate(head) if l.strip().startswith("Basis vectors")][0]
    masses = {}
    for l in head[i0 + 4:i0 + 4 + ntyp]:
        s = l.split()
        masses[int(s[0])] = float(s[-1])
    tau, ityp = [], []
    for l in head[i0 + 4 + ntyp:i0 + 4 + ntyp + nat]:
        s = l.split()
        ityp.append(int(s[1]))
        tau.append([float(x) for x in s[2:5]])
    blocks = []
    for blk in t.split("Dynamical  Matrix in cartesian axes")[1:]:
        blk = blk.split("Diagonalizing")[0].split("Dielectric Tensor")[0]
        q = [float(x) for x in re.search(r"q = \(\s*(\S+)\s+(\S+)\s+(\S+)\s*\)", blk).groups()]
        C = np.zeros((3 * nat, 3 * nat), dtype=complex)
        lines = blk.splitlines()
        k = 0
        while k < len(lines):
            s = lines[k].split()
            if len(s) == 2 and all(x.isdigit() for x in s):
                a, b = int(s[0]) - 1, int(s[1]) - 1
                for r in range(3):
                    v = [float(x) for x in lines[k + 1 + r].split()]
                    for c in range(3):
                        C[3 * a + r, 3 * b + c] = v[2 * c] + 1j * v[2 * c + 1]
                k += 4
            else:
                k += 1
        blocks.append((np.array(q), C))
    freqs = [float(m.group(1)) for m in re.finditer(r"freq \(\s*\d+\) =\s*\S+ \[THz\] =\s*(\S+) \[cm-1\]", t)]
    return alat, np.array(tau), ityp, masses, blocks, freqs


def signed_freqs(C, mvec):
    D = C / np.sqrt(np.outer(mvec, mvec))
    w2 = np.linalg.eigvalsh((D + D.conj().T) / 2)
    return np.sign(w2) * np.sqrt(np.abs(w2)) * RY_TO_CM


def alm_sites(path):
    t = Path(path).read_text()
    blk = t.split("&position", 1)[1].split("/", 1)[0]
    return [(int(l.split()[0]), np.array([float(x) for x in l.split()[1:4]])) for l in blk.strip().splitlines()]


def main():
    prefix, nstar, ref_in, out = sys.argv[1], int(sys.argv[2]), sys.argv[3], Path(sys.argv[4])
    nconf = int(sys.argv[5]) if len(sys.argv) > 5 else 12
    amp = float(sys.argv[6]) if len(sys.argv) > 6 else 0.02
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 20261119
    allq, ref_freqs = [], []
    for n in range(1, nstar + 1):
        alat, tau, ityp, masses, blocks, freqs = read_dyn(f"{prefix}{n}")
        allq += blocks
        ref_freqs.append((blocks[0][0], freqs))
    nat = len(tau)
    mvec = np.repeat([masses[t] for t in ityp], 3)
    print(f"alat = {alat} bohr; {nat} atoms; {len(allq)} q points: {[list(q) for q, _ in allq]}")
    maxim = max(np.abs(C.imag).max() for _, C in allq)
    print(f"largest |Im C(q)| over all q: {maxim:.3e} Ry/bohr^2")

    # supercell sites (ALM order): species index, fractional coordinates of the 2x2x2 cell
    sites = alm_sites(ref_in)
    spec_of_type = {1: 1, 2: 2, 3: 3}  # QE type index -> ALM species index (Sr, Ti, O in both)
    sc = []  # (kappa, L) per ALM site
    for sk, f in sites:
        x = 2.0 * f  # alat units
        hit = None
        for kap in range(nat):
            if spec_of_type[ityp[kap]] != sk:
                continue
            d = x - tau[kap]
            if np.allclose(d, np.round(d), atol=1e-8):
                hit = (kap, np.round(d).astype(int) % 2)
        if hit is None:
            raise SystemExit(f"site {f} not matched")
        sc.append(hit)
    N = len(sc)
    Phi = np.zeros((3 * N, 3 * N), dtype=complex)
    for i, (ki, Li) in enumerate(sc):
        for j, (kj, Lj) in enumerate(sc):
            dR = (Lj - Li).astype(float)
            acc = np.zeros((3, 3), dtype=complex)
            for q, C in allq:
                acc += C[3 * ki:3 * ki + 3, 3 * kj:3 * kj + 3] * np.exp(2j * np.pi * q.dot(dR))
            Phi[3 * i:3 * i + 3, 3 * j:3 * j + 3] = acc / len(allq)
    print(f"largest |Im Phi|: {np.abs(Phi.imag).max():.3e}")
    Phi = Phi.real
    print(f"largest |Phi - Phi^T|: {np.abs(Phi - Phi.T).max():.3e} Ry/bohr^2; largest |Phi|: {np.abs(Phi).max():.6f}")
    asr = np.array([[Phi[3 * i + a, b::3].sum() for b in range(3)] for i in range(N) for a in range(3)])
    print(f"translational sum rule before correction: largest |sum_j Phi_ij| = {np.abs(asr).max():.3e} Ry/bohr^2")

    # back-transform check: frequencies at the irreducible q from the supercell matrix
    msc = np.repeat([masses[ityp[k]] for k, _ in sc], 3)
    for q, freqs in ref_freqs:
        Cq = np.zeros((3 * nat, 3 * nat), dtype=complex)
        for i, (ki, Li) in enumerate(sc):
            if Li.any():
                continue  # rows from cell 0
            for j, (kj, Lj) in enumerate(sc):
                Cq[3 * ki:3 * ki + 3, 3 * kj:3 * kj + 3] += Phi[3 * i:3 * i + 3, 3 * j:3 * j + 3] * np.exp(-2j * np.pi * q.dot(Lj.astype(float)))
        w = signed_freqs(Cq, mvec)
        print(f"q = {q}: max |freq(back-transform) - freq(ph.x)| = {np.abs(np.sort(w) - np.sort(freqs)).max():.4f} cm^-1")

    rng = np.random.default_rng(seed)
    lines = []
    for c in range(nconf):
        u = rng.uniform(-amp, amp, size=3 * N)
        F = -Phi @ u
        lines.append(f"# synthetic configuration {c + 1}: harmonic forces from DFPT q2k6 supercell FC2")
        for i in range(N):
            lines.append("%20.12f %20.12f %20.12f %22.14e %22.14e %22.14e" % (*u[3 * i:3 * i + 3], *F[3 * i:3 * i + 3]))
    out.write_text("\n".join(lines) + "\n")
    np.save(out.with_suffix(".phi.npy"), Phi)
    print(f"wrote {out} ({nconf} configurations, |u| <= {amp} bohr, seed {seed}) and {out.with_suffix('.phi.npy')}")


if __name__ == "__main__":
    main()
