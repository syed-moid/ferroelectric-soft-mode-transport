#!/usr/bin/env python3
"""Verify an ALM FC2 XML against QE DFPT frequencies at the commensurate q of a 2x2x2 supercell.

Reads Data.ForceConstants.HARMONIC.FC2 (pair1 = "prim_atom xyz", pair2 = "super_atom xyz cell"; values already
divided by the image multiplicity) and Data.Structure of the XML. At q commensurate with the supercell, all
periodic images of a pair carry the same phase, so summing over the image cells gives the periodic supercell
matrix exactly. D(q) is built without any non-analytic correction and diagonalised with the QE masses.

Usage: b1_h1_fc2_verify.py FC2.xml DYN_PREFIX NSTAR REF_ALM_IN
(primitive atoms of the XML are identified through the representative supercell atoms in its <Structure>)
"""
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from b1_h1_fc2_conversion import RY_TO_CM, alm_sites, read_dyn, signed_freqs  # noqa: E402


def main():
    xml_path, prefix, nstar, ref_in = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
    root = ET.parse(xml_path).getroot()
    # supercell positions from the XML (fractional)
    pos = [np.array([float(x) for x in p.text.split()]) for p in root.find("Structure").find("Position")]
    sites = alm_sites(ref_in)
    for k, (s, f) in enumerate(sites):
        if not np.allclose(pos[k], f, atol=1e-8):
            raise SystemExit(f"XML atom {k + 1} at {pos[k]} differs from the reference site {f}")
    # primitive atoms of the XML: first supercell atom of each primitive index
    trans = root.find("Symmetry").find("Translations")
    prim_rep = {}
    for t in trans:
        if t.attrib["tran"] == "1":
            prim_rep[int(t.attrib["atom"])] = int(t.text)
    alat, tau, ityp, masses, blocks, _ = read_dyn(f"{prefix}1")
    nat = len(tau)
    mvec = np.repeat([masses[t] for t in ityp], 3)
    sc = []
    for sk, f in sites:
        x = 2.0 * f
        for kap in range(nat):
            d = x - tau[kap]
            if ityp[kap] == sk and np.allclose(d, np.round(d), atol=1e-8):
                sc.append((kap, np.round(d).astype(int) % 2))
    N = len(sc)
    # periodic supercell matrix rows of the primitive atoms (sum over image cells)
    rows = {}
    for el in root.find("ForceConstants").find("HARMONIC"):
        p, a = map(int, el.attrib["pair1"].split())
        j, b, _cell = map(int, el.attrib["pair2"].split())
        rows.setdefault((p, a), np.zeros(3 * N))[3 * (j - 1) + (b - 1)] += float(el.text)
    asr = max(abs(v[c::3].sum()) for v in rows.values() for c in range(3))
    print(f"{xml_path}: {len(rows)} rows; translational sum rule: largest |sum_j Phi_ij| = {asr:.3e} Ry/bohr^2")
    worst = 0.0
    for n in range(1, nstar + 1):
        _, _, _, _, blks, freqs = read_dyn(f"{prefix}{n}")
        q = blks[0][0]
        Cq = np.zeros((3 * nat, 3 * nat), dtype=complex)
        for (p, a), v in rows.items():
            s = prim_rep[p] - 1
            ki, Li = sc[s]
            for j in range(N):
                kj, Lj = sc[j]
                ph = np.exp(-2j * np.pi * q.dot((Lj - Li).astype(float)))
                Cq[3 * ki + a - 1, 3 * kj:3 * kj + 3] += v[3 * j:3 * j + 3] * ph
        w = np.sort(signed_freqs(Cq, mvec))
        ref = np.sort(np.array(freqs))
        d = np.abs(w - ref)
        # at Gamma, exclude the three acoustic modes (smallest |freq| in the reference); their change is the
        # translational sum rule imposed by ALM
        keep = np.ones(len(d), bool)
        if n == 1:
            keep[np.argsort(np.abs(ref))[:3]] = False
        dopt = d[keep].max()
        worst = max(worst, dopt)
        print(f"q = {q}: XML  {np.array2string(w, precision=3, max_line_width=250)}")
        print(f"            ph.x {np.array2string(ref, precision=3, max_line_width=250)}")
        print(f"            max |diff| = {d.max():.4f} cm^-1; excluding the three Gamma acoustic modes: {dopt:.4f}")
    print(f"largest difference over Gamma optical, X, M, R: {worst:.4f} cm^-1 (criterion 0.5 cm^-1): {'PASS' if worst <= 0.5 else 'FAIL'}")


if __name__ == "__main__":
    main()
