#!/usr/bin/env python3
"""Descriptive inventory of a QE displacement-force set (DISP.in.NNN / DISP.out.NNN pairs).

For every configuration: parses the supercell and the crystal positions of DISP.in, and the total forces
of the last 'Forces acting on atoms' block of DISP.out. Displacements are taken relative to the ideal
perovskite sites of a 2x2x2 supercell of the primitive cell in CELL_IN, with the minimum-image convention.

Usage: b0_inventory.py SET_DIR CELL_IN LABEL
Prints a summary to stdout (Angstrom, eV/Angstrom).
"""
import math
import re
import sys
from pathlib import Path

RY_BOHR_TO_EV_A = 13.605693122994 / 0.529177210903


def read_cell_in(path):
    t = Path(path).read_text()
    cell = [list(map(float, l.split())) for l in
            re.search(r"CELL_PARAMETERS.*?\n((?:\s*[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n){3})", t).group(1).strip().splitlines()]
    block = t.split("ATOMIC_POSITIONS", 1)[1].splitlines()[1:]
    prim = [(l.split()[0], tuple(map(float, l.split()[1:4]))) for l in block if len(l.split()) >= 4]
    return cell, prim


def ideal_sites(prim, n=2):
    sites = []
    for sp, f in prim:
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    sites.append((sp, ((f[0] + i) / n, (f[1] + j) / n, (f[2] + k) / n)))
    return sites


def read_in(path):
    t = Path(path).read_text()
    a = re.search(r"CELL_PARAMETERS.*?\n((?:\s*[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n){3})", t).group(1)
    cell = [list(map(float, l.split())) for l in a.strip().splitlines()]
    pos = []
    for l in t.split("ATOMIC_POSITIONS", 1)[1].splitlines()[1:]:
        s = l.split()
        if len(s) >= 4 and s[0].isalpha():
            pos.append((s[0], tuple(map(float, s[1:4]))))
    return cell, pos


def read_out(path):
    t = Path(path).read_text()
    if "JOB DONE" not in t:
        raise SystemExit(f"{path}: no JOB DONE")
    if "convergence has been achieved" not in t:
        raise SystemExit(f"{path}: SCF not converged")
    blk = t.rsplit("Forces acting on atoms", 1)[1]
    f = [tuple(float(x) * RY_BOHR_TO_EV_A for x in m.groups())
         for m in re.finditer(r"atom\s+\d+\s+type\s+\d+\s+force =\s+(\S+)\s+(\S+)\s+(\S+)", blk.split("Total force")[0])]
    return f


def main():
    d, cell_in, label = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
    pcell, prim = read_cell_in(cell_in)
    sites = ideal_sites(prim)
    ins = sorted(d.glob("DISP.in.*"))
    stats = {}
    maxnet, nconf, natoms, Fsq, order_ok = 0.0, 0, None, 0.0, True
    for pin in ins:
        pout = d / pin.name.replace("DISP.in", "DISP.out")
        cell, pos = read_in(pin)
        L = cell[0][0]
        forces = read_out(pout)
        natoms = len(pos)
        if len(forces) != natoms:
            raise SystemExit(f"{pout}: {len(forces)} forces for {natoms} atoms")
        for i, (sp, f) in enumerate(pos):
            # nearest ideal site of the same species (minimum image)
            best = None
            for sp2, g in sites:
                if sp2 != sp:
                    continue
                dv = [(f[c] - g[c]) - round(f[c] - g[c]) for c in range(3)]
                r = math.sqrt(sum((x * L) ** 2 for x in dv))
                if best is None or r < best[0]:
                    best = (r, dv)
            if sites[i][0] != sp:
                order_ok = False
            stats.setdefault(sp, []).append(best[0])
        for c in range(3):
            maxnet = max(maxnet, abs(sum(fr[c] for fr in forces)))
        Fsq += sum(x * x for fr in forces for x in fr)
        nconf += 1
    print(f"[{label}] configurations: {nconf}; atoms per configuration: {natoms}; supercell edge: {L:.9f} A "
          f"(a = {L / 2:.9f} A)")
    print(f"[{label}] species order of DISP.in matches Sr8 Ti8 O24 blocks: {order_ok}")
    for sp, v in stats.items():
        rms = math.sqrt(sum(x * x for x in v) / len(v))
        print(f"[{label}] {sp:2s}: n = {len(v)}, rms |u| = {rms:.4f} A, max |u| = {max(v):.4f} A")
    print(f"[{label}] largest |net force| component over all configurations: {maxnet:.3e} eV/A")
    print(f"[{label}] rms force component: {math.sqrt(Fsq / (nconf * natoms * 3)):.4f} eV/A")


if __name__ == "__main__":
    main()
