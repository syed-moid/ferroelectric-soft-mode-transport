#!/usr/bin/env python3
"""Convert QE displacement-force pairs (DISP.in.NNN / DISP.out.NNN) to an ALM DFSET.

- Positions: the crystal coordinates of DISP.in.
- Ideal sites: the &position list of an ALM input (REFERENCE_ALM_IN, fractional coordinates of the 2x2x2
  supercell with species index 1=Sr, 2=Ti, 3=O). Each atom is mapped to the nearest ideal site of the
  same species (minimum image). The map must be a bijection and the same for every configuration.
- Displacements u = (r - r_ideal) with the minimum-image convention, in bohr (cubic cell: u = df * L).
- Forces: the last 'Forces acting on atoms' block of DISP.out, in Ry/bohr as printed by pw.x (no conversion).
- Output rows follow the ideal-site order of REFERENCE_ALM_IN: "ux uy uz fx fy fz" per atom, NAT rows per
  configuration, configurations in file-number order. A header comment line precedes each configuration.

Usage: b1_convert_dfset.py SET_DIR REFERENCE_ALM_IN OUT_DFSET [--check-atom N]
"""
import re
import sys
from pathlib import Path

BOHR_A = 0.529177210903
SPECIES = {"Sr": 1, "Ti": 2, "O": 3}


def alm_sites(path):
    t = Path(path).read_text()
    blk = t.split("&position", 1)[1].split("/", 1)[0]
    return [(int(l.split()[0]), tuple(float(x) for x in l.split()[1:4])) for l in blk.strip().splitlines()]


def read_in(path):
    t = Path(path).read_text()
    a = re.search(r"CELL_PARAMETERS.*?\n((?:\s*[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n){3})", t).group(1)
    cell = [list(map(float, l.split())) for l in a.strip().splitlines()]
    for i in range(3):
        for j in range(3):
            if i != j and abs(cell[i][j]) > 1e-12:
                raise SystemExit(f"{path}: cell is not diagonal")
    if not (abs(cell[0][0] - cell[1][1]) < 1e-12 and abs(cell[0][0] - cell[2][2]) < 1e-12):
        raise SystemExit(f"{path}: cell is not cubic")
    if "ATOMIC_POSITIONS (crystal)" not in t and "ATOMIC_POSITIONS crystal" not in t:
        raise SystemExit(f"{path}: positions are not in crystal units")
    pos = []
    for l in t.split("ATOMIC_POSITIONS", 1)[1].splitlines()[1:]:
        s = l.split()
        if len(s) >= 4 and s[0] in SPECIES:
            pos.append((SPECIES[s[0]], tuple(float(x) for x in s[1:4])))
    return cell[0][0], pos


def read_out(path):
    t = Path(path).read_text()
    if "JOB DONE" not in t or "convergence has been achieved" not in t:
        raise SystemExit(f"{path}: incomplete or unconverged")
    blk = t.rsplit("Forces acting on atoms", 1)[1].split("Total force")[0]
    if "Ry/au" not in t.rsplit("Forces acting on atoms", 1)[1][:60]:
        raise SystemExit(f"{path}: force units are not Ry/au")
    return [tuple(float(x) for x in m.groups()) for m in
            re.finditer(r"atom\s+\d+\s+type\s+\d+\s+force =\s+(\S+)\s+(\S+)\s+(\S+)", blk)]


def mic(d):
    return d - round(d)


def main():
    set_dir, ref_in, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
    check_atom = int(sys.argv[sys.argv.index("--check-atom") + 1]) if "--check-atom" in sys.argv else None
    sites = alm_sites(ref_in)
    nat = len(sites)
    files = sorted(set_dir.glob("DISP.in.*"), key=lambda p: int(p.suffix[1:]))
    perm0, lines, umax, netmax = None, [], 0.0, 0.0
    for pin in files:
        L_A, pos = read_in(pin)
        L = L_A / BOHR_A
        forces = read_out(set_dir / pin.name.replace("DISP.in", "DISP.out"))
        if len(pos) != nat or len(forces) != nat:
            raise SystemExit(f"{pin}: {len(pos)} positions, {len(forces)} forces, expected {nat}")
        perm = []  # perm[k] = dataset atom index for ALM site k
        for k, (sk, fk) in enumerate(sites):
            best = min(((sum((mic(f[c] - fk[c]) * L) ** 2 for c in range(3)), i)
                        for i, (s, f) in enumerate(pos) if s == sk))
            perm.append(best[1])
        if len(set(perm)) != nat:
            raise SystemExit(f"{pin}: atom-to-site map is not a bijection")
        if perm0 is None:
            perm0 = perm
        elif perm != perm0:
            raise SystemExit(f"{pin}: atom order differs from the first configuration")
        lines.append(f"# {pin.name}  L = {L:.10f} bohr")
        for k, (sk, fk) in enumerate(sites):
            i = perm[k]
            u = [mic(pos[i][1][c] - fk[c]) * L for c in range(3)]
            umax = max(umax, max(abs(x) for x in u))
            f = forces[i]
            lines.append("%20.12f %20.12f %20.12f %20.12e %20.12e %20.12e" % (*u, *f))
            if check_atom is not None and pin == files[0] and k == check_atom - 1:
                print(f"check: {pin.name} ALM site {k + 1} <- dataset atom {i + 1}; frac {pos[i][1]} ideal {fk}; "
                      f"L = {L_A} A = {L:.10f} bohr; u = {u}; F = {f}")
        for c in range(3):
            netmax = max(netmax, abs(sum(fr[c] for fr in forces)))
    out.write_text("\n".join(lines) + "\n")
    print(f"{set_dir}: {len(files)} configurations -> {out}")
    print(f"atom map (ALM site -> dataset atom, 1-based): {[p + 1 for p in perm0]}")
    print(f"largest |u| component: {umax:.6f} bohr ({umax * BOHR_A:.4f} A); half cell edge {L / 2:.4f} bohr")
    print(f"largest |net force| component: {netmax:.3e} Ry/bohr")


if __name__ == "__main__":
    main()
