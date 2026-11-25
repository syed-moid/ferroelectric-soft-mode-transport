#!/usr/bin/env python3
"""Forces predicted by an ALM force-constant set (.fcs, all orders) on a DFSET; E_F and RMSE (G-B1).

The '**FCn' tables of the .fcs list, for each irreducible parameter (value in the '# FCn_k' header), the
Cartesian elements Phi(i; j2 <= ... <= jn) with coefficients; i is a representative atom, j2..jn are sorted
supercell indices. Elements for the other atoms follow from the 8 lattice translations of the 2x2x2 supercell.
Force: F_i = -sum_n sum_{j2..jn ordered} Phi(i,j2..jn) u_j2..u_jn / (n-1)!; a sorted tuple stands for
(n-1)!/prod(m_k!) ordered tuples, so its weight is 1/prod(m_k!) (m_k = multiplicities of repeated indices).
This is the row of ALM's sensing matrix A times x, evaluated on the supercell indices as ALM does.

E_F = 100 sqrt(sum (F_fit - F_DFT)^2 / sum F_DFT^2) over complete structures, atoms, components (B1 criteria).

Usage: b1_force_eval.py FCS ALM_INPUT DFSET NDATA [--save PRED.npy] [--xml] [--orders 2,3,4]
  ALM_INPUT: the alm input of the fit (cell and &position are read from it).
  --xml: FCS is an ALM/anphon XML instead of a .fcs file. Its HARMONIC and ANHARMn blocks list, for
        each primitive atom (pair1, mapped to the supercell by translation 1 of the XML), the elements with the other
        indices once per sorted tuple, split over periodic images; the image index is dropped and the entries are
        used exactly as the .fcs elements (same weights, same 8 translations). Validated against the .fcs route on
        the distributed fit (docs/v4/ANALYSIS_STEPS.md, step P4).
  --orders: restrict the evaluation to these orders (validation only; default all orders in the file).
"""
import argparse
import re
from collections import Counter
from math import factorial
from pathlib import Path

import numpy as np

RY_BOHR_TO_MEV_A = 13.605693122994e3 / 0.529177210903  # Ry/bohr -> meV/Angstrom


def read_alm_input(p):
    txt = Path(p).read_text()
    cell = re.search(r"&cell\s*\n(.*?)/", txt, re.S).group(1).split("\n")
    cell = [l.split() for l in cell if l.split() and not l.strip().startswith("#")]
    a = float(cell[0][0])
    lat = a * np.array([[float(x) for x in r] for r in cell[1:4]])
    pos = re.search(r"&position\s*\n(.*?)/", txt, re.S).group(1).split("\n")
    pos = [l.split() for l in pos if l.split()]
    kd = np.array([int(r[0]) for r in pos])
    frac = np.array([[float(x) for x in r[1:4]] for r in pos])
    return lat, kd, frac


def translation_maps(frac):
    shifts = [np.array(s) * 0.5 for s in np.ndindex(2, 2, 2)]
    maps = []
    for s in shifts:
        f = (frac + s) % 1.0
        d = f[:, None, :] - frac[None, :, :]
        d -= np.round(d)
        m = np.argmin(np.abs(d).sum(-1), axis=1)
        assert np.abs(d[np.arange(len(frac)), m]).max() < 1e-6
        assert len(set(m)) == len(frac)
        maps.append(m)
    return np.array(maps)  # (8, nat)


def parse_fcs(p):
    """Return {order: (idx (nterm, n) of 3*atom+xyz, 0-based; value (nterm,))} for the representative elements."""
    out = {}
    lines = Path(p).read_text().splitlines()
    i = next(k for k, l in enumerate(lines) if l.strip().startswith("**FC"))
    order = None
    idx, val = [], []
    cur = None
    xyz = {"x": 0, "y": 1, "z": 2}
    while i < len(lines):
        l = lines[i]
        if l.strip().startswith("**FC"):
            if order is not None:
                out[order] = (np.array(idx, dtype=np.int64), np.array(val))
            order = int(l.strip()[4:].split()[0])
            idx, val = [], []
        elif l.strip().startswith("#"):
            s = l.split()
            cur = float(s[3])
        elif l.strip():
            s = l.split()
            c = float(s[1])
            idx.append([3 * (int(t[:-1]) - 1) + xyz[t[-1]] for t in s[2:]])
            val.append(c * cur)
        i += 1
    if order is not None:
        out[order] = (np.array(idx, dtype=np.int64), np.array(val))
    return out


def parse_xml(p):
    """Return {order: (idx, value)} from the HARMONIC / ANHARMn blocks of an ALM XML (same layout as parse_fcs)."""
    p2s, out, order = {}, {}, None
    pat = re.compile(r'pair\d="([^"]+)"')
    with open(p) as f:
        for l in f:
            if "<map " in l and 'tran="1"' in l:
                m = re.search(r'atom="(\d+)">(\d+)<', l)
                p2s[int(m.group(1))] = int(m.group(2))
            elif "<HARMONIC>" in l:
                order = 2; out[order] = ([], [])
            elif re.search(r"<ANHARM(\d)>", l):
                order = int(re.search(r"<ANHARM(\d)>", l).group(1)); out[order] = ([], [])
            elif re.search(r"</(HARMONIC|ANHARM\d)>", l):
                order = None
            elif order is not None and "<FC" in l:
                ps = [tuple(int(x) for x in t.split()) for t in pat.findall(l)]
                first = 3 * (p2s[ps[0][0]] - 1) + ps[0][1] - 1
                rest = sorted(3 * (a[0] - 1) + a[1] - 1 for a in ps[1:])
                out[order][0].append([first] + rest)
                out[order][1].append(float(re.search(r">([^<]+)<", l).group(1)))
    return {n: (np.array(i, dtype=np.int64), np.array(v)) for n, (i, v) in out.items()}


def expand(fcs, tmaps):
    """Apply the 8 translations; return {order: (first (N,), rest (N, n-1), weight*value (N,))}."""
    nat = tmaps.shape[1]
    full = {}
    for n, (idx, val) in fcs.items():
        firsts = set(idx[:, 0] // 3)
        orbits = {tuple(sorted(tmaps[:, a])) for a in firsts}
        assert len(orbits) == len(firsts), "representative atoms not in distinct translation orbits"
        w = np.array([1.0 / np.prod([factorial(m) for m in Counter(r[1:]).values()]) for r in map(tuple, idx)])
        F, R, V = [], [], []
        for t in tmaps:
            cmap = (3 * t[:, None] + np.arange(3)[None, :]).ravel()  # component map
            tidx = cmap[idx]
            F.append(tidx[:, 0])
            R.append(tidx[:, 1:])
            V.append(val * w)  # [(n-1)!/prod m_k!] ordered tuples / (n-1)! = 1/prod m_k!
        full[n] = (np.concatenate(F), np.concatenate(R), np.concatenate(V))
        assert nat * 3 > full[n][0].max()
    return full


def predict(full, u):
    """u: (3 nat,) displacement vector -> force (3 nat,)."""
    f = np.zeros_like(u)
    for n, (first, rest, v) in full.items():
        prod = v * np.prod(u[rest], axis=1) if rest.shape[1] else v
        np.add.at(f, first, -prod)
    return f


def read_dfset(p, nat, ndata):
    rows = np.array([[float(x) for x in l.split()] for l in Path(p).read_text().splitlines()
                     if l.split() and not l.lstrip().startswith("#")])
    assert rows.shape == (nat * ndata, 6), rows.shape
    return rows[:, :3].reshape(ndata, nat * 3), rows[:, 3:].reshape(ndata, nat * 3)


def metrics(pred, ref, kd):
    d = pred - ref
    out = {"E_F_pct": 100 * np.sqrt((d ** 2).sum() / (ref ** 2).sum()),
           "RMSE_meV_A": np.sqrt((d ** 2).mean()) * RY_BOHR_TO_MEV_A,
           "rms_F_meV_A": np.sqrt((ref ** 2).mean()) * RY_BOHR_TO_MEV_A}
    comp_kd = np.repeat(kd, 3)
    for k in sorted(set(kd)):
        m = comp_kd == k
        out[f"E_F_kd{k}_pct"] = 100 * np.sqrt((d[:, m] ** 2).sum() / (ref[:, m] ** 2).sum())
        out[f"RMSE_kd{k}_meV_A"] = np.sqrt((d[:, m] ** 2).mean()) * RY_BOHR_TO_MEV_A
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fcs")
    ap.add_argument("alm_input")
    ap.add_argument("dfset")
    ap.add_argument("ndata", type=int)
    ap.add_argument("--save")
    ap.add_argument("--xml", action="store_true")
    ap.add_argument("--orders")
    a = ap.parse_args()
    lat, kd, frac = read_alm_input(a.alm_input)
    tmaps = translation_maps(frac)
    fcs = parse_xml(a.fcs) if a.xml else parse_fcs(a.fcs)
    if a.orders:
        keep = {int(x) for x in a.orders.split(",")}
        fcs = {n: v for n, v in fcs.items() if n in keep}
    full = expand(fcs, tmaps)
    for n, (f, r, v) in sorted(full.items()):
        print(f"order {n}: {len(v)} Cartesian elements (8 translations)")
    U, F = read_dfset(a.dfset, len(kd), a.ndata)
    P = np.array([predict(full, u) for u in U])
    if a.save:
        np.save(a.save, P)
    for k, v in metrics(P, F, kd).items():
        print(f"{k}: {v:.5f}")


if __name__ == "__main__":
    main()
