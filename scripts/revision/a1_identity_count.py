"""Four-phonon audit support: independent count of the quartets the 4ph audit build excludes from diagram (c), and a
classification of the 'near-miss' quartets it keeps (TOL_DEG <= min |resonant denominator| < 1e-6 Ry).
Reproduces selfenergy_c's loop (ik1, ik2 over the full mesh; k3 = q - k1 - k2; all branch triples) from the
MODE = phonons eigenvalues of the same FC2 and mesh (evec_pm_q12_T300_m<n>.mesh.evec; eigenvalues are omega^2 in Ry^2).
Usage: python a1_identity_count.py <mesh.evec> <qx qy qz> <branch 1-based> [...more targets]
Note: the .mesh.evec eigenvalues carry 7 significant digits, i.e. omega to ~5e-8 relative (~1e-10 Ry at
1e-3 Ry): adequate for the 1e-9 Ry tolerance only marginally; near-miss classification is resolved above ~1e-10 Ry.
"""
import sys
from pathlib import Path

import numpy as np

from fekappa.evec_io import read_evec

TOL, NEAR = 1e-9, 1e-6
xk, ev2, _ = read_evec(sys.argv[1])
om = np.sign(ev2) * np.sqrt(np.abs(ev2))  # the .mesh.evec eigenvalues are omega^2 (Ry^2)
nk, ns = om.shape
n = round(nk ** (1 / 3))
idx = {tuple(np.round((x * n) % n).astype(int)): i for i, x in enumerate(xk)}


def kindex(frac):
    return idx[tuple(np.round((np.asarray(frac) * n) % n).astype(int) % n)]


args = sys.argv[2:]
for t in range(0, len(args), 4):
    q = [float(v) for v in args[t:t + 3]]; s = int(args[t + 3]) - 1
    kq = kindex(q); w = om[kq, s]
    ex = chk = near = 0
    cls = {}
    for ik1 in range(nk):
        k3 = [kindex(np.array(q) - xk[ik1] - xk[ik2]) for ik2 in range(nk)]
        w1 = om[ik1][:, None, None, None]            # [s1, 1, 1, 1]
        w2 = om[:, None, :, None].transpose(1, 0, 2, 3)  # [1, ik2, s2, 1]
        w3 = om[k3][None, :, None, :]                # [1, ik2, 1, s3]
        d = np.minimum(np.minimum(np.abs(w - w1 - w2 + w3), np.abs(w + w1 - w2 - w3)), np.abs(w - w1 + w2 - w3))
        ik2a = np.arange(nk)[None, :, None, None]
        k3a = np.array(k3)[None, :, None, None]
        leg1 = (ik1 == kq) & (np.abs(w1 - w) < TOL) & (np.abs(w2 - w3) < TOL)
        leg2 = (ik2a == kq) & (np.abs(w2 - w) < TOL) & (np.abs(w1 - w3) < TOL)
        leg3 = (k3a == kq) & (np.abs(w3 - w) < TOL) & (np.abs(w1 - w2) < TOL)
        e = leg1 | leg2 | leg3
        ex += int(e.sum()); chk += int((d < TOL).sum())
        nm = (~e) & (d >= TOL) & (d < NEAR)
        near += int(nm.sum())
        if nm.any():
            # structure of near misses: does a leg sit at k = q, and is its frequency within NEAR of w?
            on_q = ((ik1 == kq) & (np.abs(w1 - w) < NEAR)) | ((ik2a == kq) & (np.abs(w2 - w) < NEAR)) | \
                   ((k3a == kq) & (np.abs(w3 - w) < NEAR))
            on_q = np.broadcast_to(on_q, d.shape)
            for lab, m in (("leg at q with |w_i - w| < 1e-6 (identity-like, split degeneracy)", nm & on_q),
                           ("no leg at q (other exact-resonance family)", nm & ~on_q)):
                cls[lab] = cls.get(lab, 0) + int(m.sum())
    print(f"target q = {q}, branch {s + 1}, w = {w:.6e} Ry ({w * 109737.316:.2f} cm-1), mesh {n}^3: "
          f"excluded {ex}, exact-zero denominators {chk}, near-misses {near}: {cls}")
