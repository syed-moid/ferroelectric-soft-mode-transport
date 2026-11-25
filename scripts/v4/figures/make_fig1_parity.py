#!/usr/bin/env python3
"""Fig. 1 (also Fig. S1): predicted against first-principles force components on the 128 held-out configurations of
model B (40 atoms, 3 components each). Three panels: H1 (harmonic terms fixed), H2 (harmonic terms fitted) and model A
applied to the same configurations. Each panel states E_F and the RMSE computed here from the stored predictions
(same definitions as scripts/v4/b1_force_eval.py).

Inputs: the test DFSET of model B (displacements in bohr, forces in Ry/bohr; ALM format) and the stored predictions
pred_h1_test.npy, pred_h2_test.npy, pred_modelA_a345_test.npy (written by b1_force_eval.py --save; shape 128 x 120,
Ry/bohr). The predictions are derived from the model B dataset (CC BY-NC 4.0).

Usage: make_fig1_parity.py DFSET_TEST PRED_H1 PRED_H2 PRED_MODEL_A OUT_PNG
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from b1_force_eval import RY_BOHR_TO_MEV_A, read_dfset  # noqa: E402

NAT, NDATA = 40, 128
EV_A = RY_BOHR_TO_MEV_A / 1000.0
INK, INK2, GRID, BLUE = "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6"
PANELS = (("H1: FC2 fixed (DFPT k6)", 1), ("H2: orders 2–6 fitted together", 2),
          ("Model A (a345) on model B test set\n(descriptive, different DFT setup)", 3))


def main():
    dfset, p1, p2, pa, out = sys.argv[1:6]
    _, ref = read_dfset(dfset, NAT, NDATA)
    preds = [np.load(p) for p in (p1, p2, pa)]
    for p in preds:
        assert p.shape == ref.shape, (p.shape, ref.shape)
    x = ref.ravel() * EV_A
    lim = 1.05 * max(np.abs(x).max(), max(np.abs(p).max() for p in preds) * EV_A)
    fig, axs = plt.subplots(1, 3, figsize=(10.2, 4.2), dpi=150, sharey=True)
    for ax, (title, _), p in zip(axs, PANELS, preds):
        d = p - ref
        ef = 100 * np.sqrt((d ** 2).sum() / (ref ** 2).sum())
        rmse = np.sqrt((d ** 2).mean()) * RY_BOHR_TO_MEV_A
        ax.plot([-lim, lim], [-lim, lim], color="#9a9996", lw=1, zorder=1)
        ax.scatter(x, p.ravel() * EV_A, s=1.5, color=BLUE, alpha=0.35, lw=0, zorder=2)
        ax.text(0.04, 0.96, f"E_F = {ef:.2f} %\nRMSE = {rmse:.1f} meV/Å", transform=ax.transAxes, va="top",
                fontsize=9, color=INK)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=9.5, color=INK)
        ax.set_xlabel("DFT force component (eV/Å)")
        ax.grid(color=GRID, lw=0.6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axs[0].set_ylabel("predicted force component (eV/Å)")
    fig.suptitle(f"Model B test set: {NDATA} structures × {NAT} atoms × 3 components", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    print("written:", out)


if __name__ == "__main__":
    main()
