#!/usr/bin/env python3
"""Rotational-family modes below 105 cm^-1 and the family's frequency range, 300 K and 600 K, both models. Labels and frequencies only: the frozen label tables on the full 16^3 mesh (61 440 modes).

Usage: rot_below_105.py LABELS_DIR
  LABELS_DIR holds labels_{A,B}_m16.csv (300 K) and T600/labels_{A,B}_m16_T600.csv (frozen label tables).
"""
import sys
from pathlib import Path

import pandas as pd

F = Path(sys.argv[1])
TAB = {("A", 300): F / "labels_A_m16.csv", ("B", 300): F / "labels_B_m16.csv",
       ("A", 600): F / "T600/labels_A_m16_T600.csv", ("B", 600): F / "T600/labels_B_m16_T600.csv"}
print("| model | T (K) | definition | modes | rot modes | rot below 105 cm⁻¹ | rot Ω range (cm⁻¹) | pol Ω range (cm⁻¹) | source |")
print("|---|---|---|---|---|---|---|---|---|")
for (m, t), p in TAB.items():
    L = pd.read_csv(p)
    for dn in ("primary_c0.5", "primary_c0.4", "primary_c0.6", "secondary_c0.5"):
        r, q = L[L[dn] == "rot"], L[L[dn] == "pol"]
        print(f"| {m} | {t} | {dn} | {len(L)} | {len(r)} | {int((r.Omega_cm1 < 105).sum())} | "
              f"{r.Omega_cm1.min():.2f}–{r.Omega_cm1.max():.2f} | {q.Omega_cm1.min():.2f}–{q.Omega_cm1.max():.2f} | "
              f"`{p.name}` |")
