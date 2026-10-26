#!/usr/bin/env python3
"""Evaluate one diagram-(c) quartet by hand and compare it with the same quartet printed by the trace build
(alamode_4ph_trace: ALAMODE c5b015c with one added print in Selfenergy::selfenergy_c; the sum is unchanged).

Usage: trace_one_quartet.py <rta.log of the trace run>

The printed line carries the four frequencies (Ry, hbar = 1), epsilon (Ry), |V4|^2 and V4 (code vertex), the
three occupations, the four frequency sums S_j, the prefactor and the Ry -> cm^-1 factor. Everything else is
recomputed here from the frequencies alone and compared with the printed values.
"""
import math
import re
import sys

line = next(l for l in open(sys.argv[1], errors="replace") if "TRACE no exclusion" in l)
tok = line.split("TRACE no exclusion;")[1].split()
v = {}
i = 0
while i < len(tok):
    if tok[i + 1] == "=":
        key = tok[i]
        if key in ("S0", "S1", "S2", "S3"):
            v[key] = complex(float(tok[i + 2]), float(tok[i + 3])); i += 4
        else:
            v[key] = float(tok[i + 2]); i += 3
    else:
        i += 1

w, w1, w2, w3, eps, T = v["omega"], v["omega1"], v["omega2"], v["omega3"], v["epsilon"], v["T"]
to_cm = v["Ry_to_kayser"]
kB = v["kB_over_Ry"]
nk = 64  # 4^3 trace run
rows = []


def check(name, hand, code, tol=1e-12):
    rel = abs(hand - code) / max(abs(code), 1e-300)
    rows.append((name, hand, code, rel))
    assert rel < tol, (name, hand, code, rel)


# occupations n(omega) = 1/(exp(omega/(kB T)) - 1), omega in Ry
n = [1.0 / math.expm1(x / (kB * T)) for x in (w1, w2, w3)]
for j, x in enumerate(n, 1):
    check(f"n{j}", x, v[f"n{j}"])

# frequency sums: S_j = L(Delta_j), L(D) = 1/(w + i eps - D) - 1/(w + i eps + D)
wt = complex(w, eps)
L = lambda D: 1 / (wt - D) - 1 / (wt + D)  # noqa: E731
D = [w1 + w2 + w3, w1 + w2 - w3, -w1 + w2 + w3, w1 - w2 + w3]
for j in range(4):
    check(f"S{j} = L(Delta{j})", L(D[j]), v[f"S{j}"])
    # closed form of Im L
    imL = -eps / ((w - D[j]) ** 2 + eps ** 2) + eps / ((w + D[j]) ** 2 + eps ** 2)
    check(f"Im L(Delta{j}) closed form", imL, L(D[j]).imag)

# thermal coefficients as written in the code, and their product forms
n1, n2, n3 = n
C = [n1 * n2 + n2 * n3 + n3 * n1 + n1 + n2 + n3 + 1,
     n3 * n1 + n2 * n3 + n3 - n1 * n2,
     n1 * n2 + n3 * n1 + n1 - n2 * n3,
     n2 * n3 + n1 * n2 + n2 - n3 * n1]
x = [wi / (kB * T) for wi in (w1, w2, w3)]
prod = [None,
        n1 * n2 * n3 * (math.exp(x[0] + x[1]) - math.exp(x[2])),
        n1 * n2 * n3 * (math.exp(x[1] + x[2]) - math.exp(x[0])),
        n1 * n2 * n3 * (math.exp(x[0] + x[2]) - math.exp(x[1]))]
for j in (1, 2, 3):
    check(f"C{j} product form", prod[j], C[j], 1e-9)
signs = [C[j] * L(D[j]).imag for j in range(4)]
assert all(s <= 0 for s in signs), signs

# vertex and prefactor
check("|V4|^2", v["V4re"] ** 2 + v["V4im"] ** 2, v["v4_tmp"])
check("prefactor -1/(96 N_k^2)", -1.0 / (96 * nk ** 2), v["factor"])
check("epsilon (cm^-1)", eps * to_cm, 5.0, 1e-12)

term = v["v4_tmp"] * sum(C[j] * L(D[j]) for j in range(4))
contrib = v["factor"] * term  # this quartet's contribution to the code quantity Sigma_code (Ry)

print("quartet: ik1 ik2 ik3 =", int(v["ik1"]), int(v["ik2"]), int(v["ik3"]),
      " is1 is2 is3 =", int(v["is1"]), int(v["is2"]), int(v["is3"]), " T =", T, "K")
print("frequencies (cm^-1): omega = %.6f  omega1..3 = %.6f %.6f %.6f" % tuple(f * to_cm for f in (w, w1, w2, w3)))
print("C_j * Im L(Delta_j) (each <= 0):", ["%.3e" % s for s in signs])
print("|V4|^2 (code vertex, Ry^2 units) = %.6e" % v["v4_tmp"])
print("contribution to Sigma_code: Re = %.6e Ry, Im = %.6e Ry  ->  Im = %.6e cm^-1" %
      (contrib.real, contrib.imag, contrib.imag * to_cm))
assert contrib.imag >= 0
print("Im Sigma_code >= 0 for this quartet; the code prints 2 * Im Sigma_code (cm^-1) as '2*Gamma4' (FWHM),")
print("so Gamma4 (HWHM) = + Im Sigma_code, i.e. Sigma_code = -Sigma^R.")
print()
print("%-28s %-26s %-26s %s" % ("check", "by hand", "printed by code", "rel. diff"))
for name, h, c, r in rows:
    print("%-28s %-26s %-26s %.1e" % (name, f"{h:.15g}", f"{c:.15g}", r))
