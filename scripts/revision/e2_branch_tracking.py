"""E2: eigenvector-character sector assignment vs frequency-rule assignment, and the
coherence pair-sector blocks recomputed under each (8x8x8, distributed tutorial dfc2).

Inputs
  eigenvectors: $FEKAPPA_RUNS/evec8/ev_<T>K.mesh.evec  (MODE=phonons, same
                FC2XML = dfc2(T) and NONANALYTIC=3 as the RTA runs; full 8^3 mesh)
  coherence:    data/raw/alamode_sto/STO_RTA_scph_<T>K.kc_elem
  q list:       data/raw/alamode_sto/STO_RTA_scph_<T>K.result (#KPOINT block)

Sector schemes
  freq     : as submitted -- acoustic = branch 1-3; soft = branch>=4 & omega<175 cm^-1; hard = rest.
  tracked  : labels fixed by `freq` at the reference temperature T_ref, then carried to
             neighbouring T (T_ref -> up and -> down, 50 K steps) by maximum-overlap matching
             (Hungarian on |<e_i(T)|e_j(T')>|^2) at each irreducible q.  The label follows the
             eigenvector, not the frequency ordering.
  frac     : fractional character c_j(S) = sum_{i in S(T_ref)} |<e_i(T_ref)|e_j(T)>|^2 (no chaining);
             a pair (a,b) contributes kc_ab * c_a(S) * c_b(S') to block {S,S'}.

Outputs
  data/processed/revision/e2_coherence_blocks_tracked.csv   (T, scheme, T_ref, six blocks, closure)
  data/processed/revision/e2_label_changes.csv               (modes whose tracked label != freq label)
"""

from fekappa.provenance import external_dir
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fekappa.kappa_decomposition import classify_sector, read_kc_elem  # noqa: E402

RAW = ROOT / "data/raw/alamode_sto"
EV = external_dir("FEKAPPA_RUNS", "the pre-queue run directory (Zenodo prequeue_runs)") / "evec8"
OUT = ROOT / "data/processed/revision"
TEMPS = list(range(50, 1001, 50))  # 0 K excluded: SCPH unconverged and kappa_C = 0
NK, NS = 512, 15
RY_TO_CM1 = 109737.31568160
BLOCKS = ["A/A", "A/hard", "A/soft", "hard/hard", "hard/soft", "soft/soft"]
SHORT = {"acoustic": "A", "soft_optic": "soft", "hard_optic": "hard"}


def read_evec(path):
    xk, om, ev = [], [], []
    cur = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("## kpoint"):
                xk.append([float(v) for v in line.split(":")[1].split()])
            elif line.startswith("### mode"):
                om.append(float(line.split(":")[1]))
                cur = []
                ev.append(cur)
            elif cur is not None and line.strip() and not line.startswith("#"):
                re_, im_ = line.split()
                cur.append(complex(float(re_), float(im_)))
    xk = np.array(xk)
    om = np.array(om).reshape(len(xk), NS)
    ev = np.array(ev).reshape(len(xk), NS, NS)  # [k, mode, component]
    return xk, om, ev


def irred_q(result_path):
    lines = Path(result_path).read_text().splitlines()
    i = lines.index("#KPOINT")
    n = int(lines[i + 2])
    return np.array([[float(v) for v in l.split(":")[1].split()[:3]] for l in lines[i + 3:i + 3 + n]])


def kmap(xk_mesh, q_irr):
    out = []
    for q in q_irr:
        d = np.abs(((xk_mesh - q) + 0.5) % 1.0 - 0.5).sum(axis=1)
        j = int(np.argmin(d))
        assert d[j] < 1e-6, (q, d[j])
        out.append(j)
    return np.array(out)


# ---- load eigenvectors at irreducible q for every T --------------------------------------
q_irr = irred_q(RAW / "STO_RTA_scph_300K.result")
E, W = {}, {}
for t in TEMPS:
    xk, om, ev = read_evec(EV / f"ev_{t}K.mesh.evec")
    km = kmap(xk, q_irr)
    E[t], W[t] = ev[km], om[km]

# frequency sanity check vs the RTA .result of the same T (units of the evec file detected)
res300 = pd.read_csv(RAW / "STO_RTA_scph_300K.kc_elem", sep=r"\s+", comment="#", header=None,
                     names=["T", "i", "j", "ib", "jb", "ik", "w1", "w2", "re", "im"])
diag = res300[(res300.i == 1) & (res300.j == 1) & (res300.ib == res300.jb)]
w_res = np.zeros((35, NS))
w_res[diag.ik - 1, diag.ib - 1] = diag.w1
w_ev = W[300]
scale = RY_TO_CM1 if np.nanmax(np.abs(w_ev)) < 10 else 1.0
w_ev_cm = np.sign(w_ev) * np.sqrt(np.abs(w_ev)) * RY_TO_CM1 if scale != 1.0 and np.nanmax(np.abs(w_ev)) < 0.1 else w_ev * scale
mask = w_res > 1.0
rel = np.max(np.abs(w_ev_cm[mask] - w_res[mask]) / w_res[mask])
print(f"evec-vs-kc_elem frequency check at 300 K: max rel dev = {rel:.2e}")
assert rel < 1e-4, "eigenvector run does not reproduce the RTA frequencies"
for t in TEMPS:
    W[t] = (np.sign(W[t]) * np.sqrt(np.abs(W[t])) * RY_TO_CM1) if scale != 1.0 and np.nanmax(np.abs(W[t])) < 0.1 else W[t] * scale


def freq_labels(t):
    br = np.tile(np.arange(1, NS + 1), (35, 1))
    return np.vectorize(SHORT.get)(classify_sector(br.ravel(), W[t].ravel())).reshape(35, NS)


def overlap(a, b):  # |<a_i|b_j>|^2 per q -> [q, i, j]
    return np.abs(np.einsum("qic,qjc->qij", a.conj(), b)) ** 2


def tracked_labels(t_ref):
    lab = {t_ref: freq_labels(t_ref)}
    for seq in (TEMPS[TEMPS.index(t_ref):], TEMPS[:TEMPS.index(t_ref) + 1][::-1]):
        for t0, t1 in zip(seq[:-1], seq[1:]):
            O = overlap(E[t0], E[t1])
            new = np.empty((35, NS), dtype=object)
            for q in range(35):
                r, c = linear_sum_assignment(-O[q])
                new[q, c] = lab[t0][q, r]
            lab[t1] = new
    return lab


def frac_character(t_ref, t):
    ref = freq_labels(t_ref)
    O = overlap(E[t_ref], E[t])  # [q, i(ref), j(t)]
    return {S: np.einsum("qij,qi->qj", O, (ref == S).astype(float)) for S in ("A", "soft", "hard")}


def blocks_from_labels(kc, lab):
    s1 = lab[kc.ik.to_numpy() - 1, kc.ibranch.to_numpy() - 1]
    s2 = lab[kc.ik.to_numpy() - 1, kc.jbranch.to_numpy() - 1]
    key = ["/".join(sorted(p)) for p in zip(s1, s2)]
    sums = kc.assign(block=key).groupby("block")["kc_real"].sum() / 3.0 / NK
    return {b: float(sums.get(b, 0.0)) for b in BLOCKS}


def blocks_from_frac(kc, ch):
    i, a, b = kc.ik.to_numpy() - 1, kc.ibranch.to_numpy() - 1, kc.jbranch.to_numpy() - 1
    v = kc.kc_real.to_numpy() / 3.0 / NK
    out = dict.fromkeys(BLOCKS, 0.0)
    for S in ("A", "soft", "hard"):
        for S2 in ("A", "soft", "hard"):
            out["/".join(sorted((S, S2)))] += float((v * ch[S][i, a] * ch[S2][i, b]).sum())
    return out


rows, changes = [], []
tracked = {tr: tracked_labels(tr) for tr in (100, 300)}
for t in TEMPS:
    kc = read_kc_elem(RAW / f"STO_RTA_scph_{t}K.kc_elem")
    kc = kc[(kc.i == kc.j) & (kc.ibranch != kc.jbranch)]
    total = kc.kc_real.sum() / 3.0 / NK
    fl = freq_labels(t)
    rows.append({"T_K": t, "scheme": "freq", "T_ref": np.nan, **blocks_from_labels(kc, fl), "pairsum": total})
    for tr, lab in tracked.items():
        rows.append({"T_K": t, "scheme": "tracked", "T_ref": tr, **blocks_from_labels(kc, lab[t]), "pairsum": total})
        diff = np.argwhere(lab[t] != fl)
        for q, s in diff:
            changes.append({"T_K": t, "T_ref": tr, "q_index": q + 1, "branch": s + 1,
                            "omega_cm1": W[t][q, s], "freq_label": fl[q, s], "tracked_label": lab[t][q, s]})
        rows.append({"T_K": t, "scheme": "frac", "T_ref": tr, **blocks_from_frac(kc, frac_character(tr, t)), "pairsum": total})

df = pd.DataFrame(rows)
df["soft_involved"] = df["A/soft"] + df["hard/soft"] + df["soft/soft"]
df["acsoft_plus_softsoft"] = df["A/soft"] + df["soft/soft"]
df.to_csv(OUT / "e2_coherence_blocks_tracked.csv", index=False, float_format="%.6g")
pd.DataFrame(changes).to_csv(OUT / "e2_label_changes.csv", index=False, float_format="%.5g")

pd.set_option("display.width", 220)
piv = df.pivot_table(index="T_K", columns=["scheme", "T_ref"], values="acsoft_plus_softsoft", dropna=False)
print("A/soft + soft/soft block (W/mK):")
print(piv.round(4).to_string())
piv2 = df.pivot_table(index="T_K", columns=["scheme", "T_ref"], values="soft_involved", dropna=False)
print("all soft-involved (A/soft + hard/soft + soft/soft):")
print(piv2.round(4).to_string())
ch = pd.DataFrame(changes)
if len(ch):
    print(ch.groupby(["T_ref", "T_K"]).size().unstack(0).fillna(0).astype(int).to_string())


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"e2_coherence_blocks_tracked.csv", f"e2_label_changes.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/raw/alamode_sto (distributed example outputs; IMPORT_NOTE.md)', '$FEKAPPA_RUNS/evec8 (zenodo:prequeue_runs.tar.gz)'])
