"""Pair counts and per-pair contributions of the coherence blocks under both sector rules.

A pair is an ordered branch pair (s, s'), s != s', at one irreducible q point of the transport mesh (165 points on 16^3,
15 x 14 = 210 ordered pairs per point, 34650 in all), as listed in anphon's .kc_elem (diagonal Cartesian components
summed). The block of a pair is the unordered pair of sector labels, so a cross block (e.g. hard/soft) counts both orders.
Sector sizes are counted per irreducible q (unweighted, as the pairs are).
Star-weighted columns: each irreducible record stands for mult(q) full-zone pairs, mult(q) the
size of the star of q (mode table `multiplicity`; sum = N_k), so n_pairs_w = sum_q mult(q) n_block(q) and the full-zone
total is 210 N_k. kappa_block and sum_V already carry the star sum, so per full-zone pair they are divided by n_pairs_w.
  freq rule: branches 1-3 acoustic; other modes soft if Omega < 175 cm^-1, else hard (sizes vary with q)
  char rule: projection onto the Gamma subspaces, degenerate subspaces labelled as a whole (prod_analysis.py)
Inputs: e1_prod_block_factors.csv (e1_production.py), data/processed/sto_production/mode_table_T.csv.gz
Outputs: data/processed/revision/e1_pair_counts.csv, e1_sector_sizes.csv
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PJ = ROOT / "data/processed/revision"
bf = pd.read_csv(PJ / "e1_prod_block_factors.csv", comment="#")
modes = pd.read_csv(ROOT / "data/processed/sto_production/mode_table_T.csv.gz",
                    usecols=["q_index", "branch", "T_K", "mesh", "multiplicity", "sector_freq", "sector_char"])


def block_counts(d, rule):
    """ordered pairs s != s' per block at each q (from the group sizes), unweighted and star-weighted"""
    n = d.groupby("q_index")[f"sector_{rule}"].value_counts().unstack(fill_value=0)
    mult = d.groupby("q_index").multiplicity.first().loc[n.index]
    g = {S: n.get(S, pd.Series(0, index=n.index)) for S in ("A", "soft", "hard")}
    out = {}
    for a in ("A", "hard", "soft"):
        for b in ("A", "hard", "soft"):
            if a > b:
                continue
            c = g[a] * (g[a] - 1) if a == b else 2 * g[a] * g[b]
            out["/".join(sorted((a, b)))] = (int(c.sum()), int((c * mult).sum()))
    return out

sizes = []
for (T, mesh), d in modes[modes.mesh.isin([12, 16])].groupby(["T_K", "mesh"]):
    for rule in ("freq", "char"):
        n = d.groupby("q_index")[f"sector_{rule}"].value_counts().unstack(fill_value=0)
        sizes.append({"T_K": T, "mesh": mesh, "labels": rule, "n_q": len(n),
                      **{f"mean_n_{S}": float(n.get(S, pd.Series(0, index=n.index)).mean()) for S in ("A", "soft", "hard")}})
sizes = pd.DataFrame(sizes)

pc = bf[["mesh", "T_K", "labels", "block", "n_pairs", "kappa_block", "share_of_kC", "sum_V", "frac_pairs_dw_lt_G"]].copy()
W = {}
for (T, mesh), d in modes[modes.mesh.isin([12, 16])].groupby(["T_K", "mesh"]):
    assert int(d.groupby("q_index").multiplicity.first().sum()) == mesh ** 3
    for rule in ("freq", "char"):
        for blk, (n0, nw) in block_counts(d, rule).items():
            W[(mesh, T, rule, blk)] = (n0, nw)
key = list(zip(pc.mesh, pc.T_K, pc.labels, pc.block))
assert all(W[k][0] == n for k, n in zip(key, pc.n_pairs)), "group-size pair counts differ from the .kc_elem counts"
pc["n_pairs_w"] = [W[k][1] for k in key]
tot = pc.groupby(["mesh", "T_K", "labels"]).n_pairs.transform("sum")
pc["pair_frac"] = pc.n_pairs / tot
pc["kappa_per_pair"] = pc.kappa_block / pc.n_pairs
pc["per_pair_rel"] = pc.share_of_kC / pc.pair_frac   # contribution per pair relative to the all-pair average
pc["pair_frac_w"] = pc.n_pairs_w / pc.groupby(["mesh", "T_K", "labels"]).n_pairs_w.transform("sum")
pc["kappa_per_pair_w"] = pc.kappa_block / pc.n_pairs_w
pc["V_per_pair_w"] = pc.sum_V / pc.n_pairs_w
pc["per_pair_rel_w"] = pc.share_of_kC / pc.pair_frac_w
pc = pc[["mesh", "T_K", "labels", "block", "n_pairs", "pair_frac", "kappa_block", "share_of_kC", "kappa_per_pair",
         "per_pair_rel", "n_pairs_w", "pair_frac_w", "kappa_per_pair_w", "V_per_pair_w", "per_pair_rel_w",
         "frac_pairs_dw_lt_G"]]
for f, d in (("e1_pair_counts.csv", pc), ("e1_sector_sizes.csv", sizes)):
    with open(PJ / f, "w") as fh:
        fh.write("# scripts/e1_pair_counts.py; ordered pairs s != s' at the irreducible q of the mesh (n_pairs, per irreducible "
                 "record) and star-weighted (n_pairs_w, per full-zone pair); sector rules as in prod_analysis.py\n")
        d.to_csv(fh, index=False, float_format="%.5g")
pd.set_option("display.width", 220)
print(sizes[(sizes.mesh == 16) & sizes.T_K.isin([300, 1000])].round(2).to_string(index=False))
print(pc[(pc.mesh == 16) & pc.T_K.isin([300, 1000])].round(4).to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"e1_pair_counts.csv", f"e1_sector_sizes.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/processed/revision/e1_prod_block_factors.csv', 'data/processed/sto_production/mode_table_T.csv.gz'])
