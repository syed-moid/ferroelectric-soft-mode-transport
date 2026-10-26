"""Derived tables and prose numbers of Supplementary Information v4, from production and revision tables only.

Inputs: data/processed/sto_production/*.csv, mode_table_T.csv.gz; data/processed/revision/*.csv (revision tables);
data/provenance/manifest_final.json (four-phonon job metrics); data/raw/alamode_sto/revision/{pm_q8,pm_q8_mix0p1,pm_q12}
(SCPH inputs/logs, read for run parameters and the mixing-parameter check); data/literature/STO/{tadano_2015,
vogt_1995}. Output: data/processed/si_v4/*.csv (one provenance line each) and data/processed/si_v4/si_numbers_v4.csv
(placeholder, value, description, source) for numbers quoted in the SI prose that are not table cells.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from fekappa.provenance import external_dir, stamp

INPUTS = ["data/processed/sto_production/", "data/processed/revision/", "data/provenance/manifest_final.json",
          "data/literature/STO/", "tools/alamode_patches/", "zenodo:pm_q12_rta.tar.gz, scph_q12.tar.gz, pm_q8.tar.gz, "
          "pm_q12_4ph_audit.tar.gz (data/raw/alamode_sto/revision/)", "$FEKAPPA_TOOLCHAIN/sto_example (STO_anharm.xml, BORN)"]

ROOT = Path(__file__).resolve().parents[2]
PROD = ROOT / "data/processed/sto_production"
PJ = ROOT / "data/processed/revision"
RAW = ROOT / "data/raw/alamode_sto/revision"
OUT = ROOT / "data/processed/si_v4"
OUT.mkdir(parents=True, exist_ok=True)
rd = lambda f: pd.read_csv(f, comment="#")
MEV = 8.0655  # cm^-1 per meV
man = json.loads((ROOT / "data/provenance/manifest_final.json").read_text())["jobs"]
km, eta, sec, coh = (rd(PROD / f) for f in ("kappa_mesh_T.csv", "eta_validity_T.csv", "kappa_sectors_T.csv",
                                             "coherence_blocks_T.csv"))
d2 = rd(PJ / "d2_summary.csv")
mt = rd(PROD / "mode_table_T.csv.gz")
NUM = []


def num(name, value, desc, src):
    NUM.append({"placeholder": name, "value": value, "description": desc, "source": src})


def write(name, df, prov):
    with open(OUT / name, "w") as f:
        f.write(f"# {prov}\n")
        df.to_csv(f, index=False)
    stamp(OUT / name, __file__, INPUTS)


def namelist(path):
    """key = value pairs of an ALAMODE input file (first occurrence wins)."""
    out = {}
    for line in Path(path).read_text().splitlines():
        for part in line.split(";"):
            if "=" in part:
                k, v = (s.strip() for s in part.split("=", 1))
                out.setdefault(k, " ".join(v.split()))
    return out


# ------------------------------------------------------------------ S1 run parameters
scph = namelist(RAW / "pm_q12/scph_pm_q12/scph.in")
rta = namelist(RAW / "pm_q12/rta_pm_q12_T300_m16/rta.in")
a4 = namelist(next((RAW / "pm_q12_4ph_audit").glob("m4ph_audit_q12_T300_m6_e5*/*.in")))
log = (RAW / "pm_q12/scph_pm_q12/scph.log").read_text()
tol = re.search(r"TOL_SCPH\s*=\s*(\S+)", log)
build = man["rta_pm_q12_T300_m16"]["build"]
abuild = man["m4ph_audit_q12_T300_m6_e5"]["build"] if "m4ph_audit_q12_T300_m6_e5" in man else ""
sha = lambda b: re.search(r"anphon-sha256:([0-9a-f]+)", b).group(1) if b else ""
# build manifest: every anphon build that ran a finished job, with its source patch and the
# checksums of the shared inputs. The job manifest records the first 12 hex digits of each binary's sha256; the full
# hashes of the node builds are added to the release manifest when the node is next reachable.
import hashlib  # noqa: E402
_h = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()
TOOL = external_dir("FEKAPPA_TOOLCHAIN", "the ALAMODE toolchain directory (sto_example/ holds the force-constant inputs)")
PATCH = ROOT / "tools/alamode_patches"
HOSTS = {man["rta_pm_q12_T300_m16"]["host"]: "HOME"}
VARIANTS = {None: ("production (unmodified c5b015c)", None),
            "alamode_4ph_audit": ("four-phonon rule: enlarged exclusion, tolerance 1e-9 Ry", "4ph_audit.diff"),
            "alamode_4ph_audit_tol1e6": ("tolerance sensitivity, 1e-6 Ry", "4ph_audit_tol1e6.diff"),
            "alamode_4ph_trace": ("one-quartet debug trace (TRACE.md)", "4ph_trace.diff"),
            "alamode_4ph_v1omit": ("validation (1): omitted-sum accumulator", "4ph_v1omit.diff"),
            "alamode_4ph_v2strict": ("validation (2): strict-index exclusion", "4ph_v2strict.diff"),
            "alamode_4ph_v3dump": ("validation (3): full-precision dump", "4ph_v3dump.diff")}
_bm = {}
for jid, j in man.items():
    b = j.get("build") or ""
    if j.get("status") != "done" or "anphon-sha256" not in b:
        continue
    tree = (j.get("params") or {}).get("anphon_tree")
    k = (tree, sha(b), HOSTS.get(j.get("host"), "NODE"))
    _bm.setdefault(k, []).append(jid)
bm = [{"build": VARIANTS[t][0], "tree": t or "alamode", "anphon_sha256_12": h, "host": host, "jobs": len(v),
       "patch": VARIANTS[t][1] or "", "patch_sha256": _h(PATCH / VARIANTS[t][1]) if VARIANTS[t][1] else ""}
      for (t, h, host), v in sorted(_bm.items(), key=lambda x: (x[0][0] or "", x[0][2]))]
bm += [{"build": "source", "tree": "ALAMODE", "anphon_sha256_12": "", "host": "", "jobs": "",
        "patch": "commit c5b015caaf72c1a86a5250b6e34a86c0e09c7986", "patch_sha256": ""},
       {"build": "input", "tree": "STO_anharm.xml (force constants)", "anphon_sha256_12": "", "host": "", "jobs": "",
        "patch": "", "patch_sha256": _h(TOOL / "sto_example/STO_anharm.xml")},
       {"build": "input", "tree": "BORN (Born charges, dielectric tensor)", "anphon_sha256_12": "", "host": "",
        "jobs": "", "patch": "", "patch_sha256": _h(TOOL / "sto_example/BORN")}]
write("build_manifest.csv", pd.DataFrame(bm), "anphon builds of every finished job (job manifest build ids), source patches "
      "(tools/alamode_patches) and input checksums; per-job input/output sha256 in data/provenance/manifest_final.json")
meshes = km.groupby("T_K").mesh.apply(lambda s: ", ".join(f"{m}³" for m in sorted(s))).to_dict()
p = [
    ("code", "ALAMODE (anphon), commit c5b015c", "all", "build id"),
    ("source commit", "c5b015caaf72c1a86a5250b6e34a86c0e09c7986", "all", "toolchain"),
    *[(f"anphon binary, {r['build']} (sha256, first 12)", f"{r['anphon_sha256_12']} ({r['host']}, {r['jobs']} jobs)",
       "SCPH, RTA" if r["tree"] == "alamode" else "diagram (c)", "job manifest")
      for r in bm if r["anphon_sha256_12"]],
    ("checksums", "force constants STO_anharm.xml sha256 " + bm[-2]["patch_sha256"][:12] + "…, BORN "
     + bm[-1]["patch_sha256"][:12] + "…; patches, full hashes and per-job input/output sha256 in the build manifest of "
     "the release", "all", "build_manifest.csv"),
    ("force constants", f"{scph['FCSXML']} (to fourth order)", "all", "scph.in"),
    ("lattice parameter", f"{rta.get('7.363', '7.363')} bohr (cubic, fixed)", "all", "rta.in &cell"),
    ("non-analytic correction", f"NONANALYTIC = {scph['NONANALYTIC']} (BORN)", "all", "scph.in"),
    ("SCPH temperatures", f"{scph['TMIN']}–{scph['TMAX']} K in steps of {scph['DT']} K", "SCPH", "scph.in"),
    ("SCPH grid q₁ (KMESH_SCPH)", "12³ (production); 8³, 16³ (checks)", "SCPH", "scph.in"),
    ("SCPH interpolation grid (KMESH_INTERPOLATE)", scph["KMESH_INTERPOLATE"].replace(" ", "×"), "SCPH", "scph.in"),
    ("polarization mixing (SELF_OFFDIAG)", scph["SELF_OFFDIAG"], "SCPH", "scph.in"),
    ("SCPH mixing parameter (MIXALPHA)", scph["MIXALPHA"], "SCPH", "scph.in"),
    ("SCPH maximum iterations (MAXITER)", scph["MAXITER"], "SCPH", "scph.in"),
    ("SCPH tolerance (TOL_SCPH)", tol.group(1) if tol else "1e-10 (default)", "SCPH", "scph.log"),
    ("transport", "Wigner RTA: particle-like κ_P and coherence κ_C (KAPPA_COHERENT = "
                  f"{rta['KAPPA_COHERENT']})", "RTA", "rta.in"),
    ("Brillouin-zone integration", "tetrahedron (ISMEAR = -1, default)", "RTA", "rta.in (keyword absent)"),
    ("isotope scattering", "none (ISOTOPE absent, default off)", "RTA", "rta.in (keyword absent)"),
    ("transport meshes", "8³, 12³, 16³ at every T; 20³ at 300, 400, 500, 650, 800, 1000 K; "
                         "15³, 17³, 24³ at 300 K", "RTA", "kappa_mesh_T.csv"),
    ("four-phonon self-energy", f"QUARTIC = {a4['QUARTIC']}, SELF_ENERGY = {a4['SELF_ENERGY']} (diagram (c)); "
                                "Lorentzian half-width ε = 2, 5, 10 cm⁻¹", "diagram (c)", "audit input"),
    ("four-phonon meshes", "4³ (gate), 6³, 8³ (rule); external q on the mesh", "diagram (c)", "audit input"),
    ("identity tolerance", "|Δω| < 1e-9 Ry (sensitivity variant 1e-6 Ry)", "diagram (c)", "audit patch"),
]
# draft 2: settings inherited from the generation of the force constants (distributed example), separated from the settings
# chosen or tested here
p = [["inherited: " + r[0], *r[1:]] if r[0] in ("force constants", "lattice parameter") else r for r in p]
p = [["inherited: DFT functional and code", "PBEsol, plane-wave cutoff 550 eV (example documentation, not in the "
      "distributed input files)", "force constants", "ALAMODE documentation sto_scph"],
     ["inherited: supercell", "2×2×2 (40 atoms), a = 14.726 bohr", "force constants", "example opt.in &cell"],
     ["inherited: fit", "harmonic: least squares; cubic to sixth order: elastic net (L1 penalty 3.4548 × 10⁻⁶) on 40 "
      "displacement–force sets from molecular dynamics and random displacements; NBODY = 2 3 3 2 2", "force constants", "example opt.in"],
     ["inherited: interaction cutoffs", "harmonic and cubic: none (all pairs in the supercell); fourth to sixth order: "
      "12 bohr", "force constants", "example opt.in &cutoff"],
     ["inherited: Born charges and dielectric tensor", "BORN file of the example (computation settings not documented)",
      "non-analytic term", "example BORN"]] + p
write("s1_parameters.csv", pd.DataFrame(p, columns=["parameter", "value", "applies_to", "source"]),
      "run parameters, read from the archived inputs and logs (pm_q12)")
assert meshes[300].startswith("8³")

rpc = rd(PJ / "r_point_check.csv")
conv = (rpc[rpc.converged].groupby(["PM", "q1"]).T_K.agg(["min", "max", "count"]).reset_index())
nconv = rpc[~rpc.converged].groupby(["PM", "q1"]).T_K.apply(lambda s: ", ".join(str(int(t)) for t in sorted(s)))
conv["not_converged_K"] = [nconv.get((a, b), "") for a, b in zip(conv.PM, conv.q1)]
write("s1_scph_convergence.csv", conv.rename(columns={"min": "T_min_converged", "max": "T_max_converged",
                                                      "count": "n_converged"}),
      "temperatures with 'convergence achieved' in scph.log, per SCPH run (r_point_check.csv)")

# ------------------------------------------------------------------ S3 pair diagnostic
r_med = d2.set_index(["T_K", "mesh"]).kw_median_r
kcp = d2.set_index(["T_K", "mesh"]).kappa_C
pairs = [(300, 16, 20), (400, 16, 20), (300, 20, 24), (500, 16, 20), (650, 16, 20), (800, 16, 20), (1000, 16, 20)]
rows = []
for T, a, b in pairs:
    rows.append({"T_K": T, "pair": f"{a}³→{b}³", "r_coarser": round(r_med[(T, a)], 2),
                 "dkC_pct": round(100 * (kcp[(T, b)] / kcp[(T, a)] - 1), 1),
                 "registered": {(300, 20, 24): "prediction 1 (failed)", (1000, 16, 20): "prediction 2 (held)"}
                 .get((T, a, b), "")})
s3 = pd.DataFrame(rows).sort_values("r_coarser")
write("s3_pairs.csv", s3, "κ_C change between meshes against r of the coarser mesh (d2_summary.csv; pair sums)")
num("dkC_1000_solver", round(100 * (k := km.set_index(["T_K", "mesh"]).kappa_C)[(1000, 20)] / k[(1000, 16)] - 100, 2),
    "Δκ_C(16³→20³, 1000 K) from the solver values", "kappa_mesh_T.csv")
num("dkC_1000_pairsum", round(100 * (kcp[(1000, 20)] / kcp[(1000, 16)] - 1), 2),
    "Δκ_C(16³→20³, 1000 K) from the pair sums (registered basis)", "d2_summary.csv")
num("kC_1000_16_pairsum", round(kcp[(1000, 16)], 5), "κ_C(16³, 1000 K), pair sum", "d2_summary.csv")
num("kC_1000_20_pairsum", round(kcp[(1000, 20)], 5), "κ_C(20³, 1000 K), pair sum", "d2_summary.csv")
num("dkC_300_20_24", round(100 * (kcp[(300, 24)] / kcp[(300, 20)] - 1), 1), "Δκ_C(20³→24³, 300 K)", "d2_summary.csv")
num("r_300_24", round(r_med[(300, 24)], 2), "κ_C-weighted median r at 24³, 300 K", "d2_summary.csv")
num("r_1000_20", round(r_med[(1000, 20)], 2), "κ_C-weighted median r at 20³, 1000 K", "d2_summary.csv")

# ------------------------------------------------------------------ S4 250 K boundary
n3 = rd(PJ / "n3_q1_16_check.csv")
n3 = n3[n3.T_K.isin([250, 300, 500, 800])].copy()   # 500 and 800 K added in the revision
LAB = {"omega_Gamma": "Ω Γ-TO1 (cm⁻¹)", "omega_X": "Ω X (cm⁻¹)", "omega_M": "Ω M (cm⁻¹)", "omega_R": "Ω R (cm⁻¹)",
       "kappa_P": "κ_P, 12³ mesh (W m⁻¹ K⁻¹)", "kappa_C": "κ_C, 12³ mesh (W m⁻¹ K⁻¹)"}
n3 = n3[n3.quantity.isin(LAB)].assign(quantity=lambda d: d.quantity.map(LAB))
lim = {"κ_P, 12³ mesh (W m⁻¹ K⁻¹)": 3.0, "Ω R (cm⁻¹)": 5.0, "Ω M (cm⁻¹)": 5.0}
n3["D0_limit_pct"] = n3.quantity.map(lim)
n3["D0_met"] = [("" if np.isnan(l) else ("yes" if abs(s) < l else "no")) for s, l in zip(n3.shift_pct, n3.D0_limit_pct)]
n3["shift_pct"] = n3.shift_pct.round(2)
write("s4_boundary.csv", n3, "SCPH q₁ = 12³ vs 16³ at 250, 300, 500 and 800 K, transport on 12³; D0 rule (n3_q1_16_check.csv)")

# ------------------------------------------------------------------ S5 labelling
c = coh[coh.T_K >= 300].copy()
c["hh_share_pct"] = 100 * c.kC_hard_hard / c.kappa_C_pairsum
hh = c.pivot_table(index=["T_K", "mesh"], columns="labels", values="hh_share_pct").reset_index()
hh24 = d2[d2.mesh == 24][["T_K", "mesh", "hard_hard_share"]].rename(columns={"hard_hard_share": "freq"})
hh24["freq"] *= 100
hh = pd.concat([hh, hh24.assign(char=np.nan)]).sort_values(["T_K", "mesh"])
hh = hh.rename(columns={"char": "hh_char_pct", "freq": "hh_freq_pct"}).round(1)
write("s5_hh_share.csv", hh[["T_K", "mesh", "hh_char_pct", "hh_freq_pct"]],
      "hard/hard share of κ_C under both labelling rules (coherence_blocks_T.csv; 24³ frequency rule from d2_summary.csv)")
s = sec[(sec.mesh == 16) & (sec.T_K >= 300) & (sec.T_K % 100 == 0)].copy()
for col, n in (("kappa_P_A", "A"), ("kappa_P_soft", "soft"), ("kappa_P_hard", "hard")):
    s[f"{n}_pct"] = (100 * s[col] / s.kappa_P).round(1)
write("s5_particle_sectors.csv", s[["T_K", "labels", "A_pct", "soft_pct", "hard_pct"]],
      "particle sectors, share of κ_P on 16³ under both rules (kappa_sectors_T.csv)")
R = mt[(mt.q_frac_x.abs() == 0.5) & (mt.q_frac_y.abs() == 0.5) & (mt.q_frac_z.abs() == 0.5) & (mt.branch == 1)
       & (mt.mesh == 16) & mt.T_K.isin([300, 650, 1000])]
write("s5_r_mode_label.csv", R[["T_K", "omega_cm1", "c_A", "c_soft", "c_hard", "sector_char", "sector_freq"]].round(3),
      "R-point lowest mode on 16³: eigenvector character weights and the label under each rule (mode_table_T.csv.gz)")
ff = rd(PJ / "e1_prod_factor_freezing.csv")
ff = ff[(ff.mesh == 16) & (ff.labels == "char") & (ff.block == "hard/hard") & (ff.T_ref == 300)
        & ff.T_K.isin([300, 500, 650, 800, 1000])]
ffp = ff.pivot_table(index="T_K", columns="frozen", values="kappa_block").reset_index().round(3)
write("s5_factor_freezing.csv", ffp[["T_K", "none", "H", "L", "V"]],
      "hard/hard κ_C (16³, character rule) with one factor of each pair held at its 300 K value (e1_prod_factor_freezing.csv)")
bf = rd(PJ / "e1_prod_block_factors.csv")
bf = bf[(bf.mesh == 16) & bf.T_K.isin([300, 1000])]
pcn = rd(PJ / "e1_pair_counts.csv")[["mesh", "T_K", "labels", "block", "pair_frac", "per_pair_rel", "n_pairs_w",
                                     "pair_frac_w", "per_pair_rel_w"]]   # star-weighted columns
bf = bf.merge(pcn, on=["mesh", "T_K", "labels", "block"], how="left")
bf["labels"] = pd.Categorical(bf.labels, ["char", "freq"])
bf = bf.sort_values(["labels", "T_K", "block"])
write("s5_block_factors.csv", bf[["labels", "T_K", "block", "kappa_block", "share_of_kC", "n_pairs", "pair_frac", "per_pair_rel",
                                  "n_pairs_w", "pair_frac_w", "per_pair_rel_w", "kw_mean_abs_dw_cm1", "kw_mean_G1pG2_cm1", "kw_mean_H", "sum_V", "frac_pairs_dw_lt_G"]].round(4),
      "coherence blocks on 16³, both rules: size, pair counts and κ-weighted pair factors (e1_prod_block_factors.csv, "
      "e1_pair_counts.csv)")
sz = rd(PJ / "e1_sector_sizes.csv")
write("s5_sector_sizes.csv", sz[(sz.mesh == 16) & sz.T_K.isin([300, 1000])].round(2),
      "mean number of modes per irreducible q in each group, 16³, both rules (e1_sector_sizes.csv)")

# ------------------------------------------------------------------ S6 damping under both rules
rows = []
for T in (300, 650, 1000):
    for m in sorted(mt[mt.T_K == T].mesh.unique()):
        x = mt[(mt.T_K == T) & (mt.mesh == m) & (mt.omega_cm1 > 1)]
        for rule in ("char", "freq"):
            for sname in ("A", "soft", "hard"):
                y = x[x[f"sector_{rule}"] == sname]
                if y.empty:
                    continue
                w = y.multiplicity
                rows.append({"T_K": T, "mesh": m, "rule": rule, "sector": sname,
                             "eta_max": round(y.eta.max(), 3),
                             "eta_median_mesh": round(float(np.median(np.repeat(y.eta, w))), 4),
                             "frac_eta_gt_0p05_mesh_pct": round(100 * w[y.eta > 0.05].sum() / w.sum(), 2)})
write("s6_sector_eta.csv", pd.DataFrame(rows), "damping ratio η = Γ/Ω per sector under both rules (mode_table_T.csv.gz)")
vg = rd(ROOT / "data/literature/STO/vogt_1995/gamma_T.csv"); vg = vg[vg.phase_mode == "cubic_F1u"]
vo = rd(ROOT / "data/literature/STO/vogt_1995/omega_T.csv"); vo = vo[vo.phase_mode == "cubic_F1u"]
x16 = mt[(mt.T_K == 300) & (mt.mesh == 16)]
to1 = x16[(x16.q_frac_x.abs() < 1e-6) & (x16.q_frac_y.abs() < 1e-6) & (x16.q_frac_z.abs() < 1e-6) & (x16.branch == 4)].iloc[0]
g_exp = float(np.interp(300, vg.T_K, vg.gamma_meV)) * MEV
o_exp = float(np.interp(300, vo.T_K, vo.Omega_meV)) * MEV
soft = (x16.sector_char == "soft").to_numpy()
rows = []
for conv_, f in (("full width (γ = 2Γ), adopted", 0.5), ("half width (γ = Γ), sensitivity", 1.0)):
    r = g_exp * f / to1.gamma_hwhm_cm1
    rows.append({"reading": conv_, "gamma_exp_cm1": round(g_exp, 2), "Gamma_exp_hwhm_cm1": round(g_exp * f, 2),
                 "Gamma_calc_hwhm_cm1": round(to1.gamma_hwhm_cm1, 3), "ratio": round(r, 1),
                 "eta_max_rescaled": round(float(np.nanmax(np.where(soft, r * x16.eta, x16.eta))), 2),
                 "eta_max_rescaled_with_Omega": round(float(np.nanmax(np.where(soft, r * x16.eta * to1.omega_cm1 / o_exp,
                                                                                x16.eta))), 2)})
write("s6_width.csv", pd.DataFrame(rows),
      "soft-sector linewidth rescaled to Vogt's zone-centre width at 300 K, 16³ (vogt_1995; mode_table_T.csv.gz)")
num("Omega_vogt_300", round(o_exp, 1), "Vogt 1995 cubic F1u frequency at 300 K (interpolated), cm⁻¹", "vogt_1995/omega_T.csv")
num("Omega_to1_300", round(to1.omega_cm1, 1), "our Γ-TO1 frequency, 300 K, 16³, cm⁻¹", "mode_table_T.csv.gz")

# ------------------------------------------------------------------ S7 four-phonon
a1 = rd(PJ / "a1_per_mesh.csv")
rows = []
MODEKEY = {"Gamma-TO1": "0,0,0,4", "M": "0.5,0.5,0,1", "R": "0.5,0.5,0.5,1"}
for mode, key in MODEKEY.items():
    un = {e: man[f"m4ph_pm_q12_T300_m4_e{e}"]["metrics"]["modes"][key]["Gamma4_hwhm"][0] for e in (2, 5, 10)}
    e = np.array([2, 5, 10.0]); y = np.array([un[2], un[5], un[10]])
    (a, B), *_ = np.linalg.lstsq(np.c_[1 / e, np.ones(3)], y, rcond=None)
    for eps in (2, 5, 10):
        pat = float(a1[(a1["mode"] == mode) & (a1.mesh == 4) & (a1.eps == eps)].Gamma4.iloc[0])
        rows.append({"mode": mode, "eps": eps, "G4_unpatched": round(un[eps], 2), "G4_patched": pat,
                     "removed": round(un[eps] - pat, 2), "fit_a_over_eps": round(a / eps, 2), "fit_B": round(B, 2),
                     "removed_over_fit_pct": round(100 * (un[eps] - pat) / (a / eps), 1)})
write("s7_gate_4cubed.csv", pd.DataFrame(rows),
      "4³, 300 K: diagram (c) Γ₄ (HWHM, cm⁻¹) of the unmodified and the audit build; fit Γ₄ = a/ε + B to the unmodified values "
      "(manifest jobs m4ph_pm_q12_T300_m4_e*, a1_per_mesh.csv)")
g = pd.DataFrame(rows)
num("g4_span_unpatched", f"{min(g.groupby('mode').G4_unpatched.max() / g.groupby('mode').G4_unpatched.min()):.1f}–"
    f"{max(g.groupby('mode').G4_unpatched.max() / g.groupby('mode').G4_unpatched.min()):.1f}",
    "max/min of unmodified Γ₄ over ε = 2–10 cm⁻¹ at 4³", "s7_gate_4cubed.csv")
sp = g.groupby("mode").G4_patched.agg(lambda s: 100 * (s.max() - s.min()) / s.mean()).round(0)
num("g4_span_patched", "; ".join(f"{m}: {sp[m]:.0f}%" for m in ("Gamma-TO1", "R", "M")),
    "range/mean of the audit-build Γ₄ over ε at 4³", "s7_gate_4cubed.csv")
un6 = man["m4ph_pm_q12_T300_m6_e5_timing"]["metrics"]["modes"]["0,0,0,4"]["Gamma4_hwhm"][0]
num("g4_unpatched_6_to1", round(un6, 3), "Γ₄ of Γ-TO1, unmodified code, 6³, ε = 5, 300 K", "manifest m4ph_pm_q12_T300_m6_e5_timing")
num("g4_removed_6_to1", round(un6 - float(a1[(a1['mode'] == 'Gamma-TO1') & (a1.mesh == 6) & (a1.eps == 5)].Gamma4.iloc[0]), 3),
    "excluded-set contribution at 6³ (Γ-TO1, ε = 5)", "manifest; a1_per_mesh.csv")
rows = []
for mode, key in MODEKEY.items():
    for eps in (2, 5):
        j = man.get(f"m4ph_audittol6_q12_T300_m6_e{eps}")
        if not j:
            continue
        g6 = j["metrics"]["modes"][key]["Gamma4_hwhm"][0]
        base = a1[(a1["mode"] == mode) & (a1.mesh == 6) & (a1.eps == eps)].iloc[0]
        # both columns from the full printed widths (2Γ, 6 digits); a1_per_mesh.csv is rounded to 3 decimals, which
        # produced a spurious +0.01% for M at ε = 5 in SI draft 3
        g9 = man[f"m4ph_audit_q12_T300_m6_e{eps}"]["metrics"]["modes"][key]["Gamma4_hwhm"][0]
        excl = [int(re.search(r"excluded = (\d+)", s).group(1)) for s in j["metrics"].get("audit_4ph_c", [])]
        order = [",".join(l.split()) for l in (RAW / "pm_q12_4ph_audit_tol1e6" / f"m4ph_audittol6_q12_T300_m6_e{eps}" /
                                               "modes.txt").read_text().splitlines()[1:]]  # audit lines follow modes.txt
        rows.append({"mode": mode, "eps": eps, "G4_tol1e9": f"{g9:.6g}", "G4_tol1e6": f"{g6:.6g}",
                     "change_pct": f"{100 * (g6 / g9 - 1):+.3f}", "excluded_tol1e9": int(base.excluded),
                     "excluded_tol1e6": excl[order.index(key)]})
write("s7_tolerance.csv", pd.DataFrame(rows),
      "6³, 300 K: identity tolerance 1e-9 Ry (rule) vs 1e-6 Ry (sensitivity) (manifest m4ph_audittol6_*, a1_per_mesh.csv)")
# Table S18 (draft 2): one formula for every row. Gamma3_same = three-phonon width on the four-phonon mesh at
# the same settings; Gamma3_prod = production 16^3 width; eta sums use Gamma3_prod (eta_34_same also given).
cs = rd(PJ / "eta4_margin_cases.csv")
cs = cs[cs.group.str.startswith(("soft", "hard optic"))].copy()
a1all = pd.concat([a1.assign(T_K=300), rd(PJ / "a1_per_mesh_T1000.csv").assign(T_K=1000)])
prod = a1all.drop_duplicates(["mode", "T_K"]).set_index(["mode", "T_K"]).Gamma3_prod16
m16 = mt[(mt.T_K == 1000) & (mt.mesh == 16)]
qs = np.sort(np.abs(m16[["q_frac_x", "q_frac_y", "q_frac_z"]].to_numpy()), axis=1)
off = m16[np.isclose(qs, [0.25, 0.5, 0.5]).all(axis=1) & (m16.branch == 13)].iloc[0]
prod[("(1/4,1/2,1/2) b13", 1000)] = off.gamma_hwhm_cm1
cs["Gamma3_prod16"] = [prod[(m, T)] for m, T in zip(cs["mode"], cs.T_K)]
cs["eta3_prod"] = cs.Gamma3_prod16 / cs.omega_cm1
cs["eta4"] = cs.Gamma4 / cs.omega_cm1
cs["eta34_prod"] = cs.eta3_prod + cs.eta4
cs["eta34_same"] = (cs.Gamma3 + cs.Gamma4) / cs.omega_cm1
cs["rule"] = np.where(cs.group.str.contains("descriptive"), "descriptive", "decision rule")
s18 = cs[["T_K", "mode", "rule", "mesh", "eps", "omega_cm1", "Gamma3", "Gamma3_prod16", "Gamma4", "eta3_prod", "eta4",
          "eta34_prod", "eta34_same"]].sort_values(["T_K", "rule", "mode", "mesh", "eps"])
write("s7_widths.csv", s18.round({"omega_cm1": 2, "Gamma3": 3, "Gamma3_prod16": 3, "Gamma4": 3, "eta3_prod": 4, "eta4": 4,
                                  "eta34_prod": 4, "eta34_same": 4}),
      "four-phonon widths, one formula: η₃ = Γ₃(16³ production)/Ω, η₄ = Γ₄/Ω, η₃ + η₄; η₃₊₄ at the same settings uses Γ₃ "
      "on the four-phonon mesh (eta4_margin_cases.csv, a1_per_mesh*.csv, mode_table_T.csv.gz)")
mg = s18[s18.rule == "decision rule"].groupby(["T_K", "mode"]).agg(omega=("omega_cm1", "first"),
                                                                  G3p=("Gamma3_prod16", "first"), G4max=("Gamma4", "max"))
mg["allowance_half"] = 0.5 * (0.5 * mg.omega - mg.G3p)
mg["holds"] = mg.G4max <= mg.allowance_half
mg["eta34_max"] = (mg.G3p + mg.G4max) / mg.omega
write("s7_margin_inequality.csv", mg.reset_index().round(4),
      "inequality Γ₄,max ≤ ½(0.5Ω − Γ₃) with Γ₃ the 16³ production width; Γ₄,max over the 6³ (ε = 2, 5, 10) and 8³ (ε = 5) "
      "runs; η₃ + η₄,max (s7_widths.csv)")
num("off_sym_mode", f"{off.omega_cm1:.1f}; {off.gamma_hwhm_cm1:.2f}", "(1/4,1/2,1/2) b13, 16³, 1000 K: Ω; Γ₃ (cm⁻¹)",
    "mode_table_T.csv.gz")
n4 = rd(PJ / "n4_results.csv")
n4s = n4[n4.eps == 5]
rows = []
for (smp, cls), y in n4s.groupby(["sample", "cls"], sort=False):
    r = {"sample": f"{smp}³", "class_300K": cls}
    for T in (300, 800):
        z = y[y["T"] == T]
        r[f"G4_over_G3_{T}"] = f"{z.r.min():.2f}–{z.r.max():.2f}"
        r[f"lifetime_factor_{T}"] = f"{z.life.min():.2f}–{z.life.max():.2f}"
    rows.append(r)
write("s7_n4_classes.csv", pd.DataFrame(rows).sort_values("sample", ascending=False),
      "heat-carrier samples, ε = 5 cm⁻¹, two modes per row: Γ₄/Γ₃ and 1/(1 + Γ₄/Γ₃) at the same settings (n4_results.csv)")
num("n4_ratio_300", f"{n4[n4['T'] == 300].r.min():.2f}–{n4[n4['T'] == 300].r.max():.2f}", "Γ₄/Γ₃, all twelve modes, all ε, 300 K",
    "n4_results.csv")
num("n4_ratio_800", f"{n4[n4['T'] == 800].r.min():.2f}–{n4[n4['T'] == 800].r.max():.2f}", "Γ₄/Γ₃, all twelve modes, 800 K",
    "n4_results.csv")
sh = lambda T: (100 * (1 - n4[n4['T'] == T].life)).agg(["min", "max"]).round(0)
num("n4_shorten_300", f"{sh(300)['min']:.0f}–{sh(300)['max']:.0f}", "lifetime shortening 1 − 1/(1+Γ₄/Γ₃), %, 300 K", "n4_results.csv")
num("n4_shorten_800", f"{sh(800)['min']:.0f}–{sh(800)['max']:.0f}", "lifetime shortening, %, 800 K", "n4_results.csv")
e6 = n4[(n4["sample"] == 6) & (n4["T"] == 300)].groupby(["q", "b"]).G4.agg(lambda s: s.max() / s.min())
num("n4_eps_span", f"{e6.min():.2f}–{e6.max():.2f}", "6³ sample, 300 K: max/min Γ₄ over ε = 2, 5, 10 per mode", "n4_results.csv")
sel = rd(PJ / "n4_selection.csv")
for smp in ("8^3", "6^3"):
    num(f"n4_share_{smp[0]}", round(sel[sel["sample"] == smp].kappa_share_pct.sum(), 2),
        f"share of κ_P (300 K) carried by the {smp} sample, %", "n4_selection.csv")

# ------------------------------------------------------------------ S7b R point
tad = rd(ROOT / "data/literature/STO/tadano_2015/scph_convergence_300K.csv")
tad = tad[tad.method == "SCPH"].assign(q1=lambda d: d.branch.str.extract(r"_(\d+)x").astype(int)[0])
rows = []
for (pm, q1), y in rpc[rpc.T_K == 300].groupby(["PM", "q1"]):
    y = y.iloc[0]
    r = {"q1": q1, "PM": pm}
    for lab, tl in (("Gamma", "Gamma"), ("R", "R"), ("M", "M")):
        t = tad[(tad.q_label == tl) & (tad.q1 == q1)]
        tv = float(t.omega_meV.iloc[0]) * MEV if len(t) else np.nan   # unrounded
        r[f"{lab}_ours"] = round(float(y[lab]), 2)
        r[f"{lab}_tadano"] = tv
        r[f"{lab}_diff_pct"] = 100 * (y[lab] / tv - 1) if tv == tv else np.nan
    rows.append(r)
s7b = pd.DataFrame(rows).sort_values(["PM", "q1"], ascending=[False, True])
# every benchmark percentage of main-text Sec. 3.1 and SI S7b/Table S27 comes from these unrounded differences
_on, _off = s7b[s7b.PM == "on"].dropna(), s7b[s7b.PM == "off"].dropna()
_R = _on.set_index("q1").R_diff_pct
_p1 = lambda v: f"{v:.1f}"
num("bench_R2", f"{_R[2]:.1f}; {round(abs(_R[2]))}", "R at q1 = 2^3 with mixing: difference from Tadano 2015 Table I, %; "
    "same rounded to an integer ('about')", "s7b_tadano.csv")
num("bench_R_q1", f"{_R[4]:.1f}; {_R[8]:.1f}; {_R[12]:.1f}", "R difference at q1 = 4^3, 8^3, 12^3 with mixing, %", "s7b_tadano.csv")
num("bench_R_ge8", f"{round(abs(_R[_R.index >= 8]).min())}–{round(abs(_R[_R.index >= 8]).max())}",
    "R softer than Tadano 2015 at q1 >= 8^3 with mixing, %, integer range", "s7b_tadano.csv")
_G, _M = abs(_on.Gamma_diff_pct).max(), abs(_on.M_diff_pct).max()
num("bench_G_max", f"{_G:.1f}; {round(2 * _G) / 2:g}", "largest |Gamma difference| at any q1 with mixing, %; same to the "
    "nearest 0.5 ('about')", "s7b_tadano.csv")
num("bench_G_8", f"{abs(_on.set_index('q1').Gamma_diff_pct[8]):.1f}", "|Gamma difference| at q1 = 8^3 with mixing, %",
    "s7b_tadano.csv")
num("bench_M_max", f"{_M:.1f}; {round(_M)}", "largest |M difference| at any q1 with mixing, %; same to an integer",
    "s7b_tadano.csv")
num("bench_G_off", f"{_off.Gamma_diff_pct.min():.1f}–{_off.Gamma_diff_pct.max():.1f}",
    "Gamma difference without mixing, %, range over q1", "s7b_tadano.csv")
for c in s7b.columns:
    if c.endswith("_tadano") or c.endswith("_diff_pct"):
        s7b[c] = s7b[c].round(4)
write("s7b_tadano.csv", s7b, "lowest SCPH frequencies at 300 K (cm⁻¹) against Tadano & Tsuneyuki 2015 Table I "
      "(meV × 8.0655, not rounded; differences from the unrounded values) (r_point_check.csv; "
      "data/literature/STO/tadano_2015/scph_convergence_300K.csv)")
fq = pd.read_csv(ROOT / "data/literature/STO/fauque_2022/dho_fits_T.csv")   # tracked; history in its README.md
fq = fq[(fq["mode"] == "AFD") & (fq["method"] == "experimental") & fq.T_K.isin([300, 600, 1000])]
ours = rpc[(rpc.PM == "on") & (rpc.q1 == 12)].set_index("T_K").R
ex = pd.DataFrame({"T_K": fq.T_K.to_numpy(), "R_exp_meV": fq.omega_meV.to_numpy(),
                   "R_exp_cm_read": (fq.omega_meV * MEV).round(1).to_numpy(), "R_exp_unc_cm": 5.0,
                   "R_ours_cm": ours.loc[fq.T_K].round(1).to_numpy()})
ex["R_deficit_pct"] = (100 * (ex.R_ours_cm / ex.R_exp_cm_read - 1)).round(1)
write("s7b_experiment.csv", ex, "R-point frequency, experiment read from Fauqué et al. Fig. 1(d) (uncertainty = marker "
      "half-size 0.62 meV) against ours (q₁ = 12³, PM) (data/literature/STO/fauque_2022/dho_fits_T.csv; r_point_check.csv)")
rows = []
for q1 in (8, 12, 16):
    y = rpc[(rpc.PM == "on") & (rpc.q1 == q1) & rpc.converged & rpc.T_K.between(300, 600)] if q1 != 8 else None
    if q1 == 8:
        y = rpc[(rpc.PM == "on") & (rpc.q1 == 8) & rpc.converged & rpc.T_K.between(300, 600)]
    A, b = np.polyfit(y.T_K, y.R ** 2, 1)
    rows.append({"q1": q1, "n_T": len(y), "A_cm2_per_K": round(A, 3), "Tc_K": round(-b / A)})
write("s7b_tc.csv", pd.DataFrame(rows), "mean-field fit Ω_R² = A(T − T_c) over the converged 300–600 K points (r_point_check.csv)")


def scph_R(d):
    f = next(d.glob("*.scph_bands"))
    head = [float(x) for x in f.read_text().splitlines()[1].lstrip("#").split()]
    a = np.loadtxt(f)
    it = {float(m.group(1)): None for m in re.finditer(r"Temperature = (\S+) K", (d / "scph.log").read_text())}
    return {T: a[(a[:, 0] == T)][np.argmin(np.abs(a[a[:, 0] == T][:, 1] - head[4]))][2] for T in sorted(set(a[:, 0]))}


r2, r1 = scph_R(RAW / "pm_q8/scph_pm_q8"), scph_R(RAW / "pm_q8_mix0p1/scph_pm_q8_mix0p1")
mix = pd.DataFrame([{"T_K": T, "R_alpha0p2": round(r2[T], 4), "R_alpha0p1": round(r1[T], 4),
                     "abs_diff_cm1": abs(round(r2[T] - r1[T], 4))} for T in sorted(r1) if 300 <= T <= 1000])
write("s7b_mixalpha.csv", mix, "R-point frequency (cm⁻¹), SCPH q₁ = 8³, PM, MIXALPHA 0.2 vs 0.1 (scph_bands of both runs)")
num("mix_maxdiff", f"{mix.abs_diff_cm1.max():.4f}", "largest |Ω_R(α=0.2) − Ω_R(α=0.1)| over 300–1000 K, cm⁻¹", "s7b_mixalpha.csv")
eta_R = x16[(x16.q_frac_x.abs() == 0.5) & (x16.q_frac_y.abs() == 0.5) & (x16.q_frac_z.abs() == 0.5) & (x16.branch == 1)].iloc[0]
tad12R = round(float(tad[(tad.q_label == "R") & (tad.q1 == 12)].omega_meV.iloc[0]) * MEV)
num("eta_R_300_at_tadano", round(eta_R.gamma_hwhm_cm1 / tad12R, 2),
    "η_R(300 K, 16³) with the same Γ at Tadano's converged Ω_R (direct factor only)", "mode_table_T.csv.gz; tadano_2015")
num("eta_R_300", round(eta_R.eta, 3), "η_R(300 K, 16³)", "mode_table_T.csv.gz")

# ------------------------------------------------------------------ S11 literature table
lit = (PJ / "prior_work_comparison.md").read_text()  # the comparison rows of SI Sec. S10, from the full texts
lines = [l for l in lit.splitlines() if l.startswith("|")]
cells = [[c.strip() for c in l.strip("|").split("|")] for l in lines]
hdr, body = cells[0], [c for c in cells[2:]]
CLEAN = [(r"\*\*|\*", ""), (r"`[^`]*`", ""), (r"\((?:D5|D7|N\d|A\d)[^)]*\)", ""), (r";\s*\)", ")"), (r"\s{2,}", " "),
         (r"\(\s*\)", ""), (r"\s*\(§\d+\)", ""), (r"\s+([,.;])", r"\1")]


def clean(t):
    for a, b in CLEAN:
        t = re.sub(a, b, t)
    return t.strip(" ;")


keep = ["Citation", "Material / phase", "FC source, functional", "Renormalization", "Transport", "3ph / 4ph", "Meshes",
        "T range", "κ(300 K)", "Mode-resolved Γ / η", "Coherence by branch pair", "Convergence of coherent term"]
rows = [[clean(c) for c in r] for r in body if r[0].strip("* ") in keep]
lt = pd.DataFrame(rows, columns=["row"] + [clean(h) for h in hdr[1:]])
# draft 2: Zhao et al. 2021 column from the full text (zhao_values.csv), and the
# "This work" column restated in the revised wording; "NOT-STATED" -> "not reported"
ZHAO = {"Citation": "PRB 104, 224304 (2021), doi:10.1103/PhysRevB.104.224304",
        "Material / phase": "16 cubic oxide and fluoride perovskites, including SrTiO₃, cubic",
        "FC source, functional": "VASP PAW, GGA-PBE, PBE-optimized a = 3.94 Å (expt 3.91 Å); 2×2×2 supercell; cubic and "
                                 "quartic IFCs by compressive sensing from 300 K AIMD snapshots (p. 3)",
        "Renormalization": "SCP (ALAMODE), q₁ = 8³, off-diagonal self-energy included, mixing 0.1; no cubic frequency "
                           "shift, no quasiharmonic expansion (pp. 3, 12)",
        "Transport": "particle BTE only (ShengBTE/FourPhonon); 3ph iterative, 4ph in SMRTA; no coherence term (p. 3)",
        "3ph / 4ph": "3ph and 3,4ph κ for every material",
        "Meshes": "12³ q-mesh (BTE); Gaussian smearing, scalebroad 0.05 (p. 3)",
        "T range": "200–1000 K for SrTiO₃ (p. 4)",
        "κ(300 K)": "graph only: κ_SCP^3ph ≈ 9.4 ± 0.1, κ_SCP^3,4ph ≈ 7.0 ± 0.1 W m⁻¹ K⁻¹ (read from figure, Fig. 2(c); "
                    "± is the reading uncertainty); ratio 0.74–0.75, i.e. four-phonon scattering lowers κ by about 25% "
                    "(24–27%, read from figure). Experimental lattice: 11.6 / 8.6 (read from figure, Fig. 11(a))",
        "Mode-resolved Γ / η": "not reported for SrTiO₃",
        "Coherence by branch pair": "none",
        "Convergence of coherent term": "—"}
THIS = {"FC source, functional": "force constants distributed with the ALAMODE SrTiO₃ example (to 4th order); checked "
                                 "harmonic frequencies agree with the published rounded values; identity with the 2015 "
                                 "production set not documented",
        "Meshes": "8³, 12³, 16³ at every T; 20³ at 300, 400, 500, 650, 800, 1000 K; 15³, 17³, 24³ at 300 K",
        "κ(300 K)": "κ_P = 8.89 (16³ reference), 8.89–9.21 across tested meshes; κ_C = 0.442–0.541 (15³–24³, observed "
                    "range); paired totals 9.399–9.651 W m⁻¹ K⁻¹",
        "Mode-resolved Γ / η": "η₃ = Γ/Ω for every mode on 12³, 16³, 20³, 300–1000 K; largest sampled η₃ = 0.184 (R-point "
                               "AFD mode, 300 K)",
        "Coherence by branch pair": "pair-resolved; six contributions by reference character under two rules",
        "Convergence of coherent term": "6 meshes at 300 K; 16³→20³ change ≤ 2.3% at 650, 800, 1000 K; 7–14% at "
                                        "300–500 K; mesh diagnostic r (empirical threshold)"}
lt.insert(3, "Zhao et al. 2021", lt.row.map(ZHAO).fillna("—"))
for k, v in THIS.items():
    lt.loc[lt.row == k, "This work"] = v
lt = lt.replace(r"NOT-STATED", "not reported", regex=True)
write("s11_literature.csv", lt, "comparison with prior SrTiO₃ transport studies, from the full texts (prior_work_comparison.md)")

num("tadano_R12", tad12R, "Tadano Table I, R25 at q₁ = 12³, 300 K, cm⁻¹ (rounded)", "tadano_2015/scph_convergence_300K.csv")
# ------------------------------------------------------------------ further prose numbers
dq = rd(PJ / "d_q1_freq_shift.csv").set_index(["T_K", "point"])
num("q1_8_12_R_300", round(dq.loc[(300, "R")].shift_pct, 1), "Ω_R(300 K), q₁ 8³ → 12³, %", "d_q1_freq_shift.csv")
num("R250_q8", round(dq.loc[(250, "R")].omega_q1_8_cm1, 1), "Ω_R(250 K), q₁ = 8³, cm⁻¹", "d_q1_freq_shift.csv")
n3a = rd(PJ / "n3_q1_16_check.csv").set_index(["T_K", "quantity"])
num("R250_q12", round(n3a.loc[(250, "omega_R")].q1_12, 1), "Ω_R(250 K), q₁ = 12³, cm⁻¹", "n3_q1_16_check.csv")
num("R250_q16", round(n3a.loc[(250, "omega_R")].q1_16, 1), "Ω_R(250 K), q₁ = 16³, cm⁻¹", "n3_q1_16_check.csv")
h = rd(PJ / "a2_harmonic_xcheck.csv")
h = h[(h.alamode_fc2 == "STO222") & (h.branch == 1)].set_index("point").omega_alamode_cm1
th = tad[tad.method == "SCPH"].iloc[:0]
tharm = rd(ROOT / "data/literature/STO/tadano_2015/scph_convergence_300K.csv")
tharm = tharm[tharm.method == "harmonic"].set_index("q_label").omega_meV
for lab, pt in (("G", "Gamma"), ("R", "R"), ("M", "M")):
    num(f"harm_{lab}", round(abs(float(h[pt])), 1), f"harmonic frequency at {pt}, lowest branch, our force constants (cm⁻¹; "
        "imaginary shown as magnitude)", "a2_harmonic_xcheck.csv")
    num(f"harm_{lab}_tad", round(abs(float(tharm[pt])) * MEV), f"harmonic frequency at {pt}, Tadano Table I (meV × 8.0655, "
        "rounded)", "tadano_2015/scph_convergence_300K.csv")
sp = mt[(mt.mesh == 16) & (mt.branch == 7) & np.isclose(mt.q_frac_x.abs(), 3 / 16) & np.isclose(mt.q_frac_y.abs(), 0.5)
        & np.isclose(mt.q_frac_z.abs(), 0.5)]
for T in (900, 950, 1000):
    num(f"spike_G{T}", round(float(sp[sp.T_K == T].gamma_hwhm_cm1.iloc[0]), 1),
        f"Γ of q = (3/16, 1/2, 1/2) branch 7 on 16³ at {T} K, cm⁻¹", "mode_table_T.csv.gz")
# strict-set diagnostic for R at 300 K: the inequality diagnostic evaluated with the
# strict-index Γ₄ on 6³ instead of the enlarged-set Γ₄,max; contains the residual finite-mesh resonances of S7.7
_v1 = rd(PJ / "v1_validation.csv")
_dec = rd(PJ / "a1_decision.csv").set_index("mode").loc["R"]
_pr = a1[(a1["mode"] == "R") & (a1.mesh == 6) & (a1.eps == 5)].iloc[0]
_om, _g3 = float(_pr.omega_cm1), float(_pr.Gamma3_prod16)
_st = {e: float(_v1[(_v1["mode"] == "R") & (_v1.mesh == 6) & (_v1.eps_cm1 == e)].Gamma4_strict.iloc[0]) for e in (2, 5)}
_allow = 0.5 * _om - _g3
num("strict_diag_R", f"{(_g3 + float(_dec.Gamma4_max)) / _om:.3f}; {(_g3 + _st[2]) / _om:.3f}; {(_g3 + _st[5]) / _om:.3f}; "
    f"{100 * _st[2] / _allow:.0f}; {_allow / 2:.2f}; {_st[2]:.3f}; {_st[5]:.3f}",
    "R, 300 K: (Γ₃ + Γ₄)/Ω with the enlarged-set Γ₄,max; with the strict-set Γ₄ on 6³ at ε = 2 and 5 cm⁻¹; strict Γ₄ "
    "(ε = 2) as % of the allowance 0.5Ω − Γ₃; half the allowance (cm⁻¹); strict Γ₄ at ε = 2 and 5 (cm⁻¹)",
    "v1_validation.csv; a1_decision.csv; a1_per_mesh.csv")
_e20 = eta[(eta.mesh == 20) & (eta.T_K == 800)].iloc[0]
num("eta20_800_omega", f"{float(_e20.eta_max_omega_cm1):.1f}", "frequency (cm⁻¹) of the largest-η mode on 20³ at 800 K "
    f"(q index {int(_e20.eta_max_q_index)}, branch {int(_e20.eta_max_branch)}, {_e20.eta_max_sector_char})", "eta_validity_T.csv")
pd.DataFrame(NUM).to_csv(OUT / "si_numbers_v4.csv", index=False)
stamp(OUT / "si_numbers_v4.csv", __file__, INPUTS)
print(f"wrote {len(list(OUT.glob('*.csv')))} files to {OUT.relative_to(ROOT)}; {len(NUM)} prose numbers")
