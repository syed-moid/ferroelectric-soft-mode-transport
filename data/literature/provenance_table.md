# Literature provenance table — Turn 4A

Updated as CSVs are filled. `Status` values: `extracted` (filled and verified),
`pre-filled` (filled using known values, no PDF digitization),
`pending` (template only, awaiting user digitization), `n/a` (paper not
obtained or quantity absent from source).

| Quantity | Material | Source | Type | Used as | CSV path | Status |
| --- | --- | --- | --- | --- | --- | --- |
| Ω_soft(T) (2 anchors) | BTO | Ehsan 2021 | First-principles SCP | calibration anchor | `BTO/ehsan_2021/anchors.csv` | pre-filled |
| Ω_soft(T) (T-curve) | BTO | Ehsan 2021 | First-principles SCP | cross-validation | `BTO/ehsan_2021/scp_omega_T.csv` | n/a (paper has no Ω(T) curve; see extraction log) |
| F_cubic(T), F_tetra(T) | BTO | Ehsan 2021 Fig. 5 | First-principles | independent F(T) crossing for T_c | `BTO/ehsan_2021/free_energy_T.csv` | extracted |
| structural a, c, V_cell (cubic + tetra) | BTO | Ehsan 2021 Table I | First-principles PBEsol | cross-check `a = 4.00 Å` lock | `BTO/ehsan_2021/structural_properties.csv` | extracted |
| cubic phonon dispersion (0 K + 500 K) | BTO | Ehsan 2021 Fig. 2 | First-principles SCP | independent Λ² estimate at q=0 calibration condition | `BTO/ehsan_2021/cubic_dispersion.csv` | extracted |
| Born Z\*, ε_∞ (cubic + tetra) | BTO | Ehsan 2021 Table II | First-principles | LO-TO splitting, polar character | `BTO/ehsan_2021/born_charges_epsilon_infty.csv` | extracted |
| ℱ(T), four phases | BTO | Tadano 2022 Fig. 4 | First-principles | T_c via free-energy crossings | `BTO/tadano_2022/free_energy_T.csv` | extracted |
| Ω_soft(T) at +1% lattice | BTO | Tadano 2022 Fig. 15 | First-principles SCP | cross-validation only (lattice-flag) | `BTO/tadano_2022/scp_dispersion_plus1pct.csv` | extracted |
| T_c (three transitions) | BTO | Tadano 2022 Table I | First-principles SCP | cross-validation | `BTO/tadano_2022/transition_temperatures.csv` | extracted |
| total κ(T) + experimental | BTO | Tadano 2022 Fig. 7 | First-principles + experimental | scale-overlay only (not validation) | `BTO/tadano_2022/kappa_total_T.csv` | extracted (28 rows) |
| Λ² | BTO | Shirane 1970 (prior fit) | Neutron INS | input | `BTO/shirane_1970/Lambda_sq.csv` | pre-filled |
| ω(q) along Σ (e‖[110]) and Δ (e‖[100]) at 295 K | BTO | Shirane 1970 Figs. 3 & 6 | Neutron INS | (re)fit Λ² with two-direction comparison | `BTO/shirane_1970/cubic_dispersion.csv` | extracted (21 points; tetragonal-phase caveat) |
| TA + soft-TO dispersion along [100] (+ low-lying dispersion) | BTO | Harada-Axe-Shirane 1971 Figs. 5 + 8 | Neutron INS | **cubic-phase** Λ² (removes tetragonal caveat) | `BTO/harada_axe_shirane_1971/cubic_dispersion.csv` | verified (Wave 2; cubic 423/503 K; soft-TO ζ≤0.3) |
| Λ²_cubic (Fig. 5 (Ω_TO)²=(Ω_TO0)²+Aq² fit) | BTO | Harada-Axe-Shirane 1971 Fig. 5 / p.160 | Neutron INS | **soft [100] Λ²=972±20% CENTRAL** (kernel); anisotropic (stiff [110]/[111]=4750, ~5×≈Shirane 5100) | `BTO/harada_axe_shirane_1971/Lambda_sq.csv` | verified; NB 06 Λ² 10%-gate **resolved by sign-off**: keep paper-text 972±20%; NB-06 3-pt digitized refit (818, q=0/0.313/0.478 Å⁻¹@423K) corroborates within ±20% envelope [778,1166] (movement +15.8% but consistent, not a discrepancy — cf. 1988/1995 freq cross-check) |
| ω_TO, ω_TA, Γ_TO, Γ_TA, coupling vs q | BTO | Harada-Axe-Shirane 1971 Table I | Neutron INS | cubic soft-TO damping Γ; **overdamped η≈2.3** (BTO methodological note) | `BTO/harada_axe_shirane_1971/mode_parameters.csv` | verified (Wave 2; phase→cubic + T °C→K corrected) |
| soft-TO anisotropy near q=0 (Ω²(angle)) | BTO | Harada-Axe-Shirane 1971 Fig. 7 | Neutron INS | directional Λ²_∥/Λ²_⊥ anisotropy ([111]–[100]–[110] arc) | `BTO/harada_axe_shirane_1971/anisotropy.csv` | verified (Wave 2; signed angles −54.74° to +45°) |
| Ω(T) for TO1 + TO2 (SCPH±PM, TDEP, exp) | STO | Tadano 2015 Fig. 5 | First-principles SCPH + INS exp | **primary STO benchmark** | `STO/tadano_2015/scph_omega_T.csv` | extracted (30 rows) |
| ω(q) at 300 K (SCPH + harm + exp) | STO | Tadano 2015 Fig. 4 | First-principles | constrain Λ²_STO | `STO/tadano_2015/scph_dispersion.csv` | extracted (14 rows) |
| ω(300 K, q₁-grid) convergence | STO | Tadano 2015 Table I | First-principles SCPH | convergence cross-check at 300 K | `STO/tadano_2015/scph_convergence_300K.csv` | extracted (21 rows) |
| total κ(T) + Muta/Popuri exp | STO | Tadano 2015 Fig. 7 | First-principles + experimental | scale-overlay only (not validation) | `STO/tadano_2015/kappa_total_T.csv` | pending (optional) |
| Ω_TO1(T) (MLFF+SSCHA shifted, INS, √-fit) | STO | Verdi 2023 Fig. 4(a) | SSCHA + INS exp | independent STO benchmark (shift_meV flag attached) | `STO/verdi_2023/sscha_omega_T.csv` | extracted (28 rows) |
| ω(q) at 0/300 K (harm + DFT_SSCHA + MLFF_SSCHA + exp) | STO | Verdi 2023 Figs. 3 + 5 | First-principles + SSCHA + exp | independent Λ²_STO + 300 K cross-check | `STO/verdi_2023/dispersion.csv` | extracted (22 rows) |
| 1/ε(T) (SSCHA + Curie-Weiss) | STO | Verdi 2023 Fig. 4(b) | SSCHA + classical fit | LST consistency check on Ω_TO1(T) | `STO/verdi_2023/inverse_dielectric_T.csv` | extracted (26 rows) |
| Ω_soft(q=0, T) + ω(q) dispersion | STO | Yamada-Shirane 1969 Table II (+ Figs. 2/3) | Neutron INS | **gold-standard exp STO Ω(T)** 4.5–293 K | `STO/yamada_shirane_1969/omega_T.csv` | extracted (Table II pre-filled: 19 q=0 points + 4.5 K dispersion) |
| Γ_soft(T) at q=0 | STO | Yamada-Shirane 1969 Fig. 6 | Neutron INS | **experimental STO Γ(T)** (enables eta_fit, NB 08); underdamped η≈0.09 | `STO/yamada_shirane_1969/gamma_T.csv` | verified (Wave 2; 0.24→0.84 meV 20–250 K, matches Fig. 6) |
| DHO ω(T) at Γ, R, M, 4-1400 K + DHO Γ(T) if available | STO | Fauqué 2022 Fig. 1(d), 2(c), 3(d) | INS DHO fits + TDEP | **central STO Γ(T) source if DHO Γ reported**; Ω(T) extension to 1400 K | `STO/fauque_2022/dho_fits_T.csv` | extracted (28 rows; Γ/R/M; q_point_rlu re-quoted) |
| TDEP dispersion at 300 K (M-Γ-R) | STO | Fauqué 2022 Fig. 1(c) | TDEP | independent Λ²_STO at 300 K | `STO/fauque_2022/tdep_dispersion_300K.csv` | extracted (14 rows) |
| Ω_Eu(T) for ST¹⁶O + ST¹⁸O + experimental | STO | Schmidt-Spaldin 2025 Fig. 3 top-left | SSCHA+MLIP + INS exp [Ref 50] | third independent SSCHA Ω(T) for STO; isotope cross-check | `STO/schmidt_2025/sscha_omega_T.csv` | extracted (24 rows) |
| SSCHA energy vs polarization (0 K, 16O + 18O) | STO | Schmidt-Spaldin 2025 Fig. 2 | SSCHA + RPA-MLIP | anharmonic-potential illustration §3.3 | `STO/schmidt_2025/energy_polarization.csv` | extracted (25 rows) |
| Ω²_TO(T): SSCHA + Shirane exp | KTO | Meier 2022 Fig. 3(a) | SSCHA+MLIP + Shirane INS exp | **primary KTaO₃ benchmark** (Section 4 third material) | `KTO/meier_2022/omega_T.csv` | extracted (18 rows) |
| ε(T): calc + Fujishita exp + Barrett fit | KTO | Meier 2022 Fig. 3(b) | SSCHA+MLIP + Fujishita exp | LST cross-check on KTaO₃ | `KTO/meier_2022/inverse_dielectric_T.csv` | extracted (18 rows) |
| T-dependent phonon spectrum | KTO | Meier 2022 Fig. 2 | SSCHA+MLIP | multi-T Λ²_KTO sampling | `KTO/meier_2022/dispersion_T.csv` | extracted (16 rows) |
| Ω_TO1(T) (MLFF+SSCHA + INS + √-fit) | KTO | Ranalli 2023 Fig. 4(a) | MLFF+SSCHA + INS exp | **2nd KTaO₃ source**; 2-source Ω(T) triangulation vs Meier (exp anchors agree) | `KTO/ranalli_2023/omega_T.csv` | extracted+verified (30 rows; native cm⁻¹→meV; MLFF shift_meV=−4.30) |
| 1/ε(T) + Curie-Weiss | KTO | Ranalli 2023 Fig. 4(b) | MLFF+SSCHA + classical fit | KTaO₃ LST cross-check (1/ε vs Ω²) + Barrett saturation | `KTO/ranalli_2023/inverse_dielectric_T.csv` | extracted+verified (22 rows; NB 05b Fig. F2) |
| Ω(T) exp q=0 | KTO | Vogt 1995 Fig. 5 | Hyper-Raman (exp) | KTaO₃ Ω(T) tie-breaker (Meier vs Ranalli); η_KTO numerator-pair | `KTO/vogt_1995/omega_T.csv` | dense WPD (39 pts), figure-faithful; 1995 top 294 K=10.15 (NOT 10.48); 5 K=2.44 plateau; +1988 Table III pt ±0.37 |
| γ(T) exp q=0 | KTO | Vogt 1995 Fig. 5 | Hyper-Raman (exp) | **PRIMARY η_KTO** → eta_fit (case 'a' resolved); plotted γ=FWHM | `KTO/vogt_1995/gamma_T.csv` | dense WPD (39 pts; raw, never scaled — no slope artifact), verified vs Fig.5 (boundary 7→9.7 cm⁻¹; endpoint 2.57 meV=20.7 cm⁻¹); η_KTO=(γ/2)/Ω≈0.127; 296/298 K double-pick dropped; unc 0.12/0.17 stat⊕sys |
| Ω(T) exp q=0 (phase-split) | STO | Vogt 1995 Fig. 6 | Hyper-Raman (exp) | 6th STO Ω(T) source; cubic_F1u for kernel (F1u top 299 K=10.96, η denom) | `STO/vogt_1995/omega_T.csv` | dense WPD (47 pts); A2u/Eu split, 3 sub-105K merged pts → transition_merged, cubic_F1u≥110K; unc 0.12/0.17 (60K) |
| γ(T) exp q=0 (phase-split) | STO | Vogt 1995 Fig. 6 | Hyper-Raman (exp) | **PRIMARY η_STO** → eta_fit (true q=0, central peak separated) | `STO/vogt_1995/gamma_T.csv` | dense WPD (36 pts); **η_STO=0.109 central, ±0.02 reading unc** (dense marker-center vs eye-read ~0.125; not corrected toward round); single shared-γ branch (Flag-2); YS-vs-Vogt ~2.4× (was 3.3×); unc 0.12/0.17 (60K) |
| Vogt γ/Ω uncertainty convention | STO+KTO | Vogt 1995 p.8051 | — | ±1 cm⁻¹≈0.12 meV (stat, T≥60 K); 0.17 meV=√(0.12²+0.12²) stat⊕sys for T<60 K (defect-induced floor; "below about 60 K"). 1988 Table III pt keeps ±0.37. | all 4 Vogt files | locked (recorded; applied consistently) |
| AOM/Barrett fit params | KTO | Vogt 1995 Fig. 5 caption | metadata | §3.3 Barrett/T₁ discussion only — **NOT data** | `KTO/vogt_1995/model_params.md` | metadata (no loader; no rows reconstructed) |
| Akimov 2000 STO film E-field Raman | STO | Akimov 2000 | thin-film + E-field Raman | **not used** (thin-film, bulk reference already in Tadano Fig. 5) | — | n/a (reviewed and skipped) |
| Maity 2025 AFD R-point soft mode + central peak | STO | Maity 2025 (PRL) | Neutron INS + Cp + κ | **not used** (AFD R₂₅ octahedral-rotation mode, not the FE Γ soft mode; same class skipped in Tadano 2015 Fig. 6). Fig. 1(e) κ(T) optional scale-overlay only | — | n/a (reviewed and skipped; optional κ overlay pending user decision) |
| Ω(T), Γ(T) | BTO | Old hybrid CSVs | Phenomenological reconstruction | internal consistency check only | `data/raw/BTO/phonon_modes_*K.csv` | already extracted (Turn 1) |
