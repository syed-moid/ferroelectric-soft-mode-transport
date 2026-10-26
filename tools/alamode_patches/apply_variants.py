#!/usr/bin/env python3
"""Write the four-phonon validation variants of ALAMODE c5b015c anphon/selfenergy.cpp (diagram (c) only).

Each variant lives in its own clean copy of the c5b015c tree under the toolchain directory and is built
with its own build directory; the production tree and the earlier audit trees are never modified.

    alamode_4ph_trace     unmodified sum; one added print of a single fixed quartet (sign/convention trace)
    alamode_4ph_v1omit    enlarged (frequency) exclusion, as in the audit build, plus an accumulator of the
                          omitted terms
    alamode_4ph_v2strict  strict-index exclusion (leg at k = q carries the external branch index; the other
                          two legs are the same branch at k and -k), plus accumulators of the omitted terms
                          and of the enlarged-minus-strict terms (kept in this build)
    alamode_4ph_v3dump    enlarged exclusion (as the audit build) plus a per-rank dump of every excluded
                          quartet and every kept near-miss, frequencies at full double precision

Usage: apply_variants.py <toolchain dir>   (each tree must be a clean c5b015c checkout)
All variants print one line per external mode starting with "AUDIT 4ph(c) identity set:", which the queue
validator requires of non-production builds; the fields before the first ';' keep the audit build's meaning
(enlarged set, 1e-9 Ry).
"""
import sys
from pathlib import Path

HEAD = """void Selfenergy::selfenergy_c(const unsigned int N,
                              const double *T,
                              const double omega,
                              const unsigned int knum,
                              const unsigned int snum,
                              const KpointMeshUniform *kmesh_in,
                              const double *const *eval_in,
                              const std::complex<double> *const *const *evec_in,
                              std::complex<double> *ret) const
{"""

# ---------------------------------------------------------------------------------------------- trace build
TRACE_OLD = """                            ret_mpi[i] += v4_tmp
                                          * ((n12 + n23 + n31 + n1 + n2 + n3 + 1.0) * omega_sum[0]
                                             + (n31 + n23 + n3 - n12) * omega_sum[1]
                                             + (n12 + n31 + n1 - n23) * omega_sum[2]
                                             + (n23 + n12 + n2 - n31) * omega_sum[3]);
"""
TRACE_NEW = TRACE_OLD + """
                            // TRACE BUILD (not upstream ALAMODE): print one fixed quartet; the sum is unchanged.
                            if (ik1 == 5 && ik2 == 22 && is1 == 4 && is2 == 9 && is3 == 13) {
                                const auto v4c = anharmonic_core->V4(arr_quartic);
                                std::cout << std::setprecision(17) << std::scientific
                                          << "  AUDIT 4ph(c) identity set: TRACE no exclusion; T = " << T_tmp
                                          << " knum = " << knum << " snum = " << snum
                                          << " ik1 = " << ik1 << " ik2 = " << ik2 << " ik3 = " << ik3
                                          << " is1 = " << is1 << " is2 = " << is2 << " is3 = " << is3
                                          << " omega = " << omega << " omega1 = " << omega1
                                          << " omega2 = " << omega2 << " omega3 = " << omega3
                                          << " epsilon = " << epsilon
                                          << " V4re = " << v4c.real() << " V4im = " << v4c.imag()
                                          << " v4_tmp = " << v4_tmp
                                          << " n1 = " << n1 << " n2 = " << n2 << " n3 = " << n3
                                          << " S0 = " << omega_sum[0].real() << " " << omega_sum[0].imag()
                                          << " S1 = " << omega_sum[1].real() << " " << omega_sum[1].imag()
                                          << " S2 = " << omega_sum[2].real() << " " << omega_sum[2].imag()
                                          << " S3 = " << omega_sum[3].real() << " " << omega_sum[3].imag()
                                          << " factor = " << -1.0 / (std::pow(static_cast<double>(nk), 2) * std::pow(2.0, 5) * 3.0)
                                          << " Ry_to_kayser = " << Hz_to_kayser / time_ry
                                          << " kB_over_Ry = " << k_Boltzmann / Ryd
                                          << std::endl;
                            }
"""

# ---------------------------------------------------------------------------------------------- validation builds
VAR = HEAD + r"""
    /*
    Diagram (c), VALIDATION BUILD @NAME@ (not upstream ALAMODE).
    Same loop, vertex, thermal coefficients, frequency sums and prefactor as c5b015c. Predicates (1e-9 Ry):
      enlarged : one leg i has k_i = q and omega_i = omega (external mode or an exactly degenerate partner)
                 and the other two legs (at k and -k) have equal frequencies  [the audit build's set]
      strict   : one leg i has k_i = q and the external branch index, and the other two legs are the same
                 branch index at k and -k
    Excluded here: @EXCL@. Accumulators (all with the same prefactor, reduced over ranks):
      ret   : kept terms (the printed Gamma4 = 2 Im ret, FWHM, as upstream)
      omit  : the excluded terms
      extra : enlarged-but-not-strict terms (kept in the strict build; zero otherwise)
    */

    unsigned int i;
    unsigned int arr_quartic[4];

    const auto nk = kmesh_in->nk;
    const auto xk = kmesh_in->xk;
    double xk_tmp[3];

    std::complex<double> omega_sum[4];
    std::complex<double> *ret_mpi, *omit_mpi, *extra_mpi, *term;

    allocate(ret_mpi, N);
    allocate(omit_mpi, N);
    allocate(extra_mpi, N);
    allocate(term, N);

    std::complex<double> omega_shift = omega + im * epsilon;

    for (i = 0; i < N; ++i) {
        ret_mpi[i] = std::complex<double>(0.0, 0.0);
        omit_mpi[i] = std::complex<double>(0.0, 0.0);
        extra_mpi[i] = std::complex<double>(0.0, 0.0);
    }

    arr_quartic[0] = ns * kmesh_in->kindex_minus_xk[knum] + snum;

    const double TOL_DEG = 1.0e-9; // Ry
    const bool strict_build = @STRICT@;
    // cnt: 0 enlarged; 1 enlarged with the external branch on the leg; 2 exact-zero resonant denominators;
    //      3 near-misses kept (TOL_DEG <= min denominator < 1e-6 Ry, not enlarged); 4 strict; 5 strict, not enlarged
    double cnt[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
@DUMPOPEN@
    for (unsigned int ik1 = mympi->my_rank; ik1 < nk; ik1 += mympi->nprocs) {
        for (unsigned int ik2 = 0; ik2 < nk; ++ik2) {

            xk_tmp[0] = xk[knum][0] - xk[ik1][0] - xk[ik2][0];
            xk_tmp[1] = xk[knum][1] - xk[ik1][1] - xk[ik2][1];
            xk_tmp[2] = xk[knum][2] - xk[ik1][2] - xk[ik2][2];

            const auto ik3 = kmesh_in->get_knum(xk_tmp);

            for (unsigned int is1 = 0; is1 < ns; ++is1) {

                arr_quartic[1] = ns * ik1 + is1;
                double omega1 = eval_in[ik1][is1];

                for (unsigned int is2 = 0; is2 < ns; ++is2) {

                    arr_quartic[2] = ns * ik2 + is2;
                    double omega2 = eval_in[ik2][is2];

                    for (unsigned int is3 = 0; is3 < ns; ++is3) {

                        arr_quartic[3] = ns * ik3 + is3;
                        double omega3 = eval_in[ik3][is3];

                        const double d1 = std::abs(omega - omega1 - omega2 + omega3);
                        const double d2 = std::abs(omega + omega1 - omega2 - omega3);
                        const double d3 = std::abs(omega - omega1 + omega2 - omega3);
                        const double dmin = std::min({d1, d2, d3});
                        if (dmin < TOL_DEG) cnt[2] += 1.0;
                        const bool e1 = ik1 == knum && std::abs(omega1 - omega) < TOL_DEG
                                        && std::abs(omega2 - omega3) < TOL_DEG;
                        const bool e2 = ik2 == knum && std::abs(omega2 - omega) < TOL_DEG
                                        && std::abs(omega1 - omega3) < TOL_DEG;
                        const bool e3 = ik3 == knum && std::abs(omega3 - omega) < TOL_DEG
                                        && std::abs(omega1 - omega2) < TOL_DEG;
                        const bool enlarged = e1 || e2 || e3;
                        const bool s1 = ik1 == knum && is1 == snum && is2 == is3;
                        const bool s2 = ik2 == knum && is2 == snum && is1 == is3;
                        const bool s3 = ik3 == knum && is3 == snum && is1 == is2;
                        const bool strict = s1 || s2 || s3;
                        if (enlarged) cnt[0] += 1.0;
                        if ((e1 && is1 == snum) || (e2 && is2 == snum) || (e3 && is3 == snum)) cnt[1] += 1.0;
                        if (strict) cnt[4] += 1.0;
                        if (strict && !enlarged) cnt[5] += 1.0;
                        if (!enlarged && dmin >= TOL_DEG && dmin < 1.0e-6) cnt[3] += 1.0;
@DUMP@
                        const bool excluded = strict_build ? strict : enlarged;

                        double v4_tmp = std::norm(anharmonic_core->V4(arr_quartic));

                        omega_sum[0]
                                = 1.0 / (omega_shift - omega1 - omega2 - omega3)
                                  - 1.0 / (omega_shift + omega1 + omega2 + omega3);
                        omega_sum[1]
                                = 1.0 / (omega_shift - omega1 - omega2 + omega3)
                                  - 1.0 / (omega_shift + omega1 + omega2 - omega3);
                        omega_sum[2]
                                = 1.0 / (omega_shift + omega1 - omega2 - omega3)
                                  - 1.0 / (omega_shift - omega1 + omega2 + omega3);
                        omega_sum[3]
                                = 1.0 / (omega_shift - omega1 + omega2 - omega3)
                                  - 1.0 / (omega_shift + omega1 - omega2 + omega3);

                        for (i = 0; i < N; ++i) {
                            double T_tmp = T[i];

                            double n1 = thermodynamics->fB(omega1, T_tmp);
                            double n2 = thermodynamics->fB(omega2, T_tmp);
                            double n3 = thermodynamics->fB(omega3, T_tmp);

                            double n12 = n1 * n2;
                            double n23 = n2 * n3;
                            double n31 = n3 * n1;

                            term[i] = v4_tmp
                                      * ((n12 + n23 + n31 + n1 + n2 + n3 + 1.0) * omega_sum[0]
                                         + (n31 + n23 + n3 - n12) * omega_sum[1]
                                         + (n12 + n31 + n1 - n23) * omega_sum[2]
                                         + (n23 + n12 + n2 - n31) * omega_sum[3]);
                        }
                        if (excluded) {
                            for (i = 0; i < N; ++i) omit_mpi[i] += term[i];
                            continue;
                        }
                        if (enlarged) { // strict build only: enlarged-but-not-strict, kept
                            for (i = 0; i < N; ++i) extra_mpi[i] += term[i];
                        }
                        for (i = 0; i < N; ++i) ret_mpi[i] += term[i];
                    }
                }
            }
        }
    }
@DUMPCLOSE@
    double factor = -1.0 / (std::pow(static_cast<double>(nk), 2) * std::pow(2.0, 5) * 3.0);
    for (i = 0; i < N; ++i) {
        ret_mpi[i] *= factor;
        omit_mpi[i] *= factor;
        extra_mpi[i] *= factor;
    }

    mpi_reduce_complex(N, ret_mpi, ret);

    std::complex<double> *omit, *extra;
    allocate(omit, N);
    allocate(extra, N);
    mpi_reduce_complex(N, omit_mpi, omit);
    mpi_reduce_complex(N, extra_mpi, extra);

    {
        std::complex<double> cnt_c[6], cnt_sum[6];
        for (i = 0; i < 6; ++i) cnt_c[i] = std::complex<double>(cnt[i], 0.0);
        mpi_reduce_complex(6, cnt_c, cnt_sum);
        if (mympi->my_rank == 0) {
            const double to_cm = Hz_to_kayser / time_ry; // Ry -> cm^-1, as Writes::in_kayser
            std::cout << "  AUDIT 4ph(c) identity set: nk = " << nk
                      << ", excluded = " << static_cast<long long>(cnt_sum[0].real() + 0.5)
                      << " (external mode itself on the leg: " << static_cast<long long>(cnt_sum[1].real() + 0.5)
                      << "), denominator check = " << static_cast<long long>(cnt_sum[2].real() + 0.5)
                      << ", near-misses kept (< 1e-6 Ry) = " << static_cast<long long>(cnt_sum[3].real() + 0.5)
                      << "; VARIANT @NAME@: strict = " << static_cast<long long>(cnt_sum[4].real() + 0.5)
                      << ", strict-not-enlarged = " << static_cast<long long>(cnt_sum[5].real() + 0.5)
                      << std::setprecision(12) << std::scientific;
            for (i = 0; i < N; ++i) {
                std::cout << "; T = " << T[i]
                          << " ImKept_cm1 = " << ret[i].imag() * to_cm
                          << " ImOmitted_cm1 = " << omit[i].imag() * to_cm
                          << " ImEnlargedNotStrict_cm1 = " << extra[i].imag() * to_cm
                          << " ReKept_cm1 = " << ret[i].real() * to_cm
                          << " ReOmitted_cm1 = " << omit[i].real() * to_cm;
            }
            std::cout << std::defaultfloat << std::endl;
        }
    }

    deallocate(ret_mpi);
    deallocate(omit_mpi);
    deallocate(extra_mpi);
    deallocate(term);
    deallocate(omit);
    deallocate(extra);
}

"""

DUMPOPEN = r"""    std::ofstream ofs_dump;
    {
        const std::string fname = "AUDITDUMP_k" + std::to_string(knum) + "_s" + std::to_string(snum)
                                  + "_r" + std::to_string(mympi->my_rank) + ".txt";
        ofs_dump.open(fname.c_str(), std::ios::out);
        ofs_dump << std::setprecision(17) << std::scientific;
        if (mympi->my_rank == 0) {
            ofs_dump << "# kpoints ik x y z (fractional)\n";
            for (unsigned int ik = 0; ik < nk; ++ik)
                ofs_dump << "K " << ik << " " << xk[ik][0] << " " << xk[ik][1] << " " << xk[ik][2] << "\n";
        }
        ofs_dump << "# kind ik1 ik2 ik3 is1 is2 is3 omega omega1 omega2 omega3 d1 d2 d3 e1 e2 e3 s1 s2 s3 (Ry)\n";
    }
"""
DUMP = r"""                        if (enlarged || strict || dmin < 1.0e-6) {
                            ofs_dump << (enlarged ? "X " : "N ") << ik1 << " " << ik2 << " " << ik3 << " "
                                     << is1 << " " << is2 << " " << is3 << " " << omega << " " << omega1 << " "
                                     << omega2 << " " << omega3 << " " << d1 << " " << d2 << " " << d3 << " "
                                     << e1 << " " << e2 << " " << e3 << " " << s1 << " " << s2 << " " << s3 << "\n";
                        }
"""
DUMPCLOSE = "    ofs_dump.close();\n"

VARIANTS = {
    "alamode_4ph_v1omit": dict(NAME="omitted-sum", STRICT="false", EXCL="the enlarged set", dump=False),
    "alamode_4ph_v2strict": dict(NAME="strict-index", STRICT="true", EXCL="the strict set only", dump=False),
    "alamode_4ph_v3dump": dict(NAME="full-precision-dump", STRICT="false", EXCL="the enlarged set", dump=True),
}


def function_span(src):
    a = src.index(HEAD)
    b = src.index("void Selfenergy::selfenergy_d(")
    return a, b


def main(tool):
    tool = Path(tool)
    f = tool / "alamode_4ph_trace/anphon/selfenergy.cpp"
    s = f.read_text()
    assert s.count(TRACE_OLD) == 1 and "TRACE BUILD" not in s
    s = s.replace(TRACE_OLD, TRACE_NEW)
    f.write_text(s.replace('#include "mpi_common.h"', '#include <iomanip>\n#include "mpi_common.h"', 1))
    for tree, v in VARIANTS.items():
        f = tool / tree / "anphon/selfenergy.cpp"
        s = f.read_text()
        assert "VALIDATION BUILD" not in s
        a, b = function_span(s)
        body = VAR.replace("@NAME@", v["NAME"]).replace("@STRICT@", v["STRICT"]).replace("@EXCL@", v["EXCL"])
        body = body.replace("@DUMPOPEN@", DUMPOPEN if v["dump"] else "")
        body = body.replace("@DUMP@", DUMP if v["dump"] else "")
        body = body.replace("@DUMPCLOSE@", DUMPCLOSE if v["dump"] else "")
        s = s[:a] + body + s[b:]
        if "#include <algorithm>" not in s:
            s = s.replace('#include "mpi_common.h"', '#include <algorithm>\n#include <fstream>\n#include <iomanip>\n#include <string>\n#include "mpi_common.h"', 1)
        f.write_text(s)
        print("patched", tree)


if __name__ == "__main__":
    main(sys.argv[1])
