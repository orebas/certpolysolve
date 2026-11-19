// The compensated certificate's loop (compensated_kernel.hpp), at one point and at four points at a time.
//
// The scheme (CompensatedCertificate's evaluation, compensated_certificate.cpp, states it with its bounds). With S_k
// the exact partial sums for the coefficients hi + lo (S_n = a_n, S_k = S_(k+1) x + a_k) and s_k the rounded ones,
// every step is split exactly:
//     s_(k+1) x + a_k = s_k + t_k,  t_k = sigma + h1 - h2 + h5 + lo_re + i (sigma' + h3 + h4 + h6 + lo_im)
// (Cameron and Graillat 2024, Theorems 2.4 and 2.5), up to 2^-1074 for each of the four TwoProd below the normal
// range. Hence S_k - s_k = sum_(j >= k) t_j x^(j - k), and p(x) = s_0 + E(x), E(x) = sum_j t_j x^j, which the
// correction c evaluates by Horner's scheme on the rounded t_j. The accumulators, nonnegative and rounded to nearest:
//   terms        = sum_j T_j X^j,             T_j >= |t_j| (the moduli of t_j's terms, by |Re| + |Im|)
//   terms_slope  = sum_j j T_j X^(j - 1)
//   values       = sum_(k < n) (|Re s_(k+1)| + |Im s_(k+1)|) X^k
//   errors       = sum_k error_k X^k
//   errors_slope = sum_k k error_k X^(k - 1)
//
// avx2_fma is `reference` written a second time, each scalar operation as the AVX2 instruction that makes it at four
// points: a sum, a difference, a product and a fused product round each lane as the scalar operation rounds
// (IEEE 754, to nearest), |a| clears the sign bit and -a flips it, as std::abs and the unary minus do. Neither
// function reorders, fuses or drops an operation (the library's -ffp-contract=off). tests/test_kernel_variants.cpp
// compares the two on every field.
#include "certify/attempts/compensated_kernel.hpp"
#include "numeric/arithmetic/machine_arithmetic.hpp"
#include <cmath>
#ifdef CERTPOLYSOLVE_HAVE_AVX2_FMA_KERNELS
#include <immintrin.h>
#endif

namespace certpolysolve::detail::compensated_kernel {
Sums reference(const Term *terms_, uint64_t n, double xr, double xi, double X) noexcept {
    const Term &top = terms_[n];
    double sr = top.re_hi, si = top.im_hi;
    double cr = top.re_lo, ci = top.im_lo; // t_n = the leading coefficient's low word, exactly
    double dr = 0, di = 0;                 // the derivative recurrence on the rounded s
    double terms = std::abs(top.re_lo) + std::abs(top.im_lo), terms_slope = 0, values = 0;
    double errors = top.error, errors_slope = 0;
    for (uint64_t k = n; k-- > 0;) {
        const Term &a = terms_[k];
        // d = d x + s_(k+1)
        {
            const double pr = dr * xr - di * xi, pi = dr * xi + di * xr;
            dr = pr + sr;
            di = pi + si;
        }
        values = values * X + (std::abs(sr) + std::abs(si));
        // s x as w + x' + y + z exactly (TwoProdCmplx), then w + hi = s' + sigma (TwoSumCmplx).
        const double g1 = sr * xr, h1 = std::fma(sr, xr, -g1);
        const double g2 = si * xi, h2 = std::fma(si, xi, -g2);
        const double g3 = sr * xi, h3 = std::fma(sr, xi, -g3);
        const double g4 = si * xr, h4 = std::fma(si, xr, -g4);
        const auto w_re = machine::DoubleDouble::sum(g1, -g2);
        const auto w_im = machine::DoubleDouble::sum(g3, g4);
        const auto s_re = machine::DoubleDouble::sum(w_re.hi, a.re_hi);
        const auto s_im = machine::DoubleDouble::sum(w_im.hi, a.im_hi);
        const double tr = (((h1 - h2) + w_re.lo) + s_re.lo) + a.re_lo;
        const double ti = (((h3 + h4) + w_im.lo) + s_im.lo) + a.im_lo;
        const double T = (((((((((std::abs(h1) + std::abs(h2)) + std::abs(w_re.lo)) + std::abs(s_re.lo)) +
                               std::abs(a.re_lo)) + std::abs(h3)) + std::abs(h4)) + std::abs(w_im.lo)) +
                           std::abs(s_im.lo)) + std::abs(a.im_lo));
        // c = c x + t
        {
            const double pr = cr * xr - ci * xi, pi = cr * xi + ci * xr;
            cr = pr + tr;
            ci = pi + ti;
        }
        terms_slope = terms_slope * X + terms;
        terms = terms * X + T;
        errors_slope = errors_slope * X + errors;
        errors = errors * X + a.error;
        sr = s_re.hi;
        si = s_im.hi;
    }
    return {sr, si, cr, ci, dr, di, terms, terms_slope, values, errors, errors_slope};
}

double majorant(const double *coefficients, uint64_t top, uint64_t bottom, double x) noexcept {
    double sum = coefficients[top];
    for (uint64_t k = top; k-- > bottom;)
        sum = sum * x + coefficients[k];
    return sum;
}

#ifdef CERTPOLYSOLVE_HAVE_AVX2_FMA_KERNELS
// `majorant` at sixteen points: each step's product and sum as the AVX2 instructions that make them at four points,
// in four registers; no fused product (the function is built without the FMA instruction).
__attribute__((target("avx2"))) void majorants_avx2(const double *coefficients, uint64_t top, uint64_t bottom,
                                                    const double *x, double *out) noexcept {
    constexpr size_t registers = 4, lanes = 4;
    static_assert(majorant_group == registers * lanes, "four AVX2 registers of four binary64 numbers");
    __m256d at[registers], sum[registers];
    for (size_t r = 0; r < registers; ++r) {
        at[r] = _mm256_loadu_pd(x + lanes * r);
        sum[r] = _mm256_broadcast_sd(coefficients + top);
    }
    for (uint64_t k = top; k-- > bottom;) {
        const __m256d c = _mm256_broadcast_sd(coefficients + k);
        for (size_t r = 0; r < registers; ++r)
            sum[r] = _mm256_add_pd(_mm256_mul_pd(sum[r], at[r]), c);
    }
    for (size_t r = 0; r < registers; ++r)
        _mm256_storeu_pd(out + lanes * r, sum[r]);
}

namespace {
#define CERTPOLYSOLVE_AVX2_FMA __attribute__((target("avx2,fma")))
#define CERTPOLYSOLVE_AVX2_FMA_INLINE __attribute__((always_inline, target("avx2,fma"))) inline
using Four = __m256d;
// |a|, as std::abs: the sign bit cleared.
CERTPOLYSOLVE_AVX2_FMA_INLINE Four absolute(Four a) { return _mm256_andnot_pd(_mm256_set1_pd(-0.0), a); }
// -a, as the unary minus: the sign bit flipped.
CERTPOLYSOLVE_AVX2_FMA_INLINE Four negated(Four a) { return _mm256_xor_pd(a, _mm256_set1_pd(-0.0)); }
// machine::DoubleDouble::sum(a, b): s = a + b, v = s - a, {s, (a - (s - v)) + (b - v)}.
CERTPOLYSOLVE_AVX2_FMA_INLINE void two_sum(Four a, Four b, Four &hi, Four &lo) {
    const Four s = _mm256_add_pd(a, b), v = _mm256_sub_pd(s, a);
    hi = s;
    lo = _mm256_add_pd(_mm256_sub_pd(a, _mm256_sub_pd(s, v)), _mm256_sub_pd(b, v));
}
// a b - c d and a b + c d, each product rounded, then the difference or the sum.
CERTPOLYSOLVE_AVX2_FMA_INLINE Four products_difference(Four a, Four b, Four c, Four d) {
    return _mm256_sub_pd(_mm256_mul_pd(a, b), _mm256_mul_pd(c, d));
}
CERTPOLYSOLVE_AVX2_FMA_INLINE Four products_sum(Four a, Four b, Four c, Four d) {
    return _mm256_add_pd(_mm256_mul_pd(a, b), _mm256_mul_pd(c, d));
}
} // namespace

CERTPOLYSOLVE_AVX2_FMA void avx2_fma(const Term *terms_, uint64_t n, const double *xr4, const double *xi4,
                                 const double *X4, Sums *out) noexcept {
    static_assert(group == 4, "one AVX2 register holds four binary64 numbers");
    const Four xr = _mm256_loadu_pd(xr4), xi = _mm256_loadu_pd(xi4), X = _mm256_loadu_pd(X4);
    const Term &top = terms_[n];
    Four sr = _mm256_set1_pd(top.re_hi), si = _mm256_set1_pd(top.im_hi);
    Four cr = _mm256_set1_pd(top.re_lo), ci = _mm256_set1_pd(top.im_lo);
    Four dr = _mm256_setzero_pd(), di = _mm256_setzero_pd();
    Four terms = _mm256_set1_pd(std::abs(top.re_lo) + std::abs(top.im_lo));
    Four terms_slope = _mm256_setzero_pd(), values = _mm256_setzero_pd();
    Four errors = _mm256_set1_pd(top.error), errors_slope = _mm256_setzero_pd();
    for (uint64_t k = n; k-- > 0;) {
        const Term &a = terms_[k];
        const Four a_re_hi = _mm256_broadcast_sd(&a.re_hi), a_im_hi = _mm256_broadcast_sd(&a.im_hi);
        const Four a_re_lo = _mm256_broadcast_sd(&a.re_lo), a_im_lo = _mm256_broadcast_sd(&a.im_lo);
        const Four a_error = _mm256_broadcast_sd(&a.error);
        // d = d x + s_(k+1)
        {
            const Four pr = products_difference(dr, xr, di, xi), pi = products_sum(dr, xi, di, xr);
            dr = _mm256_add_pd(pr, sr);
            di = _mm256_add_pd(pi, si);
        }
        values = _mm256_add_pd(_mm256_mul_pd(values, X), _mm256_add_pd(absolute(sr), absolute(si)));
        // TwoProdCmplx: std::fma(a, b, -g) is the fused a b - g.
        const Four g1 = _mm256_mul_pd(sr, xr), h1 = _mm256_fmsub_pd(sr, xr, g1);
        const Four g2 = _mm256_mul_pd(si, xi), h2 = _mm256_fmsub_pd(si, xi, g2);
        const Four g3 = _mm256_mul_pd(sr, xi), h3 = _mm256_fmsub_pd(sr, xi, g3);
        const Four g4 = _mm256_mul_pd(si, xr), h4 = _mm256_fmsub_pd(si, xr, g4);
        Four w_re_hi, w_re_lo, w_im_hi, w_im_lo, s_re_hi, s_re_lo, s_im_hi, s_im_lo;
        two_sum(g1, negated(g2), w_re_hi, w_re_lo);
        two_sum(g3, g4, w_im_hi, w_im_lo);
        two_sum(w_re_hi, a_re_hi, s_re_hi, s_re_lo);
        two_sum(w_im_hi, a_im_hi, s_im_hi, s_im_lo);
        const Four tr =
            _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_sub_pd(h1, h2), w_re_lo), s_re_lo), a_re_lo);
        const Four ti =
            _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(h3, h4), w_im_lo), s_im_lo), a_im_lo);
        Four T = _mm256_add_pd(absolute(h1), absolute(h2));
        T = _mm256_add_pd(T, absolute(w_re_lo));
        T = _mm256_add_pd(T, absolute(s_re_lo));
        T = _mm256_add_pd(T, absolute(a_re_lo));
        T = _mm256_add_pd(T, absolute(h3));
        T = _mm256_add_pd(T, absolute(h4));
        T = _mm256_add_pd(T, absolute(w_im_lo));
        T = _mm256_add_pd(T, absolute(s_im_lo));
        T = _mm256_add_pd(T, absolute(a_im_lo));
        // c = c x + t
        {
            const Four pr = products_difference(cr, xr, ci, xi), pi = products_sum(cr, xi, ci, xr);
            cr = _mm256_add_pd(pr, tr);
            ci = _mm256_add_pd(pi, ti);
        }
        terms_slope = _mm256_add_pd(_mm256_mul_pd(terms_slope, X), terms);
        terms = _mm256_add_pd(_mm256_mul_pd(terms, X), T);
        errors_slope = _mm256_add_pd(_mm256_mul_pd(errors_slope, X), errors);
        errors = _mm256_add_pd(_mm256_mul_pd(errors, X), a_error);
        sr = s_re_hi;
        si = s_im_hi;
    }
    constexpr size_t fields = 11;
    static_assert(sizeof(Sums) == fields * sizeof(double), "Sums is its eleven numbers");
    alignas(32) double lanes[fields][group];
    const Four results[fields] = {sr, si, cr, ci, dr, di, terms, terms_slope, values, errors, errors_slope};
    for (size_t f = 0; f < fields; ++f)
        _mm256_store_pd(lanes[f], results[f]);
    for (size_t m = 0; m < group; ++m)
        out[m] = {lanes[0][m], lanes[1][m], lanes[2][m], lanes[3][m], lanes[4][m], lanes[5][m],
                  lanes[6][m], lanes[7][m], lanes[8][m], lanes[9][m], lanes[10][m]};
}
#undef CERTPOLYSOLVE_AVX2_FMA
#undef CERTPOLYSOLVE_AVX2_FMA_INLINE
#endif
} // namespace certpolysolve::detail::compensated_kernel
