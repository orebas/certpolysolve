// The direct interaction sum's terms, four at a time (interaction_kernel.hpp).
//
// The reference, a term (machine_candidates.cpp with machine::Complex's operators, machine_arithmetic.hpp): the
// difference b = z + (-z_j), part by part; then 1/b by Smith's division with the numerator a = (1, 0),
//     |b.im| < |b.re|:  r = b.im / b.re,  d = b.re + b.im r,  term = ((a.re + a.im r) / d, (a.im - a.re r) / d)
//     otherwise:        r = b.re / b.im,  d = b.im + b.re r,  term = ((a.re r + a.im) / d, (a.im r - a.re) / d)
// and sum = sum + term, part by part. Here p is the part of b the case divides by and q the other, so r = q / p and
// d = p + q r in both cases; a.im r = 0 r and a.re r = 1 r are made once, as each case makes them; the numerators of
// both cases are then formed and each lane takes its own. The numerators' sums and differences are exact (a zero or
// a one with a product by zero or by one), so forming the case a lane does not take raises no flag; every other
// operation is one the reference makes at that term, rounded as the reference's (IEEE 754, to nearest), and no
// operation is reordered or fused (the library's -ffp-contract=off; the function is built without the FMA
// instruction). The comparisons are quiet. A group with a zero difference is left to the caller before anything of
// it is divided. tests/test_kernel_variants.cpp compares the sums and the flags with the reference's.
#include "numeric/candidates/interaction_kernel.hpp"
#ifdef CERTPOLYSOLVE_HAVE_AVX2_FMA_KERNELS
#include <immintrin.h>

namespace certpolysolve::detail::machine::interaction_kernel {
__attribute__((target("avx2"))) size_t add_groups(const Complex<double> *points, size_t first, size_t last,
                                                  double z_re, double z_im, double &sum_re, double &sum_im) noexcept {
    static_assert(group == 4, "one AVX2 register holds four binary64 numbers");
    static_assert(sizeof(Complex<double>) == 2 * sizeof(double), "a point is its two parts, the real one first");
    using Four = __m256d;
    const Four sign = _mm256_set1_pd(-0.0), zero = _mm256_setzero_pd(), one = _mm256_set1_pd(1.0);
    const Four zr = _mm256_set1_pd(z_re), zi = _mm256_set1_pd(z_im);
    double sr = sum_re, si = sum_im;
    size_t j = first;
    for (; j + group <= last; j += group) {
        // Two loads of two points each; the unpacking leaves the neighbours in the lane order j, j + 2, j + 1, j + 3.
        const Four low = _mm256_loadu_pd(&points[j].re), high = _mm256_loadu_pd(&points[j + 2].re);
        const Four pr = _mm256_unpacklo_pd(low, high), pi = _mm256_unpackhi_pd(low, high);
        // b = z + (-z_j)
        const Four br = _mm256_add_pd(zr, _mm256_xor_pd(pr, sign)), bi = _mm256_add_pd(zi, _mm256_xor_pd(pi, sign));
        const Four coincident =
            _mm256_and_pd(_mm256_cmp_pd(br, zero, _CMP_EQ_OQ), _mm256_cmp_pd(bi, zero, _CMP_EQ_OQ));
        if (_mm256_movemask_pd(coincident))
            break;
        // |b.im| < |b.re|: the case that divides by the real part.
        const Four by_real = _mm256_cmp_pd(_mm256_andnot_pd(sign, bi), _mm256_andnot_pd(sign, br), _CMP_LT_OQ);
        const Four p = _mm256_blendv_pd(bi, br, by_real), q = _mm256_blendv_pd(br, bi, by_real);
        const Four r = _mm256_div_pd(q, p);
        const Four d = _mm256_add_pd(p, _mm256_mul_pd(q, r));
        const Four zero_r = _mm256_mul_pd(zero, r), one_r = _mm256_mul_pd(one, r); // a.im r and a.re r
        const Four numerator_re = _mm256_blendv_pd(_mm256_add_pd(one_r, zero), _mm256_add_pd(one, zero_r), by_real);
        const Four numerator_im = _mm256_blendv_pd(_mm256_sub_pd(zero_r, one), _mm256_sub_pd(zero, one_r), by_real);
        alignas(32) double term_re[group], term_im[group];
        _mm256_store_pd(term_re, _mm256_div_pd(numerator_re, d));
        _mm256_store_pd(term_im, _mm256_div_pd(numerator_im, d));
        // sum = sum + term, in the neighbours' order: lanes 0, 2, 1, 3.
        sr = sr + term_re[0], si = si + term_im[0];
        sr = sr + term_re[2], si = si + term_im[2];
        sr = sr + term_re[1], si = si + term_im[1];
        sr = sr + term_re[3], si = si + term_im[3];
    }
    sum_re = sr;
    sum_im = si;
    return j;
}
} // namespace certpolysolve::detail::machine::interaction_kernel
#endif
