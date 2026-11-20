// The machine iteration's dense Horner loop at four points at a time (horner_kernel.hpp).
//
// The reference is PreparedJets<double, 1>::evaluate's dense loop (candidate_jets.hpp), with machine::Complex's
// operators (machine_arithmetic.hpp) written out: a product a b is (a.re b.re - a.im b.im, a.re b.im + a.im b.re),
// each of its four products and its two sums rounded, and a sum adds the parts. A step is
//     derivative = derivative z + value;  value = value z + coefficient;  absolute_sum = absolute_sum radius + modulus
// in that order, the derivative from the value before its own step. Each operation here is the AVX2 instruction that
// makes the scalar one at four points, each lane rounded as the scalar operation rounds (IEEE 754, to nearest); no
// operation is reordered or fused (the function is built without the FMA instruction, and the library with
// -ffp-contract=off). tests/test_kernel_variants.cpp compares the two on every field.
#include "numeric/candidates/horner_kernel.hpp"
#ifdef CERTPOLYSOLVE_HAVE_AVX2_FMA_KERNELS
#include <immintrin.h>

namespace certpolysolve::detail::machine::horner_kernel {
__attribute__((target("avx2"))) void avx2(const Complex<double> *coefficients, const double *magnitudes,
                                          size_t count, const double *zr, const double *zi, const double *radius,
                                          Sums *out) noexcept {
    static_assert(group == 4, "one AVX2 register holds four binary64 numbers");
    using Four = __m256d;
    const Four xr = _mm256_loadu_pd(zr), xi = _mm256_loadu_pd(zi), R = _mm256_loadu_pd(radius);
    Four vr = _mm256_broadcast_sd(&coefficients[count - 1].re), vi = _mm256_broadcast_sd(&coefficients[count - 1].im);
    Four dr = _mm256_setzero_pd(), di = _mm256_setzero_pd();
    Four sum = _mm256_broadcast_sd(&magnitudes[count - 1]);
    for (size_t j = count - 1; j-- > 0;) {
        // derivative = derivative z + value
        const Four nr = _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(dr, xr), _mm256_mul_pd(di, xi)), vr);
        const Four ni = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(dr, xi), _mm256_mul_pd(di, xr)), vi);
        dr = nr;
        di = ni;
        // value = value z + coefficient
        const Four ur = _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(vr, xr), _mm256_mul_pd(vi, xi)),
                                      _mm256_broadcast_sd(&coefficients[j].re));
        const Four ui = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(vr, xi), _mm256_mul_pd(vi, xr)),
                                      _mm256_broadcast_sd(&coefficients[j].im));
        vr = ur;
        vi = ui;
        // absolute_sum = absolute_sum radius + modulus
        sum = _mm256_add_pd(_mm256_mul_pd(sum, R), _mm256_broadcast_sd(&magnitudes[j]));
    }
    constexpr size_t fields = 5;
    static_assert(sizeof(Sums) == fields * sizeof(double), "Sums is its five numbers");
    alignas(32) double lanes[fields][group];
    _mm256_store_pd(lanes[0], vr);
    _mm256_store_pd(lanes[1], vi);
    _mm256_store_pd(lanes[2], dr);
    _mm256_store_pd(lanes[3], di);
    _mm256_store_pd(lanes[4], sum);
    for (size_t m = 0; m < group; ++m)
        out[m] = {lanes[0][m], lanes[1][m], lanes[2][m], lanes[3][m], lanes[4][m]};
}
} // namespace certpolysolve::detail::machine::horner_kernel
#endif
