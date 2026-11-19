#pragma once
// The loop of the compensated certificate's evaluation (compensated_certificate.cpp: CompensatedCertificate's
// evaluations have the scheme, its bounds and what happens before and after the loop): the compensated Horner
// recurrences for the value, its correction and the slope, and the accumulators of their bounds. `reference` runs it
// at one point; `avx2_fma` at four points at a time, each point's operations the reference's in the reference's order
// (platform/kernel_variant.hpp), so that it returns the reference's bits.
//
// The certificate's other loops are real Horner recurrences on nonnegative coefficients, the majorants of the
// coefficients' errors, of p' and of p'' (may_certify, second_derivative_bound): `majorant` at one point,
// `majorants_avx2` at sixteen points at a time, again each point's operations the reference's in its order.
//
// The functions take plain arrays, and the file that builds the variants includes no project header beyond this one
// and the arithmetic the reference uses: their instructions stay inside the functions that have them.
#include "platform/kernel_variant.hpp"
#include <cstddef>
#include <cstdint>

namespace certpolysolve::detail::compensated_kernel {
// A coefficient as the certificate holds it: a_k 2^-scale = (re_hi + re_lo) + i (im_hi + im_lo) + delta,
// |delta| <= error, for every a_k in its ball.
struct Term {
    double re_hi = 0, im_hi = 0, re_lo = 0, im_lo = 0, error = 0;
};
// What the loop leaves at a point x, under the names of CompensatedCertificate's evaluation: the rounded value
// s = sr + i si, the correction c = cr + i ci, the slope d = dr + i di, and the accumulators of the bounds.
struct Sums {
    double sr = 0, si = 0, cr = 0, ci = 0, dr = 0, di = 0;
    double terms = 0, terms_slope = 0, values = 0, errors = 0, errors_slope = 0;
};
// The loop at x = xr + i xi with X >= |x|, on the coefficients terms[0 .. degree] (terms[degree] the leading one).
Sums reference(const Term *terms, uint64_t degree, double xr, double xi, double X) noexcept;
// The points one call of avx2_fma takes.
inline constexpr size_t group = 4;
#ifdef CERTPOLYSOLVE_HAVE_AVX2_FMA_KERNELS
// The loop at the points xr[m] + i xi[m] with X[m], m < group, into out[m]. To be called only where
// kernel_variant() is KernelVariant::Avx2Fma: the function is built for AVX2 and FMA.
void avx2_fma(const Term *terms, uint64_t degree, const double *xr, const double *xi, const double *X,
              Sums *out) noexcept;
#endif
// A majorant's Horner recurrence at x: the sum c_top, then sum = sum x + c_k for k = top - 1, ..., bottom, c_k
// being coefficients[k].
double majorant(const double *coefficients, uint64_t top, uint64_t bottom, double x) noexcept;
// The points one call of majorants_avx2 takes: four registers of four, since one recurrence is a chain of dependent
// operations and the processor runs four such chains side by side at no cost to any of them.
inline constexpr size_t majorant_group = 16;
#ifdef CERTPOLYSOLVE_HAVE_AVX2_FMA_KERNELS
// The same recurrence at the points x[m], m < majorant_group, into out[m]. To be called only where kernel_variant()
// is KernelVariant::Avx2Fma: the function is built for AVX2.
void majorants_avx2(const double *coefficients, uint64_t top, uint64_t bottom, const double *x,
                    double *out) noexcept;
#endif
} // namespace certpolysolve::detail::compensated_kernel
