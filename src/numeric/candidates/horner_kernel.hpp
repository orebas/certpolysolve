#pragma once
// The dense loop of the machine iteration's evaluation in binary64 (candidate_jets.hpp: PreparedJets<double, 1>'s
// evaluate is the reference), at four points at a time: the value, the derivative and the condition estimate's sum by
// Horner's scheme, each point's operations the reference's in the reference's order (platform/kernel_variant.hpp), so
// that it returns the reference's bits and raises the reference's floating-point flags.
#include "numeric/arithmetic/machine_arithmetic.hpp"
#include "platform/kernel_variant.hpp"
#include <cstddef>

namespace certpolysolve::detail::machine::horner_kernel {
// The points one call takes.
inline constexpr size_t group = 4;
// What the loop leaves at a point: Jet<double>'s value, derivative and absolute_sum.
struct Sums {
    double value_re = 0, value_im = 0, derivative_re = 0, derivative_im = 0, absolute_sum = 0;
};
#ifdef CERTPOLYSOLVE_HAVE_AVX2_FMA_KERNELS
// At the points zr[m] + i zi[m], m < group, with radius[m] the reference's norm of the point: the loop over
// coefficients[count - 1 .. 0] (count >= 1) and their moduli magnitudes, into out[m]. To be called only where
// kernel_variant() is KernelVariant::Avx2Fma: the function is built for AVX2.
void avx2(const Complex<double> *coefficients, const double *magnitudes, size_t count, const double *zr,
          const double *zi, const double *radius, Sums *out) noexcept;
#endif
} // namespace certpolysolve::detail::machine::horner_kernel
