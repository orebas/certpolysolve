#pragma once
// The terms of the machine iteration's direct interaction sum in binary64, S = sum_{j != i} 1/(z_i - z_j)
// (machine_candidates.cpp; aberth_step.hpp has the term, machine::interaction_term), made four at a time: the
// differences, Smith's division and its three quotients run on four neighbours in AVX2 instructions, and the four
// terms are then added to the one accumulator in the neighbours' order, as the reference adds them. Every operation
// of a term is the reference's, and so is the order of the sum: the same bits and the same floating-point flags
// (platform/kernel_variant.hpp).
#include "numeric/arithmetic/machine_arithmetic.hpp"
#include "platform/kernel_variant.hpp"
#include <cstddef>

namespace certpolysolve::detail::machine::interaction_kernel {
// The neighbours one step takes.
inline constexpr size_t group = 4;
#ifdef CERTPOLYSOLVE_HAVE_AVX2_FMA_KERNELS
// Adds to (sum_re, sum_im) the terms 1/(z - points[j]) for j = first, first + 1, ..., a group at a time while a whole
// group lies below `last`, and returns the j it stopped at: `last` less a remainder below a group, or the first j of
// a group that holds a point equal to z (nothing of that group is computed or added: the caller's own loop then
// meets the coincident point as the reference does). To be called only where kernel_variant() is
// KernelVariant::Avx2Fma: the function is built for AVX2.
size_t add_groups(const Complex<double> *points, size_t first, size_t last, double z_re, double z_im, double &sum_re,
                  double &sum_im) noexcept;
#endif
} // namespace certpolysolve::detail::machine::interaction_kernel
