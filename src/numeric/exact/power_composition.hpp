#pragma once
#include "numeric/exact/exact_complex.hpp"
#include "input/coefficients/polynomial_input.hpp"

namespace certpolysolve::detail {
// Replace p by q only when the exact identity p(z)=q(z^k), k>1, holds
// and q(0)!=0. The latter avoids the ramified zero preimage.
ulong reduce_power_composition(IntegerPolynomial &polynomial);
ulong reduce_power_composition(GaussianPolynomial &polynomial);
// Read exact stored coefficients without expanding rational denominators or
// large exponents. Return 1 when the identity or nonzero constant is unknown.
ulong coefficient_power(std::span<const Coefficient> coefficients);

// Owned coefficient-index reduction. An absent result means no exact nontrivial
// identity was found; ordinary inputs allocate nothing. The coefficients retain
// their compressed representation, independently of denominator-expansion limits.
struct ExactPowerComposition {
    DensePolynomial polynomial;
    ulong power = 1;
};
std::optional<ExactPowerComposition> prepare_power_composition(std::span<const Coefficient>);

} // namespace certpolysolve::detail
