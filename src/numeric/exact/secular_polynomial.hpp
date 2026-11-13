#pragma once
// Exact numerator conversion, bounded by the caller's exact-work and storage limits.
#include "input/coefficients/polynomial_input.hpp"
#include "input/secular/secular.hpp"
#include "core/query.hpp"
#include <string>

namespace certpolysolve::detail {
bool secular_polynomial(DensePolynomial &out, const cps_secular &source, const cps_options &options,
                        std::string &diagnostic);
} // namespace certpolysolve::detail
