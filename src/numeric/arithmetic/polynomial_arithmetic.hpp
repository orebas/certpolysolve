#pragma once
// Exact polynomial arithmetic on the input object (the 0.5 spec's C39 and C40), so that a caller can build up
// products, polynomials from sets of roots and recurrences of their own, sparse input included. Sums, differences,
// products, scaling by a number, powers, derivatives and products of linear factors, over the Gaussian rationals,
// on the object's nonzero terms: a polynomial of few terms stays few terms whatever its degree. A product (and each
// step of a power) is formed term by term when that pairs no more terms than the dense result has coefficients, and
// by FLINT's dense multiplication otherwise; products of linear factors are dense. The operands' coefficients must
// be exact (std::invalid_argument otherwise); a result's degree beyond max_exponent (2^63 - 1), or a dense product
// beyond what FLINT can address, is std::out_of_range.
// History: docs/code-history-20261021.md, polynomial_arithmetic.hpp 1.
#include "input/coefficients/polynomial_input.hpp"
#include "core/number.hpp"
#include <stdexcept>
#include <cstdint>
#include <span>

namespace certpolysolve::detail {
// The input object's coefficients from x^0 to its degree, the form the routes that read every coefficient take
// (polynomial_arithmetic.cpp). ExpansionRefused for a degree beyond what a vector can address.
DensePolynomial dense(const cps_polynomial &p);
cps_polynomial polynomial_add(const cps_polynomial &a, const cps_polynomial &b);
cps_polynomial polynomial_subtract(const cps_polynomial &a, const cps_polynomial &b);
cps_polynomial polynomial_multiply(const cps_polynomial &a, const cps_polynomial &b);
cps_polynomial polynomial_scale(const cps_polynomial &a, const Coefficient &factor);
cps_polynomial polynomial_power(const cps_polynomial &a, uint64_t exponent);
cps_polynomial polynomial_derivative(const cps_polynomial &a);
// The product of (x - roots[k])^multiplicities[k]; an empty list gives the polynomial 1.
cps_polynomial polynomial_from_roots(std::span<const Coefficient> roots, std::span<const uint64_t> multiplicities);
} // namespace certpolysolve::detail
