#pragma once
// The polynomial as a solve receives it: DensePolynomial (every coefficient from x^0 up), the public input object
// cps_polynomial (its nonzero terms by exponent), the largest supported exponent, and the input-limit exception.
// Expansion and exact polynomial operations are declared in polynomial_arithmetic.hpp.
#include <certpolysolve/certpolysolve.h>
#include "core/number.hpp"
#include <cstdint>
#include <map>
#include <stdexcept>
#include <vector>

namespace certpolysolve::detail {
// A polynomial as the solver reads it: every coefficient from x^0 up, dense. A solve expands the input object's terms
// (cps_polynomial below) into one; the routes build their own for factors, compositions and secular numerators.
struct DensePolynomial {
    std::vector<Coefficient> coefficients;
};
} // namespace certpolysolve::detail

// The public input object (the 0.5 spec's C40): its nonzero terms by exponent, so any exponent can be set and a
// polynomial of few terms stays small whatever its degree. An exact zero is no term; a ball whose value is zero is
// one. Several threads may read one at once (the spec's T3), so a read never fills a cache in it.
struct cps_polynomial {
    std::map<uint64_t, certpolysolve::detail::Coefficient> terms;
};
namespace certpolysolve::detail {
// The largest exponent a polynomial can hold: its degree is an int64_t (cps_polynomial_degree, cps_result_degree).
inline constexpr uint64_t max_exponent = static_cast<uint64_t>(INT64_MAX);
// A degree beyond what a coefficient list can address: the input's limit, not a numerical failure, so the session
// passes it to the caller as it is (CPS_OUT_OF_RANGE), wherever the expansion happens.
struct ExpansionRefused : std::out_of_range {
    using std::out_of_range::out_of_range;
};
} // namespace certpolysolve::detail
