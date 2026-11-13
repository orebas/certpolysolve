#pragma once
// Auto's direct route for exact polynomials of degree 1 to 4 (small_degree.cpp, small_degree_forms.cpp).
#include "input/coefficients/polynomial_input.hpp"
#include "core/query.hpp"
#include <optional>

namespace certpolysolve::detail {
// An attempt on an exact polynomial of degree 1 to 4 from the closed forms, in ball arithmetic at doubling precision,
// with the multiplicities decided exactly. An all-roots query under the accuracy or approximation goal only: no
// real-only query, no region, no isolation or cluster-radius goal. nullopt means no numerical attempt ran (including
// inapplicable input, refused exact preparation, or no remaining work). Otherwise the result carries ALL spent work
// and the last complete disjoint cover, if any. Only goal_met means success. A caller taking the general route after
// a stopped attempt must deduct that work and carry its statistics; a spent budget or cancellation returns the stop.
std::optional<cps_result> direct_small_degree(const DensePolynomial &input, const cps_options &options);
} // namespace certpolysolve::detail
