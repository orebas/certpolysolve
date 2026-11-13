#pragma once
// Stateless real-query contracts, independent of route/controller state.
#include "input/coefficients/polynomial_input.hpp"
#include "core/query.hpp"
#include "certify/real_region.hpp"
#include <memory>

namespace certpolysolve::detail::real_query {
// Whether `o` is a query this route answers on `input`: real_only under CPS_METHOD_AUTO on exact real
// coefficients (the caller checks that Auto routes by the selector).
bool applies(const cps_options &o, const DensePolynomial &input);
bool applies(const cps_options &o, const cps_polynomial &input);
// Whether every coefficient is exact and real: what this route needs, and what the all-real seed pass reads
// (real_line_seeds.hpp), which also seeds Auto's queries for all roots of such input (session.cpp).
bool exact_input(const DensePolynomial &input);
bool exact_input(const cps_polynomial &input);
// The nested complex query. Every real-only goal asks for individual roots: the complex query's accuracy goal
// does. Its bits are the caller's, except under the isolation goal, which asks no width: there the library's
// default goal, so that isolation costs on these routes what the default accuracy query costs. A lower goal saves
// only the last refinement and starts the precision ladder on other rungs, which can cost more than it saves (the
// D1 note, Part 1). The region is the route's to apply: the complex query asks for every root. History:
// docs/code-history-20261021.md, real_only_route.hpp 6.
cps_options complex_options(cps_options o);
// The query's region as the route applies it (real_region.hpp), or null when none is set. Its refusal() names an
// end the exact storage policy refuses; the caller reports it and builds no route.
std::unique_ptr<RealRegionFilter> region(const cps_options &o, const DensePolynomial &input);
std::unique_ptr<RealRegionFilter> region(const cps_options &o, const cps_polynomial &input);

} // namespace certpolysolve::detail::real_query
