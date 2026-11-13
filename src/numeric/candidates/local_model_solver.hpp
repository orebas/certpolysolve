#pragma once
// A synchronous capability for proposing roots of an owned local polynomial.
// The controller supplies the implementation; a numerical method cannot name or
// construct the controller's state. The caller supplies each attempt's options
// and owns cancellation state until the call returns. No callback or borrow is
// retained by the implementation. Results concern only this local polynomial:
// the original source must independently certify any hint used by a method.
#include "core/query.hpp"
#include "input/coefficients/polynomial_input.hpp"

namespace certpolysolve::detail {
// A function reference makes a capability mandatory and excludes captured state
// and null dispatch. Calls happen once per local-model attempt, outside kernels.
using LocalModelSolver = cps_result (&)(DensePolynomial &&, const cps_options &);
} // namespace certpolysolve::detail
