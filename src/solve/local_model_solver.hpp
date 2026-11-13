#pragma once
#include "numeric/candidates/local_model_solver.hpp"

namespace certpolysolve::detail {
// Maintained secular recipe with geometric local restarts disabled, so a local
// solve cannot recursively request another local solve. Options are the caller's
// bounded local_model_options; counters and precision are returned even on stop.
cps_result solve_local_model(DensePolynomial &&, const cps_options &);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
cps_result solve_local_model_numerical(DensePolynomial &&, const cps_options &);
#endif
// The maintained inner recipe, including its fallback, with private per-action
// prices. No outer trigger, service callback or trial-local policy is inherited.
cps_result solve_local_model_priced(DensePolynomial &&, const cps_options &);
} // namespace certpolysolve::detail
