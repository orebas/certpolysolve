#pragma once
// Ordinary engine orchestration and handoffs, independent of RealOnlyRoute and
// SessionState. Its input borrow lasts until destruction; a session owns both in
// declaration order, and a fresh complex query owns both on its stack.
#include "solve/auto_selection.hpp"
#include "core/session_statistics.hpp"
#include "input/coefficients/polynomial_input.hpp"
#include <memory>

namespace certpolysolve::detail {
struct QueryInput;
class SolveState;
class PreparedPiecewiseState;
class ComplexQueryState {
    QueryInput &source_;
    std::unique_ptr<SolveState> global;
    std::unique_ptr<PreparedPiecewiseState> piecewise;
    std::unique_ptr<DensePolynomial> fallback;
    // What the engines before the current one counted (CarriedStatistics).
    CarriedStatistics carried;
    std::string handoff;
    AutoSelection selection;
    bool advanced_ = false;
  public:
    explicit ComplexQueryState(QueryInput &source);
    ~ComplexQueryState();
    ComplexQueryState(const ComplexQueryState &) = delete;
    ComplexQueryState &operator=(const ComplexQueryState &) = delete;
    // Caller has validated/resolved next and owns its cancellation latch.
    // D5 also handles the legacy real/native-secular cases outside RealOnlyRoute.
    cps_result advance(cps_options &next, uint64_t degree);
};
// A fresh all-complex solve, also used for real-route refinement. Isolation
// requires exact coefficients. No route/session state is reused; budgets are
// the caller's remaining allowances.
cps_result solve_complex_query(const DensePolynomial &, const cps_options &);
} // namespace certpolysolve::detail
