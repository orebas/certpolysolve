#pragma once
// The solve's entry points and states: solve() on a coefficient list or the borrowed input object, SolveState (the
// global driver, solver.cpp), compare_numbers, maintained_method_policy (policies.cpp), the start of a precision stop's
// note, and SessionState (a call's route: the selector, the hand-offs, session.cpp). Part of internal.hpp.
#include <certpolysolve/certpolysolve.h>
#include "core/engine_policy.hpp"
#include "solve/policies.hpp"
#include "core/number.hpp"
#include "input/coefficients/polynomial_input.hpp"
#include "numeric/exact/prepared_structure.hpp"
#include "core/query.hpp"
#include <memory>
#include <limits>
#include <optional>

namespace certpolysolve::detail {
cps_result solve(const DensePolynomial &polynomial, const cps_options &options);
cps_result solve(DensePolynomial &&polynomial, const cps_options &options);
// The input object's terms, borrowed for the call (C40's second step): the session expands them only when it builds
// a route that reads every coefficient, so a real-only query the sparse recursion answers never forms the list.
cps_result solve(const cps_polynomial &borrowed, const cps_options &options);
struct RealLineSeedHint; // real_line_routing.hpp
class SolveState {
    struct Impl;
    std::unique_ptr<Impl> impl_;

  public:
    explicit SolveState(const DensePolynomial &, EnginePolicy = production_policy());
    explicit SolveState(DensePolynomial &&, EnginePolicy = production_policy());
    explicit SolveState(const cps_secular &, EnginePolicy = production_policy());
    // Program input (docs/slp-input-spec-20261015.md): the structured route on a ProgramSource, in Aberth mode
    // under the classical recipe's policy and in secular mode under the secular recipe's (session.cpp). The
    // program must have an output (std::invalid_argument otherwise); its degree is proved at the first advance.
    explicit SolveState(const cps_program &, EnginePolicy);
    ~SolveState();
    cps_result advance(const cps_options *);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    struct ResearchCheckpoint {
        bool prepared = false, warm_eligible = false;
        uint64_t factors = 0, ready_factors = 0;
    };
    struct ResearchSubdivisionTurn {
        bool suspended = false, idle = false;
        double pending_work = 0;
    };
    // Coefficient-only private driver. Each turn is one existing safe action;
    // these are not resumable public solves. Settings/cancellation remain owned.
    void research_start(const cps_options &);
    // Install before the first private turn. The policy must explicitly
    // enable research_local_packet_meter. data is borrowed until this state
    // is finished/destroyed; callbacks are synchronous and never forked.
    void research_set_local_packet_service(ResearchLocalPacketService, void *);
    bool research_step(uint64_t remaining_work = UINT64_MAX, uint64_t cell_limit = UINT64_MAX,
                       double subdivision_allowance = std::numeric_limits<double>::infinity());
    const cps_result &research_result() const;
    const cps_result &research_finish(cps_stop_reason = CPS_STOP_WORK_BUDGET);
    uint64_t research_live_cells() const;
    double research_native_work() const;
    // Includes the current outer action and nested search units whose public
    // progress merge has not happened yet. Safe to read from a service hook.
    uint64_t research_inflight_work_units() const;
    bool research_local_service_stopped() const;
    uint64_t research_engine_bits() const;
    bool research_handoff_requested() const;
    ResearchCheckpoint research_checkpoint() const;
    ResearchSubdivisionTurn research_subdivision_turn() const;
    bool research_warm_handoff(uint64_t fallback_bits = 0);
    // Copy only exact prepared factors and their certified covers. Numerical
    // candidates remain owned by the parent; no algebra or numerical solve is
    // repeated. The caller owns the observation callback for the fork's life.
    // A bounded M4 query must retain its factor-owned search instead: forking
    // that query is refused, so its synchronous budget borrow cannot escape.
    std::unique_ptr<SolveState> research_warm_fork(uint64_t fallback_bits,
        void (*trace)(void *, const TraceEvent &) = nullptr, void *trace_data = nullptr) const;
#endif
    // The all-real seed pass's outcome on the input, computed before routing (session.cpp): the factor whose integer
    // polynomial equals the hint's takes its seeds instead of running the pass again. Before the first advance.
    void hint_real_line_seeds(RealLineSeedHint &&);
    // Inspect committed preparation without running or repeating algebra.
    // Call explicitly when needed; normal solves pay no observation cost.
    std::optional<PreparedStructure> prepared_structure() const;
};
// The exact order of two numbers (-1, 0, 1): equal representations, then their exact rationals within the
// options' exact storage policy, then ball enclosures at doubling precision up to options.max_work_bits; nothing
// where all three fail (two distinct numbers agreeing to the precision cap, or equal ones too large to expand).
std::optional<int> compare_numbers(const Number &a, const Number &b, const cps_options &options);
// Auto's route: the selector (auto_selection.hpp), adopted
// after the selector qualification (docs/selector-campaign-20260928.md).
// D5 remains the selector's answer for the queries the explicit recipes do
// not support (AutoRoute::D5, docs/d5-route-note-20260930.md). The switch
// that let a session skip the selector went in W8's S2.
// History: docs/code-history-20261021.md, internal.hpp 15.
// The start of what a precision stop's diagnostics say about the cap (solver.cpp).
inline constexpr const char *precision_cap_note = "stopped at the precision cap of ";
class SessionState {
    struct Impl;
    std::unique_ptr<Impl> impl_;
  public:
    explicit SessionState(const DensePolynomial &);
    explicit SessionState(DensePolynomial &&);
    // Borrows the input object's terms for the one advance (cps_solve's call): the caller keeps the object alive and
    // unchanged during it, as the C API requires; the session expands the terms only for a route that reads every
    // coefficient.
    explicit SessionState(const cps_polynomial &borrowed);
    explicit SessionState(const cps_secular &);
    explicit SessionState(const cps_program &);
    ~SessionState();
    cps_result advance(const cps_options *);
};
} // namespace certpolysolve::detail
