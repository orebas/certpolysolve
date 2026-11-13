#pragma once
// What one call asks and what it answers: the cancellation callback and token, cps_options (the goal, the budgets,
// the region, the cancellation and its checkpoint), cps_result (the cover, the stop, the statistics and C14's
// facts), and the three helpers every route reads the options through: the call's work units, the goal resolved
// for the input kind, and the goal an engine applies. Part of internal.hpp.
#include <certpolysolve/certpolysolve.h>
#include "core/cover.hpp"
#include "core/engine_statistics.hpp"
#include "core/heuristics.hpp"
#include "core/number.hpp"
#include "platform/parallel.hpp"
#include "core/resource_limits.hpp"
#include "core/heuristic_budget.hpp"
#include "core/research_scheduling_metrics.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace certpolysolve::detail {
// A cancellation poll: nonzero cancels. Internal: the public API has no callback (C4); the token, the time limit
// and certpolysolve_testing.h's counter are polled through cps_options::cancelled, and the tests set one.
// History: docs/code-history-20261021.md, internal.hpp 1.
using CancelCallback = int (*)(void *);
} // namespace certpolysolve::detail
// A cancellation token: the flag the library owns (cps_cancel_token_create); options that name it share it.
struct cps_cancel_token {
    std::shared_ptr<std::atomic<bool>> flag = std::make_shared<std::atomic<bool>>(false);
};
struct cps_options {
    // CPS_GOAL_DEFAULT until a goal is set: every solve resolves it for its input kind first (with_resolved_goal
    // below), so no engine reads it.
    cps_goal goal = CPS_GOAL_DEFAULT;
    cps_method method = CPS_METHOD_AUTO;
    // max_iterations: the sweep budget of one call (0: none), counted by the candidate statistics' sweeps of the
    // solve state that owns the call. A state is advanced once and counts no sweep before its advance, so its count
    // is the call's (a hand-off and a nested solve are fresh states, with their part of the budget).
    uint64_t bits = CPS_DEFAULT_BITS, max_work_bits = CPS_DEFAULT_MAX_WORK_BITS, max_iterations = 0;
    uint64_t exact_budget = UINT64_MAX;
    uint64_t exact_operations = UINT64_MAX, temporary_bytes = UINT64_MAX;
    bool allow_exponent_expansion = false;
    uint64_t work_units = CPS_DEFAULT_WORK_UNITS, live_cells = CPS_DEFAULT_LIVE_CELLS;
    // The work budget was left at its default (cps_options_set_work_budget never called, or with 0):
    // a session then allows at least heuristics::default_work_units_per_root
    // units a root (SessionState, the input's degree), never fewer than
    // work_units.
    bool work_units_default = true;
    // The goal is a stage of the caller's question (SecularRefinement::step stages the library's default goal first; the
    // real-only route answers from a complex cover, real_query::complex_options): the certificate attempt's cluster start
    // is left out there (certificate_attempt.cpp), since what it certifies at the stage sends the caller's question down
    // a path the stage's iteration would not.
    bool staged_goal = false;
    certpolysolve::detail::Number absolute_tolerance;
    certpolysolve::detail::Number cluster_radius;
    // The literals as the caller last set them, trimmed, for the getters (C17): "0" when nothing was set; the
    // region's ends empty where a side is unbounded.
    std::string absolute_tolerance_text = "0", cluster_radius_text = "0";
    std::array<std::string, 4> region_text;
    // cps_options_set_region's closed rectangle, each side absent where unbounded (none: the whole plane).
    // Real-only queries use its trace on the real axis (real_region.hpp).
    struct Region {
        std::optional<certpolysolve::detail::Number> real_lower, real_upper, imag_lower, imag_upper;
        bool set() const { return real_lower || real_upper || imag_lower || imag_upper; }
    };
    Region region;
    bool real_only = false;
    // The all-real seed pass (real_line_seeds.hpp, docs/real-line-seed-pass-design-20261020.md): an exact real factor
    // starts from the pass's seeds where it gives any. Set for Auto's queries on exact real coefficients: by the
    // real-only route on its nested complex query (real_only_route.cpp, complex_query) and by the session for a
    // query for all roots (session.cpp); forced methods, programs and secular input leave it off.
    bool real_line_seeds = false;
    // A factor of the program route's product lift, solved by a nested solve state (solver_program.cpp): its finish
    // proves no region real, so the lift unites, merges and judges the factors' covers as published before the proof
    // existed, and the outer solve's finish proves realness on the cover the lift publishes. Set by the lift alone.
    bool program_lift_factor = false;
    // Internal synchronous borrow. A local model uses the outer query's
    // remaining heuristic allowance; no child can renew it. Never cancellation.
    certpolysolve::detail::HeuristicBudget *heuristic_budget = nullptr;
    certpolysolve::detail::CancelCallback cancel = nullptr;
    void *cancel_data = nullptr;
    // The cancellation token's flag (cps_cancel_token; shared, so that options outlive the caller's handle), the time
    // limit in seconds (0: none), and the deadline a solve call fixes from it at its entry (c_api.cpp, timed_options):
    // every solve the call runs copies these options and so keeps the call's deadline instead of starting its own.
    std::shared_ptr<const std::atomic<bool>> cancel_flag;
    double time_limit = 0;
    std::optional<std::chrono::steady_clock::time_point> deadline;
    // certpolysolve_testing.h's cps_testing_options_stop_at_poll: 0, counting off; k, cancelled from the call's k-th poll
    // on. The solve call installs its counter as the callback (c_api.cpp, timed_options), so it is polled exactly
    // where a caller's callback is, and the internal wrappers and clears of the callback apply to it alike.
    uint64_t testing_stop_at_poll = 0;
    // cps_options_set_workers: the threads a solve call may run its split loops on (W9, parallel.hpp); 1, none started.
    uint32_t workers = 1;
    // A checkpoint's question: the token first (one atomic load), the clock only when a limit is set, the caller's
    // callback last. With neither a token nor a limit this is the callback alone.
    // History: docs/code-history-20261021.md, internal.hpp 3.
    bool cancelled() const {
        if (cancel_flag && cancel_flag->load(std::memory_order_relaxed)) return true;
        if (deadline && std::chrono::steady_clock::now() >= *deadline) return true;
        return cancel && cancel(cancel_data);
    }
    // The token and the clock alone: what a block of a split loop polls on a worker (W9, parallel.hpp). The callback,
    // a latch of the solve's thread, stays with that thread's own polls (cancelled()).
    bool interrupted() const {
        if (cancel_flag && cancel_flag->load(std::memory_order_relaxed)) return true;
        return deadline && std::chrono::steady_clock::now() >= *deadline;
    }
    uint64_t exact_storage_limit() const {
        return std::min(exact_budget, certpolysolve::detail::resources::multiply(temporary_bytes, 8));
    }
    // The growth that writing one compressed number out as a rational may add to the storage the input owns, under
    // the compact-exponent rule (compact_expansion_fits, docs/resource-policy.md): one precision cap. The numerical
    // solve holds each coefficient at no more than max_work_bits bits, so exact preparation may write each number
    // out with as much again. With allow_exponent_expansion the rule does not apply.
    uint64_t exact_expansion_allowance() const { return max_work_bits; }
    bool exact_effort_fits(uint64_t estimate) const { return estimate != UINT64_MAX && estimate <= exact_operations; }
    // The sweeps the call may still make under max_iterations. `statistics` is the calling state's own.
    uint64_t sweeps_left(const certpolysolve::detail::CandidateStatistics &statistics) const {
        if (!max_iterations)
            return UINT64_MAX;
        const uint64_t used = statistics.candidate_sweeps;
        return used < max_iterations ? max_iterations - used : 0;
    }
    // The sweeps a packet of the given length may make: its length, or what
    // the budget leaves. Zero means the budget is spent and no packet runs.
    unsigned packet_sweeps(unsigned length, const certpolysolve::detail::CandidateStatistics &statistics) const {
        return static_cast<unsigned>(std::min<uint64_t>(length, sweeps_left(statistics)));
    }
};
struct cps_result {
    int64_t degree = -1;
    cps_termination termination = CPS_PRECISION_LIMIT;
    bool complete = false, goal_met = false;
    cps_goal goal = CPS_ROOT_ACCURACY;
    uint64_t work_bits = 0;
    std::vector<certpolysolve::detail::Root> roots;
    std::string diagnostics;
    certpolysolve::detail::StageStatistics statistics{};
    certpolysolve::detail::SearchStatistics statistics_v2{};
    certpolysolve::detail::CandidateStatistics statistics_v3{};
    certpolysolve::detail::SearchProgress progress{};
    // A synchronous local solver may return its already-priced action totals.
    // The parent consumes that charge once instead of repricing aggregate
    // counters at the final precision. No borrowed observer escapes the child.
    bool research_priced_work_available = false;
    certpolysolve::detail::ResearchSchedulingMetrics research_scheduling_metrics{};
    std::vector<cps_component_progress> components;
    // Set by the C entry once the call returns (c_api.cpp): the call's wall time on its monotonic clock
    // (cps_statistics.total_ns), and the cancellation polls it made while the testing counter was on.
    uint64_t total_ns = 0, testing_polls = 0;
    // The jobs each split loop ran on the call's workers (W9, parallel.hpp; certpolysolve_testing.h).
    certpolysolve::detail::SplitJobs testing_split_jobs{};
    // The v1 counters of the engines that handed over before the answering one (session_statistics.hpp; every
    // field but exact_factors): the 0.5 totals add them, the old v1 getter does not.
    certpolysolve::detail::StageStatistics carried_v1{};
    // C14's facts (certpolysolve.h): set where the diagnostics state them, by the engine that answers.
    uint64_t sufficient_work_bits = 0;
    cps_squarefree squarefree = CPS_SQUAREFREE_UNKNOWN;
    cps_handoff handoff = CPS_HANDOFF_NONE;
    // The regions' ends as the API's text, made once when the result is published (certpolysolve::detail::publish_text,
    // which every C entry's finish_call runs after the regions take their published order).
    std::vector<std::array<certpolysolve::detail::EndpointText, 4>> exported;
};

namespace certpolysolve::detail {
// The text of every region's ends (cps_result::exported), from the binary ends; text is made nowhere else.
void publish_text(cps_result &result);
// A result as a C entry hands it out, its text made: for the tests that publish a solve state's answer.
cps_result *published_result(cps_result result);
// The work units one call of a session may spend: the caller's budget, or, when
// it was left at its default, heuristics::default_work_units_per_root a root of
// the input, never fewer than the fixed default (cps_options::work_units_default).
inline uint64_t call_work_units(const cps_options &options, uint64_t degree) {
    if (!options.work_units_default)
        return options.work_units;
    return std::max(options.work_units, resources::multiply(heuristics::default_work_units_per_root, degree));
}
// The goal a solve applies (C13): CPS_GOAL_DEFAULT is CPS_ROOT_APPROXIMATION (counted clusters) for a program, the
// default the program spec gives (docs/slp-input-spec-20261015.md, section 2), and CPS_ROOT_ACCURACY for coefficient
// and secular input. Resolved on every call, by the session before it validates the call and by SolveState as it
// configures, so a call with null options and a call with untouched options read the same.
inline cps_options with_resolved_goal(cps_options options, bool program) {
    if (options.goal == CPS_GOAL_DEFAULT) options.goal = program ? CPS_ROOT_APPROXIMATION : CPS_ROOT_ACCURACY;
    return options;
}
// The goal an engine applies. Every solve resolves CPS_GOAL_DEFAULT for its input kind before an engine runs
// (with_resolved_goal), so this is the options' own goal there; an engine driven directly with the default still
// unresolved (the tests that construct an engine's options) reads it as CPS_ROOT_ACCURACY, the coefficient default.
// History: docs/code-history-20261021.md, internal.hpp 14.
inline cps_goal engine_goal(const cps_options &options) {
    return options.goal == CPS_GOAL_DEFAULT ? CPS_ROOT_ACCURACY : options.goal;
}
} // namespace certpolysolve::detail
