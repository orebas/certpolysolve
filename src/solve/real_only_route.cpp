#include "solve/prepared_piecewise.hpp"
// The real-only route: RealOnlyRoute answers a real-only Auto query on exact real polynomial input, under its goal
// and region, by the first of three routes that answers: the sparse recursion, the axis charts, the complex route.
//
// The constructors decide the attempts (Impl::decide), only under a work budget of at least
// heuristics::default_work_units_per_root a root and no region off the real axis: the sparse recursion on two terms
// or more where sparse_selector_units_per_multiplication x sparse_recursion_multiplications is below the estimate
// of the route select_auto_route takes for complex_query without the square-free proof (piecewise_operations, or
// initial_sweep_allowance global sweeps), and the axis charts (axis_only, at most axis_attempt_rounds rounds) where
// that route is piecewise. Impl::advance answers a region off the axis empty, without work; runs each attempt once,
// a stopped one's work and statistics carried (CarriedStatistics) and taken from the call's budget; then a fresh
// solve_complex_query under the selector (the accuracy goal; the default bits under the isolation goal).
// Each answer passes answer_from: real_roots keeps the roots proved real (CPS_REAL), the answer complete only if
// the route's met its goal, every kept root meets the caller's (meets_goal) and no other root's box meets the axis;
// RealRegionFilter::keep applies the region. While a kept root of even multiplicity has an end of the region in its
// interval, or (complex route) a box meeting the axis is not proved real, the route is solved afresh at twice the
// bits until the precision cap, the work budget, the sweep budget (complex route) or a cancellation stops it.
//
// State: RealOnlyRoute::Impl, for its session's one call: the input (input, or the borrowed terms that dense_input
// expands once), degree, zero_input, region, the attempts (sparse_attempt, axis_attempt, axis_config) and carried.
// Entry points: the constructors and advance, called by session.cpp. Stateless query predicates and
// the nested options belong to real_query.cpp; the complex engine route belongs to complex_query.cpp.
// Not here: the sparse recursion (SparseRealAttempt, sparse_real_attempt.hpp), the axis mode (piecewise_roots.cpp),
// the region's exact tests (real_region.hpp), the realness proofs (certification.cpp and the engines).
// Design notes: docs/api.md (section 7), docs/real-interval-contract-note-20261015.md (the isolation goal, regions),
// docs/b6-sparse-real-roots-note-20261015.md, docs/c1-axis-charts-note-20261015.md,
// docs/session-state-design-note-20261017.md (fresh solves).
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "solve/complex_query.hpp"
#include "solve/real_query.hpp"
#include "solve/real_only_route.hpp"
#include "solve/auto_selection.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include <algorithm>
#include <optional>
#include <string>

namespace certpolysolve::detail {

struct RealOnlyRoute::Impl {
    // The input's degree and whether it is the zero polynomial, as the session recorded them.
    uint64_t degree = 0;
    bool zero_input = false;
    // The input, kept for the routes' fresh solves: each solve of a route, a refinement's included, builds the
    // route again from it (docs/session-state-design-note-20261017.md: no solve state is advanced twice). History:
    // docs/code-history-20261021.md, real_only_route.cpp 1.
    std::unique_ptr<DensePolynomial> input;
    // Or the input object's terms, borrowed for the call (C40's second step), from which `input` is formed when a
    // route that reads every coefficient is first reached: the sparse recursion reads the terms.
    const cps_polynomial *terms = nullptr;
    const DensePolynomial &dense_input() {
        if (!input) input = std::make_unique<DensePolynomial>(dense(*terms));
        return *input;
    }
    // The attempts this query makes (decided at construction from the first query's budget and the selector's
    // route), each dropped after its stop, and the axis charts' configuration.
    bool sparse_attempt = false, axis_attempt = false;
    std::optional<PiecewiseRootConfig> axis_config;
    std::unique_ptr<RealRegionFilter> region;
    // What the attempts that stopped had counted (CarriedStatistics).
    CarriedStatistics carried;

    // The real roots in the complex route's result: its regions proved real.
    // A region on the axis not proved real (its reflection meets another
    // region) leaves the answer incomplete; smaller regions, at a higher
    // accuracy, settle it: `unresolved_out` reports such a region in an answer
    // that is otherwise complete, and answer_from refines until it is proved
    // one way or the other, the cap or the budget. History:
    // docs/code-history-20261021.md, real_only_route.cpp 2.
    // The working precision of a route's answer, at which its roots are tested and its signs start (RealRegionFilter).
    static slong answer_precision(const cps_result &result) {
        return static_cast<slong>(std::max<uint64_t>(64, result.work_bits));
    }
    static cps_result real_roots(cps_result result, const cps_options &query, bool *unresolved_out = nullptr) {
        std::vector<Root> real;
        bool unresolved = false;
        for (auto &root : result.roots) {
            if (root.realness == CPS_REAL)
                real.push_back(root);
            else if (arb_contains_zero(acb_imagref(root.box.value())))
                unresolved = true;
        }
        if (unresolved_out) *unresolved_out = unresolved && result.complete && result.goal_met;
        const auto precision = answer_precision(result);
        const bool answered = result.complete && result.goal_met && !unresolved &&
                              std::all_of(real.begin(), real.end(),
                                          [&](const Root &root) { return meets_goal(root, query, precision); });
        result.roots = std::move(real);
        result.goal = engine_goal(query);
        result.complete = result.goal_met = answered;
        if (!answered && result.termination == CPS_SUCCESS) {
            result.termination = CPS_PRECISION_LIMIT;
            result.progress.stop_reason = CPS_STOP_PRECISION;
            result.diagnostics = "a region on the real axis is not proved real at this accuracy; " + result.diagnostics;
        }
        return result;
    }
    // One route's answer: its roots proved real (real_roots), and under a region only those in it. A kept root's
    // interval that holds an end of the region that is not its root is placed by the input's signs where the root has
    // odd multiplicity (RealRegionFilter::place); while one of even multiplicity still holds such an end, the route
    // is solved afresh at twice the bits (the accuracy goal), within what is left of the call's work budget, until
    // the root lies on one side or a stop: on the complex route and the axis charts this refines every root, which is
    // the filter's cost (the D1 note). History: docs/code-history-20261021.md, real_only_route.cpp 3.
    // The call's sweep budget (cps_options_set_max_iterations) covers every advance of a route whose sweeps are under
    // it, the refinements' too: `route_sweeps` is that route's cumulative count of sweeps (the complex route's;
    // null for the attempts, whose work is not under the budget), and a refinement gets what the call has left.
    // History: docs/code-history-20261021.md, real_only_route.cpp 4.
    // With `refine_realness` (the complex route, whose answer is the final one; an attempt's unresolved region
    // sends the query on to the next route instead) a region on the axis not proved real is refined the same way.
    template <class Advance>
    cps_result answer_from(Advance &&advance, cps_options query, const cps_options &next, uint64_t *route_sweeps = nullptr,
                          bool refine_realness = false) {
        const uint64_t origin = route_sweeps ? *route_sweeps : 0;
        bool unresolved = false;
        // Each solve of the route is fresh, its counters from zero: its sweeps add to the route's, and the
        // statistics of the solves a refinement supersedes are carried into the answer (counters count the whole
        // solve, docs/api.md).
        CarriedStatistics superseded;
        uint64_t solves = 0;
        const auto answered = [&](cps_result raw) {
            ++solves;
            if (route_sweeps) *route_sweeps = resources::add(*route_sweeps, raw.statistics_v3.candidate_sweeps);
            return real_roots(std::move(raw), next, &unresolved);
        };
        auto result = answered(advance(query));
        if (!region && !refine_realness) return result;
        const uint64_t budget = query.work_units;
        uint64_t units = result.progress.work_units_this_call;
        for (;;) {
            bool undecided = false;
            if (region) {
                const auto outcome = region->keep(result.roots, next, answer_precision(result));
                if (outcome.refusal) {
                    result.complete = result.goal_met = false;
                    result.termination = CPS_UNSUPPORTED_QUERY;
                    result.progress.stop_reason = CPS_STOP_INPUT;
                    result.diagnostics = std::string(outcome.refusal) + "; " + result.diagnostics;
                    break;
                }
                undecided = outcome.undecided && result.goal_met;
            }
            const bool realness = refine_realness && unresolved;
            if (!undecided && !realness) break;
            const bool capped = query.bits >= query.max_work_bits;
            const uint64_t swept = route_sweeps ? *route_sweeps - origin : 0;
            const bool sweeps_spent = route_sweeps && next.max_iterations && swept >= next.max_iterations;
            if (capped || sweeps_spent || next.cancelled() || units >= budget) {
                result.complete = result.goal_met = false;
                result.termination = next.cancelled() ? CPS_CANCELLED : capped ? CPS_PRECISION_LIMIT : CPS_WORK_LIMIT;
                result.progress.stop_reason = next.cancelled() ? CPS_STOP_CANCELLED
                                              : capped         ? CPS_STOP_PRECISION
                                              : sweeps_spent   ? CPS_STOP_ITERATION_BUDGET
                                                               : CPS_STOP_WORK_BUDGET;
                if (undecided)
                    result.diagnostics = "a root's interval still holds an end of the region; " + result.diagnostics;
                break;
            }
            query.goal = CPS_ROOT_ACCURACY;
            query.bits = std::min(resources::multiply(query.bits, 2), query.max_work_bits);
            query.work_units = budget - units;
            if (route_sweeps && next.max_iterations) query.max_iterations = next.max_iterations - swept;
            superseded.carry(result);
            result = answered(advance(query));
            units = resources::add(units, result.progress.work_units_this_call);
        }
        result.progress.work_units_this_call = units;
        superseded.add_to(result);
        if (solves > 1)
            result.diagnostics =
                "refined by " + std::to_string(solves) + " fresh solves at doubled accuracy; " + result.diagnostics;
        return result;
    }
    // A region that misses the real axis holds no real root: an empty, complete answer, without work.
    cps_result empty_answer(const cps_options &next) const {
        cps_result result;
        result.degree = static_cast<int64_t>(degree);
        result.goal = engine_goal(next);
        result.complete = result.goal_met = true;
        result.termination = CPS_SUCCESS;
        result.progress.stop_reason = CPS_STOP_NONE;
        result.diagnostics = "the region misses the real axis: no real root lies in it";
        return result;
    }
    cps_result advance(const cps_options &next) {
        // (The zero polynomial's routes report it: every number is its root.)
        if (region && !region->meets_axis() && !zero_input) return empty_answer(next);
        cps_options query = real_query::complex_options(next);
        // The complex route's sweeps over its solves (answer_from: the solve's sweep budget).
        uint64_t complex_sweeps = 0;
        uint64_t spent = 0;
        std::string stopped;
        if (sparse_attempt) {
            // Under the isolation goal the recursion publishes its isolating intervals without refining them.
            cps_options sparse_query = query;
            if (next.goal == CPS_ROOT_ISOLATION) sparse_query.goal = CPS_ROOT_ISOLATION;
            auto attempt = answer_from(
                [&](const cps_options &q) {
                    SparseRealAttempt fresh = terms ? SparseRealAttempt(*terms) : SparseRealAttempt(*input);
                    return fresh.advance(q);
                },
                sparse_query, next);
            if (attempt.goal_met || next.cancelled()) {
                attempt.diagnostics = "real roots by the sparse recursion; " + attempt.diagnostics;
                return attempt;
            }
            // A stop of any kind, or an answer a tighter goal withdrew: its work and statistics are the session's,
            // and the other routes answer within what is left of the call's budget.
            spent += attempt.progress.work_units_this_call;
            carried.carry(attempt);
            sparse_attempt = false;
            stopped += "the sparse recursion stopped; ";
            query.work_units = next.work_units > spent ? next.work_units - spent : 0;
        }
        if (axis_attempt) {
            auto attempt = answer_from(
                [&](const cps_options &q) {
                    PreparedPiecewiseState fresh(dense_input(), *axis_config);
                    return fresh.advance(q);
                },
                query, next);
            if (attempt.goal_met || next.cancelled()) {
                attempt.progress.work_units_this_call += spent;
                carried.add_to(attempt);
                attempt.diagnostics = stopped + "real roots by the axis charts; " + attempt.diagnostics;
                return attempt;
            }
            // The attempt runs once, within the call: a stop of any kind, the call's spent budget included, ends it
            // (the prepared state's direct solve of a small factor starts again at every call, so an attempt kept
            // across calls of one unit never moved: test_real_refinement's close dyadic pair). Its work is the
            // session's: the complex route answers within what is left of the call's budget.
            spent += attempt.progress.work_units_this_call;
            carried.carry(attempt);
            axis_attempt = false;
            stopped += "the axis charts' attempt stopped; ";
            // spent is cumulative across both attempts; subtract it from
            // the original allowance, not the axis attempt's remainder.
            query.work_units = next.work_units > spent ? next.work_units - spent : 0;
        }
        auto result = answer_from(
            [&](cps_options q) {
                return solve_complex_query(dense_input(), q);
            },
            query, next, &complex_sweeps, true);
        result.progress.work_units_this_call += spent;
        carried.add_to(result);
        result.diagnostics = stopped + "real roots by the complex route; " + result.diagnostics;
        return result;
    }
    // The attempts this query makes, from the selector's route for the complex query (`select`) and the input's
    // number of nonzero terms (`count_support`), each asked only when the budget and the region allow an attempt.
    // The attempts, only under an ample work budget: at least the default's heuristics::default_work_units_per_root
    // a root (about four times the most a measured solve used; a default budget is raised to it, and the
    // benchmarks' unlimited one exceeds it). A tighter budget may hold just what the complex route needs, and an
    // attempt that spends it starves the route. History: docs/code-history-20261021.md, real_only_route.cpp 5. And
    // only where the selector sends the complex query to the piecewise method, its square-free proof aside (the
    // attempt's own preparation makes it): the axis charts are that method on the axis alone, and that group holds
    // the attempt's gain; where the selector chooses the global recipes, the attempt costs time whether it answers
    // or not. History: docs/code-history-20261021.md, real_only_route.cpp 6. And not for a region that misses the
    // real axis, which advance answers empty.
    template <class Select, class Support> void decide(const cps_options &first, Select &&select, Support &&count_support) {
        if (first.work_units >= resources::multiply(heuristics::default_work_units_per_root, degree) &&
            (!region || region->meets_axis())) {
            const auto complex_selection = select();
            // The sparse recursion, where its prediction is below the estimate of the route the selector would take
            // (heuristics.hpp: 800 of the selector's units a counted multiplication, the measured crossovers against
            // the axis charts).
            const uint64_t support = count_support();
            // The route's estimate, an exact count of operations; the comparison with the recursion's prediction (a
            // measured ratio) in binary64.
            const uint64_t route_estimate =
                complex_selection.route == AutoRoute::Piecewise
                    ? complex_selection.piecewise_operations
                    : resources::multiply(resources::multiply(heuristics::initial_sweep_allowance, degree),
                                          resources::add(degree, support));
            if (support >= 2 && heuristics::sparse_selector_units_per_multiplication *
                                        heuristics::sparse_recursion_multiplications(degree, support) <
                                    static_cast<double>(route_estimate))
                sparse_attempt = true;
            if (complex_selection.route == AutoRoute::Piecewise) {
                auto config = maintained_piecewise_config();
                config.axis_only = true;
                config.axis_rounds = heuristics::axis_attempt_rounds;
                axis_attempt = true;
                axis_config = std::move(config);
            }
        }
    }
};

RealOnlyRoute::RealOnlyRoute(DensePolynomial input, const cps_options &first, uint64_t degree, bool zero_input,
                             std::unique_ptr<RealRegionFilter> region)
    : impl_(std::make_unique<Impl>()) {
    auto &state = *impl_;
    state.degree = degree;
    state.zero_input = zero_input;
    state.region = std::move(region);
    state.input = std::make_unique<DensePolynomial>(std::move(input));
    state.decide(
        first, [&] { return select_auto_route(*state.input, real_query::complex_options(first), false); },
        [&] {
            uint64_t support = 0;
            for (const auto &c : state.input->coefficients) support += !c.zero();
            return support;
        });
}

RealOnlyRoute::RealOnlyRoute(const cps_polynomial &borrowed, const cps_options &first, uint64_t degree, bool zero_input,
                             std::unique_ptr<RealRegionFilter> region)
    : impl_(std::make_unique<Impl>()) {
    auto &state = *impl_;
    state.degree = degree;
    state.zero_input = zero_input;
    state.region = std::move(region);
    state.terms = &borrowed;
    state.decide(
        first,
        [&] {
            // From the terms where the selector needs no coefficient list (an admission its structure refuses: a
            // power composite's), else from the list, formed here once and kept for the routes.
            if (auto selection = select_auto_route(borrowed, real_query::complex_options(first), false)) return *selection;
            return select_auto_route(state.dense_input(), real_query::complex_options(first), false);
        },
        [&] { return static_cast<uint64_t>(borrowed.terms.size()); });
}

RealOnlyRoute::~RealOnlyRoute() = default;

cps_result RealOnlyRoute::advance(const cps_options &next) { return impl_->advance(next); }

} // namespace certpolysolve::detail
