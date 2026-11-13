#include "solve/complex_query.hpp"
#include "solve/query_execution.hpp"
#include "solve/solve_state.hpp"
#include "solve/prepared_piecewise.hpp"
#include "solve/policies.hpp"
#include "solve/real_line_routing.hpp"
#include "methods/direct/small_degree.hpp"

namespace certpolysolve::detail {
ComplexQueryState::ComplexQueryState(QueryInput &source) : source_(source) {}
ComplexQueryState::~ComplexQueryState() = default;
cps_result ComplexQueryState::advance(cps_options &next, uint64_t degree) {
    if (advanced_) throw std::logic_error("a complex query state is advanced once");
    advanced_ = true;
    auto &input = source_.input;
    auto &terms = source_.terms;
    auto &program = source_.program;
    auto &secular = source_.secular;
    const auto polynomial = [&] { return source_.polynomial(); };
    const auto dense_input = [&]() -> DensePolynomial & { return source_.dense_input(); };
    if (next.goal == CPS_ROOT_ISOLATION && !next.real_only) {
        // Validated exact coefficient input. This goal requests no width and
        // gets exclusive subdivision service, independently of the numerical selector.
        cps_options isolation = next;
        isolation.bits = 0; // no accuracy target; the driver starts at its 64-bit minimum
        isolation.absolute_tolerance = Number();
        global = std::make_unique<SolveState>(std::move(dense_input()), standalone_subdivision_policy());
        input.reset();
        auto result = global->advance(&isolation);
        result.diagnostics = "auto-route=complex-isolation (native subdivision); " + result.diagnostics;
        return result;
    }
    // The selector reads only the input and the query; the chosen
    // engine performs the one exact preparation.
    bool all_real_seeds = false;
    std::optional<RealLineSeedHint> all_real_hint;
    // Degree 1 to 4: the closed forms on exact input, before the selector (small_degree.cpp). A stopped
    // attempt retains its work and any proved cover. Fallback gets only the unspent allowance.
    std::optional<cps_result> direct_attempt;
    if (next.method == CPS_METHOD_AUTO && polynomial() && degree >= 1 && degree <= 4) {
        direct_attempt = direct_small_degree(dense_input(), next);
        if (direct_attempt) {
            const uint64_t spent = direct_attempt->progress.work_units_this_call;
            if (direct_attempt->goal_met || direct_attempt->termination == CPS_CANCELLED ||
                spent >= next.work_units) {
                direct_attempt->diagnostics = "auto-route=direct (closed form); " + direct_attempt->diagnostics;
                return std::move(*direct_attempt);
            }
            next.work_units -= spent;
        }
    }
    if (next.method == CPS_METHOD_AUTO && polynomial()) {
        // From the terms where the selection needs no coefficient list (a real-only query, a power
        // composite), else from the list.
        std::optional<AutoSelection> by_terms;
        if (!input) by_terms = select_auto_route(*terms, next);
        selection = by_terms ? *by_terms : select_auto_route(dense_input(), next);
        // Exact real coefficients: a factor whose roots are all real starts from the all-real seed pass's
        // seeds instead of the circles, on the secular route and on piecewise's hand-off to it, as the
        // real-only route's nested query does (docs/real-line-seed-pass-design-20261020.md). Read before
        // the input moves into the route below. Forced methods keep the circles.
        all_real_seeds = input ? real_query::exact_input(*input) : real_query::exact_input(*terms);
        // Where the selector sends such an input to piecewise, a query for all roots takes the secular route
        // instead when the pass seeds every root of the input's primitive polynomial: on all-real inputs
        // piecewise ran to its deadline where the seeded secular route finishes (docs/allroots-seed-routing-
        // note-20261020.md). Dense input only: a term form is not expanded for it. A pass that seeds fewer
        // or stops is kept for piecewise's hand-off to secular (below), whose start would run it again on
        // the same polynomial (dev/docs/slowdowns-20261029.md).
        if (all_real_seeds && input && selection.route == AutoRoute::Piecewise &&
            (next.goal == CPS_ROOT_ACCURACY || next.goal == CPS_ROOT_APPROXIMATION)) {
            all_real_hint = real_line_check(input->coefficients, next);
            if (all_real_hint && seeds_every_root(*all_real_hint)) {
                selection.route = AutoRoute::Secular;
                selection.reason = "every-root-seeded-before-piecewise";
            }
        }
    }
    if (next.method == CPS_METHOD_PIECEWISE || selection.route == AutoRoute::Piecewise) {
        auto config = maintained_piecewise_config();
        config.certificate_attempt = selection.attempt_first;
        config.certificate_attempt_fallback = selection.piecewise_operations;
        if (selection.route == AutoRoute::Piecewise)
            config.regime_exit_estimate = selection.piecewise_operations;
        else if (next.method == CPS_METHOD_PIECEWISE) {
            // The caller's piecewise leaves as Auto's route does where the selector would choose
            // piecewise (its estimate), and after the global route's minimum budget in binary64
            // operations where it would not (PiecewiseRootConfig::caller_method; chosen among
            // measured exit rules). The selector reads only the input and the query, not the method.
            // History: docs/code-history-20261021.md, session.cpp 8.
            config.caller_method = true;
            auto as_auto = next;
            as_auto.method = CPS_METHOD_AUTO;
            const auto would = select_auto_route(dense_input(), as_auto);
            if (would.route == AutoRoute::Piecewise)
                config.regime_exit_estimate = would.piecewise_operations;
        }
        // Piecewise hands its non-final stops to secular (below).
        fallback = std::make_unique<DensePolynomial>(dense_input());
        piecewise = std::make_unique<PreparedPiecewiseState>(std::move(dense_input()), std::move(config));
    } else if (program) {
        // Auto and classical take the classical recipe (refinement Aberth: the direct packet on the
        // program's evaluations), secular the secular recipe. No selector: Auto's coefficient
        // routing is untouched.
        global = std::make_unique<SolveState>(
            *program, maintained_method_policy(next.method == CPS_METHOD_SECULAR ? CPS_METHOD_SECULAR
                                                                                : CPS_METHOD_CLASSICAL));
    } else {
        const auto method = selection.route == AutoRoute::Secular ? CPS_METHOD_SECULAR : next.method;
        // Enable the complete global policy for direct exact coefficient Auto
        // queries. Staged real/program wrappers and piecewise's handoff below
        // retain their own recipes and separately stated whole-query accounts.
        const bool complete = next.method == CPS_METHOD_AUTO && polynomial() && source_.exact_coefficients() &&
                              !next.staged_goal && !next.program_lift_factor;
        const auto policy = complete ? complete_global_policy(method) : maintained_method_policy(method);
        global = polynomial() ? std::make_unique<SolveState>(std::move(dense_input()), policy) :
                                std::make_unique<SolveState>(*secular, policy);
        if (all_real_hint) global->hint_real_line_seeds(std::move(*all_real_hint));
    }
    input.reset();
    secular.reset();
    program.reset();
    cps_options routed = next;
    if (all_real_seeds) routed.real_line_seeds = true;
    auto result = piecewise ? piecewise->advance(next) : global->advance(&routed);
    // Piecewise hands off every stop that is not the caller's: the
    // regime exit (it spent more than a global solve's predicted cost),
    // its precision budget, or a storage or range refusal, which is
    // permanent for this state. Secular reaches the complete method within the
    // same budgets. A spent work budget or a cancellation is the
    // caller's stop and stays with piecewise.
    // History: docs/code-history-20261021.md, session.cpp 9.
    // Own the snapshot: the handoff below destroys the piecewise state before
    // using its refusal flags to describe the new route.
    const auto piecewise_stats = piecewise ? piecewise->statistics() : PiecewiseRootStatistics{};
    const bool hand_off = piecewise && fallback && !result.goal_met && !next.cancelled() &&
                          (piecewise_stats.regime_exits || piecewise_stats.storage_refusals ||
                           piecewise_stats.range_refusals ||
                           result.progress.stop_reason == CPS_STOP_PRECISION);
    if (hand_off) {
        // Its certificates do not transfer; the maintained secular
        // recipe solves the original input. The call's work budget is
        // charged once: secular gets what piecewise left.
        const uint64_t spent = result.progress.work_units_this_call;
        carried.carry(result);
        piecewise.reset();
        global = std::make_unique<SolveState>(std::move(*fallback),
                                              maintained_method_policy(CPS_METHOD_SECULAR));
        // The routing check's pass, where it ran: the factor that equals its polynomial starts from its
        // outcome instead of running it again (SolveState::Impl, the factors' preparation).
        if (all_real_hint) global->hint_real_line_seeds(std::move(*all_real_hint));
        fallback.reset();
        routed.method = CPS_METHOD_SECULAR;
        routed.work_units = routed.work_units > spent ? routed.work_units - spent : 0;
        handoff = piecewise_stats.regime_exits      ? "piecewise regime exit to secular; "
                  : piecewise_stats.storage_refusals ? "piecewise storage refusal to secular; "
                  : piecewise_stats.range_refusals   ? "piecewise range refusal to secular; "
                                                     : "piecewise precision budget to secular; ";
        result = global->advance(&routed);
        result.progress.work_units_this_call += spent;
        carried.add_to(result);
        result.diagnostics = handoff + result.diagnostics;
        result.handoff = piecewise_stats.regime_exits      ? CPS_HANDOFF_REGIME_EXIT
                         : piecewise_stats.storage_refusals ? CPS_HANDOFF_STORAGE
                         : piecewise_stats.range_refusals   ? CPS_HANDOFF_RANGE
                                                            : CPS_HANDOFF_PRECISION;
    }
    if (next.method != CPS_METHOD_AUTO) {
        const char *name = next.method == CPS_METHOD_CLASSICAL ? "classical" :
                           next.method == CPS_METHOD_SECULAR ? "secular" : "piecewise";
        result.diagnostics = std::string("method=") + name + "; " + result.diagnostics;
    } else if (*selection.reason) { // the selector's reason; empty where it did not run
        result.diagnostics = std::string("auto-route=") + detail::name(selection.route) + " (" +
                             selection.reason + "); " + result.diagnostics;
    }
    if (direct_attempt) {
        CarriedStatistics direct_statistics;
        direct_statistics.carry(*direct_attempt);
        direct_statistics.add_to(result);
        const uint64_t spent = direct_attempt->progress.work_units_this_call;
        result.progress.work_units_this_call = resources::add(result.progress.work_units_this_call, spent);
        result.diagnostics = "closed-form attempt spent " + std::to_string(spent) + " work units; " +
                             result.diagnostics;
    }
    return result;
}
cps_result solve_complex_query(const DensePolynomial &input, const cps_options &options) {
    if (options.real_only || options.region.set())
        throw std::invalid_argument("complex query requires all-complex roots without a region");
    QueryInput source(input);
    ComplexQueryState state(source);
    QueryExecution query(source, &options);
    if (query.refusal) return std::move(*query.refusal);
    return query.run([&] { return state.advance(query.next, query.degree); });
}
} // namespace certpolysolve::detail
