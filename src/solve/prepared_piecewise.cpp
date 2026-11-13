#include "certify/attempts/newton_cover.hpp"
#include "core/engine_statistics.hpp"
#include "solve/prepared_piecewise.hpp"
// The piecewise method on a whole exact input (PreparedPiecewiseState; solve_prepared_piecewise_roots is one advance
// of a fresh state): an exact preparation, the input or its square-free factors solved, their covers mapped back.
//
// Impl::advance spends its first work unit on prepare: prepare_squarefree_input (exact_preprocessing.cpp) gives the
// zero order and a proof that the nonzero part is square-free, or its exact square-free factors, each q(x^power) to a
// multiplicity. One direct PiecewiseRootState takes the whole input if that part is proved square-free with no power
// reduction (and is not of degree 1 or 2) or the preparation is unavailable; otherwise a factor of degree 3 or more
// gets a state, one of degree 1 or 2 the linear or quadratic formula in ball arithmetic (small_squarefree_factor).
// Under config.certificate_attempt a direct state is preceded once by certificate_attempt (machine_newton_cover: a
// binary64 Aberth packet, then a certificate at each point), whose cover is the answer if it meets the goal.
//
// The factor loop advances each factor short of the goal, round-robin from cursor, and maps the covers back: the exact
// zero root; at a budget, cap or regime-exit stop, an unfinished factor's cover completed (complete_cover_at_stop); a
// complete cover of q lifted to the roots of q(z^power) (lift_power_cover); counts times the multiplicity. When every
// factor's cover is complete and mapped, overlapping regions are merged (merge_overlapping) and the cover is
// published if it meets the goal or prefer_cover ranks it ahead of the answer, as the global driver publishes;
// otherwise the last complete cover stays, the input's square before any (publish_input_square). When
// every factor meets its goal and the mapped cover does not, only the factor whose region is the wider of the first
// meeting pair gets a tighter goal (its separation bits, separation.hpp), or a power factor whose lifted regions miss
// the goal 32 more bits, and it is advanced again. PreparedPiecewiseState::advance then adds the preparation's
// square-free verdict and, on a complete cover of real input, realness by reflection (prove_real_by_reflection).
//
// State: PreparedPiecewiseState::Impl keeps preparation, direct, states, answers, cursor, separations, attempted,
// statistics (the direct state's, or the factors' summed by PiecewiseRootStatistics::operator+= in piecewise_statistics.cpp).
// Entry points: PreparedPiecewiseState (piecewise_roots.hpp), built by session.cpp (forced piecewise, the selector's
// piecewise route) and real_only_route.cpp (the axis attempt, config.axis_only, which stops at once on a factored
// input); solve_prepared_piecewise_roots and prepared_structure are called only from tests/ and benchmarks/.
// Not here: PiecewiseRootState's charts, certificates, regime exit and axis mode, and machine_newton_cover
// (piecewise_roots.cpp); the preparation (exact_preprocessing.cpp, power_composition.cpp); the cover algebra
// (certification.cpp); the selector (auto_selection.cpp); the hand-off to secular (session.cpp).
// Design notes: docs/piecewise-preprocessing-progress-findings-20260920.md (factors, the exclusion test),
// docs/preparation-handoff-findings-20260925.md (power lifting), docs/work-queue-20260930.md (item 3.2: a complete
// cover at every stop), docs/c1-axis-charts-note-20261015.md (why a factored input ends the axis attempt).
#include "certify/power_cover.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include "core/heuristics.hpp"
#include "numeric/exact/exact_preprocessing.hpp"
#include "certify/factor_cover.hpp"
#include "numeric/exact/power_composition.hpp"
#include "certify/separation.hpp"
#include "numeric/evaluate/structured_source.hpp"
#include "core/telemetry.hpp"
#include "core/statistics_accumulation.hpp"
#include <algorithm>
#include <stdexcept>

namespace certpolysolve::detail {
namespace {
cps_result small_squarefree_factor(const DensePolynomial &input, const cps_options &options,
                                  PiecewiseRootStatistics &statistics) {
    StageTimer timer(statistics.local_solve_ns);
    cps_result result;
    result.degree = static_cast<int64_t>(input.coefficients.size()) - 1;
    result.goal = engine_goal(options);
    const auto requested = std::max<uint64_t>(64, options.bits + 32);
    uint64_t precision = std::min(requested, options.max_work_bits);
    while (precision && !options.cancelled() && result.progress.work_units_this_call < options.work_units) {
        ++result.progress.work_units_this_call;
        result.work_bits = precision;
        ComplexBall a, b, c, discriminant, radical, product, q, alternative;
        input.coefficients[0].enclose(c.value(), static_cast<slong>(precision));
        input.coefficients[1].enclose(b.value(), static_cast<slong>(precision));
        std::vector<Root> roots(static_cast<size_t>(result.degree));
        if (result.degree == 1) {
            acb_div(roots[0].box.value(), c.value(), b.value(), static_cast<slong>(precision));
            acb_neg(roots[0].box.value(), roots[0].box.value());
        } else {
            input.coefficients[2].enclose(a.value(), static_cast<slong>(precision));
            acb_mul(discriminant.value(), b.value(), b.value(), static_cast<slong>(precision));
            acb_mul(product.value(), a.value(), c.value(), static_cast<slong>(precision));
            acb_mul_2exp_si(product.value(), product.value(), 2);
            acb_sub(discriminant.value(), discriminant.value(), product.value(), static_cast<slong>(precision));
            if (arb_is_negative(acb_realref(discriminant.value()))) {
                // Either square-root branch gives the same pair of roots.
                // Rotate a left-half-plane discriminant away from the cut so
                // a tiny imaginary uncertainty does not enclose both signs.
                acb_neg(radical.value(), discriminant.value());
                acb_sqrt(radical.value(), radical.value(), static_cast<slong>(precision));
                acb_mul_onei(radical.value(), radical.value());
            } else {
                acb_sqrt(radical.value(), discriminant.value(), static_cast<slong>(precision));
            }
            // Choose a safely large q = -(b ± sqrt(D))/2. Vieta's product
            // gives the other root as c/q, avoiding catastrophic subtraction
            // when the two root magnitudes differ greatly. All values enclose
            // their exact expressions; an unresolved divisor requests a retry.
            acb_add(q.value(), b.value(), radical.value(), static_cast<slong>(precision));
            acb_sub(alternative.value(), b.value(), radical.value(), static_cast<slong>(precision));
            Float q_lower, alternative_lower;
            acb_get_abs_lbound_arf(q_lower.value(), q.value(), static_cast<slong>(precision));
            acb_get_abs_lbound_arf(alternative_lower.value(), alternative.value(), static_cast<slong>(precision));
            if (arf_cmp(alternative_lower.value(), q_lower.value()) > 0)
                acb_swap(q.value(), alternative.value());
            acb_neg(q.value(), q.value());
            acb_mul_2exp_si(q.value(), q.value(), -1);
            acb_div(roots[0].box.value(), q.value(), a.value(), static_cast<slong>(precision));
            acb_div(roots[1].box.value(), c.value(), q.value(), static_cast<slong>(precision));
        }
        bool valid = true;
        for (auto &root : roots) {
            if (!acb_is_finite(root.box.value())) { valid = false; break; }
            root.count = root.multiplicity = 1;
            root.certificate = "outward linear/quadratic formula on exact square-free factor";
            if (arb_is_zero(acb_imagref(root.box.value()))) root.realness = CPS_REAL;
            else if (!arb_contains_zero(acb_imagref(root.box.value()))) root.realness = CPS_NONREAL;
            root.publish(static_cast<slong>(precision));
            valid = valid && meets_goal(root, options, static_cast<slong>(precision));
        }
        if (valid && disjoint(roots)) {
            result.roots = std::move(roots);
            result.complete = result.goal_met = true;
            result.termination = CPS_SUCCESS;
            return result;
        }
        if (precision == options.max_work_bits) break;
        precision = precision > options.max_work_bits / 2 ? options.max_work_bits : 2 * precision;
    }
    result.termination = options.cancelled() ? CPS_CANCELLED :
        result.progress.work_units_this_call >= options.work_units ? CPS_WORK_LIMIT : CPS_PRECISION_LIMIT;
    return result;
}
} // namespace

struct PreparedPiecewiseState::Impl {
    DensePolynomial input;
    PiecewiseRootConfig config;
    SquarefreePreparation preparation;
    std::unique_ptr<PiecewiseRootState> direct;
    std::vector<std::unique_ptr<PiecewiseRootState>> states;
    std::vector<cps_result> answers;
    PiecewiseRootStatistics statistics, leaf_statistics;
    cps_result result;
    uint64_t preparation_ns = 0, total_work = 0;
    // Each factor's separation bits (separation.hpp: its goal is the caller's with these added), kept across calls.
    std::vector<uint64_t> separations;
    size_t cursor = 0;
    bool prepared = false;
    // The certificate attempt before piecewise (config.certificate_attempt):
    // made at most once; attempt_cover says that `result` is its cover; its
    // observations are kept for the results piecewise returns after a refusal.
    bool attempted = false, attempt_cover = false;
    certpolysolve::detail::CandidateStatistics attempt_statistics{};

    Impl(DensePolynomial p, PiecewiseRootConfig c) : input(std::move(p)), config(std::move(c)) {
        while (!input.coefficients.empty() && input.coefficients.back().zero()) input.coefficients.pop_back();
        result.degree = static_cast<int64_t>(input.coefficients.size()) - 1;
    }
    bool prepare(const cps_options &options) {
        SquarefreePreparation next;
        {
            StageTimer timer(preparation_ns);
            next = prepare_squarefree_input(input, options, config.power_composition);
        }
        if (next.cancelled || options.cancelled()) return false;
        // A modular/degree proof establishes square-freeness of the nonzero
        // part. An elementary part then needs only the same certified formulas
        // used for elementary factors, not a collection of local charts.
        const size_t nonzero_length = input.coefficients.size() - static_cast<size_t>(next.zero_order);
        if (next.nonzero_squarefree && nonzero_length >= 2 && nonzero_length <= 3) {
            ExactInputFactor factor;
            factor.polynomial.coefficients.assign(input.coefficients.begin() +
                static_cast<ptrdiff_t>(next.zero_order), input.coefficients.end());
            next.factors.push_back(std::move(factor));
            next.nonzero_squarefree = false;
        }
        if (config.regime_exit) {
            // A global solve's minimum budget on the whole input, in complex
            // coefficient operations (the selector's unit, auto_selection.cpp): an exact count, saturating.
            const auto degree = static_cast<uint64_t>(std::max<int64_t>(result.degree, 0));
            const auto support = static_cast<uint64_t>(std::count_if(input.coefficients.begin(),
                input.coefficients.end(), [](const Coefficient &c) { return !c.zero(); }));
            config.regime_exit_units = resources::multiply(
                resources::multiply(heuristics::initial_sweep_allowance, degree), resources::add(degree, support));
        }
        std::vector<std::unique_ptr<PiecewiseRootState>> next_states;
        std::unique_ptr<PiecewiseRootState> next_direct;
        if (next.nonzero_squarefree || !next.available()) {
            next_direct = std::make_unique<PiecewiseRootState>(input, config);
        } else {
            for (const auto &factor : next.factors) {
                if (factor.polynomial.coefficients.size() <= 3) next_states.push_back(nullptr);
                else {
                    auto factor_config = config;
                    factor_config.diagnostic_source_id = next_states.size() + 1;
                    factor_config.diagnostic_multiplicity = factor.multiplicity;
                    next_states.push_back(std::make_unique<PiecewiseRootState>(factor.polynomial,
                                                                             std::move(factor_config)));
                }
            }
        }
        std::vector<cps_result> next_answers(next_states.size());
        preparation = std::move(next);
        states = std::move(next_states);
        answers = std::move(next_answers);
        direct = std::move(next_direct);
        prepared = true;
        return true;
    }
    // The certificate attempt the global recipes begin with
    // (machine_newton_cover: a binary64 Aberth packet, then a certificate at
    // every point, with further packets on the points it refuses), on the
    // squarefree nonzero part of the input, before any piecewise model is
    // built. For an input it finishes it costs one packet and two evaluation
    // passes; the selector asks for it when that is well below piecewise's
    // own estimate (auto_selection.cpp). Only the binary64 packet: outside
    // its range the scaled packet costs more than piecewise does, so such an
    // input is refused at once. A refusal leaves nothing behind but its
    // statistics.
    // History: docs/code-history-20261021.md, prepared_piecewise.cpp 1.
    //
    // On success `result` is the complete cover: the attempt's regions (each
    // certified from the original coefficients to hold one root, n of them
    // pairwise disjoint), the exact zero root, and for real coefficients the
    // realness proof the global driver gives the same cover
    // (prove_real_by_reflection).
    bool certificate_attempt(const cps_options &options) {
        attempted = true;
        const auto zeros = static_cast<size_t>(preparation.zero_order);
        if (input.coefficients.size() < zeros + 2)
            return false;
        const auto first = input.coefficients.begin() + static_cast<std::ptrdiff_t>(zeros);
        const slong degree = static_cast<slong>(input.coefficients.size() - zeros) - 1;
        const slong support = static_cast<slong>(
            std::count_if(first, input.coefficients.end(), [](const Coefficient &c) { return !c.zero(); }));
        const bool real = std::all_of(first, input.coefficients.end(), [](const Coefficient &c) { return c.is_real(); });
        // The source the secular recipe gives its attempt: the exact
        // coefficients enclosed at the precision asked for, evaluated term by
        // term when at most a quarter of them are nonzero.
        PolynomialSourceOptions source_options;
        source_options.sparse_coefficients = support > 2 && support * 4 <= degree;
        source_options.input = InputSemantics::Exact;
        PolynomialSource source(
            [this, zeros](acb_poly_t out, slong bits) {
                ComplexBall coefficient;
                for (size_t j = zeros; j < input.coefficients.size(); ++j) {
                    input.coefficients[j].enclose(coefficient.value(), bits);
                    acb_poly_set_coeff_acb(out, static_cast<slong>(j - zeros), coefficient.value());
                }
            },
            degree, real, source_options);
        // The global recipes' first precision (SolveState: goal bits + 32, at least 64).
        const slong precision =
            static_cast<slong>(std::min(options.max_work_bits, std::max<uint64_t>(64, options.bits + 32)));
        std::vector<Root> cover;
        slong certified_at = precision;
        ++statistics.certificate_attempts;
        // The route a failed attempt leaves the solve to is piecewise: the
        // attempt may spend what the selector predicts it costs.
        const bool certified = machine_newton_cover(cover, source, nullptr, options, precision, true, attempt_statistics,
                                                    true, nullptr, nullptr, nullptr, &certified_at, false, nullptr, {},
                                                    config.certificate_attempt_fallback);
        if (!certified || options.cancelled())
            return false;
        const slong bits = std::max(precision, certified_at);
        if (zeros) cover.push_back(exact_zero_point(preparation.zero_order));
        if (!disjoint(cover))
            return false;
        if (real)
            prove_real_by_reflection(cover, bits);
        if (!std::all_of(cover.begin(), cover.end(), [&](const Root &root) { return meets_goal(root, options, bits); }))
            return false;
        ++statistics.certificate_attempt_finishes;
        cps_result finished;
        finished.degree = result.degree;
        finished.goal = engine_goal(options);
        finished.roots = std::move(cover);
        finished.complete = finished.goal_met = true;
        finished.termination = CPS_SUCCESS;
        finished.progress.stop_reason = CPS_STOP_NONE;
        finished.work_bits = static_cast<uint64_t>(bits);
        finished.progress.max_outer_bits = static_cast<uint64_t>(bits);
        finished.components = {{0, static_cast<uint64_t>(bits), 0, static_cast<uint64_t>(result.degree),
                                CPS_COMPONENT_FINISHED}};
        finished.statistics_v3 = attempt_statistics;
        finished.statistics.preprocessing_ns = preparation_ns;
        finished.progress.preprocessing_limited = !preparation.available();
        finished.diagnostics = "certificate attempt before piecewise: complete cover; " + preparation.status;
        result = std::move(finished);
        statistics.preprocessing_ns = preparation_ns;
        statistics.preprocessing_limited = !preparation.available();
        attempt_cover = true;
        return true;
    }
    // Cauchy's global square on the input, whose count is its degree (every coefficient exact), narrowed by
    // complete_cover_at_stop's ladder where its disk tests certify a smaller square: the cover published when a
    // stop finds no complete assembly of the factors (before the preparation, or with a factor unfinished). Never
    // above the caller's precision cap.
    bool publish_input_square(const cps_options &options) {
        const slong bits =
            static_cast<slong>(std::min<uint64_t>(std::max<uint64_t>(64, options.bits), options.max_work_bits));
        ComplexPolynomial polynomial;
        ComplexBall coefficient;
        for (size_t j = 0; j < input.coefficients.size(); ++j) {
            input.coefficients[j].enclose(coefficient.value(), bits);
            acb_poly_set_coeff_acb(polynomial.value(), static_cast<slong>(j), coefficient.value());
        }
        std::vector<Root> cover;
        if (!complete_cover_at_stop(cover, polynomial.value(), bits, &options) || cover.empty()) return false;
        result.roots = std::move(cover);
        result.complete = true;
        result.work_bits = std::max<uint64_t>(result.work_bits, static_cast<uint64_t>(bits));
        return true;
    }
    cps_result advance(const cps_options &options) {
        if (options.real_only || (engine_goal(options) != CPS_ROOT_APPROXIMATION && engine_goal(options) != CPS_ROOT_ACCURACY) ||
            !std::all_of(input.coefficients.begin(), input.coefficients.end(),
                         [](const Coefficient &c) { return c.exact(); })) {
            cps_result unsupported;
            unsupported.degree = result.degree;
            unsupported.goal = engine_goal(options);
            unsupported.termination = CPS_UNSUPPORTED_QUERY;
            return unsupported;
        }
        result.goal = engine_goal(options);
        result.progress.work_units_this_call = 0;
        if (prepared && (!options.work_units || options.cancelled())) {
            // Checking a saved cover is allowed, but a tighter query with no
            // work allowance or a cancelled request must retain its power images.
            result.goal_met = result.complete && std::all_of(result.roots.begin(), result.roots.end(),
                [&](const Root &root) { return meets_goal(root, options,
                    static_cast<slong>(std::max<uint64_t>(64, result.work_bits))); });
            result.termination = result.goal_met ? CPS_SUCCESS : options.cancelled() ? CPS_CANCELLED : CPS_WORK_LIMIT;
            result.progress.stop_reason = result.goal_met ? CPS_STOP_NONE :
                result.termination == CPS_CANCELLED ? CPS_STOP_CANCELLED : CPS_STOP_WORK_BUDGET;
            return result;
        }
        uint64_t work = 0;
        if (!prepared) {
            if (!options.cancelled() && options.work_units) {
                ++work;
                prepare(options);
            }
            if (!prepared) {
                total_work = resources::add(total_work, work);
                result.termination = options.cancelled() ? CPS_CANCELLED : CPS_WORK_LIMIT;
                result.progress.stop_reason = options.cancelled() ? CPS_STOP_CANCELLED : CPS_STOP_WORK_BUDGET;
                result.progress.work_units_this_call = work;
                result.statistics_v2.work_units = total_work;
                result.statistics.preprocessing_ns = preparation_ns;
                result.diagnostics = "piecewise preparation pending";
                // The stop's cover: the input's square (publish_input_square). Not with zero work, which returns
                // a snapshot without mathematical progress (docs/api.md). History:
                // docs/code-history-20261021.md, prepared_piecewise.cpp 2.
                if (!options.cancelled() && options.work_units && publish_input_square(options))
                    result.diagnostics += " with a complete cover";
                return result;
            }
        }
        if (config.axis_only && !direct) {
            // The axis charts answer an input whose nonzero part is square-free (the direct state). The factor loop
            // below counts a factor's answer complete when it covers all of the factor's roots, which an axis answer
            // (the real ones) never does: the loop would raise the goal to the precision cap, refining every real
            // root at each step, and stop. An input with repeated or power factors ends the attempt at once, a stop
            // of the attempt itself. History: docs/code-history-20261021.md, prepared_piecewise.cpp 3.
            total_work = resources::add(total_work, work);
            result.termination = CPS_PRECISION_LIMIT;
            result.progress.stop_reason = CPS_STOP_PRECISION;
            result.progress.work_units_this_call = work;
            result.statistics_v2.work_units = total_work;
            result.statistics.preprocessing_ns = preparation_ns;
            result.diagnostics = "the axis charts take a square-free input; " + preparation.status;
            return result;
        }
        if (attempt_cover) {
            // The attempt's cover answers a later call with the same goal; a
            // tighter one is piecewise's.
            result.goal = engine_goal(options);
            result.goal_met = std::all_of(result.roots.begin(), result.roots.end(), [&](const Root &root) {
                return meets_goal(root, options, static_cast<slong>(std::max<uint64_t>(64, result.work_bits)));
            });
            if (result.goal_met) {
                result.termination = CPS_SUCCESS;
                result.progress.stop_reason = CPS_STOP_NONE;
                result.progress.work_units_this_call = work;
                total_work = resources::add(total_work, work);
                result.statistics_v2.work_units = total_work;
                return result;
            }
            attempt_cover = false;
        }
        if (direct && config.certificate_attempt && !attempted && work < options.work_units && !options.cancelled()) {
            ++work;
            if (certificate_attempt(options)) {
                total_work = resources::add(total_work, work);
                result.progress.work_units_this_call = work;
                result.statistics_v2.work_units = total_work;
                return result;
            }
        }
        cps_options query = options;
        query.work_units -= work;
        if (direct) {
            result = direct->advance(query);
            if (attempted)
                // Observed attempt work belongs to the solve whatever finishes it.
                statistics::add_ordinary(result.statistics_v3, attempt_statistics, statistics::candidate_fields);
            result.progress.work_units_this_call += work;
            total_work = resources::add(total_work, result.progress.work_units_this_call);
            result.statistics_v2.work_units = total_work;
            statistics = direct->statistics();
            // The attempt is this state's, not the direct state's.
            statistics.certificate_attempts = attempted ? 1 : 0;
            statistics.preprocessing_ns = preparation_ns;
            statistics.preprocessing_limited = !preparation.available();
            result.statistics.preprocessing_ns = preparation_ns;
            result.progress.preprocessing_limited = statistics.preprocessing_limited != 0;
            result.diagnostics += "; " + preparation.status;
            return result;
        }
        if (separations.size() != states.size()) separations.assign(states.size(), 0);
        bool exited = false;
        for (;;) {
            bool all_complete = true;
            cps_termination termination = CPS_PRECISION_LIMIT;
            const size_t first = cursor;
            for (size_t step = 0; step < states.size(); ++step) {
                const size_t i = (first + step) % states.size();
                auto &answer = answers[i];
                auto goal = separated_options(options, separations[i]);
                answer.goal_met = answer.complete && std::all_of(answer.roots.begin(), answer.roots.end(),
                    [&](const Root &root) { return meets_goal(root, goal, static_cast<slong>(
                        std::max<uint64_t>(64, answer.work_bits))); });
                if (!answer.goal_met && !exited && work < options.work_units && !options.cancelled()) {
                    goal.work_units = options.work_units - work;
                    answer = states[i] ? states[i]->advance(goal) :
                        small_squarefree_factor(preparation.factors[i].polynomial, goal, leaf_statistics);
                    if (answer.progress.work_units_this_call > goal.work_units)
                        throw std::logic_error("factor exceeded its work allowance");
                    work += answer.progress.work_units_this_call;
                    cursor = (i + 1) % states.size();
                    exited = states[i] && states[i]->statistics().regime_exits;
                }
                result.work_bits = std::max(result.work_bits, answer.work_bits);
                all_complete = all_complete && answer.complete && answer.goal_met;
                if (!answer.goal_met) termination = answer.termination;
            }
            std::vector<Root> zero_point;
            std::vector<std::vector<Root>> mapped(answers.size());
            std::vector<bool> mapped_complete(answers.size(), false);
            if (preparation.zero_order) zero_point.push_back(exact_zero_point(preparation.zero_order));
            // A factor's goal reaching the cap is a stop of the loop.
            const bool at_cap = std::any_of(separations.begin(), separations.end(), [&](uint64_t separation) {
                return separated_options(options, separation).bits >= options.max_work_bits;
            });
            // A stop short of the goal (not a cancellation) returns a complete
            // cover: each factor's partial cover is completed on the factor
            // (complete_cover_at_stop), as the factor's own state does for a
            // call that reaches it; here it also covers the factors this call
            // never reached.
            const bool stopping_short = !options.cancelled() &&
                (work >= options.work_units || at_cap || exited);
            std::vector<std::vector<Root>> completed(answers.size());
            for (size_t i = 0; i < answers.size(); ++i) {
                const auto &factor = preparation.factors[i];
                const auto multiplicity = factor.multiplicity;
                std::vector<Root> lifted;
                const auto *cover = &answers[i].roots;
                bool complete_here = answers[i].complete;
                if (!answers[i].complete && stopping_short) {
                    StageTimer timer(leaf_statistics.certification_ns);
                    const slong bits = static_cast<slong>(std::max<uint64_t>(64, answers[i].work_bits ? answers[i].work_bits :
                        separated_options(options, separations[i]).bits));
                    ComplexPolynomial polynomial;
                    ComplexBall coefficient;
                    for (size_t j = 0; j < factor.polynomial.coefficients.size(); ++j) {
                        factor.polynomial.coefficients[j].enclose(coefficient.value(), bits);
                        acb_poly_set_coeff_acb(polynomial.value(), static_cast<slong>(j), coefficient.value());
                    }
                    completed[i] = answers[i].roots;
                    if (complete_cover_at_stop(completed[i], polynomial.value(), bits, &options)) {
                        uint64_t counted = 0;
                        for (const auto &root : completed[i]) counted += root.count;
                        if (counted == static_cast<uint64_t>(acb_poly_degree(polynomial.value()))) {
                            cover = &completed[i];
                            complete_here = true;
                        }
                    }
                }
                if (factor.power > 1) {
                    StageTimer timer(leaf_statistics.certification_ns);
                    // Exact reduced roots can satisfy every later query without
                    // another child iteration. Lifting still needs arithmetic
                    // appropriate to the new query, independently of that child.
                    const auto lifting_bits = std::min(options.max_work_bits,
                        std::max(result.work_bits, resources::add(separated_options(options, separations[i]).bits, 32)));
                    // Inflated images of a partial cover could contain other,
                    // unaccounted roots of this same factor. Only lift complete
                    // factor covers; their full set of branches accounts for
                    // every root, and the overlap ledger below checks the boxes.
                    if (!complete_here || options.cancelled() ||
                        !lift_power_cover(lifted, *cover, factor.power,
                            static_cast<slong>(lifting_bits))) continue;
                    result.work_bits = std::max(result.work_bits, lifting_bits);
                    cover = &lifted;
                }
                mapped_complete[i] = complete_here;
                for (auto root : *cover) {
                    root.count *= multiplicity;
                    if (root.multiplicity) root.multiplicity *= multiplicity;
                    root.certificate += "; exact square-free factor, exponent " + std::to_string(multiplicity);
                    mapped[i].push_back(std::move(root));
                }
            }
            // The mapped covers of coprime factors are merged where they overlap (merge_overlapping): a region may
            // then hold roots of two factors as one counted cluster, which meets the goal if its size does, as the
            // global driver's aggregate merges (W8's S8, rule 1). A complete assembly (every factor's cover complete
            // and mapped) may be published (C4, below), its counts summing to the degree; one with a factor unfinished (a
            // cancellation, a refused lift, a stop whose completion failed) publishes nothing, so the last complete
            // cover stays the answer, the input's square if no assembly has completed (rule 10).
            if (std::all_of(mapped_complete.begin(), mapped_complete.end(), [](bool value) { return value; })) {
                std::vector<Root> proposed = zero_point;
                for (const auto &cover : mapped) proposed.insert(proposed.end(), cover.begin(), cover.end());
                const auto bits = static_cast<slong>(std::max<uint64_t>(64, result.work_bits));
                assemble_factor_covers(proposed, static_cast<uint64_t>(result.degree), bits);
                // Published as the global driver's aggregate publishes (W8's S8, C4): when it meets the goal, when no
                // complete answer stands yet, or when prefer_cover ranks it ahead.
                const bool attained = std::all_of(proposed.begin(), proposed.end(),
                                                  [&](const Root &root) { return meets_goal(root, options, bits); });
                publish_assembly(result, proposed, attained, bits);
            } else if (!result.complete) {
                publish_input_square(options);
            }
            result.goal_met = result.complete && std::all_of(result.roots.begin(), result.roots.end(),
                [&](const Root &root) { return meets_goal(root, options,
                    static_cast<slong>(std::max<uint64_t>(64, result.work_bits))); });
            result.termination = result.goal_met ? CPS_SUCCESS : options.cancelled() ? CPS_CANCELLED :
                work >= options.work_units ? CPS_WORK_LIMIT : termination;
            if (result.goal_met || !all_complete || options.cancelled() || work >= options.work_units) break;
            // Every factor meets its goal and the mapped cover does not: two factors' covers meet, or a lifted region
            // misses the goal. Rule 2 (W8's S8, the global driver's rule, separation.hpp): only the factor whose
            // region is the wider of the first meeting pair gets a tighter goal (the exact zero point, first in the
            // order, is never the wider); with no two covers meeting, a power factor whose lifted regions miss the
            // goal gets 32 more bits; with neither, nothing more can be done. The tightened factor's persistent
            // solver is advanced again, the others keep their covers.
            std::vector<const std::vector<Root> *> covers;
            if (!zero_point.empty()) covers.push_back(&zero_point);
            for (const auto &cover : mapped) covers.push_back(&cover);
            const auto chosen = factor_to_separate(covers);
            if (chosen) {
                const size_t k = *chosen - (zero_point.empty() ? 0 : 1);
                if (separations[k] >= resources::precision_max - options.bits) break;
                separations[k] = next_separation(options.bits, separations[k]);
                continue;
            }
            const auto missing = [&](size_t k) {
                return preparation.factors[k].power > 1 && std::any_of(mapped[k].begin(), mapped[k].end(),
                    [&](const Root &root) { return !meets_goal(root, options,
                        static_cast<slong>(std::max<uint64_t>(64, result.work_bits))); });
            };
            size_t k = 0;
            while (k < mapped.size() && !missing(k)) ++k;
            if (k == mapped.size() || separations[k] >= resources::precision_max - options.bits) break;
            separations[k] = std::min<uint64_t>(resources::precision_max - options.bits, separations[k] + 32);
        }
        statistics = {};
        for (const auto &state : states) if (state) statistics += state->statistics();
        statistics += leaf_statistics;
        statistics.discovery_count = config.durable_ownership ? preparation.zero_order : 0;
        statistics.final_accurate_count = 0;
        for (size_t i = 0; i < states.size(); ++i) if (config.durable_ownership) {
            uint64_t discovered = states[i] ? states[i]->statistics().discovery_count : 0;
            if (!states[i]) for (const auto &root : answers[i].roots) discovered += root.count;
            statistics.discovery_count += discovered * preparation.factors[i].multiplicity *
                                          preparation.factors[i].power;
        }
        for (const auto &root : result.roots)
            if (meets_goal(root, options, static_cast<slong>(std::max<uint64_t>(64, result.work_bits))))
                statistics.final_accurate_count += root.count;
        statistics.preprocessing_ns = preparation_ns;
        statistics.prepared_factors = states.size();
        result.statistics.preprocessing_ns = preparation_ns;
        result.statistics.exact_factors = states.size();
        total_work = resources::add(total_work, work);
        result.statistics_v2.work_units = total_work;
        result.progress.work_units_this_call = work;
        result.progress.max_outer_bits = result.work_bits;
        result.progress.pending_components = result.goal_met ? 0 : states.size();
        result.progress.stop_reason = result.goal_met ? CPS_STOP_NONE : result.termination == CPS_CANCELLED ?
            CPS_STOP_CANCELLED : result.termination == CPS_WORK_LIMIT ? CPS_STOP_WORK_BUDGET : CPS_STOP_PRECISION;
        result.diagnostics = (result.complete ? "complete" : "partial") +
            std::string(" piecewise cover after ") + preparation.status;
        return result;
    }
};
PreparedPiecewiseState::PreparedPiecewiseState(DensePolynomial input, PiecewiseRootConfig config)
    : impl_(std::make_unique<Impl>(std::move(input), std::move(config))) {}
PreparedPiecewiseState::~PreparedPiecewiseState() = default;
cps_result PreparedPiecewiseState::advance(const cps_options &options) {
    auto result = impl_->advance(options);
    // The square-free preparation's verdict (C14, cps_result_squarefree): an exact zero root of order two or more,
    // or a factor of exact multiplicity above one, proves a repeated root; the nonzero part proved square-free (by a
    // modular test or as the factors' product, each in x^power with a nonzero constant) with at most a simple zero
    // root proves every root simple; an unavailable preparation proves nothing.
    if (impl_->prepared) {
        const auto &p = impl_->preparation;
        const bool repeated = p.zero_order > 1 || std::any_of(p.factors.begin(), p.factors.end(),
                                                              [](const ExactInputFactor &f) { return f.multiplicity > 1; });
        result.squarefree = !p.available() ? CPS_SQUAREFREE_UNKNOWN
                            : repeated     ? CPS_SQUAREFREE_DISPROVED
                                           : CPS_SQUAREFREE_PROVED;
    }
    // A complete cover of a real polynomial proves the realness of every
    // region that holds one distinct root and meets no other region's
    // reflection (prove_real_by_reflection), as the global driver's
    // aggregate() and the certificate attempt's assembly do.
    // History: docs/code-history-20261021.md, prepared_piecewise.cpp 4.
    if (result.complete && std::all_of(impl_->input.coefficients.begin(), impl_->input.coefficients.end(),
                                       [](const Coefficient &c) { return c.is_real(); }))
        prove_real_by_reflection(result.roots, static_cast<slong>(std::max<uint64_t>(64, result.work_bits)));
    return result;
}
const PiecewiseRootStatistics &PreparedPiecewiseState::statistics() const { return impl_->statistics; }
std::optional<PreparedStructure> PreparedPiecewiseState::prepared_structure() const {
    if (!impl_->prepared || impl_->result.degree < 0) return std::nullopt;
    PreparedStructure out;
    out.input_degree = impl_->result.degree;
    if (impl_->direct) {
        out.append({input_structure(impl_->input.coefficients), 1, 1,
                    impl_->preparation.nonzero_squarefree, PreparedRepresentation::Compressed});
    } else {
        out.zero_order = impl_->preparation.zero_order;
        for (const auto &factor : impl_->preparation.factors)
            out.append({input_structure(factor.polynomial.coefficients), factor.multiplicity, factor.power,
                        true, PreparedRepresentation::Compressed});
    }
    return out;
}
cps_result solve_prepared_piecewise_roots(const DensePolynomial &input, const cps_options &options,
                                        PiecewiseRootStatistics &statistics, PiecewiseRootConfig config) {
    PreparedPiecewiseState state(input, std::move(config));
    auto result = state.advance(options);
    statistics = state.statistics();
    return result;
}
} // namespace certpolysolve::detail
