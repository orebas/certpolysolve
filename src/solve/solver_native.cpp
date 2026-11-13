#include "solve/local_model_solver.hpp"
// SolveState::Impl's native secular route (native_advance; the class: solve_state_impl.hpp): SecularRefinement on a
// SecularSource seeded at the poles, falling back to the coefficients (secular_polynomial) at its precision stop.
#include "certify/factor_cover.hpp"
#include "solve/solve_state_impl.hpp"
#include "methods/secular/secular_refinement.hpp"
#include "numeric/evaluate/structured_source.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

namespace certpolysolve::detail {

bool SolveState::Impl::native_advance(const cps_options &o) {
    auto update = [&] {
        if (!native_secular)
            return;
        result.work_bits = native_secular->precision();
        result.progress.max_outer_bits = std::max(result.progress.max_outer_bits, result.work_bits);
        std::vector<Root> cover;
        if (native_secular->cover(cover)) {
            const auto prec = std::max<uint64_t>(64, result.work_bits);
            publish_assembly(result, cover, goal(cover, prec), static_cast<slong>(prec));
            result.goal_met = goal(result.roots, prec);
        }
    };
    if (result.complete)
        result.goal_met = goal(result.roots, std::max<uint64_t>(64, result.work_bits));
    bool iteration_stop = false;
    while (!result.goal_met && !o.cancelled() && result.progress.work_units_this_call < o.work_units) {
        if (!o.sweeps_left(result.statistics_v3)) {
            iteration_stop = true;
            break;
        }
        ++result.statistics_v2.work_units;
        ++result.progress.work_units_this_call;
        if (!native_secular) {
            TraceScope trace(&policy, "native_secular_prepare", initial_precision,
                             static_cast<uint64_t>(result.degree), &result.statistics_v3,
                             &result.statistics_v2);
            auto source = std::make_unique<SecularSource>(*secular_input);
            slong seed_precision = static_cast<slong>(initial_precision);
            for (;;) {
                source->enclose_terms(seed_precision);
                bool separate = true;
                for (size_t i = 0; i < source->poles().size() && separate; ++i)
                    for (size_t j = 0; j < i; ++j)
                        if (acb_overlaps(source->poles()[i].value(), source->poles()[j].value())) {
                            separate = false;
                            break;
                        }
                if (separate || seed_precision == static_cast<slong>(o.max_work_bits) || o.cancelled())
                    break;
                seed_precision = std::min<slong>(static_cast<slong>(o.max_work_bits), resources::precision_product(seed_precision, 2));
            }
            RootVector seeds(result.degree);
            RealBall perturbation;
            ComplexBall direction;
            for (slong i = 0; i < result.degree; ++i) {
                acb_get_mid(seeds.data() + i, source->poles()[static_cast<size_t>(i)].value());
                if (result.degree > 1) {
                    // Small, differently oriented pole perturbations break
                    // real-axis and lattice symmetries without changing the
                    // represented equation. Regeneration checks every node.
                    acb_abs(perturbation.value(), seeds.data() + i, 64);
                    arb_add_ui(perturbation.value(), perturbation.value(), 1, 64);
                    arb_mul_2exp_si(perturbation.value(), perturbation.value(),
                                    -std::min<slong>(seed_precision - 12, 48));
                    ComplexBall difference;
                    Float separation;
                    for (slong j = 0; j < result.degree; ++j)
                        if (i != j) {
                            acb_sub(difference.value(), source->poles()[static_cast<size_t>(i)].value(),
                                    source->poles()[static_cast<size_t>(j)].value(), seed_precision);
                            acb_get_abs_lbound_arf(separation.value(), difference.value(), 64);
                            arf_mul_2exp_si(separation.value(), separation.value(), -4);
                            if (arf_cmp(separation.value(), arb_midref(perturbation.value())) < 0)
                                arb_set_arf(perturbation.value(), separation.value());
                        }
                    const double angle = 2.39996322972865332 * static_cast<double>(i + 1);
                    acb_set_d_d(direction.value(), std::cos(angle), std::sin(angle));
                    acb_mul_arb(direction.value(), direction.value(), perturbation.value(), seed_precision);
                    acb_add(seeds.data() + i, seeds.data() + i, direction.value(), seed_precision);
                }
                acb_get_mid(seeds.data() + i, seeds.data() + i);
            }
            RootVector block_seeds(result.degree);
            if (secular_block_seeds(block_seeds.data(), *source, seed_precision, o))
                _acb_vec_set(seeds.data(), block_seeds.data(), result.degree);
            if (o.cancelled())
                break;
            native_secular = std::make_unique<SecularRefinement>(std::move(source), seeds.data(), 1,
                                                                 seed_precision, policy, solve_local_model);
        }
        const auto reason = native_secular->step(o, result.statistics_v2, result.statistics_v3);
        update();
        if (reason == CPS_STOP_CANCELLED)
            break;
        if (reason == CPS_STOP_ITERATION_BUDGET) {
            iteration_stop = !result.goal_met;
            break;
        }
        if (!result.goal_met &&
            (reason == CPS_STOP_PRECISION ||
             (result.degree <= heuristics::native_secular_fallback_degree &&
              native_secular->rounds() >= heuristics::native_secular_fallback_rounds &&
              native_secular->precision() >= heuristics::native_secular_fallback_bits))) {
            native_fallback = true;
            return false;
        }
    }
    // A met goal is a success whatever the last checkpoint polls (I2; C25).
    const bool cancelled = !result.goal_met && o.cancelled();
    result.progress.stop_reason = result.goal_met ? CPS_STOP_NONE
                                  : cancelled     ? CPS_STOP_CANCELLED
                                  : iteration_stop ? CPS_STOP_ITERATION_BUDGET
                                                   : CPS_STOP_WORK_BUDGET;
    result.termination = result.goal_met ? CPS_SUCCESS : cancelled ? CPS_CANCELLED : CPS_WORK_LIMIT;
    result.progress.pending_components = result.goal_met  ? 0
                                         : native_secular ? native_secular->active()
                                                          : 1;
    result.progress.min_pending_bits = result.progress.max_pending_bits =
        result.goal_met ? 0 : result.work_bits;
    result.components.clear();
    for (size_t i = 0; i < result.roots.size(); ++i)
        result.components.push_back(
            {0, result.work_bits, 0, result.roots[i].count,
             meets_goal(result.roots[i], o, static_cast<slong>(std::max<uint64_t>(64, result.work_bits)))
                 ? CPS_COMPONENT_FINISHED
                 : CPS_COMPONENT_ACTIVE});
    result.diagnostics = "native secular representation; checked interpolation inclusion cover; "
                         "coefficient expansion deferred";
    sort_roots(result.roots);
    return true;
}

} // namespace certpolysolve::detail
