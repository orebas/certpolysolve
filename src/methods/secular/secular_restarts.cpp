// The secular refinement's restarts: the seeds of a component (component_seeds, place_undetermined_hints by saddle
// points), the groups of nodes (node_groups), a group's local model (group_model, locate_mean, dominant_count), and
// the restarts themselves: local_component_restart, regroup_restart, geometric_restart, certificate_reseed. Part of
// the secular refinement (secular_refinement.cpp), split from it verbatim in W8's S5 (docs/w8-plan-20261021.md);
// state: the SecularRefinement members these members change. Entry points: SecularRefinement::step.
#include "methods/secular/secular_refinement.hpp"
#include "methods/secular/secular_internal.hpp"
#include "methods/secular/secular_geometry.hpp"
#include "numeric/evaluate/disk_blocked_polynomial_evaluator.hpp"
#include "numeric/evaluate/disk_sparse_polynomial_evaluator.hpp"
#include "core/dovetail.hpp"
#include "core/service_meter.hpp"
#include "numeric/evaluate/limb_horner_evaluator.hpp"
#include "numeric/evaluate/node_products.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <bit>
#include <span>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>

namespace certpolysolve::detail {
// The hints of a local model that the model does not determine, placed by
// their scale alone. The model is the component's Taylor polynomial truncated
// at its count, with coefficient balls, so its value at a point is known up to
// those balls and the dropped terms (`dropped`: a bound on them at the point,
// from their magnitudes, at most sum |a_j| |x|^j, or from the source's value
// there; component_seeds has both). Two hints are told apart by the model
// only if that enclosure excludes zero at their midpoint: otherwise an
// admissible polynomial has a root there, and the model does not say how the
// two split.
//
// It matters because of what the iteration then does with them. For two real
// roots closer than the model resolves, the truncated real model has a
// complex-conjugate pair, and a Jacobi sweep of a real source keeps a
// conjugate pair conjugate: the pair settles on the bisector of the two roots
// and never reaches them.
//
// The midpoint is where two hints with no other root near them would meet. With
// a third root near, they meet elsewhere: two simple roots of the model merge
// under a perturbation exactly at a critical point of the model, where the
// perturbation is the model's value there with its sign changed. So the pair is
// also tested at the model's critical point between them (saddle_between): the
// smallest perturbation that makes them one double root is admissible if the
// enclosure contains zero there. History:
// docs/code-history-20261021.md, secular_refinement.cpp 4.
//
// Each group of mutually undetermined hints goes on the circle of its own
// spread about its centroid at the placement angles (candidates.hpp), as an
// undetermined model's hints go on its polygon's circles: distinct points at
// the scale the evidence gives, with no symmetry the roots need not have.
// Returns the number of hints placed.
namespace {
// Newton's steps for a critical point, and the bits of the two hints' half-distance it is located to. The
// model's value is stationary at the critical point, so an error of 2^-s of the half-distance changes it by
// about 2^-2s of its size: 2^-16 here, far inside what the dropped terms' bound is known to. Newton's
// iteration doubles the bits from the midpoint, which has one where the critical point lies between the hints
// at all: one bit to eight in three steps, and one more for a start that has none.
inline constexpr unsigned saddle_steps = 4;
inline constexpr slong saddle_bits = 8;
// The model's critical point between two of its roots: Newton's iteration on its derivative from their
// midpoint `point`, which it replaces. False, with `point` unchanged, if the iteration leaves the disk on the
// two roots as a diameter or does not settle within its steps: then no critical point is known between them.
bool saddle_between(acb_t point, const acb_poly_t first, const acb_t one, const acb_t other, slong precision) {
    ComplexBall middle, at, slope, curvature, step, offset;
    Magnitude half, moved, reach;
    acb_set(middle.value(), point);
    acb_sub(offset.value(), one, other, precision);
    acb_get_mag(half.value(), offset.value());
    mag_mul_2exp_si(half.value(), half.value(), -1);
    acb_set(at.value(), point);
    for (unsigned k = 0; k < saddle_steps; ++k) {
        acb_poly_evaluate2(slope.value(), curvature.value(), first, at.value(), precision);
        if (acb_contains_zero(curvature.value()))
            return false;
        acb_div(step.value(), slope.value(), curvature.value(), precision);
        acb_sub(at.value(), at.value(), step.value(), precision);
        acb_get_mid(at.value(), at.value());
        if (!acb_is_finite(at.value()))
            return false;
        acb_sub(offset.value(), at.value(), middle.value(), precision);
        acb_get_mag(reach.value(), offset.value());
        if (mag_cmp(reach.value(), half.value()) >= 0)
            return false;
        acb_get_mag(moved.value(), step.value());
        mag_mul_2exp_si(moved.value(), moved.value(), saddle_bits);
        if (mag_cmp(moved.value(), half.value()) <= 0) {
            acb_set(point, at.value());
            return true;
        }
    }
    return false;
}
} // namespace
slong place_undetermined_hints(acb_ptr hints, slong count, const acb_poly_t model,
                               const std::function<void(mag_t, const acb_t)> &dropped, slong precision) {
    if (count < 2)
        return 0;
    ComplexPolynomial first;
    acb_poly_derivative(first.value(), model, precision);
    std::vector<slong> group(static_cast<size_t>(count));
    std::iota(group.begin(), group.end(), slong{0});
    const auto find = [&](slong i) {
        while (group[static_cast<size_t>(i)] != i)
            i = group[static_cast<size_t>(i)] = group[static_cast<size_t>(group[static_cast<size_t>(i)])];
        return i;
    };
    ComplexBall difference, middle, value;
    Magnitude distance, nearest, tau;
    bool any = false;
    for (slong i = 0; i < count; ++i) {
        slong neighbour = -1;
        mag_inf(nearest.value());
        for (slong j = 0; j < count; ++j)
            if (j != i) {
                acb_sub(difference.value(), hints + i, hints + j, precision);
                acb_get_mag(distance.value(), difference.value());
                if (mag_cmp(distance.value(), nearest.value()) < 0) {
                    mag_set(nearest.value(), distance.value());
                    neighbour = j;
                }
            }
        if (neighbour < 0 || find(i) == find(neighbour))
            continue;
        acb_add(middle.value(), hints + i, hints + neighbour, precision);
        acb_mul_2exp_si(middle.value(), middle.value(), -1);
        acb_get_mid(middle.value(), middle.value());
        acb_poly_evaluate(value.value(), model, middle.value(), precision);
        dropped(tau.value(), middle.value());
        acb_add_error_mag(value.value(), tau.value());
        if (!acb_contains_zero(value.value()) &&
            saddle_between(middle.value(), first.value(), hints + i, hints + neighbour, precision)) {
            acb_poly_evaluate(value.value(), model, middle.value(), precision);
            dropped(tau.value(), middle.value());
            acb_add_error_mag(value.value(), tau.value());
        }
        if (acb_contains_zero(value.value())) {
            group[static_cast<size_t>(find(i))] = find(neighbour);
            any = true;
        }
    }
    if (!any)
        return 0;
    slong placed_count = 0;
    for (slong representative = 0; representative < count; ++representative) {
        std::vector<slong> members;
        for (slong i = 0; i < count; ++i)
            if (find(i) == representative)
                members.push_back(i);
        if (members.size() < 2)
            continue;
        ComplexBall centre;
        for (const auto i : members)
            acb_add(centre.value(), centre.value(), hints + i, precision);
        acb_div_ui(centre.value(), centre.value(), static_cast<ulong>(members.size()), precision);
        acb_get_mid(centre.value(), centre.value());
        Magnitude spread;
        for (const auto i : members) {
            acb_sub(difference.value(), hints + i, centre.value(), precision);
            acb_get_mag(distance.value(), difference.value());
            mag_max(spread.value(), spread.value(), distance.value());
        }
        // Coincident hints have no scale of their own; the caller's
        // separation of equal proposals moves them.
        if (mag_is_zero(spread.value()))
            continue;
        RealBall radius;
        arf_set_mag(arb_midref(radius.value()), spread.value());
        RootVector placed(static_cast<slong>(members.size()));
        circle_candidates(placed.data(), static_cast<slong>(members.size()), centre.value(), radius.value(), precision);
        for (size_t j = 0; j < members.size(); ++j)
            acb_set(hints + members[j], placed.data() + static_cast<slong>(j));
        placed_count += static_cast<slong>(members.size());
    }
    return placed_count;
}
namespace {
struct LocalPacketMeter {
    ResearchLocalPacketAccounting &owner;
    uint64_t terms = 0, interactions = 0;
    void charge(const ResearchLocalPacketProgress &progress) {
        if (progress.evaluation_terms < terms || progress.root_interactions < interactions ||
            progress.evaluation_terms > UINT64_MAX - owner.prepaid_terms)
            throw std::logic_error("local packet accounting counters moved backwards or overflowed");
        const auto new_terms = progress.evaluation_terms - terms;
        const auto new_interactions = progress.root_interactions - interactions;
        const auto bits = static_cast<int64_t>(progress.precision);
        if (owner.policy && owner.policy->research_service_meter) {
            owner.work += service_charge(owner.policy, ServiceOperation::Evaluation,
                                          static_cast<double>(new_terms) * operation_weight(bits), bits);
            owner.work += service_charge(owner.policy, ServiceOperation::Interaction,
                                          static_cast<double>(new_interactions) * operation_weight(bits), bits);
        } else
            owner.work += (static_cast<double>(new_terms) + static_cast<double>(new_interactions)) * operation_weight(bits);
        terms = progress.evaluation_terms;
        interactions = progress.root_interactions;
    }
    static ResearchServiceAction service(void *pointer, const ResearchLocalPacketProgress &progress) {
        auto &meter = *static_cast<LocalPacketMeter *>(pointer);
        meter.charge(progress);
        return meter.owner.service ? meter.owner.service(meter.owner.data, progress)
                                   : ResearchServiceAction::Continue;
    }
};
} // namespace
bool component_seeds(acb_ptr out, PolynomialSource &source, const Root &component, slong k, slong &working,
                     const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics, const EnginePolicy &policy,
                     const std::optional<ComplexBall> &center_hint, slong resolution, bool *unresolved,
                     slong counted_hint_resolution
                     , ResearchLocalPacketAccounting *accounting
                     ) {
    if (unresolved)
        *unresolved = false;
    // The models one call may compute: the first, and one for each time a deficiency raises the precision.
    constexpr unsigned attempts = 8;
    for (unsigned retry = 0; retry < attempts; ++retry) {
        ComplexBall center, step;
        ComplexPolynomial local;
        acb_get_mid(center.value(), center_hint ? center_hint->value() : component.box.value());
        if (source.direct()) {
            // A cheap exact centre inside the retained component often exposes
            // coefficient cancellation before any high-precision arithmetic.
            // It is only a local-model proposal; the original source certifies.
            ComplexBall coarse;
            constexpr slong local_center_bits = heuristics::local_center_bits;
            arf_set_round(arb_midref(acb_realref(coarse.value())), arb_midref(acb_realref(center.value())),
                          local_center_bits, ARF_RND_NEAR);
            arf_set_round(arb_midref(acb_imagref(coarse.value())), arb_midref(acb_imagref(center.value())),
                          local_center_bits, ARF_RND_NEAR);
            if (source.real() && arb_contains_zero(acb_imagref(component.box.value())))
                arb_zero(acb_imagref(coarse.value()));
            if (acb_contains(component.box.value(), coarse.value()))
                acb_set(center.value(), coarse.value());
        }
        bool finite = true;
        // p^(k-1) locates a component centre more stably than p locates its
        // individual roots. All polynomial tails remain in certification.
        // The step is Newton's on p^(k-1): the expansion's coefficient k - 1
        // over k times its coefficient k, which are that derivative's value
        // and slope at the centre (PolynomialSource::taylor_pair, one joint
        // evaluation, O(n), where a whole shift for each step costs far more).
        // The whole shift is for the model below, whose dropped terms' sizes
        // place_undetermined_hints reads. The expansion's first k + 1
        // coefficients by synthetic division can be cheaper or dearer than a
        // whole shift, and the model itself comes from the expansion's head
        // where that is the cheaper, below.
        //
        // Each step is taken at the precision its own accuracy needs. The model wants its centre within its
        // hints' radius, over 4 k, of the roots' mean: the expansion's coefficient k - 1 is k times that
        // offset and then lies under the model's Newton polygon, where it places nothing. Two things bound the
        // hints' radius. The model's noise: at 2^-N it places its hints 2^(-N / k) out, working / k bits of
        // the centre (noise_bits). And the component's box, which every hint has to lie in: the bits of its
        // smaller radius at the centre's modulus (box_bits). The centre is located to the finer of the two,
        // and the bits of 4 k (mean_bits). A Newton step doubles the centre's bits, so the three steps aim at
        // a quarter, a half and the whole of that, and none at less than the box's bits: the centre starts
        // inside the box, and a step that resolves or keeps less puts it outside again, beyond the reach of
        // Newton's iteration for the roots' mean: the model's hints then land outside the box, and the
        // component is re-seeded on a circle.
        //
        // An exact local model (PolynomialSource::exact_local_models, small Gaussian factors) has no noise:
        // it resolves its roots at every radius, so no radius bounds its hints from below and its steps take
        // the working precision.
        //
        // A step's two values need its aim, the regeneration's guard, and the bits the derivative's
        // evaluation cancels: the precision starts with the first two and doubles until the step's own radius
        // is within the aim at the centre's modulus (a centre at the origin has no such measure: the working
        // precision). The centre keeps the bits of its aim: the model's passes multiply by it. At the working
        // precision the three steps can be a large share of a restart. History:
        // docs/code-history-20261021.md, secular_refinement.cpp 5.
        constexpr unsigned centre_steps = 3;
        const slong mean_bits = static_cast<slong>(FLINT_BIT_COUNT(static_cast<ulong>(k))) + 2;
        const slong noise_bits = resources::precision_sum(working / k, mean_bits);
        const bool noiseless = source.exact_local_models();
        Magnitude box_radius;
        mag_min(box_radius.value(), arb_radref(acb_realref(component.box.value())),
                arb_radref(acb_imagref(component.box.value())));
        ComplexBall lower, upper;
        slong pair_precision = 0;
        for (unsigned attempt = 0; attempt < centre_steps && finite; ++attempt) {
            if (options.cancelled())
                return false;
            Magnitude size, bound;
            acb_get_mag_lower(bound.value(), center.value());
            slong aim = working;
            if (!noiseless && !mag_is_zero(bound.value()) && !mag_is_zero(box_radius.value())) {
                const double box_log2 = std::ceil(mag_get_d_log2_approx(bound.value()) -
                                                  mag_get_d_log2_approx(box_radius.value()));
                // A box wider than the centre's modulus bounds nothing; one finer than the working precision
                // asks for all of it.
                const slong box_bits =
                    box_log2 >= static_cast<double>(working)
                        ? working
                        : resources::precision_sum(static_cast<slong>(std::max(box_log2, 0.0)), mean_bits);
                const slong centre_bits = std::max(noise_bits, box_bits);
                aim = std::clamp<slong>(std::max(centre_bits >> (centre_steps - 1 - attempt), box_bits), 1, working);
            }
            mag_mul_2exp_si(bound.value(), bound.value(), -aim);
            pair_precision = mag_is_zero(bound.value())
                                 ? working
                                 : std::max(pair_precision,
                                            std::min(working, resources::whole_limbs(resources::precision_sum(
                                                                  aim, regeneration_guard_bits))));
            for (;;) {
                finite = source.taylor_pair(lower.value(), upper.value(), center.value(), pair_precision, k) &&
                         !acb_contains_zero(upper.value());
                if (finite) {
                    acb_div(step.value(), lower.value(), upper.value(), pair_precision);
                    acb_div_ui(step.value(), step.value(), static_cast<ulong>(k), pair_precision);
                    finite = acb_is_finite(step.value());
                }
                if (finite) {
                    mag_add(size.value(), arb_radref(acb_realref(step.value())), arb_radref(acb_imagref(step.value())));
                    if (mag_cmp(size.value(), bound.value()) <= 0)
                        break;
                }
                if (pair_precision >= working)
                    break;
                pair_precision = std::min(working, resources::precision_product(pair_precision, 2));
            }
            if (finite) {
                const slong kept = std::min(working, resources::precision_sum(aim, regeneration_guard_bits));
                acb_sub(center.value(), center.value(), step.value(), working);
                acb_get_mid(center.value(), center.value());
                arf_set_round(arb_midref(acb_realref(center.value())), arb_midref(acb_realref(center.value())), kept,
                              ARF_RND_NEAR);
                arf_set_round(arb_midref(acb_imagref(center.value())), arb_midref(acb_imagref(center.value())), kept,
                              ARF_RND_NEAR);
                finite = acb_is_finite(center.value());
            }
        }
        // The model: the expansion's first k + 1 coefficients, whose whole shift can dominate a solve's work.
        // Where the head's passes are the fewer (heuristics::restart_model_from_head) the head alone is
        // computed, and the dropped terms at a midpoint come from the source's value there (headed, below).
        // History: docs/code-history-20261021.md, secular_refinement.cpp 6.
        const auto head_passes = source.taylor_head_passes(working, k + 1);
        const auto value_passes = source.taylor_head_passes(working, 1);
        const bool headed = finite && head_passes && value_passes &&
                            heuristics::restart_model_from_head(source.degree() + 1, working, *head_passes,
                                                                *head_passes +
                                                                    static_cast<double>(k) * *value_passes) &&
                            source.taylor_head(local.value(), center.value(), working, k + 1);
        if (finite && (headed || source.taylor(local.value(), center.value(), working))) {
            // The dropped terms' magnitudes, for place_undetermined_hints.
            std::vector<Magnitude> tail;
            for (slong j = k + 1; j < acb_poly_length(local.value()); ++j) {
                tail.emplace_back();
                acb_get_mag(tail.back().value(), local.value()->coeffs + j);
            }
            acb_poly_truncate(local.value(), k + 1);
            if (acb_poly_degree(local.value()) == k && !acb_contains_zero(local.value()->coeffs + k)) {
                slong zeros = 0;
                while (zeros < k && acb_is_zero(local.value()->coeffs + zeros))
                    ++zeros;
                ComplexPolynomial reduced;
                acb_poly_shift_right(reduced.value(), local.value(), zeros);
                RootVector initial(k - zeros);
                // The model's roots are hints only where this precision
                // determines them. Where a vertex of the model's Newton
                // polygon is a coefficient with no correct bit, the roots of
                // the midpoint polynomial are roots of rounding noise, and
                // zero midpoints can put two hints on the centre itself, a
                // pair the iteration does not separate. Such a model gives
                // scales only, and the hints go on its polygon's circles:
                // distinct points, as far from the centre as the roots can be
                // on this evidence. History:
                // docs/code-history-20261021.md, secular_refinement.cpp 7.
                bool determined = true;
                UndeterminedReach reach;
                const bool initialized =
                    zeros < k && scale_candidates(initial.data(), reduced.value(), working, determined, &reach);
                // A coarse round can precede a much finer caller request.
                // Keep scale hints only when they fit the final counted
                // placement target. They remain proposals: the original
                // source's regeneration and certificates still decide success.
                const slong hint_resolution = std::max(resolution, counted_hint_resolution);
                bool within_resolution = initialized && resolution > 0 && counted_hint_resolution > 0;
                if (within_resolution) {
                    Magnitude floor, radius;
                    acb_get_mag_lower(floor.value(), center.value());
                    mag_mul_2exp_si(floor.value(), floor.value(), -hint_resolution);
                    within_resolution = !mag_is_zero(floor.value());
                    for (slong j = 0; j < k - zeros && within_resolution; ++j) {
                        acb_get_mag(radius.value(), initial[j]);
                        within_resolution = mag_cmp(radius.value(), floor.value()) <= 0;
                    }
                }
                const bool use_scale_hints = engine_goal(options) == CPS_ROOT_APPROXIMATION && within_resolution;
                // Such hints are as far from the centre as the unresolved coefficients' bounds allow: with the
                // model's noise at 2^-N of its leading coefficient, up to 2^(-N / k). The members are then
                // iterated at `resolution` bits, where a node is settled once its correction is below
                // 2^(2 - resolution) of its modulus (model_stop). Hints within 2^-resolution of the centre's
                // modulus are settled as they stand, and the next precision's restart places them again. Hints
                // farther out the iteration has to carry in, and k nodes close on a cluster of k roots
                // linearly. The model costs one expansion; a sweep costs an evaluation at every node. So a
                // model whose noise places hints beyond that resolution is computed again with the precision
                // that brings them inside it: the largest circle an unresolved coefficient bounds
                // (scale_candidates' reach), where every bit of the excess takes as many bits of noise as the
                // circle has roots (one more for the two logarithms' rounding). A circle between resolved
                // coefficients is the roots' own: five roots of which four stand 2^-300 from the fifth keep
                // that circle at every precision, and raising the model for it gains nothing. History:
                // docs/code-history-20261021.md, secular_refinement.cpp 8.
                if (initialized && !determined && resolution > 0 && retry + 1 < attempts &&
                    static_cast<uint64_t>(working) < options.max_work_bits && reach.roots > 0) {
                    Magnitude floor;
                    acb_get_mag_lower(floor.value(), center.value());
                    mag_mul_2exp_si(floor.value(), floor.value(), -resolution);
                    const double excess = mag_is_zero(floor.value())
                                              ? 0
                                              : std::ceil(reach.radius_log2 - mag_get_d_log2_approx(floor.value()));
                    if (std::isfinite(excess) && excess > 0) {
                        // In whole limbs: the arithmetic costs the same, and the precisions the source
                        // prepares its coefficients at stay few.
                        const slong more = resources::precision_product(reach.roots, static_cast<slong>(excess) + 1);
                        working = std::min<slong>(static_cast<slong>(options.max_work_bits),
                                                  resources::whole_limbs(resources::precision_sum(working, more)));
                        if (policy.trace) {
                            TraceEvent event{"local_model_raised", "precision", static_cast<uint64_t>(working),
                                             static_cast<uint64_t>(k), &statistics, nullptr};
                            event.engine = "secular";
                            policy.trace(policy.trace_data, event);
                        }
                        continue;
                    }
                }
                for (slong j = 0; j < zeros; ++j)
                    acb_zero(out + j);
                // A unit-circle start needs thousands of iterations for a
                // tiny cluster even when its local quadratic is well resolved.
                // Use coefficient-scale radii for this low-degree polynomial.
                // A bounded native local packet, not a recursive global solve.
                // Even converged local-model roots remain proposals; the
                // original source and its retained cover certify every result.
                certpolysolve::detail::CandidateStatistics local_statistics{};
                uint64_t prepaid_local_terms = 0;
                TraceScope trace(&policy, "local_candidates", static_cast<uint64_t>(working),
                                 static_cast<uint64_t>(k), &statistics);
                CandidateProgress progress;
                // Hints this model does not determine: all of them (scales only), or those placed below.
                slong placed = 0;
                if (zeros == k)
                    progress = CandidateProgress{0, true, true, false};
                else if (!determined || use_scale_hints) {
                    if (initialized)
                        _acb_vec_set(out + zeros, initial.data(), k - zeros);
                    progress.finite = initialized;
                } else if (const unsigned sweeps =
                               options.packet_sweeps(heuristics::local_restart_sweeps, statistics);
                           !sweeps) {
                    // The call's iteration budget is spent: this component
                    // gets no local seeds and takes the caller's placement by
                    // scale. (Stopping the restart here to resume it would
                    // repeat the components before this one on every call.)
                    progress = CandidateProgress{};
                } else {
                    if (accounting) {
                        LocalPacketMeter meter{*accounting};
                        // The old call used a transient workspace with a
                        // default policy. Keep those exact defaults; only
                        // attach this packet's synchronous service hook.
                        CandidateWorkspace workspace;
                        workspace.research_set_local_packet_service(LocalPacketMeter::service, &meter);
                        progress = candidate_packet(out + zeros, reduced.value(), initialized ? initial.data() : nullptr,
                                                    working, sweeps, CandidateEngine::Aberth, options, local_statistics,
                                                    {}, false, &workspace);
                        // Cancellation can leave a partially completed sweep.
                        // Charge that remainder at the same actual precision,
                        // without creating a service boundary or another unit.
                        meter.charge({local_statistics.evaluation_terms, local_statistics.root_interactions,
                                      local_statistics.candidate_sweeps, static_cast<uint64_t>(working), 0});
                        prepaid_local_terms = meter.terms;
                    } else
                    progress = candidate_packet(out + zeros, reduced.value(), initialized ? initial.data() : nullptr,
                                                working, sweeps, CandidateEngine::Aberth, options, local_statistics);
                    if (progress.finite && !progress.cancelled) {
                        // The dropped terms at a midpoint x of two hints, over x^zeros: a bound on what the
                        // source's value there differs from the model's by. From the whole shift it is the
                        // sum of their magnitudes, sum |a_j| |x|^(j - zeros) over j > k. From the head it is
                        // the difference itself: the source's value at centre + x (one more pass, on the
                        // exact sum) less the model's, enclosed. Either way the model's value with that
                        // radius encloses the source's value over x^zeros; the second is the tighter where
                        // the dropped terms cancel, and no bound where the evaluation fails.
                        const slong degree = acb_poly_degree(reduced.value());
                        const std::function<void(mag_t, const acb_t)> summed = [&](mag_t tau, const acb_t middle) {
                            Magnitude rho, power, term;
                            acb_get_mag(rho.value(), middle);
                            mag_pow_ui(power.value(), rho.value(), static_cast<ulong>(degree + 1));
                            mag_zero(tau);
                            for (const auto &coefficient : tail) {
                                mag_mul(term.value(), coefficient.value(), power.value());
                                mag_add(tau, tau, term.value());
                                mag_mul(power.value(), power.value(), rho.value());
                            }
                        };
                        //
                        // The difference is the evaluation's own error where the source's value cancels to
                        // below its precision, as next to a cluster of k roots: then every midpoint's bound
                        // holds zero, the hints are grouped and the component is taken for unresolved, and its
                        // next model starts at k times the precision. The dropped terms' magnitudes summed term
                        // by term over the source's terms (taylor_tail_bound) are a bound with no evaluation in
                        // it; the smaller of the two is taken, and where the magnitudes alone keep zero out the
                        // evaluation is not made. History:
                        // docs/code-history-20261021.md, secular_refinement.cpp 9.
                        const std::function<void(mag_t, const acb_t)> measured = [&](mag_t tau, const acb_t middle) {
                            ComplexBall point, value, modelled, divisor;
                            ComplexPolynomial at;
                            Magnitude bound, power;
                            mag_inf(tau);
                            acb_poly_evaluate(modelled.value(), reduced.value(), middle, working);
                            if (source.taylor_tail_bound(bound.value(), center.value(), middle, working, k + 1)) {
                                if (zeros) {
                                    acb_get_mag_lower(power.value(), middle);
                                    mag_pow_ui_lower(power.value(), power.value(), static_cast<ulong>(zeros));
                                    mag_div(bound.value(), bound.value(), power.value());
                                }
                                if (mag_is_finite(bound.value())) {
                                    mag_set(tau, bound.value());
                                    acb_set(value.value(), modelled.value());
                                    acb_add_error_mag(value.value(), tau);
                                    if (!acb_contains_zero(value.value()))
                                        return;
                                }
                            }
                            acb_add(point.value(), center.value(), middle, ARF_PREC_EXACT);
                            if (!source.taylor_head(at.value(), point.value(), working, 1))
                                return;
                            acb_poly_get_coeff_acb(value.value(), at.value(), 0);
                            if (zeros) {
                                acb_pow_ui(divisor.value(), middle, static_cast<ulong>(zeros), working);
                                if (acb_contains_zero(divisor.value()))
                                    return;
                                acb_div(value.value(), value.value(), divisor.value(), working);
                            }
                            acb_sub(value.value(), value.value(), modelled.value(), working);
                            if (acb_is_finite(value.value())) {
                                acb_get_mag(bound.value(), value.value());
                                mag_min(tau, tau, bound.value());
                            }
                        };
                        placed = place_undetermined_hints(out + zeros, k - zeros, reduced.value(),
                                                          headed ? measured : summed, working);
                        if (placed && policy.trace) {
                            TraceEvent event{"local_hints_undetermined", "placement", static_cast<uint64_t>(working),
                                             static_cast<uint64_t>(placed), &statistics, nullptr};
                            event.engine = "secular";
                            policy.trace(policy.trace_data, event);
                        }
                    }
                }
                // candidate_ns is charged by the enclosing initialization
                // timer. Preserve its nested mathematical work without double
                // charging the containing stage.
                if (accounting && prepaid_local_terms > UINT64_MAX - accounting->prepaid_terms)
                    throw std::logic_error("local prepaid term counter overflowed");
                statistics.evaluation_ns += local_statistics.evaluation_ns;
                statistics.interaction_ns += local_statistics.interaction_ns;
                statistics.evaluations += local_statistics.evaluations;
                statistics.evaluation_terms += local_statistics.evaluation_terms;
                statistics.candidate_sweeps += local_statistics.candidate_sweeps;
                statistics.root_interactions += local_statistics.root_interactions;
                // Commit the exclusion only after the same raw terms reach
                // the parent. A throwing service/allocation before this merge
                // leaves native work charged, but cannot make the outer RAII
                // meter subtract terms its statistics never received.
                if (accounting)
                    accounting->prepaid_terms += prepaid_local_terms;
                if (progress.cancelled)
                    return false;
                // A finite packet supplies restart hints whether or not it
                // converged: the hints are proposals the original source
                // certifies. Its convergence test is relative to the working
                // precision, so raising the precision until it converges buys
                // nothing. Only a measured deficiency of the local model,
                // below, raises the precision. History:
                // docs/code-history-20261021.md, secular_refinement.cpp 10.
                const bool usable = progress.finite;
                bool inside = true;
                for (slong j = 0; j < k; ++j) {
                    acb_add(out + j, out + j, center.value(), working);
                    acb_get_mid(out + j, out + j);
                    inside = inside && acb_is_finite(out + j) && acb_contains(component.box.value(), out + j);
                }
                // A low-degree model whose roots leave the certified component
                // is not a useful initializer, even if its own solve succeeds.
                if (usable) {
                    if (unresolved)
                        *unresolved = !determined || placed > 0 || use_scale_hints;
                    return inside;
                }
            }
        }
        if (working == static_cast<slong>(options.max_work_bits) || retry + 1 == attempts ||
            source.degree() <= heuristics::small_complex_degree)
            return false;
        working = std::min<slong>(static_cast<slong>(options.max_work_bits), resources::precision_product(working, 2));
    }
    return false;
}
std::vector<size_t> SecularRefinement::restart_members(const Root &component,
                                                       std::span<const ComplexBall> points) const {
    std::vector<size_t> members;
    if (component.multiplicity != 0)
        return members;
    for (size_t i = 0; i < points.size(); ++i)
        if (!certified_.inactive[i] && acb_contains(component.box.value(), points[i].value()))
            members.push_back(i);
    if (members.size() < 2 || members.size() * multiplicity_ != component.count)
        members.clear();
    return members;
}
bool SecularRefinement::unresolved_members(const std::vector<size_t> &members) const {
    return !members.empty() && std::all_of(members.begin(), members.end(), [&](size_t i) {
        return residue_unresolved(residues_[i].value());
    });
}
// The secular cover gives node i the disk D(b_i, n |w_i|) (secular_inclusion_
// cover). Put the k members of a component on a circle of radius rho about
// their centroid c. Then |w_i| is at most about rho/k + E/(k rho^(k-1)), where
// E = rad(w_i) prod_{j in C, j != i} |b_i - b_j| is the evaluation error of
// p at the members over the exterior product, which moving the members inside
// the cluster does not change. The component's radius, about
// rho (1 + n/k) + n E / (k rho^(k-1)), is least at rho^k = n (k-1) E / (k + n).
// Every quantity is read from the current nodes and residue enclosures; the
// certificate of whatever geometry results is still checked from scratch.
bool SecularRefinement::certificate_reseed(const std::vector<size_t> &members, slong precision) {
    const size_t k = members.size(), n = nodes_.size();
    if (k < 2)
        return false;
    double error_log2 = -std::numeric_limits<double>::infinity();
    Magnitude radius, distance;
    ComplexBall difference;
    for (const auto i : members) {
        mag_add(radius.value(), arb_radref(acb_realref(residues_[i].value())),
                arb_radref(acb_imagref(residues_[i].value())));
        if (mag_is_zero(radius.value()))
            return false;
        double member_log2 = mag_get_d_log2_approx(radius.value());
        for (const auto j : members)
            if (j != i) {
                acb_sub(difference.value(), nodes_[i].value(), nodes_[j].value(), precision);
                acb_get_mag(distance.value(), difference.value());
                if (mag_is_zero(distance.value()))
                    return false;
                member_log2 += mag_get_d_log2_approx(distance.value());
            }
        error_log2 = std::max(error_log2, member_log2);
    }
    const double kd = static_cast<double>(k), nd = static_cast<double>(n);
    const double radius_log2 = (std::log2(nd * (kd - 1) / (kd + nd)) + error_log2) / kd;
    if (!std::isfinite(radius_log2))
        return false;
    ComplexBall centre;
    for (const auto i : members)
        acb_add(centre.value(), centre.value(), nodes_[i].value(), precision);
    acb_div_ui(centre.value(), centre.value(), static_cast<ulong>(k), precision);
    acb_get_mid(centre.value(), centre.value());
    RealBall scale;
    arb_one(scale.value());
    // 2^radius_log2 as an exact dyadic with a 53-bit significand.
    const double whole = std::floor(radius_log2);
    arb_set_d(scale.value(), std::exp2(radius_log2 - whole));
    arb_mul_2exp_si(scale.value(), scale.value(), static_cast<slong>(whole));
    RootVector directions(static_cast<slong>(k));
    _acb_vec_unit_roots(directions.data(), static_cast<slong>(k), static_cast<slong>(k), precision);
    for (size_t j = 0; j < k; ++j) {
        auto &proposal = proposals_[members[j]];
        acb_mul_arb(proposal.value(), directions.data() + j, scale.value(), precision);
        acb_add(proposal.value(), proposal.value(), centre.value(), precision);
        acb_get_mid(proposal.value(), proposal.value());
    }
    return true;
}
bool SecularRefinement::local_restart(size_t members) const {
    const auto *polynomial = polynomial_;
    return polynomial && members <= static_cast<size_t>(heuristics::local_component_degree);
}
// The index whose Taylor coefficient dominates on the circle of radius R
// about the shift's centre: sum_{j != c} |b_j| R^j < |b_c| R^c, upper bounds
// on the left and a lower bound on the right. It is the number of roots in
// the disk (Pellet), and (Rouche against the truncation's tail, which the sum
// bounds) the degree-c truncation has its c roots there too: the local model
// describes the disk. -1 if no index dominates at this precision. No Graeffe
// step: a dominance reached only after root squaring says nothing about the
// truncation.
static slong dominant_count(const acb_poly_t local, const mag_t radius) {
    const slong length = acb_poly_length(local);
    std::vector<Magnitude> upper(static_cast<size_t>(length)), lower(static_cast<size_t>(length));
    Magnitude power_upper, power_lower, total;
    mag_one(power_upper.value());
    mag_one(power_lower.value());
    for (slong j = 0; j < length; ++j) {
        acb_get_mag(upper[static_cast<size_t>(j)].value(), local->coeffs + j);
        mag_mul(upper[static_cast<size_t>(j)].value(), upper[static_cast<size_t>(j)].value(), power_upper.value());
        acb_get_mag_lower(lower[static_cast<size_t>(j)].value(), local->coeffs + j);
        mag_mul_lower(lower[static_cast<size_t>(j)].value(), lower[static_cast<size_t>(j)].value(), power_lower.value());
        mag_add(total.value(), total.value(), upper[static_cast<size_t>(j)].value());
        mag_mul(power_upper.value(), power_upper.value(), radius);
        mag_mul_lower(power_lower.value(), power_lower.value(), radius);
    }
    if (!mag_is_finite(total.value()))
        return -1;
    // The others' sum for each index from prefix and suffix sums (upper
    // bounds both; total - upper_j is no rigorous subtraction in magnitudes).
    std::vector<Magnitude> suffix(static_cast<size_t>(length) + 1);
    for (slong j = length; j-- > 0;)
        mag_add(suffix[static_cast<size_t>(j)].value(), suffix[static_cast<size_t>(j) + 1].value(),
                upper[static_cast<size_t>(j)].value());
    Magnitude prefix, others;
    for (slong j = 0; j < length; ++j) {
        mag_add(others.value(), prefix.value(), suffix[static_cast<size_t>(j) + 1].value());
        if (!mag_is_zero(lower[static_cast<size_t>(j)].value()) &&
            mag_cmp(others.value(), lower[static_cast<size_t>(j)].value()) < 0)
            return j;
        mag_add(prefix.value(), prefix.value(), upper[static_cast<size_t>(j)].value());
    }
    return -1;
}
std::vector<std::vector<size_t>> SecularRefinement::node_groups() const {
    std::vector<size_t> index;
    for (size_t i = 0; i < nodes_.size(); ++i)
        if (!certified_.inactive[i])
            index.push_back(i);
    const size_t n = index.size();
    std::vector<std::vector<size_t>> groups;
    if (n < 3)
        return groups;
    std::vector<double> re(n), im(n);
    for (size_t a = 0; a < n; ++a) {
        re[a] = arf_get_d(arb_midref(acb_realref(nodes_[index[a]].value())), ARF_RND_NEAR);
        im[a] = arf_get_d(arb_midref(acb_imagref(nodes_[index[a]].value())), ARF_RND_NEAR);
        if (!std::isfinite(re[a]) || !std::isfinite(im[a]))
            return groups;
    }
    const auto dist = [&](size_t a, size_t b) { return std::hypot(re[a] - re[b], im[a] - im[b]); };
    // Prim's minimum spanning tree: single linkage.
    std::vector<double> best(n, std::numeric_limits<double>::infinity());
    std::vector<size_t> parent(n, 0);
    std::vector<uint8_t> in_tree(n, 0);
    std::vector<std::tuple<double, size_t, size_t>> edges;
    best[0] = 0;
    for (size_t step = 0; step < n; ++step) {
        size_t u = n;
        for (size_t a = 0; a < n; ++a)
            if (!in_tree[a] && (u == n || best[a] < best[u]))
                u = a;
        in_tree[u] = 1;
        if (step)
            edges.emplace_back(best[u], parent[u], u);
        for (size_t a = 0; a < n; ++a)
            if (!in_tree[a]) {
                const double d = dist(u, a);
                if (d < best[a]) {
                    best[a] = d;
                    parent[a] = u;
                }
            }
    }
    std::sort(edges.begin(), edges.end());
    // Kruskal over the tree's edges in increasing order: every merged set is
    // a single-linkage cluster; it is a group if its diameter is below the
    // weight of the edge that later merges it with the rest (its distance
    // to every other node).
    std::vector<size_t> root(n);
    std::iota(root.begin(), root.end(), 0);
    const auto find = [&](size_t a) {
        while (root[a] != a) a = root[a] = root[root[a]];
        return a;
    };
    std::vector<std::vector<size_t>> members(n);
    for (size_t a = 0; a < n; ++a) members[a] = {a};
    std::vector<double> diameter(n, 0);
    for (size_t e = 0; e < edges.size(); ++e) {
        const auto [w, a, b] = edges[e];
        // Sets merged at this weight that stood apart from everything else
        // by more than their own diameter are groups.
        for (const size_t side : {find(a), find(b)})
            if (members[side].size() >= 2 && members[side].size() < n && diameter[side] < w) {
                auto set = members[side];
                for (auto &i : set) i = index[i];
                std::sort(set.begin(), set.end());
                groups.push_back(std::move(set));
            }
        const size_t ra = find(a), rb = find(b);
        if (ra == rb) continue;
        double d = std::max({diameter[ra], diameter[rb], w});
        for (const size_t x : members[ra])
            for (const size_t y : members[rb])
                d = std::max(d, dist(x, y));
        root[rb] = ra;
        members[ra].insert(members[ra].end(), members[rb].begin(), members[rb].end());
        members[rb].clear();
        diameter[ra] = d;
    }
    // The active nodes as a whole, when the others are certified: a group if
    // it stands apart from every certified node by more than its diameter
    // (the groups above never take the whole active set, which a last
    // uncertified cluster can be). History:
    // docs/code-history-20261021.md, secular_refinement.cpp 23.
    if (n < nodes_.size()) {
        const size_t whole = find(0);
        double apart = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < nodes_.size(); ++i)
            if (certified_.inactive[i]) {
                const double x = arf_get_d(arb_midref(acb_realref(nodes_[i].value())), ARF_RND_NEAR);
                const double y = arf_get_d(arb_midref(acb_imagref(nodes_[i].value())), ARF_RND_NEAR);
                for (size_t a = 0; a < n; ++a)
                    apart = std::min(apart, std::hypot(re[a] - x, im[a] - y));
            }
        if (members[whole].size() == n && diameter[whole] < apart) {
            auto set = members[whole];
            for (auto &i : set) i = index[i];
            std::sort(set.begin(), set.end());
            groups.push_back(std::move(set));
        }
    }
    return groups;
}
// Whether `count` roots about `point` are a declined cluster's: the same
// count, and the point inside the disk beyond which its other roots lie.
bool SecularRefinement::inside_declined(const acb_t point, slong count) const {
    ComplexBall difference;
    Magnitude distance;
    for (const auto &declined : declined_clusters_) {
        if (declined.count != count)
            continue;
        acb_sub(difference.value(), point, declined.mean.value(), 64);
        acb_get_mag(distance.value(), difference.value());
        if (mag_is_zero(distance.value()) || mag_get_d_log2_approx(distance.value()) < declined.outer_log2)
            return true;
    }
    return false;
}
// The local model of a group of nodes. The disk's root count K is exact (an
// isolating plateau of the model's counts about the group's centroid); the
// model is the Taylor shift at the K roots' mean, truncated at degree K +
// extra, the extra terms making the truncation error (R/D)^extra meet the
// goal, where R is the cluster's own scale (the smallest circle about the
// mean on which K dominates, Pellet) and D the distance to the other roots
// (the plateau's outer circle). On success `hints` holds the K roots of the
// truncation nearest the centre, in original coordinates. The shift's
// precision rises until the model's coefficients are determined.
namespace {
// The mean of an isolated cluster's `count` roots, from an estimate of it, to the bits a model at `working` bits
// can use: Newton's steps on p^(count-1), as component_seeds takes them for a component's centre. That derivative
// has one root by the cluster, the roots' mean up to count (cluster radius)^2 over the distance to the other
// roots, and about a point much nearer the cluster than the other roots it is linear: its value and slope there
// are the expansion's coefficients count - 1 and count (PolynomialSource::taylor_pair, one joint evaluation), and
// the step is the first over count times the second.
//
// Why the estimate does not do. The moments (root_mean) give the mean to 2^-53 of their circle times the count,
// which can be far coarser than a narrow cluster. About a point d from a cluster narrower than d the expansion is
// that of (w + d)^count: the count dominates only on circles beyond count / ln 2 times d, the coefficients below
// it carry d and not the cluster's shape, and to read the shape through them takes count bits for every bit of d
// over the cluster's radius. Short of those bits the model's roots are a ring of rounding noise, and k nodes
// close on a cluster of k roots linearly, at 2 / k bits a sweep.
//
// The bits kept are component_seeds': a model with noise 2^-N of its leading coefficient resolves nothing nearer
// its centre than 2^(-N / count), so working / count bits of the centre and those of 4 count (the offset enters
// the expansion count-fold, and is to lie under its polygon), with the regeneration's guard; and the estimate's
// own limb at least, so that no precision makes the centre coarser than the moments gave it. They are bits of
// the centre's modulus and the radius of `inside`, the disk the cluster lies in (the plateau's inner circle
// about the nodes' centroid; the mean lies in it too): both parts are rounded on that one grid, and the mean
// of a real polynomial's cluster about the real axis, whose imaginary part is the moments' rounding, comes out
// real. The expansion's passes multiply by the centre, and a centre of that length leaves coefficient
// count - 1 a number the shift determines. A centre whose own rounding could pass for the cluster's shape
// would have to leave coefficient zero above the shift's noise, which these bits exclude: a model that comes
// out determined is determined by the cluster.
//
// A step is taken only while it is larger than its own enclosure's radius: then the centre is off by about
// that step, and taking it leaves what the enclosure does not resolve. A step within its radius says that
// this precision locates the mean no better than the centre already is, and the centre stands; at a low
// precision that is the moments' estimate, untouched. (Taken regardless and rounded to the precision's own few
// bits, such a step can move a centre that was inside its cluster off it.) Newton's iteration doubles the correct
// bits, so from one bit it takes at most as many steps as the kept bits have binary digits, and one more to see
// the rounded centre repeat. History: docs/code-history-20261021.md, secular_refinement.cpp 24.
//
// False, with the estimate restored, when a step leaves the disk: Newton's iteration is then not converging on
// this cluster. `operations` takes the steps' multiplications at the working precision.
bool locate_mean(acb_t centre, PolynomialSource &source, slong count, slong working, const acb_t disk_centre,
                 const arf_t disk_radius, const cps_options &options, double &operations) {
    const slong mean_bits = static_cast<slong>(FLINT_BIT_COUNT(static_cast<ulong>(count))) + 2;
    const slong kept = std::min(working, std::max<slong>(FLINT_BITS, resources::precision_sum(
                                                                         resources::precision_sum(working / count, mean_bits),
                                                                         regeneration_guard_bits)));
    const unsigned steps = static_cast<unsigned>(FLINT_BIT_COUNT(static_cast<ulong>(kept))) + 1;
    ComplexBall estimate, lower, upper, step, next, offset;
    Float distance, scaled;
    Integer grid;
    Magnitude scale, radius, size, error;
    acb_set(estimate.value(), centre);
    // The grid: 2^-kept of the centre's modulus and the disk's radius together.
    acb_get_mag(scale.value(), centre);
    arf_get_mag(radius.value(), disk_radius);
    mag_add(scale.value(), scale.value(), radius.value());
    if (mag_is_zero(scale.value()) || !mag_is_finite(scale.value()))
        return true;
    const double scale_log2 = std::ceil(mag_get_d_log2_approx(scale.value()));
    if (!std::isfinite(scale_log2))
        return true;
    const slong ulp = static_cast<slong>(scale_log2) - kept;
    // The value and the slope of a polynomial of degree n - count + 1 at the centre.
    const double pair_operations = 2 * static_cast<double>(source.degree() - count + 2);
    for (unsigned attempt = 0; attempt < steps && !options.cancelled(); ++attempt) {
        operations += pair_operations;
        if (!source.taylor_pair(lower.value(), upper.value(), centre, working, count) ||
            acb_contains_zero(upper.value()))
            break;
        acb_div(step.value(), lower.value(), upper.value(), working);
        acb_div_ui(step.value(), step.value(), static_cast<ulong>(count), working);
        if (!acb_is_finite(step.value()))
            break;
        mag_add(error.value(), arb_radref(acb_realref(step.value())), arb_radref(acb_imagref(step.value())));
        acb_get_mid(step.value(), step.value());
        acb_get_mag(size.value(), step.value());
        if (mag_cmp(size.value(), error.value()) <= 0)
            break;
        acb_sub(next.value(), centre, step.value(), ARF_PREC_EXACT);
        for (arb_ptr part : {acb_realref(next.value()), acb_imagref(next.value())}) {
            arf_mul_2exp_si(scaled.value(), arb_midref(part), -ulp);
            arf_get_fmpz(grid.value(), scaled.value(), ARF_RND_NEAR);
            arb_set_fmpz(part, grid.value());
            arb_mul_2exp_si(part, part, ulp);
        }
        acb_sub(offset.value(), next.value(), disk_centre, working);
        acb_get_abs_ubound_arf(distance.value(), offset.value(), FLINT_BITS);
        if (arf_cmp(distance.value(), disk_radius) > 0) {
            acb_set(centre, estimate.value());
            return false;
        }
        if (acb_equal(next.value(), centre))
            break;
        acb_set(centre, next.value());
    }
    return true;
}
// The worst relative condition of the inner roots on the scaled model,
// sum_j |q_j| |u|^j / (|u| |q'(u)|) with |q'(u_i)| = |q_k| prod_{l != i}
// |u_i - u_l| over the model's own roots: the root's relative motion per
// relative coefficient perturbation. An estimate in binary64 logarithms
// (the scaled roots are at least 2^-8 apart and the terms are positive).
double local_condition_log2(const ComplexPolynomial &scaled, std::span<const ComplexBall> roots, slong inner_count) {
    const slong k = acb_poly_degree(scaled.value());
    std::vector<double> coefficient_log2(static_cast<size_t>(k + 1), -std::numeric_limits<double>::infinity());
    for (slong j = 0; j <= k; ++j) {
        Magnitude size;
        acb_get_mag(size.value(), scaled.value()->coeffs + j);
        if (!mag_is_zero(size.value()))
            coefficient_log2[static_cast<size_t>(j)] = mag_get_d_log2_approx(size.value());
    }
    std::vector<std::complex<double>> u(static_cast<size_t>(k));
    for (slong j = 0; j < k; ++j) {
        const acb_srcptr r = roots[static_cast<size_t>(j)].value();
        u[static_cast<size_t>(j)] = {arf_get_d(arb_midref(acb_realref(r)), ARF_RND_NEAR),
                                     arf_get_d(arb_midref(acb_imagref(r)), ARF_RND_NEAR)};
    }
    double worst = 0;
    for (slong i = 0; i < inner_count; ++i) {
        const double size_log2 = std::log2(std::abs(u[static_cast<size_t>(i)]));
        if (!std::isfinite(size_log2))
            return std::numeric_limits<double>::infinity();
        double top = -std::numeric_limits<double>::infinity();
        for (slong j = 0; j <= k; ++j)
            top = std::max(top, coefficient_log2[static_cast<size_t>(j)] + static_cast<double>(j) * size_log2);
        double sum = 0;
        for (slong j = 0; j <= k; ++j)
            sum += std::exp2(coefficient_log2[static_cast<size_t>(j)] + static_cast<double>(j) * size_log2 - top);
        double slope_log2 = coefficient_log2[static_cast<size_t>(k)];
        for (slong l = 0; l < k; ++l)
            if (l != i)
                slope_log2 += std::log2(std::abs(u[static_cast<size_t>(i)] - u[static_cast<size_t>(l)]));
        const double condition = top + std::log2(sum) - size_log2 - slope_log2;
        if (!std::isfinite(condition))
            return std::numeric_limits<double>::infinity();
        worst = std::max(worst, condition);
    }
    return worst;
}
} // namespace
std::optional<cps_options> local_model_options(const cps_options &options, uint64_t bits,
                                             const SearchStatistics &search, const CandidateStatistics &statistics) {
    const uint64_t work_left = options.work_units > search.work_units ? options.work_units - search.work_units : 0;
    const uint64_t sweeps_left = options.sweeps_left(statistics);
    if (!work_left || !sweeps_left)
        return std::nullopt; // zero must not turn into an unlimited child allowance
    cps_options nested;
    nested.work_units = std::min(nested.work_units, work_left);
    nested.work_units_default = false;
    nested.max_iterations = options.max_iterations ? sweeps_left : 0;
    nested.goal = CPS_ROOT_APPROXIMATION;
    nested.method = CPS_METHOD_SECULAR;
    nested.bits = bits;
    nested.max_work_bits = options.max_work_bits;
    nested.temporary_bytes = options.temporary_bytes;
    nested.cancel = options.cancel;
    nested.cancel_data = options.cancel_data;
    nested.cancel_flag = options.cancel_flag;
    nested.deadline = options.deadline;
    nested.exact_budget = 0;
    nested.heuristic_budget = options.heuristic_budget;
    return nested;
}

SecularRefinement::GroupModel SecularRefinement::group_model(const std::vector<size_t> &members, const cps_options &options,
                                                           SearchStatistics &search, CandidateStatistics &statistics) {
    GroupModel model;
    if (options.heuristic_budget) {
        if (options.heuristic_budget->exhausted()) return model;
        options.heuristic_budget->quadratic();
    }
    auto *polynomial = polynomial_;
    if (!polynomial)
        return model;
    const slong m = static_cast<slong>(members.size());
    ComplexBall center;
    for (const size_t i : members)
        acb_add(center.value(), center.value(), nodes_[i].value(), 64);
    acb_div_ui(center.value(), center.value(), static_cast<ulong>(m), 64);
    acb_get_mid(center.value(), center.value());
    Magnitude spread, distance;
    ComplexBall difference;
    for (const size_t i : members) {
        acb_sub(difference.value(), nodes_[i].value(), center.value(), 64);
        acb_get_mag(distance.value(), difference.value());
        if (mag_cmp(distance.value(), spread.value()) > 0) mag_set(spread.value(), distance.value());
    }
    if (mag_is_zero(spread.value()) || !mag_is_finite(spread.value()))
        return model;
    // A declined cluster's own nodes (as many as its roots, about it) are not
    // walked again: the walk would find the same count. History:
    // docs/code-history-20261021.md, secular_refinement.cpp 25.
    if (inside_declined(center.value(), m))
        return model;
    // The count, from the circle of twice the group's radius down. While
    // nodes still converge from far away their spread is no measure of the
    // cluster, and Pellet's test on the group's own circle can fail round
    // after round where the counts certify the cluster. History:
    // docs/code-history-20261021.md, secular_refinement.cpp 26.
    DiskCountScratch scratch;
    const auto distances = node_distances_log2(nodes_, center.value());
    auto sorted = distances;
    std::sort(sorted.begin(), sorted.end());
    const auto plateau = isolating_plateau(nodes_, residues_, center.value(), distances, sorted,
                                           mag_get_d_log2_approx(spread.value()) + 1, scratch);
    if (!plateau || plateau->count + 1 >= source_->degree() || inside_declined(center.value(), plateau->count)) {
        work_ += service_charge(&policy_, ServiceOperation::Auxiliary, scratch.operations * operation_weight(FLINT_BITS));
        return model;
    }
    // The centre: the roots' mean. The nodes' centroid is off by the nodes
    // still migrating, which can put it outside a narrow cluster, and a model
    // about a point outside the cluster sees it as a cluster again. The moments
    // give an estimate of it, to their one limb; where the model about the
    // estimate shows it to be outside its cluster, the mean is located from
    // there (locate_mean, below).
    //
    // Except about the original origin of a centred source. The goals are
    // relative to it, the uncentred coefficients are exact there, and a mean
    // is only as accurate as the plateau's circle allows (2^-53 of its radius),
    // which can exceed a cluster's width: Pellet's test then certifies the
    // cluster only on circles some K / ln 2 times that error, and the model
    // cannot resolve it. When the plateau's inner disk holds the original
    // origin, the model is taken about the origin itself. History:
    // docs/code-history-20261021.md, secular_refinement.cpp 27.
    ComplexBall mean;
    bool about_origin = false;
    if (const acb_struct *offset = policy_.centred_origin_basis ? coordinates_.exact_offset() : nullptr) {
        ComplexBall original;
        Magnitude from_origin;
        acb_add(original.value(), center.value(), offset, ARF_PREC_EXACT);
        acb_get_mag(from_origin.value(), original.value());
        Float inner;
        set_radius_log2(inner.value(), plateau->inner_log2);
        Magnitude inner_bound;
        arf_get_mag_lower(inner_bound.value(), inner.value());
        if (mag_cmp(from_origin.value(), inner_bound.value()) < 0) {
            about_origin = true;
            acb_neg(mean.value(), offset);
        }
    }
    double moment_operations = 0;
    // At most one rounded binary64 precision's worth of quadrature points
    // per global root: each point visits all nodes, a quadratic model shell.
    constexpr uint64_t points_per_root = std::bit_ceil(static_cast<uint64_t>(std::numeric_limits<double>::digits));
    const bool centred = about_origin ||
                         root_mean(mean.value(), nodes_, residues_, center.value(), *plateau, sorted, moment_operations,
                             options.heuristic_budget ? resources::multiply(points_per_root, options.heuristic_budget->degree()) : UINT64_MAX);
    // Metered at one limb (the counts' and the moments' precision).
        work_ += service_charge(&policy_, ServiceOperation::Auxiliary, (scratch.operations + moment_operations) * operation_weight(FLINT_BITS));
    if (!centred)
        return model;
    Magnitude offset;
    acb_sub(difference.value(), mean.value(), center.value(), 64);
    acb_get_mag(offset.value(), difference.value());
    acb_swap(center.value(), mean.value());
    // The other roots lie beyond the plateau's outer circle about the
    // centroid, so at D = r_out - |mean - centroid| from the mean at least.
    const double offset_log2 = mag_is_zero(offset.value()) ? -std::numeric_limits<double>::infinity()
                                                            : mag_get_d_log2_approx(offset.value());
    const double outer_log2 = plateau->outer_log2 + std::log2(1 - std::exp2(offset_log2 - plateau->outer_log2));
    if (!std::isfinite(outer_log2))
        return model;
    model.count = plateau->count;
    model.centred = true;
    model.outer_log2 = outer_log2;
    acb_set(model.center.value(), center.value());
    // The hints must reach the certificate's accuracy, not only the goal's:
    // the secular inclusion disks are n times the residues, so the nodes
    // need log2(n) bits beyond the goal.
    const double target_bits = static_cast<double>(options.bits) + std::log2(static_cast<double>(source_->degree())) + 64;
    // The hints' accuracy: the certificate needs the nodes within 2^-(goal +
    // log2 n) of |z| (the inclusion disks are n times the residues), with a
    // margin. (The moments' estimate of the centre gives |z| as well as the mean will.)
    Magnitude center_size;
    acb_get_mag(center_size.value(), center.value());
    if (const acb_struct *offset = policy_.centred_origin_basis ? coordinates_.exact_offset() : nullptr) {
        // |z| measured from the original origin, which the goals are relative to: a cluster near it
        // needs its hints to that nearer scale.
        ComplexBall original;
        Magnitude from_origin;
        acb_add(original.value(), center.value(), offset, ARF_PREC_EXACT);
        acb_get_mag(from_origin.value(), original.value());
        if (!mag_is_zero(from_origin.value()))
            mag_min(center_size.value(), center_size.value(), from_origin.value());
    }
    const double need_log2 = mag_get_d_log2_approx(center_size.value()) - static_cast<double>(options.bits) -
                             std::log2(static_cast<double>(source_->degree())) - heuristics::local_restart_hint_margin_bits;
    // The circle of that accuracy about the centre. A goal that takes a cluster within its tolerance as it
    // stands (every goal but the accuracy goal, which isolates each root) needs nothing of the cluster's shape
    // once all of the cluster lies inside this circle: any distinct points on it are hints as good as the
    // goal can use, and the members' next cover is one component within the tolerance.
    const bool shapeless_goal = engine_goal(options) != CPS_ROOT_ACCURACY && std::isfinite(need_log2);
    Magnitude goal_circle;
    if (shapeless_goal) {
        Float goal_radius;
        set_radius_log2(goal_radius.value(), need_log2);
        arf_get_mag_lower(goal_circle.value(), goal_radius.value());
    }
    bool within_goal = false;
    slong working = static_cast<slong>(std::min<uint64_t>(options.max_work_bits, std::max<uint64_t>(64, maximum_precision_)));
    ComplexPolynomial local;
    Magnitude radius;
    double shrunk_log2 = 0;
    {
        TraceScope trace(&policy_, "geometric_model_test", static_cast<uint64_t>(working), static_cast<uint64_t>(m), &statistics);
        // The shift's precision rises until the coefficients up to the
        // largest model allowed (twice the count: no more extra terms than
        // the cluster's own) carry a correct bit.
        const slong needed = std::min<slong>(source_->degree(), 2 * model.count);
        // The cluster lies within the plateau's inner circle about the nodes' centroid (`mean` holds the
        // centroid since the swap above).
        Float reach;
        set_radius_log2(reach.value(), plateau->inner_log2);
        // The centre is the moments' estimate while the model about it bears the estimate out, and the mean
        // located from it (locate_mean) once the model shows the estimate to be outside its cluster: the test
        // is at the end of this loop, which then runs once more.
        bool relocated = false, found = false;
        Magnitude trial;
        for (;;) {
            for (;;) {
                if (options.cancelled()) return model;
                if (relocated) {
                    // The mean to this precision's bits, from the last one's.
                    double pair_operations = 0;
                    locate_mean(center.value(), *polynomial, model.count, working, mean.value(), reach.value(), options,
                                pair_operations);
        work_ += service_charge(&policy_, ServiceOperation::Auxiliary, pair_operations * operation_weight(working));
                    acb_set(model.center.value(), center.value());
                }
                polynomial->research_service_meter(policy_.research_service_meter);
                {
                    TraceScope service_trace(policy_.research_service_meter ? &policy_ : nullptr,
                        "service_taylor", static_cast<uint64_t>(working), static_cast<uint64_t>(source_->degree()));
                    model.determined = polynomial->taylor(local.value(), center.value(), working);
                }
                // Metered as Subdivision's shifts are (dovetail.hpp), on the
                // shifted polynomial's block structure.
                if (policy_.research_service_meter)
                    work_ += service_charge(&policy_, polynomial->research_last_taylor_operation(),
                                            polynomial->research_last_taylor_work(), working);
                else
                    work_ += block_taylor_shift_operations(local.value(), working) * operation_weight(working);
                const bool shifted = model.determined;
                for (slong j = 0; model.determined && j <= needed; ++j)
                    model.determined = acb_is_zero(local.value()->coeffs + j) || acb_rel_accuracy_bits(local.value()->coeffs + j) > 0;
                // A cluster this precision does not resolve, all of it within the goal's circle of the
                // centre: the count dominates the expansion there (Pellet, with every coefficient this
                // precision leaves undetermined at its upper bound). More precision would resolve what such
                // a goal does not ask for. History:
                // docs/code-history-20261021.md, secular_refinement.cpp 28.
                within_goal = shapeless_goal && shifted && !model.determined && !mag_is_zero(goal_circle.value()) &&
                              dominant_count(local.value(), goal_circle.value()) == model.count;
                if (within_goal || model.determined || working >= static_cast<slong>(options.max_work_bits))
                    break;
                working = std::min<slong>(static_cast<slong>(options.max_work_bits), resources::precision_product(working, 2));
            }
            model.working = working;
            if (within_goal)
                break;
            // The cluster's own scale: circles about the mean halve from D until
            // K dominates the model and while it does (Pellet: K roots in the
            // disk, and by Rouche the degree-K truncation's K roots there too).
            // The smallest such circle is where the model's coefficient spread and
            // the extra terms are least. A dominant count below K means the circle
            // entered the cluster without K dominating: no scale.
            found = false;
            for (slong halving = 1; model.determined && halving <= static_cast<slong>(std::ceil(target_bits)); ++halving) {
                const double trial_log2 = outer_log2 - static_cast<double>(halving);
                mag_set_d(trial.value(), std::exp2(trial_log2 - std::floor(trial_log2)));
                mag_mul_2exp_si(trial.value(), trial.value(), static_cast<slong>(std::floor(trial_log2)));
                const slong dominant = dominant_count(local.value(), trial.value());
                if (dominant == model.count) {
                    found = true;
                    mag_set(radius.value(), trial.value());
                    shrunk_log2 = trial_log2 - outer_log2;
                    continue;
                }
                if (found || (dominant >= 0 && dominant < model.count))
                    break;
            }
            // Whether that scale is the cluster's or the centre's own distance from it. A model wants its
            // centre within its scale, over 4 K, of the roots' mean (component_seeds): coefficient K - 1 is K
            // times the centre's offset from the root of p^(K-1), so four times it against the scale times
            // coefficient K. About a point d from a cluster narrower than d the expansion is (w + d)^K's,
            // the count dominates only beyond d K / ln 2, and that is the scale found, about three times
            // what the test allows. An estimate inside its cluster passes by a wide margin, and its model is
            // the model: the root of p^(K-1) is itself K (cluster radius)^2 over the distance to the other
            // roots from the mean, no nearer it than the moments are in a wide cluster, and a second shift
            // about it buys nothing. History: docs/code-history-20261021.md, secular_refinement.cpp 29.
            if (relocated || about_origin || !found)
                break;
            Magnitude offset_term, scale_term, leading;
            acb_get_mag(offset_term.value(), local.value()->coeffs + model.count - 1);
            mag_mul_2exp_si(offset_term.value(), offset_term.value(), 2);
            acb_get_mag_lower(leading.value(), local.value()->coeffs + model.count);
            mag_mul_lower(scale_term.value(), radius.value(), leading.value());
            if (mag_cmp(offset_term.value(), scale_term.value()) <= 0)
                break;
            relocated = true;
        }
        if (within_goal) {
            // The hints: the cluster's count of points on the goal's circle (candidates.hpp's placement).
            if (policy_.trace) {
                TraceEvent event{"geometric_local_restart", "within_goal", static_cast<uint64_t>(working),
                                 static_cast<uint64_t>(m), &statistics, nullptr};
                event.engine = "secular";
                event.algorithm_cutoff = static_cast<uint64_t>(model.count);
                policy_.trace(policy_.trace_data, event);
            }
            StageTimer timer(statistics.candidate_ns);
            Float goal_radius;
            set_radius_log2(goal_radius.value(), need_log2);
            RealBall rho;
            arb_set_arf(rho.value(), goal_radius.value());
            RootVector placed(model.count);
            circle_candidates(placed.data(), model.count, center.value(), rho.value(), working);
            model.hints.resize(static_cast<size_t>(model.count));
            for (slong j = 0; j < model.count; ++j)
                acb_set(model.hints[static_cast<size_t>(j)].value(), placed.data() + j);
            acb_set(model.center.value(), center.value());
            model.justified = model.count >= 2;
            if (model.justified)
                maximum_precision_ = std::max(maximum_precision_, static_cast<uint64_t>(working));
            return model;
        }
        if (found)
            model.extra = std::max<slong>(1, static_cast<slong>(std::ceil(target_bits / -shrunk_log2)));
        // A model whose extra terms outnumber the cluster's own is not local.
        model.justified = found && model.count >= 2 && model.extra <= model.count &&
                          model.count + model.extra < source_->degree();
    }
    if (policy_.trace) {
        TraceEvent event{"geometric_local_restart", model.justified ? "justified" : (model.determined ? "declined" : "undetermined"),
                         static_cast<uint64_t>(working), static_cast<uint64_t>(m), &statistics, nullptr};
        event.engine = "secular";
        event.algorithm_cutoff = static_cast<uint64_t>(model.count);
        event.certified = static_cast<uint64_t>(model.extra);
        policy_.trace(policy_.trace_data, event);
    }
    if (!model.justified)
        return model;
    TraceScope restart_trace(&policy_, "component_initialization", static_cast<uint64_t>(working), static_cast<uint64_t>(model.count), &statistics);
    StageTimer timer(statistics.candidate_ns);
    // The local model is solved by the library itself, in the cluster's own
    // scale u = w / 2^scale (an exact power of two near the disk's radius):
    // there binary64 resolves the cluster, which the global iteration
    // cannot. Its answers are hints; only the global validation certifies.
    // History: docs/code-history-20261021.md, secular_refinement.cpp 30.
    const slong scale_log2 = static_cast<slong>(std::floor(mag_get_d_log2_approx(radius.value())));
    // (|u| is at most two in the cluster's scale; need_log2, above, is the hints' accuracy.)
    // About the original origin the model's own coordinates share the goals'
    // origin: its relative accuracy is theirs, with the certificate's log2 n
    // and the margin.
    //
    // Elsewhere the hints need the bits between the cluster's scale and that
    // accuracy, and the goal's own at least. A cluster narrower than the
    // accuracy has none between them: the goal's bits at the cluster's scale
    // then, which keep its roots apart. A negative difference counts as none
    // (converted to an unsigned count as it stands, it would be a precision
    // beyond every bound). History:
    // docs/code-history-20261021.md, secular_refinement.cpp 31.
    const double scale_bits = std::ceil(1 + static_cast<double>(scale_log2) - need_log2);
    const uint64_t nested_bits =
        about_origin ? static_cast<uint64_t>(std::ceil(static_cast<double>(options.bits) +
                                                       std::log2(static_cast<double>(source_->degree())) +
                                                       heuristics::local_restart_hint_margin_bits))
                     : std::max<uint64_t>(options.bits, scale_bits > 0 ? static_cast<uint64_t>(scale_bits) : 0);
    struct LocalSolution { bool ok = false; std::vector<ComplexBall> roots; ComplexPolynomial scaled; };
    const auto solve_local = [&](slong extra) {
        LocalSolution solution;
        const auto nested_options = local_model_options(options, nested_bits, search, statistics);
        if (!nested_options)
            return solution;
        const slong k = model.count + extra;
        acb_poly_set(solution.scaled.value(), local.value());
        acb_poly_truncate(solution.scaled.value(), k + 1);
        if (acb_poly_degree(solution.scaled.value()) != k || acb_contains_zero(solution.scaled.value()->coeffs + k))
            return solution;
        DensePolynomial nested;
        nested.coefficients.resize(static_cast<size_t>(k + 1));
        for (slong j = 0; j <= k; ++j) {
            acb_ptr c = solution.scaled.value()->coeffs + j;
            acb_mul_2exp_si(c, c, scale_log2 * j);
            acb_get_mid(c, c);
            for (int part = 0; part < 2; ++part) {
                arf_srcptr x = arb_midref(part ? acb_imagref(c) : acb_realref(c));
                auto &number = part ? nested.coefficients[static_cast<size_t>(j)].imag
                                    : nested.coefficients[static_cast<size_t>(j)].real;
                Integer mantissa;
                arf_get_fmpz_2exp(mantissa.value(), number.exponent.value(), x);
                fmpq_set_fmpz(number.significand.value(), mantissa.value());
            }
        }
        TraceScope trace(&policy_, "local_candidates", nested_bits, static_cast<uint64_t>(k), &statistics);
        // The controller supplies a bounded synchronous solve with recursive
        // local restarts disabled. Its result proves nothing about this source.
        const cps_result nested_result = [&] {
            auto *metrics = policy_.research_scheduling_metrics;
            std::optional<StageTimer> local_timer;
            if (metrics) {
                ++metrics->local_calls;
                local_timer.emplace(metrics->local_ns);
            }
            return local_model_solver_(std::move(nested), *nested_options);
        }();
        if (auto *metrics = policy_.research_scheduling_metrics) {
            metrics->local_goal_calls += nested_result.complete && nested_result.goal_met;
            if (nested_result.research_priced_work_available) {
                const auto &child = nested_result.research_scheduling_metrics;
                metrics->local_subdivision_ns = resources::add(metrics->local_subdivision_ns,
                                                              child.subdivision_ns);
                metrics->local_subdivision_turns = resources::add(metrics->local_subdivision_turns,
                                                                 child.subdivision_turns);
                metrics->local_subdivision_work += child.subdivision_work;
            }
        }
        // Charge the promised call totals even when the model refuses or stops.
        // Arithmetic remains charged below at the child's precision; merging
        // its diagnostic counters here would make the outer step meter count
        // those operations again. Child precision is not outer precision.
        search.work_units = resources::add(search.work_units, nested_result.progress.work_units_this_call);
        statistics.candidate_sweeps = resources::add(statistics.candidate_sweeps,
                                                     nested_result.statistics_v3.candidate_sweeps);
        // A private priced child returns both engines' actual action charges.
        // Neither the final precision nor diagnostic counters are another bill.
        if (nested_result.research_priced_work_available) {
            const auto &child = nested_result.research_scheduling_metrics;
            work_ += child.numerical_work + child.subdivision_work;
        } else {
        // The legacy work is metered at its own precision (its counters ran at the
        // nested precisions; weighed at this refinement's, the charge would
        // starve the iteration). History:
        // docs/code-history-20261021.md, secular_refinement.cpp 32.
        const auto &done = nested_result.statistics_v3;
        work_ += static_cast<double>(done.evaluation_terms + done.root_interactions +
                                     done.validation_calls * static_cast<uint64_t>(k)) *
                 operation_weight(static_cast<int64_t>(nested_result.work_bits));
        }
        if (!nested_result.complete || !nested_result.goal_met)
            return solution;
        for (const auto &root : nested_result.roots) {
            if (root.count != 1)
                return solution; // the model has a cluster at the nested goal: not this restart's case
            ComplexBall u;
            acb_get_mid(u.value(), root.box.value());
            solution.roots.push_back(std::move(u));
        }
        if (static_cast<slong>(solution.roots.size()) != k)
            return solution;
        std::sort(solution.roots.begin(), solution.roots.end(), [](const ComplexBall &a, const ComplexBall &b) {
            Magnitude x, y;
            acb_get_mag(x.value(), a.value());
            acb_get_mag(y.value(), b.value());
            return mag_cmp(x.value(), y.value()) < 0;
        });
        solution.ok = true;
        return solution;
    };
    auto solution = solve_local(model.extra);
    if (!solution.ok) {
        model.justified = false;
        return model;
    }
    // The truncation, amplified by that condition, must still meet the
    // target: a Wilkinson-like cluster (equal spacing on a segment) needs far
    // more terms than the goal alone asks for; with those alone, hints land
    // between roots. History:
    // docs/code-history-20261021.md, secular_refinement.cpp 33.
    const double condition_log2 = local_condition_log2(solution.scaled, solution.roots, model.count);
    const slong needed_extra =
        std::isfinite(condition_log2)
            ? static_cast<slong>(std::ceil((target_bits + condition_log2) / -shrunk_log2))
            : std::numeric_limits<slong>::max();
    if (needed_extra > model.extra) {
        bool determined = needed_extra <= model.count && model.count + needed_extra < source_->degree();
        for (slong j = 0; determined && j <= model.count + needed_extra; ++j)
            determined = acb_is_zero(local.value()->coeffs + j) || acb_rel_accuracy_bits(local.value()->coeffs + j) > 0;
        if (!determined) {
            model.justified = false;
            return model;
        }
        model.extra = needed_extra;
        solution = solve_local(model.extra);
        if (!solution.ok) {
            model.justified = false;
            return model;
        }
    }
    model.packet_bits = static_cast<slong>(nested_bits);
    // The inner `count` roots must lie in the certified disk (|u| at most its
    // radius over 2^scale).
    Magnitude scaled_radius;
    mag_mul_2exp_si(scaled_radius.value(), radius.value(), -scale_log2);
    Magnitude inside;
    for (slong j = 0; j < model.count; ++j) {
        acb_get_mag(inside.value(), solution.roots[static_cast<size_t>(j)].value());
        if (mag_cmp(inside.value(), scaled_radius.value()) > 0) {
            model.justified = false; // the model's inner roots left the disk
            return model;
        }
    }
    acb_set(model.center.value(), center.value());
    model.hints.resize(static_cast<size_t>(model.count));
    for (slong j = 0; j < model.count; ++j) {
        acb_mul_2exp_si(model.hints[static_cast<size_t>(j)].value(), solution.roots[static_cast<size_t>(j)].value(), scale_log2);
        acb_add(model.hints[static_cast<size_t>(j)].value(), model.hints[static_cast<size_t>(j)].value(), center.value(), working);
        acb_get_mid(model.hints[static_cast<size_t>(j)].value(), model.hints[static_cast<size_t>(j)].value());
    }
    maximum_precision_ = std::max(maximum_precision_, static_cast<uint64_t>(working));
    if (policy_.research_scheduling_metrics)
        ++policy_.research_scheduling_metrics->local_accepted_calls;
    return model;
}
// A restart component's members re-initialized from its local model (component_seeds), at the round's
// iteration precision `prec`.
SecularRefinement::LocalRestart SecularRefinement::local_component_restart(const Root &component,
                                                                           const std::vector<size_t> &members,
                                                                           slong prec, const cps_options &options,
                                                                           certpolysolve::detail::CandidateStatistics &statistics,
                                                                           certpolysolve::detail::SearchStatistics &search) {
    if (!local_restart(members.size()))
        return LocalRestart::None;
    auto *polynomial = polynomial_;
    const bool interpolated_model = conditioned_secular(*source_, &policy_);
    bool local_initialized = false;
    {
        TraceScope trace(&policy_, "component_initialization", maximum_precision_, members.size(),
                         &statistics, &search);
        StageTimer timer(statistics.candidate_ns);
        slong working =
            std::min<slong>(static_cast<slong>(options.max_work_bits),
                            std::max<slong>(prec, static_cast<slong>(maximum_precision_)));
        const slong k = static_cast<slong>(members.size());
        RootVector hints(k);
        // What the hints are placed for (component_seeds' resolution). Under the accuracy goal every
        // root of the component is to be isolated, and the members are iterated to the round's
        // resolution. Under the other goals a component whose members lie within the goal's tolerance
        // of its centre is certified as it stands: a node's inclusion disk, n times its residue, meets
        // the goal when the node is within about 2^-(goal + log2 n) of its modulus from its roots
        // (heuristics::local_restart_hint_margin_bits has the margin), and hints placed finer than
        // that buy nothing. History: docs/code-history-20261021.md, secular_refinement.cpp 34.
        const slong goal_resolution = static_cast<slong>(std::ceil(
            static_cast<double>(options.bits) + std::log2(static_cast<double>(source_->degree())) +
            heuristics::local_restart_hint_margin_bits));
        // That holds for members iterated on the secular model, which keep the round's precision (below).
        // Members iterated on the source take the model's own precision for their iteration, so no model
        // puts its hints within that: raised for the round's resolution, a model of k roots would make
        // the members' iteration about k times as precise as the round. Their models are asked for no
        // resolution and keep the regeneration's precision. History:
        // docs/code-history-20261021.md, secular_refinement.cpp 35.
        const bool model_iteration = interpolated_model || policy_.refinement == RefinementMethod::Secular;
        const slong resolution =
            !model_iteration ? 0 : engine_goal(options) == CPS_ROOT_ACCURACY ? prec : std::min(prec, goal_resolution);
        // The model's precision. A model of a k-root component with noise 2^-N of its leading coefficient
        // places its hints 2^(-N / k) out (component_seeds), so k bits for each bit of the resolution, and the
        // bits the component's evaluation cancels, put them inside it; a finer model buys nothing the members'
        // iteration can use. The regeneration's precision can be several times that. The cancelled bits are
        // what the solve's models so far had to be raised by (restart_need_): the first model starts with
        // none, and component_seeds raises it by the excess it measures.
        //
        // The resolution it is sized for is the one its hints end at: the round's under the accuracy goal,
        // whose next round restarts the component again, and the goal's where the goal accepts a component
        // within its tolerance, the round being a stage on the way to it. Sized for the round there, a model
        // can stop short of what resolves its cluster, and the members then take the later rounds one by one.
        // And the goal is the one the caller asked for: while the default goal is staged first (step), a model
        // sized for that stage leaves the cluster certified within the default tolerance only, where the
        // certificate attempt at the caller's goal cannot start from it and the refinement has to resolve it.
        // History: docs/code-history-20261021.md, secular_refinement.cpp 36.
        const slong caller_resolution = static_cast<slong>(std::ceil(
            static_cast<double>(std::max(options.bits, caller_bits_)) + std::log2(static_cast<double>(source_->degree())) +
            heuristics::local_restart_hint_margin_bits));
        //
        // Where that sizing is below the regeneration's precision the model takes it, and where it is above, the
        // model first takes the regeneration's: a component whose roots that precision already tells apart
        // needs no more, and component_seeds raises a model that does. A raise is a second model, and a model
        // is the restart's cost where it is a whole shift (a source with no dense head). What a component's
        // own models have shown of it decides where its next one starts:
        //
        // - A component whose model had to be raised and which is restarted again is still unresolved, and its
        //   next model needs what the last one did: k bits for each bit of the hints' resolution and the bits
        //   the raise showed beyond them, which a cluster keeps from one growth to the next. Its members carry
        //   those bits (restart_need_), and a component all of whose members carry them starts there.
        // - Under the accuracy goal, where the sizing is what an unresolved component's hints need: a
        //   component whose last model did not determine its hints (restart_unresolved_) and which is
        //   restarted again has roots closer than that model told apart, and starts at its sizing, with the
        //   most bits beyond k a bit that a raise of the solve has shown.
        // - Every other model starts at the lower of the regeneration's precision and its sizing, with those
        //   bits. Of a component with no such record nothing says the lower precision fails.
        //
        // The bits beyond k a bit are the solve's, not a growth's: forgotten at each growth, a cluster's next
        // model can start without them, coarser than its model of the growth before, and a cluster's models
        // must not get coarser.
        //
        // Not every later model at its sizing once any model of the solve has been raised: under a goal that
        // accepts a cluster within its tolerance the sizing is that goal's, a bound and not a need, and the
        // models of clusters the regeneration's precision resolves would take several times the bits, and
        // their members' residues with them. History: docs/code-history-20261021.md, secular_refinement.cpp 37.
        const slong sized = engine_goal(options) == CPS_ROOT_ACCURACY ? prec : caller_resolution;
        if (model_iteration) {
            slong need = 0, shown = 0;
            bool raised_before = restart_need_.size() == proposals_.size();
            bool unresolved_before =
                engine_goal(options) == CPS_ROOT_ACCURACY && restart_unresolved_.size() == proposals_.size();
            for (size_t j = 0; j < members.size(); ++j) {
                raised_before = raised_before && restart_need_[members[j]] > 0;
                unresolved_before = unresolved_before && restart_unresolved_[members[j]];
                if (raised_before)
                    need = std::max(need, restart_need_[members[j]]);
            }
            for (const auto marked : restart_need_)
                shown = std::max(shown, marked);
            const auto sizing = [&](slong a_bit, slong beyond) {
                return std::max(prec, resources::whole_limbs(resources::precision_sum(
                                          resources::precision_product(k, a_bit), beyond)));
            };
            const slong ceiling = static_cast<slong>(options.max_work_bits);
            if (raised_before)
                working = std::min(ceiling, sizing(resolution, need));
            else if (unresolved_before)
                working = std::min(ceiling, sizing(sized, shown));
            else
                working = std::min(working, sizing(sized, shown));
        }
        const slong asked = working;
        bool unresolved = false;
        if (policy_.research_local_packet_meter && research_packet_started_)
            throw std::logic_error("local packet accounting is only defined before the ordinary packet");
        ResearchLocalPacketAccounting accounting{work_, research_prepaid_local_terms_,
            policy_.research_local_packet_service, policy_.research_local_packet_data, &policy_};
        local_initialized = component_seeds(hints.data(), *polynomial, component, k, working, options, statistics,
                                            policy_, std::nullopt, resolution, &unresolved, caller_resolution
                                            , policy_.research_local_packet_meter ? &accounting : nullptr
                                            );
        if (model_iteration && local_initialized) {
            if (working > asked) {
                if (restart_need_.size() != proposals_.size())
                    restart_need_.assign(proposals_.size(), 0);
                for (const auto i : members)
                    restart_need_[i] = std::max<slong>(1, working - resources::precision_product(k, resolution));
            }
            if (restart_unresolved_.size() != proposals_.size())
                restart_unresolved_.assign(proposals_.size(), 0);
            for (const auto i : members)
                restart_unresolved_[i] = unresolved;
        }
        // A failed attempt's escalation says nothing about the
        // precision the global iteration needs.
        if (local_initialized)
            maximum_precision_ = std::max(maximum_precision_, static_cast<uint64_t>(working));
        if (options.cancelled())
            return LocalRestart::Cancelled;
        if (local_initialized)
            for (size_t j = 0; j < members.size(); ++j) {
                acb_set(proposals_[members[j]].value(), hints.data() + static_cast<slong>(j));
                // Accurate coefficient preparation does not
                // determine the coordinate update precision: a
                // member iterated on the secular model keeps the
                // round's precision, not the local model's,
                // which is the source precision. History:
                // docs/code-history-20261021.md,
                // secular_refinement.cpp 38.
                iteration_bits_[members[j]] = model_iteration ? prec : working;
            }
        // The members' next residues. At the model's precision the source's value at a hint is about its own
        // rounding (the hints lie where the model's unresolved coefficients end), so a residue adequate for the
        // members' iteration needs that precision and the accuracy asked of it (adequate_accuracy_bits), with
        // the guard the regeneration's own first precision has. The regeneration starts there, not at the
        // iteration's bits, from which its doubling can pass the need by nearly a factor of two. (Members
        // iterated on the source already start at the model's precision: their iteration's bits are the
        // model's.) History: docs/code-history-20261021.md, secular_refinement.cpp 39.
        if (model_iteration && local_initialized && source_precision_hints_.size() == proposals_.size()) {
            const bool iteration_model =
                interpolated_model || policy_.model_stopping == SecularModelPolicy::Continuation;
            const slong hint = std::min<slong>(
                static_cast<slong>(options.max_work_bits),
                resources::precision_sum(working, resources::precision_sum(
                                                      adequate_accuracy_bits(prec, options.bits, iteration_model),
                                                      regeneration_guard_bits)));
            if (restart_estimated_.size() != proposals_.size())
                restart_estimated_.assign(proposals_.size(), 0);
            for (const auto i : members) {
                source_precision_hints_[i] = std::max(source_precision_hints_[i], hint);
                restart_estimated_[i] = 1;
            }
        }
    }
    return local_initialized ? LocalRestart::Initialized : LocalRestart::None;
}
// A restart component with more members than a local model takes. The certified cover makes one component of
// everything a wide inclusion disk reaches, and a node that serves no root has such a disk: when the binary64
// stage ends with a node too many on one cluster and one too few on another, the surplus cluster's disks (n
// times a residue) can cover every root, and the cover is one component of all the nodes. Re-seeding that
// component on one circle discards every converged node with them. History:
// docs/code-history-20261021.md, secular_refinement.cpp 40.
//
// The model shows where the nodes stand. A tight group of them (node_groups) whose disk holds as many roots as
// it has nodes, by an exact count (counted_group), is a component of its own in everything but the cover: it
// is restarted from its local model as a certified component of that size is (the count's square for its box).
// A group whose disk holds fewer roots than nodes releases its farthest nodes, and one that holds more takes
// released ones, as the geometric restart redistributes; a group made whole so is restarted as well, and the
// nodes of one that is not keep their places, as do the nodes in no counted group. Smaller groups first, so a
// cluster is restarted by itself and not with its neighbours; a miscounted group is taken at its largest, so
// a cluster's own sub-pairs do not stand for it.
//
// Returns whether a group with as many roots as nodes was found: the component then has structure the circle
// would destroy, and is not re-seeded. With none (nodes collapsed into one neighbourhood: the case the circle
// is for) nothing is changed.
bool SecularRefinement::regroup_restart(const std::vector<size_t> &members, slong prec, const cps_options &options,
                                        certpolysolve::detail::CandidateStatistics &statistics, certpolysolve::detail::SearchStatistics &search) {
    auto *polynomial = polynomial_;
    if (!polynomial ||
        members.size() <= static_cast<size_t>(heuristics::local_component_degree) ||
        residues_.size() != nodes_.size())
        return false;
    const size_t n = nodes_.size();
    std::vector<uint8_t> in_component(n, 0), taken(n, 0), member(n, 0);
    for (const auto i : members)
        in_component[i] = 1;
    auto groups = node_groups();
    std::erase_if(groups, [&](const std::vector<size_t> &group) {
        return group.size() > static_cast<size_t>(heuristics::local_component_degree) ||
               !std::all_of(group.begin(), group.end(), [&](size_t i) { return in_component[i] != 0; });
    });
    std::stable_sort(groups.begin(), groups.end(), [](const auto &a, const auto &b) { return a.size() < b.size(); });
    struct Planned {
        std::vector<size_t> members;
        CountedGroup counted;
    };
    std::vector<Planned> whole, miscounted;
    DiskCountScratch scratch;
    const slong count_precision = static_cast<slong>(maximum_precision_);
    for (const auto &group : groups) {
        if (options.cancelled())
            return false;
        if (std::any_of(group.begin(), group.end(), [&](size_t i) { return taken[i] != 0; }))
            continue;
        auto counted = counted_group(group, nodes_, residues_, member, count_precision, scratch);
        if (!counted)
            continue;
        if (counted->count == static_cast<slong>(group.size())) {
            for (const auto i : group)
                taken[i] = 1;
            whole.push_back({group, std::move(*counted)});
        } else
            miscounted.push_back({group, std::move(*counted)});
    }
    // Metered at one limb, as the geometric restart's counts are.
        work_ += service_charge(&policy_, ServiceOperation::Auxiliary, scratch.operations * operation_weight(FLINT_BITS));
    if (whole.empty())
        return false;
    // The miscounted groups free of the whole ones, the largest first; then the redistribution by counts.
    std::vector<Planned> balanced;
    std::vector<size_t> released;
    std::stable_sort(miscounted.begin(), miscounted.end(),
                     [](const Planned &a, const Planned &b) { return a.members.size() > b.members.size(); });
    std::vector<Planned> short_of;
    for (auto &plan : miscounted) {
        if (std::any_of(plan.members.begin(), plan.members.end(), [&](size_t i) { return taken[i] != 0; }))
            continue;
        for (const auto i : plan.members)
            taken[i] = 1;
        const slong surplus = static_cast<slong>(plan.members.size()) - plan.counted.count;
        if (surplus < 0) {
            short_of.push_back(std::move(plan));
            continue;
        }
        // Release the nodes farthest from the centre.
        std::vector<std::pair<double, size_t>> far;
        ComplexBall difference;
        Magnitude distance;
        for (const auto i : plan.members) {
            acb_sub(difference.value(), nodes_[i].value(), plan.counted.centre.value(), count_precision);
            acb_get_mag(distance.value(), difference.value());
            far.emplace_back(mag_is_zero(distance.value()) ? -std::numeric_limits<double>::infinity()
                                                           : mag_get_d_log2_approx(distance.value()),
                             i);
        }
        std::sort(far.begin(), far.end(), std::greater<>());
        for (slong j = 0; j < surplus; ++j) {
            released.push_back(far[static_cast<size_t>(j)].second);
            plan.members.erase(std::find(plan.members.begin(), plan.members.end(), far[static_cast<size_t>(j)].second));
        }
        balanced.push_back(std::move(plan));
    }
    size_t moved = 0;
    for (auto &plan : short_of) {
        const size_t deficit = static_cast<size_t>(plan.counted.count) - plan.members.size();
        if (deficit > released.size())
            continue; // not enough released nodes for its roots: its nodes keep their places
        for (size_t j = 0; j < deficit; ++j) {
            plan.members.push_back(released.back());
            released.pop_back();
            ++moved;
        }
        balanced.push_back(std::move(plan));
    }
    size_t restarted = 0;
    const auto restart = [&](const Planned &plan) {
        if (plan.members.size() < 2)
            return true;
        Root component;
        acb_set(component.box.value(), plan.counted.centre.value());
        acb_add_error_arf(component.box.value(), arb_midref(plan.counted.radius.value()));
        component.count = static_cast<uint64_t>(plan.members.size()) * multiplicity_;
        component.multiplicity = 0;
        switch (local_component_restart(component, plan.members, prec, options, statistics, search)) {
        case LocalRestart::Cancelled:
            return false;
        case LocalRestart::Initialized:
            ++restarted;
            break;
        case LocalRestart::None:
            break;
        }
        return true;
    };
    for (const auto &plan : whole)
        if (!restart(plan))
            return false;
    for (const auto &plan : balanced)
        if (!restart(plan))
            return false;
    if (policy_.trace) {
        TraceEvent event{"secular_component_regroup", "groups", static_cast<uint64_t>(prec), members.size(), &statistics,
                         &search};
        event.engine = "secular";
        event.certified = restarted;
        event.algorithm_cutoff = moved;
        policy_.trace(policy_.trace_data, event);
    }
    return true;
}
bool SecularRefinement::geometric_restart(const cps_options &options, SearchStatistics &search,
                                          CandidateStatistics &statistics) {
    auto *polynomial = polynomial_;
    if (!polynomial)
        return false;
    // Largest groups first: a group's local model resolves its sub-groups.
    auto groups = node_groups();
    std::sort(groups.begin(), groups.end(), [](const auto &a, const auto &b) { return a.size() > b.size(); });
    std::vector<uint8_t> taken(proposals_.size(), 0);
    // Groups whose model is justified, with their members and the model; a
    // group whose disk holds fewer roots than it has nodes releases its
    // farthest nodes, one that holds more takes released or ungrouped ones.
    struct Planned { std::vector<size_t> members; GroupModel model; ComplexBall center; };
    std::vector<Planned> planned;
    std::map<std::vector<size_t>, std::pair<double, double>> observed;
    for (const auto &members : groups) {
        if (options.heuristic_budget && options.heuristic_budget->exhausted()) break;
        if (options.cancelled())
            return false;
        Magnitude spread;
        ComplexBall difference;
        for (size_t a = 1; a < members.size(); ++a) {
            acb_sub(difference.value(), nodes_[members[a]].value(), nodes_[members[0]].value(), 64);
            Magnitude d;
            acb_get_mag(d.value(), difference.value());
            if (mag_cmp(d.value(), spread.value()) > 0) mag_set(spread.value(), d.value());
        }
        if (mag_is_zero(spread.value()) || !mag_is_finite(spread.value()))
            continue;
        const double spread_log2 = mag_get_d_log2_approx(spread.value());
        // Only a crawling group: linear convergence onto a cluster, the spread
        // shrinking by a steady factor between one half and one, round after
        // round. The same members in the last two validations, each shrink
        // positive and below a bit, the second at least half the first. A
        // group converging quadratically, or whose nodes sit at separate roots
        // or are done (a constant spread), is the iteration's; qualifying it
        // (counts, moments, a Taylor shift) would only cost. History:
        // docs/code-history-20261021.md, secular_refinement.cpp 41.
        const auto seen = group_spreads_.find(members);
        const double shrink = seen == group_spreads_.end() ? std::numeric_limits<double>::quiet_NaN()
                                                           : seen->second.first - spread_log2;
        observed[members] = {spread_log2, shrink};
        if (std::any_of(members.begin(), members.end(), [&](size_t i) { return taken[i] != 0; }))
            continue;
        // A group of more than half the nodes is not a local problem: its
        // model's degree would approach the source's. History:
        // docs/code-history-20261021.md, secular_refinement.cpp 42.
        if (2 * members.size() > nodes_.size())
            continue;
        const double previous_shrink = seen == group_spreads_.end() ? std::numeric_limits<double>::quiet_NaN()
                                                                    : seen->second.second;
        if (!(shrink > 0 && shrink < 1 && previous_shrink > 0 && shrink >= previous_shrink / 2))
            continue;
        const auto attempted = geometric_attempts_.find(members);
        if (attempted != geometric_attempts_.end() && spread_log2 > attempted->second - 1)
            continue;
        // Declined attempts are pure cost: they get at most the share of the
        // work the policy gives Subdivision. History:
        // docs/code-history-20261021.md, secular_refinement.cpp 43.
        if (std::ldexp(declined_work_, static_cast<int>(policy_.dovetail_share_log2)) > work_)
            continue;
        geometric_attempts_[members] = spread_log2;
        const double before = work_;
        auto model = group_model(members, options, search, statistics);
        if (!model.justified) {
            declined_work_ += work_ - before;
            if (model.centred)
                declined_clusters_.push_back({model.center, model.outer_log2, model.count});
            continue;
        }
        for (const size_t i : members) taken[i] = 1;
        Planned plan{members, std::move(model), {}};
        acb_set(plan.center.value(), plan.model.center.value());
        planned.push_back(std::move(plan));
    }
    group_spreads_ = std::move(observed);
    if (planned.empty())
        return false;
    // Redistribution by certified counts.
    std::vector<size_t> released;
    for (auto &plan : planned) {
        const slong surplus = static_cast<slong>(plan.members.size()) - plan.model.count;
        if (surplus <= 0) continue;
        // Release the nodes farthest from the centre.
        std::vector<std::pair<double, size_t>> far;
        ComplexBall difference;
        for (const size_t i : plan.members) {
            acb_sub(difference.value(), nodes_[i].value(), plan.center.value(), 64);
            Magnitude d;
            acb_get_mag(d.value(), difference.value());
            far.emplace_back(-mag_get_d_log2_approx(d.value()), i);
        }
        std::sort(far.begin(), far.end());
        for (slong j = 0; j < surplus; ++j) {
            released.push_back(far[static_cast<size_t>(j)].second);
            plan.members.erase(std::find(plan.members.begin(), plan.members.end(), far[static_cast<size_t>(j)].second));
        }
    }
    // Ungrouped active nodes are the second pool.
    std::vector<size_t> ungrouped;
    for (size_t i = 0; i < proposals_.size(); ++i)
        if (!certified_.inactive[i] && !taken[i])
            ungrouped.push_back(i);
    bool any = false;
    for (auto &plan : planned) {
        slong deficit = plan.model.count - static_cast<slong>(plan.members.size());
        while (deficit > 0 && !released.empty()) {
            plan.members.push_back(released.back());
            released.pop_back();
            --deficit;
        }
        while (deficit > 0 && !ungrouped.empty()) {
            // The ungrouped node nearest this centre.
            size_t best = 0;
            double best_log2 = std::numeric_limits<double>::infinity();
            ComplexBall difference;
            for (size_t u = 0; u < ungrouped.size(); ++u) {
                acb_sub(difference.value(), nodes_[ungrouped[u]].value(), plan.center.value(), 64);
                Magnitude d;
                acb_get_mag(d.value(), difference.value());
                const double v = mag_get_d_log2_approx(d.value());
                if (v < best_log2) { best_log2 = v; best = u; }
            }
            plan.members.push_back(ungrouped[best]);
            ungrouped.erase(ungrouped.begin() + static_cast<ptrdiff_t>(best));
            --deficit;
        }
        if (deficit > 0)
            continue; // not enough nodes to take the model's roots: leave the group alone
        // The members iterate from the hints, so their coordinates must keep
        // the hints apart: the precision the iteration gives a goal of the
        // hints' span (hint_span_log2), its guard bits included: at a coarser
        // precision a cluster's hints can round to one point. History:
        // docs/code-history-20261021.md, secular_refinement.cpp 44.
        const double span_log2 = hint_span_log2(plan.model.hints, plan.center.value(), plan.model.working);
        const slong separation_bits = std::min<slong>(
            static_cast<slong>(options.max_work_bits),
            SecularMachineWorkspace::stage_bits(static_cast<uint64_t>(std::ceil(span_log2))));
        for (size_t j = 0; j < plan.members.size(); ++j) {
            acb_set(proposals_[plan.members[j]].value(), plan.model.hints[j].value());
            iteration_bits_[plan.members[j]] = std::max(iteration_bits_[plan.members[j]], separation_bits);
        }
        if (policy_.study_statistics) {
            ++policy_.study_statistics->candidate_restarts;
            policy_.study_statistics->candidates_replaced += plan.members.size();
        }
        if (policy_.trace) {
            TraceEvent event{"geometric_local_restart", "applied", static_cast<uint64_t>(plan.model.working),
                             static_cast<uint64_t>(plan.members.size()), &statistics, nullptr};
            event.engine = "secular";
            event.algorithm_cutoff = static_cast<uint64_t>(plan.model.count);
            event.arithmetic_bits = static_cast<uint64_t>(separation_bits);
            policy_.trace(policy_.trace_data, event);
        }
        any = true;
    }
    // Released nodes with no taker: spread them on the ungrouped region?
    // They keep their positions; the next rounds move them.
    return any;
}
} // namespace certpolysolve::detail
