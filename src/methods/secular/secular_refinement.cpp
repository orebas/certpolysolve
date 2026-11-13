// SecularRefinement takes one monic source from its seeds to a certified cover of its roots over repeated calls of
// step(); cover() exports the cover in the caller's coordinates (RootCoordinates). Its model is the secular equation
// S(z) = 1 - sum_j w_j/(z - b_j) = 0 at exact nodes b_i, residues w_i = -p(b_i)/Q'(b_i), Q = prod (z - b_j): p = Q S.
//
// The controller lives here (state, step, attempts, preliminary passes and cover retention), with
// secular_model.cpp owning model comparison, refresh and validation scheduling, and secular_growth.cpp
// owning precision growth and resumable candidate restarts. secular_packets.cpp owns packet dispatch,
// proposal adoption and packet exhaustion. Numerical helpers live in
// secular_regeneration.cpp (residues), secular_validation.cpp (coordinates and inclusion covers),
// and secular_restarts.cpp (component seeds and restarts); shared details are in secular_internal.hpp.
//
// step() in its order, each optional stage where the policy (EnginePolicy) asks: a new goal clears the certified flags
// but for kept regions that meet it; once per goal, the certificate attempt (newton_cover, by machine_newton_cover),
// whose refusal keeps what a pass over every point certified (keep_certificates) and, for a goal above the attempt's
// first precision, stages the library's default goal first (staged_); once, preliminary (an Aberth packet, then
// FLINT's inclusion test); a precision growth (iteration bits doubled except at the input's floor, unsteered_nodes)
// with the restarts, or where nothing can grow CPS_STOP_PRECISION or CPS_STOP_INPUT; regeneration and validation where
// the model changed; geometric_restart; and a packet of Aberth sweeps whose Newton corrections come from the model
// (secular_candidate_packet; SecularMachineWorkspace in binary64, alone or in mixed rounds) or from p (direct).
//
// regenerate_residues computes each w_i, with p(b_i) from DiskValues or the source, at a precision raised from its
// iteration bits plus regeneration_guard_bits until it is adequate (adequate) or the budget or the input stops it;
// when few nodes moved, an unmoved node's residue is multiplied by prod (b_i - old b_j)/(b_i - b_j) over the moved
// ones. secular_inclusion_cover merges the discs D(b_i, n|w_i|) to a fixed point: they contain the Gershgorin discs of
// diag(b) + w 1^T, whose characteristic polynomial is p, so a component holds as many roots as discs. It proves a lone
// root real by conjugation, takes kept regions, and splits a component that fails the goal by exact counts
// (split_by_exact_counts; secular_disk_count, by the argument principle on S without evaluating p).
//
// The restarts place a component's nodes again: where every member's residue lacks a correct bit, on the circle that
// minimises its certified radius (certificate_reseed); with at most heuristics::local_component_degree members, from
// its local model (component_seeds, sized by local_component_restart): Newton's iteration on p^(k-1) for the centre,
// then the roots of the degree-k Taylor model by an Aberth packet (place_undetermined_hints spreads those the model
// does not tell apart), or on its Newton polygon's circles where undetermined (scale_candidates); with more, by its
// counted node groups (regroup_restart); else on a circle. geometric_restart takes a crawling node group (node_groups)
// to group_model: its exact root count from an isolating plateau (isolating_plateau), the roots' mean (root_mean,
// locate_mean), Pellet's test on the Taylor shift (dominant_count), and the truncation solved by a nested SolveState.
//
// State: SecularRefinement, for the source's life: the points and the model (proposals_, nodes_, residues_,
// iteration_bits_), certified_ (cover, restart cover, inactive flags), kept_roots_ and kept_, the regeneration's and
// restarts' records (source_precision_hints_, restart_need_, group_spreads_, declined_clusters_), growth_stage_,
// machine_, and work_ (binary64-equivalent operations, dovetail.hpp); RootCoordinates, a centred source's exact offset.
// Entry points: SecularRefinement (step, cover, precision, work, active, rounds, last_correction_log2), built and
// stepped by solver.cpp's factor_step (a PolynomialSource, the only source the attempt, preliminary and local models
// take), native_advance (SecularSource) and program_advance (ProgramSource); RootCoordinates, by factor_step;
// secular_disk_count, by secular_geometry.cpp; the header's other functions have no caller in another src/ file.
// Not here: the attempt (machine_newton_cover, certificate_attempt.cpp), the sweeps (candidates.cpp,
// study_candidates.cpp), plateaus and means (secular_geometry.cpp), the sources (structured_source.cpp,
// program_source.cpp), the cover algebra (certification.cpp), the hand-off to Subdivision (solver.cpp).
// Design notes: docs/mathematics.md ("Persistent interpolation certificates"), docs/session-state-inventory-20261017.md
// (the members), docs/method-composition-design-20260928.md (sections 1 and 2: the attempt first; regrouping).
#include "methods/secular/secular_refinement.hpp"
#include "core/query.hpp"
#include "methods/secular/secular_internal.hpp"
#include "methods/secular/secular_geometry.hpp"
#include "numeric/evaluate/disk_blocked_polynomial_evaluator.hpp"
#include "numeric/evaluate/disk_sparse_polynomial_evaluator.hpp"
#include "core/dovetail.hpp"
#include "core/service_meter.hpp"
#include "numeric/evaluate/limb_horner_evaluator.hpp"
#include "numeric/evaluate/node_products.hpp"
#include "certify/attempts/newton_cover.hpp"
#include "core/engine_policy.hpp"
#include "core/engine_statistics.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <exception>

namespace certpolysolve::detail {
bool prepared_source_policy(const EnginePolicy *policy) {
    return policy && policy->prepared_source_refinement &&
           policy->candidates == CandidateEngine::Adaptive &&
           policy->study == StudyMethod::None && !policy->explicit_source_policy();
}
bool conditioned_secular(const MonicSource &source, const EnginePolicy *policy) {
    const auto *polynomial = dynamic_cast<const PolynomialSource *>(&source);
    // Exact translated models can resolve cancellation for regeneration.
    // Reuse that capability while moving ordinary updates to the better
    // conditioned interpolating representation. Explicit studies retain
    // their selected method and evaluation semantics.
    return prepared_source_policy(policy) &&
           polynomial && polynomial->direct() && polynomial->exact_local_models() &&
           !polynomial->local();
}
SecularRefinement::SecularRefinement(std::unique_ptr<MonicSource> source, acb_srcptr seeds,
                                     uint64_t multiplicity, slong initial_precision, EnginePolicy policy,
                                     LocalModelSolver local_model_solver, RootCoordinates coordinates)
    : source_(std::move(source)), polynomial_(dynamic_cast<PolynomialSource *>(source_.get())),
      coordinates_(std::move(coordinates)),
      proposals_(static_cast<size_t>(source_->degree())),
      iteration_bits_(proposals_.size(), initial_precision), certified_(proposals_.size()),
      multiplicity_(multiplicity), policy_(policy), local_model_solver_(local_model_solver),
      kept_roots_(proposals_.size()), kept_(proposals_.size(), 0) {
    if (polynomial_) polynomial_->research_service_meter(policy_.research_service_meter);
    for (size_t i = 0; i < proposals_.size(); ++i)
        acb_get_mid(proposals_[i].value(), seeds + i);
    missed_growths_.assign(proposals_.size(), 0);
}
size_t SecularRefinement::active() const {
    return static_cast<size_t>(std::count(certified_.inactive.begin(), certified_.inactive.end(), uint8_t(0)));
}
void SecularRefinement::snapshot(const char *stage) const {
    if (!policy_.snapshot) return;
    NumericalSnapshot event;
    event.stage = stage;
    event.factor = policy_.snapshot_factor;
    event.power = policy_.snapshot_power;
    event.multiplicity = multiplicity_;
    event.offset = coordinates_.offset();
    event.revision = rounds_;
    event.source_precision = maximum_precision_;
    event.points = proposals_;
    event.nodes = nodes_;
    event.residues = residues_;
    event.coordinate_precision = iteration_bits_;
    event.inactive = certified_.inactive;
    policy_.snapshot(std::move(event));
}
bool SecularRefinement::cover(std::vector<Root> &out) const {
    if (!certified_.available)
        return false;
    if (coordinates_.identity()) {
        out = certified_.cover;
        return true;
    }
    std::vector<Root> exported;
    exported.reserve(certified_.cover.size());
    for (const auto &root : certified_.cover)
        exported.push_back(coordinates_.original(root, static_cast<slong>(maximum_precision_)));
    if (!disjoint(exported))
        return false;
    out = std::move(exported);
    return true;
}
void SecularRefinement::CertifiedState::consider(std::vector<Root> proposed,
                                                std::vector<uint8_t> proposed_inactive,
                                                slong precision) {
    const bool completed = std::all_of(proposed_inactive.begin(), proposed_inactive.end(),
                                      [](uint8_t value) { return value != 0; });
    if (!completed && available && !prefer_cover(proposed, cover, precision))
        return;
    // The vector swaps below cannot throw: the previous state is replaced whole or not at all.
    cover.swap(proposed);
    inactive.swap(proposed_inactive);
    available = true;
}
void separate_interpolation_nodes(std::span<ComplexBall> proposals, std::span<const uint8_t> inactive,
                                  std::span<const slong> iteration_bits) {
    if (proposals.size() != inactive.size() || proposals.size() != iteration_bits.size())
        throw std::invalid_argument("interpolation candidate state sizes differ");
    // Equal proposals are found by sorting (value, then index), not by comparing every pair, which would cost
    // O(n^2) every round however few proposals are equal. Within each group of equal proposals the pairs are
    // visited in this order: the later index first, then each earlier one. A moved proposal that lands exactly
    // on another one is left there (a proposal-only matter). History:
    // docs/code-history-20261021.md, secular_refinement.cpp 17.
    const auto compare_ball = [](arb_srcptr a, arb_srcptr b) {
        if (const int middle = arf_cmp(arb_midref(a), arb_midref(b)))
            return middle;
        return mag_cmp(arb_radref(a), arb_radref(b));
    };
    std::vector<size_t> order(proposals.size());
    std::iota(order.begin(), order.end(), size_t{0});
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        const auto x = proposals[a].value(), y = proposals[b].value();
        if (const int real = compare_ball(acb_realref(x), acb_realref(y)))
            return real < 0;
        if (const int imag = compare_ball(acb_imagref(x), acb_imagref(y)))
            return imag < 0;
        return a < b;
    });
    for (size_t start = 0, end = 0; start < order.size(); start = end) {
        for (end = start + 1;
             end < order.size() && acb_equal(proposals[order[start]].value(), proposals[order[end]].value()); ++end) {
        }
        for (size_t later = start + 1; later < end; ++later)
            for (size_t earlier = start; earlier < later; ++earlier) {
                const size_t i = order[later], j = order[earlier];
                if (!acb_equal(proposals[i].value(), proposals[j].value()))
                    continue;
                // Frozen indices are owned by the retained certificate. Move
                // an active duplicate, including when its index comes first.
                const auto moving = inactive[i] ? j : i;
                if (inactive[moving])
                    continue;
                RealBall shift;
                acb_abs(shift.value(), proposals[moving].value(), 64);
                arb_add_ui(shift.value(), shift.value(), 1, 64);
                arb_get_mid_arb(shift.value(), shift.value());
                arb_mul_2exp_si(shift.value(), shift.value(), -iteration_bits[moving] + 8);
                arb_mul_ui(shift.value(), shift.value(), static_cast<ulong>(moving + 1), iteration_bits[moving]);
                arb_add(acb_imagref(proposals[moving].value()), acb_imagref(proposals[moving].value()), shift.value(),
                        iteration_bits[moving]);
                acb_get_mid(proposals[moving].value(), proposals[moving].value());
            }
    }
}
double ball_horner_loss_bits(acb_srcptr points, slong count, slong degree) {
    double worst = 0;
    Float scaled;
    for (slong i = 0; i < count; ++i) {
        const arf_struct *re = arb_midref(acb_realref(points + i)), *im = arb_midref(acb_imagref(points + i));
        if (!arf_is_finite(re) || !arf_is_finite(im) || (arf_is_zero(re) && arf_is_zero(im)))
            continue;
        // The two parts at the larger one's scale, as binary64 numbers: the ratio below needs their sizes
        // against each other only.
        const slong scale = arf_is_zero(re)   ? arf_abs_bound_lt_2exp_si(im)
                            : arf_is_zero(im) ? arf_abs_bound_lt_2exp_si(re)
                                              : std::max(arf_abs_bound_lt_2exp_si(re), arf_abs_bound_lt_2exp_si(im));
        arf_mul_2exp_si(scaled.value(), re, -scale);
        const double a = std::abs(arf_get_d(scaled.value(), ARF_RND_NEAR));
        arf_mul_2exp_si(scaled.value(), im, -scale);
        const double b = std::abs(arf_get_d(scaled.value(), ARF_RND_NEAR));
        worst = std::max(worst, static_cast<double>(degree) * std::log2((a + b) / std::hypot(a, b)));
    }
    return worst;
}
void SecularRefinement::preliminary(const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics) {
    auto *polynomial = polynomial_;
    if (!polynomial)
        return;
    const slong n = source_->degree();
    slong precision = *std::max_element(iteration_bits_.begin(), iteration_bits_.end());
    auto coefficients = polynomial->coefficient_enclosure(precision);
    RootVector initial(n), output(n);
    for (slong i = 0; i < n; ++i) acb_set(initial[i], proposals_[static_cast<size_t>(i)].value());
    auto ordinary = policy_;
    ordinary.study = StudyMethod::Aberth;
    ordinary.refinement = RefinementMethod::Aberth;
    StudyCandidates candidates(ordinary);
    TraceScope trace(&policy_, "secular_ordinary_preliminary", static_cast<uint64_t>(precision),
                     static_cast<uint64_t>(n), &statistics);
    // One bounded original-source phase (heuristics::preliminary_packet_sweeps),
    // within the call's iteration budget.
    const unsigned sweeps = options.packet_sweeps(heuristics::preliminary_packet_sweeps, statistics);
    if (!sweeps)
        return;
    auto progress = candidates.advance(output.data(), coefficients->polynomial.value(), initial.data(),
                                       precision, sweeps, options, statistics, {}, false);
    auto *counts = policy_.study_statistics.get();
    if (counts) ++counts->preliminary_packets;
    if (progress.cancelled || !progress.finite)
        return;
    for (slong i = 0; i < n; ++i)
        acb_get_mid(proposals_[static_cast<size_t>(i)].value(), output[i]);
    // The validation below evaluates by ball Horner (FLINT's). Where that loses the whole working precision at
    // some candidate, the candidate's value has no correct bit, and the validation, which completes the solve
    // only with every candidate isolated, cannot: it is not run, and the candidates enter the secular model as
    // after a validation that failed. History: docs/code-history-20261021.md, secular_refinement.cpp 18.
    if (ball_horner_loss_bits(output.data(), n, n) >= static_cast<double>(precision)) {
        if (counts) ++counts->preliminary_validations_skipped;
        return;
    }
    // Three validations at most: the packet's points; where they are isolated but not sufficient, the same points at
    // the doubled precision; where those still fall short, the points of a polish at that precision. The second, added
    // on 29 October 2026, spares the polish where the doubled precision alone suffices (cluster-irrational-32-80 with
    // its real-line seeds: 194 M to 113 M instructions; dev/docs/slowdowns-20261029.md, proposal 4).
    for (unsigned attempt = 0; attempt < 3; ++attempt) {
        bool isolated = false, sufficient = true;
        std::vector<Root> cover;
        {
            StageTimer timer(statistics.validation_ns);
            TraceScope validation(&policy_, "preliminary_original_validation", static_cast<uint64_t>(precision),
                                  static_cast<uint64_t>(n), &statistics);
            ++statistics.validation_calls;
            isolated = _acb_poly_validate_roots(output.data(), coefficients->polynomial.value()->coeffs,
                                                n + 1, precision) == n;
            if (isolated) {
                const bool real = source_->real() &&
                    acb_poly_validate_real_roots(output.data(), coefficients->polynomial.value(), precision);
                cover.resize(static_cast<size_t>(n));
                for (slong i = 0; i < n; ++i) {
                    auto &root = cover[static_cast<size_t>(i)];
                    acb_set(root.box.value(), output[i]);
                    root.count = root.multiplicity = multiplicity_;
                    root.certificate = "ordinary-preliminary-original-polynomial-inclusion";
                    if (real && arb_contains_zero(acb_imagref(root.box.value()))) {
                        root.realness = CPS_REAL;
                        arb_zero(acb_imagref(root.box.value()));
                    } else if (!arb_contains_zero(acb_imagref(root.box.value())))
                        root.realness = CPS_NONREAL;
                    root.publish(precision);
                    sufficient = sufficient && coordinates_.meets(root, options, precision);
                }
            }
        }
        maximum_precision_ = std::max(maximum_precision_, static_cast<uint64_t>(precision));
        if (isolated && sufficient) {
            std::vector<Root> exported;
            exported.reserve(cover.size());
            for (const auto &root : cover)
                exported.push_back(coordinates_.original(root, precision));
            sufficient = disjoint(cover) && disjoint(exported);
        }
        if (isolated && sufficient) {
            certified_.consider(std::move(cover), std::vector<uint8_t>(static_cast<size_t>(n), 1), precision);
            if (counts) ++counts->preliminary_completions;
            return;
        }
        // Only independently isolated candidates qualify for the doubled precision and a final ordinary polish.
        // Unresolved geometry enters the secular model with the finite preliminary points intact.
        if (attempt == 0) {
            if (!isolated || precision >= static_cast<slong>(options.max_work_bits) || options.cancelled())
                return;
            precision = std::min(static_cast<slong>(options.max_work_bits), resources::precision_product(precision, 2));
            coefficients = polynomial->coefficient_enclosure(precision);
            for (slong i = 0; i < n; ++i) acb_set(output[i], proposals_[static_cast<size_t>(i)].value());
            continue;
        }
        if (attempt == 2 || options.cancelled())
            return;
        for (slong i = 0; i < n; ++i) acb_set(initial[i], proposals_[static_cast<size_t>(i)].value());
        const unsigned polish = options.packet_sweeps(heuristics::secular_packet_sweeps, statistics);
        if (!polish)
            return; // the iteration budget is spent: the isolated points go to the model as they are
        if (counts) ++counts->ordinary_polishes;
        progress = candidate_packet(output.data(), coefficients->polynomial.value(), initial.data(),
                                    precision, polish, CandidateEngine::Aberth, options, statistics);
        if (!progress.finite || progress.cancelled)
            return;
        for (slong i = 0; i < n; ++i)
            acb_get_mid(proposals_[static_cast<size_t>(i)].value(), output[i]);
    }
}
bool SecularRefinement::newton_cover(const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics, bool *attempted) {
    if (attempted) *attempted = false;
    auto *polynomial = polynomial_;
    if (!polynomial || options.real_only ||
        (engine_goal(options) != CPS_ROOT_APPROXIMATION && engine_goal(options) != CPS_ROOT_ACCURACY))
        return false;
    // At least the goal's working precision, as the driver starts a goal
    // (goal bits plus its 32 guard bits, solver.cpp subdivision_start): after
    // a tighter goal the iteration's points may still be at the last goal's
    // precision, from which the certificate's raises would not reach the
    // new goal within the attempt's allowance.
    const slong precision = std::max<slong>(
        *std::max_element(iteration_bits_.begin(), iteration_bits_.end()),
        static_cast<slong>(std::min<uint64_t>(options.max_work_bits,
                                              static_cast<uint64_t>(resources::precision_sum(options.bits, 32)))));
    if (attempted) *attempted = static_cast<uint64_t>(precision) <= options.max_work_bits;
    std::vector<Root> checked;
    RootVector seeds(source_->degree());
    for (size_t i = 0; i < proposals_.size(); ++i) acb_set(seeds[static_cast<slong>(i)], proposals_[i].value());
    // With the handover, a failed attempt whose machine run converged
    // (corrections at binary64's floor) seeds the preliminary phase. An
    // exhausted run's points are mid-flight and are not evidence of a better
    // start. Either way they are hints; every certificate uses the source.
    const bool retain = policy_.newton_point_handover && policy_.secular_preliminary && !preliminary_done_;
    std::unique_ptr<RootVector> machine;
    bool converged = false;
    const auto hand_over = [&] {
        if (retain && converged && machine && machine->size() == static_cast<slong>(proposals_.size()))
            for (size_t i = 0; i < proposals_.size(); ++i)
                acb_get_mid(proposals_[i].value(), (*machine)[static_cast<slong>(i)]);
        return false;
    };
    slong certified_at = precision;
    PartialNewtonCover partial;
    // A retry after a goal change, with a cover held: the nodes whose region
    // already meets the new goal are settled, and the attempt works on the
    // others from the iteration's points (machine_newton_cover's `settled`).
    // Retried over every point, it would redo what the next validation does
    // where most nodes are settled: a factor whose goal the driver raises to
    // lift a power, or, after the staged default goal, a cluster's points
    // beside the kept certificates. History:
    // docs/code-history-20261021.md, secular_refinement.cpp 19.
    std::vector<uint8_t> settled;
    if (certified_.available) {
        settled.assign(proposals_.size(), 0);
        std::vector<size_t> order(proposals_.size());
        std::iota(order.begin(), order.end(), 0);
        const auto real_mid = [&](size_t i) { return arb_midref(acb_realref(proposals_[i].value())); };
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return arf_cmp(real_mid(a), real_mid(b)) < 0; });
        Float lower, upper;
        for (const auto &region : certified_.cover) {
            if (!coordinates_.meets(region, options, static_cast<slong>(maximum_precision_)))
                continue;
            arb_get_lbound_arf(lower.value(), acb_realref(region.box.value()), static_cast<slong>(maximum_precision_));
            arb_get_ubound_arf(upper.value(), acb_realref(region.box.value()), static_cast<slong>(maximum_precision_));
            auto first = std::lower_bound(order.begin(), order.end(), lower,
                                          [&](size_t i, const Float &x) { return arf_cmp(real_mid(i), x.value()) < 0; });
            for (auto it = first; it != order.end() && arf_cmp(real_mid(*it), upper.value()) <= 0; ++it)
                if (acb_contains(region.box.value(), proposals_[*it].value()))
                    settled[*it] = 1;
        }
    }
    if (!machine_newton_cover(checked, *polynomial, seeds.data(), options, precision,
                              policy_.candidate_allow_underflow, statistics, policy_.coefficient_newton_bound, &policy_,
                              retain ? &machine : nullptr, &converged, &certified_at, true, &partial, settled)) {
        if (partial.points && partial.held() && partial.pass_complete) {
            keep_certificates(partial, options, precision, statistics);
            return false;
        }
        return hand_over();
    }
    // The goal and disjointness are judged on the regions as exported: in the original coordinates, which are the
    // refinement's own where they do not move the roots (no copy is made then).
    for (auto &root : checked)
        root.count = root.multiplicity = multiplicity_;
    std::vector<Root> original;
    if (!coordinates_.identity()) {
        original.reserve(checked.size());
        for (const auto &root : checked)
            original.push_back(coordinates_.original(root, precision));
    }
    const std::vector<Root> &exported = coordinates_.identity() ? checked : original;
    const bool sufficient =
        std::all_of(exported.begin(), exported.end(), [&](const Root &root) { return meets_goal(root, options, precision); });
    if (!sufficient || !disjoint(exported))
        return hand_over();
    std::vector<ComplexBall> retained(checked.size());
    for (size_t i = 0; i < checked.size(); ++i)
        acb_get_mid(retained[i].value(), checked[i].box.value());
    // The attempt may have certified at a precision above this stage's: the
    // precision reported is the precision used.
    maximum_precision_ = std::max(maximum_precision_, static_cast<uint64_t>(std::max(precision, certified_at)));
    certified_.consider(std::move(checked), std::vector<uint8_t>(proposals_.size(), 1), precision);
    proposals_.swap(retained);
    return true;
}
// The point with the shortest coordinates well inside a region: each
// coordinate of its centre rounded to a quarter of that coordinate's radius,
// and zero where the region crosses the other axis. A kept node may sit
// anywhere in its certified region; the attempt's centre carries the
// arithmetic of its last steps (after complex Newton steps a real root's
// centre keeps a tiny imaginary part, their rounding), which in binary64
// complex arithmetic underflows (Smith's division of the differences of the
// points) and refuses every binary64 packet of the model. History:
// docs/code-history-20261021.md, secular_refinement.cpp 20.
namespace {
void short_point(acb_t out, const acb_t box) {
    for (const bool imaginary : {false, true}) {
        const arb_struct *part = imaginary ? acb_imagref(box) : acb_realref(box);
        arb_ptr target = imaginary ? acb_imagref(out) : acb_realref(out);
        Float radius, half;
        arf_set_mag(radius.value(), arb_radref(part));
        arb_zero(target);
        if (arf_is_zero(radius.value())) {
            arf_set(arb_midref(target), arb_midref(part));
            continue;
        }
        arf_mul_2exp_si(half.value(), radius.value(), -1);
        if (arf_cmpabs(arb_midref(part), half.value()) <= 0)
            continue;
        // |error| <= 2^(e_mid - bits) = 2^(e_rad - 3) < radius / 4.
        const slong bits = std::max<slong>(2, arf_abs_bound_lt_2exp_si(arb_midref(part)) -
                                                  (arf_abs_bound_lt_2exp_si(radius.value()) - 3));
        arf_set_round(arb_midref(target), arb_midref(part), bits, ARF_RND_NEAR);
    }
}
} // namespace
// A refused attempt's certificates are kept (PartialNewtonCover): a node that
// holds one stays in its region (short_point), inactive, and the validation's
// cover uses its region (secular_inclusion_cover); the iteration works on the
// other nodes, from where the attempt left them. A cluster below the
// certificate's resolution refuses the attempt while the other points hold
// certificates; without them every node would be iterated at the precision
// the cluster needs. Certificates kept earlier that still meet the goal stay;
// a new region that meets a kept one is not taken.
//
// Only an attempt that tried every point keeps its certificates (newton_cover):
// what it refused is then what its own multiprecision rounds did not certify,
// and the iteration that follows is about those points. The preliminary phase
// (a binary64 packet, a polish no more precise than the rounds', an isolation
// of every point) has nothing to add, and is not run. When the refusals'
// allowance ended the pass, most points were not tried: the iteration that
// follows is a whole one, from binary64's stage, which nodes certified far
// beyond binary64 do not fit. Its certificates are not kept. History:
// docs/code-history-20261021.md, secular_refinement.cpp 21.
void SecularRefinement::keep_certificates(PartialNewtonCover &partial, const cps_options &options, slong precision,
                                          certpolysolve::detail::CandidateStatistics &statistics) {
    const size_t n = proposals_.size();
    if (!partial.points || partial.points->size() != static_cast<slong>(n) || partial.roots.size() != n ||
        partial.certified.size() != n)
        return;
    std::vector<Root> regions;
    std::vector<size_t> index;
    for (size_t i = 0; i < n; ++i)
        if (kept_[i]) {
            regions.push_back(kept_roots_[i]);
            index.push_back(i);
        }
    const size_t earlier = regions.size();
    for (size_t i = 0; i < n; ++i)
        if (partial.certified[i] && !kept_[i]) {
            auto &root = partial.roots[i];
            root.count = root.multiplicity = multiplicity_;
            if (!coordinates_.meets(root, options, precision))
                continue;
            regions.push_back(root);
            index.push_back(i);
        }
    const auto conflict = overlapping_roots(regions);
    for (size_t k = 0; k < regions.size(); ++k)
        if (k >= earlier && !conflict[k]) {
            kept_roots_[index[k]] = std::move(regions[k]);
            kept_[index[k]] = 1;
        }
    if (!coordinates_.identity()) {
        // As the validation's cover: the exported regions must be disjoint too.
        std::vector<Root> exported;
        for (size_t i = 0; i < n; ++i)
            if (kept_[i])
                exported.push_back(coordinates_.original(kept_roots_[i], precision));
        if (!disjoint(exported))
            std::fill(kept_.begin(), kept_.end(), 0);
    }
    // The attempt's other points are where its packet and rounds left them,
    // better than initial candidates. A retry once the iteration has run (a
    // goal change) starts with a binary64 packet, which rounds refined nodes
    // that binary64 does not tell apart onto one point; the iteration then
    // keeps its own points. History:
    // docs/code-history-20261021.md, secular_refinement.cpp 22.
    const bool fresh = !preliminary_done_ && rounds_ == 0;
    for (size_t i = 0; i < n; ++i) {
        if (kept_[i])
            short_point(proposals_[i].value(), kept_roots_[i].box.value());
        else if (fresh)
            acb_get_mid(proposals_[i].value(), (*partial.points)[static_cast<slong>(i)]);
        certified_.inactive[i] = certified_.inactive[i] | kept_[i];
    }
    maximum_precision_ = std::max(maximum_precision_, static_cast<uint64_t>(std::max(precision, partial.certificate_bits)));
    preliminary_done_ = true;
    if (policy_.trace) {
        TraceEvent event{"secular_kept_certificates", "attempt-refused", maximum_precision_, active(), &statistics,
                         nullptr};
        event.engine = "secular";
        event.reason = "pass-complete";
        event.certified = static_cast<uint64_t>(std::count(kept_.begin(), kept_.end(), uint8_t(1)));
        policy_.trace(policy_.trace_data, event);
    }
}
cps_stop_reason SecularRefinement::step(const cps_options &caller, certpolysolve::detail::SearchStatistics &search,
                                       certpolysolve::detail::CandidateStatistics &statistics) {
    if (caller.cancelled())
        return CPS_STOP_CANCELLED;
    // While the default goal is staged (`staged_`), the refinement works to
    // the library's default goal and the caller's goal follows: every node
    // certified at the default goal returns the refinement to the caller's,
    // a goal change that retries the certificate attempt from those points.
    if (staged_ && certified_.available && !active())
        staged_ = false;
    cps_options default_goal;
    if (staged_) {
        default_goal = caller;
        default_goal.bits = cps_options{}.bits;
        default_goal.staged_goal = true;
    }
    const cps_options &options = staged_ ? default_goal : caller;
    caller_bits_ = caller.bits;
    // Meter this step (dovetail.hpp): the regeneration's evaluations and node
    // products, each at the precision it ran at; the other source terms
    // evaluated before the packet (restart seeds) and each validation's n
    // inclusion radii at the model's precision; the packet's terms and
    // interactions at its tier.
    RegenerationWork regeneration;
    research_packet_started_ = false;
    struct Meter {
        SecularRefinement &self;
        const certpolysolve::detail::CandidateStatistics &statistics;
        const RegenerationWork &regeneration;
        uint64_t terms, validations;
        PacketWork packet;
        uint64_t prepaid_terms;
        double exact_taylor_before;
        ~Meter() {
            if (self.polynomial_)
                self.work_ += service_charge(&self.policy_, ServiceOperation::ExactTaylor,
                    self.polynomial_->research_exact_taylor_work() - exact_taylor_before);
            const auto before_packet = packet.terms.value_or(statistics.evaluation_terms);
            const double model_weight = operation_weight(static_cast<int64_t>(self.maximum_precision_));
            self.work_ += service_charge(&self.policy_, ServiceOperation::Regeneration,
                self.policy_.research_service_meter ? regeneration.service_operations : regeneration.operations);
            if (self.policy_.research_service_meter)
                self.work_ += service_charge(&self.policy_, ServiceOperation::ShortProduct,
                    static_cast<double>(regeneration.short_product_factors), ShortNodeProducts::product_bits);
            if (self.policy_.research_service_meter)
                self.work_ += service_charge(&self.policy_, ServiceOperation::LimbHorner, regeneration.limb_horner_work);
            if (self.policy_.trace && regeneration.operations > 0) {
                TraceEvent event{"regeneration", "work", self.maximum_precision_, self.proposals_.size(),
                                 &statistics, nullptr};
                event.metered_work = regeneration.operations;
                self.policy_.trace(self.policy_.trace_data, event);
            }
            const auto prepaid = self.research_prepaid_local_terms_ - prepaid_terms;
            if (before_packet < terms || before_packet - terms < regeneration.source_terms ||
                before_packet - terms - regeneration.source_terms < prepaid)
                std::terminate(); // An accounting invariant failure must not become unsigned work credit.
            self.work_ += service_charge(&self.policy_, ServiceOperation::Auxiliary,
                static_cast<double>(before_packet - terms - regeneration.source_terms - prepaid) * model_weight,
                static_cast<int64_t>(self.maximum_precision_));
            const double validation_work = static_cast<double>(statistics.validation_calls - validations) *
                                           static_cast<double>(self.proposals_.size()) * model_weight;
            self.work_ += service_charge(&self.policy_, ServiceOperation::Auxiliary, validation_work,
                                        static_cast<int64_t>(self.maximum_precision_));
            if (packet.terms) {
                if (self.policy_.research_service_meter) {
                    if (packet.service_scaled_binary64 > packet.binary64_units) std::terminate();
                    self.work_ += service_charge(&self.policy_, ServiceOperation::Binary64,
                        static_cast<double>(packet.binary64_units - packet.service_scaled_binary64), 53);
                    self.work_ += service_charge(&self.policy_, ServiceOperation::ScaledBinary64,
                        static_cast<double>(packet.service_scaled_binary64), 53);
                    self.work_ += service_charge(&self.policy_, ServiceOperation::Evaluation, packet.service_evaluation);
                    self.work_ += service_charge(&self.policy_, ServiceOperation::Interaction, packet.service_interaction);
                    return;
                }
                const uint64_t units = statistics.evaluation_terms - *packet.terms + statistics.root_interactions -
                                       packet.interactions;
                const uint64_t binary64 = std::min(packet.binary64_units, units);
                self.work_ += static_cast<double>(binary64) * operation_weight(std::numeric_limits<double>::digits) +
                              static_cast<double>(units - binary64) * operation_weight(packet.bits);
            }
        }
    } meter{*this, statistics, regeneration, statistics.evaluation_terms, statistics.validation_calls, {}
        , research_prepaid_local_terms_, polynomial_ ? polynomial_->research_exact_taylor_work() : 0
    };
    const bool interpolated_model = conditioned_secular(*source_, &policy_);
    auto same = [](const Number &a, const Number &b) {
        return a.radix == b.radix && fmpz_equal(a.exponent.value(), b.exponent.value()) &&
               fmpq_equal(a.significand.value(), b.significand.value());
    };
    if (last_query_ != options.goal || last_goal_ != options.bits || !same(last_absolute_, options.absolute_tolerance) ||
        !same(last_radius_, options.cluster_radius)) {
        // Binary64 owns the iteration in its initial stage only
        // (SecularMachineWorkspace::stage_bits). A node the growth took past
        // the last goal's stage holds more than binary64 does, and a binary64
        // packet has nothing to add to it, yet a tighter goal raises the stage
        // past it, and the packet's write-back would round it. Binary64
        // arithmetic that reached its floor has nothing to add at a later
        // goal: the machine is retired for the refinement's life, as at a
        // spent precision budget (below). History:
        // docs/code-history-20261021.md, secular_refinement.cpp 45.
        if (last_goal_ != 0 && !machine_retired_) {
            const slong stage = SecularMachineWorkspace::stage_bits(last_goal_);
            if (std::any_of(iteration_bits_.begin(), iteration_bits_.end(), [&](slong bits) { return bits > stage; })) {
                machine_retired_ = true;
                if (policy_.trace) {
                    TraceEvent event{"secular_machine_retired", "goal", maximum_precision_, active(), &statistics,
                                     &search};
                    event.engine = "secular";
                    event.reason = "outgrown-stage";
                    policy_.trace(policy_.trace_data, event);
                }
            }
        }
        std::fill(certified_.inactive.begin(), certified_.inactive.end(), 0);
        // A goal change returns certified nodes as they are (missed_growths_ serves the validation's returns).
        std::fill(missed_growths_.begin(), missed_growths_.end(), 0);
        // Kept certificates (keep_certificates) stay while their regions
        // meet the goal, and their nodes stay where they were certified.
        for (size_t i = 0; i < kept_.size(); ++i) {
            if (kept_[i] && !coordinates_.meets(kept_roots_[i], options, static_cast<slong>(maximum_precision_)))
                kept_[i] = 0;
            certified_.inactive[i] = kept_[i];
        }
        declined_clusters_.clear();
        last_goal_ = options.bits;
        last_query_ = options.goal;
        last_absolute_ = options.absolute_tolerance;
        last_radius_ = options.cluster_radius;
        // The model was built for the last goal: its residues made adequate to it and its validation judging its
        // certificates against it. It does not answer this one, whatever its nodes: the next regeneration rebuilds
        // and validates it, which puts the certificates kept for this goal (keep_certificates) into the cover.
        // History: docs/code-history-20261021.md, secular_refinement.cpp 46.
        model_bits_.clear();
        need_regeneration_ = true;
        limit_reseeded_ = false;
        // A tighter goal (the caller's after the staged default goal, or a
        // separation round's) is first tried by the certificate attempt from
        // the current points: converged at the last goal, they are the
        // attempt's best start, where the iteration would otherwise go on at
        // the new goal's precision. History:
        // docs/code-history-20261021.md, secular_refinement.cpp 47.
        newton_attempted_ = false;
    }
    if (proposals_.empty()) {
        certified_.available = true;
        return CPS_STOP_NONE;
    }
    if (policy_.newton_point_refinement && !newton_attempted_) {
        newton_attempted_ = true;
        bool attempted = false;
        if (newton_cover(options, statistics, &attempted))
            return CPS_STOP_NONE;
        // Refused from initial candidates at a goal above the attempt's ladder
        // (its first rung, twice binary64's bits), the iteration would run at
        // the goal's precision and grow from there. The default goal's
        // iteration isolates the roots in binary64 and at low precision, and
        // the attempt certifies the caller's goal from its points. Only where
        // the goal's working precision fits the budget's, where the attempt
        // can certify it. History:
        // docs/code-history-20261021.md, secular_refinement.cpp 48.
        const uint64_t working = static_cast<uint64_t>(resources::precision_sum(options.bits, 32));
        if (attempted && !staged_once_ && rounds_ == 0 && options.bits > cps_options{}.bits &&
            working > static_cast<uint64_t>(heuristics::newton_ladder_first_rung) && working <= options.max_work_bits) {
            staged_ = staged_once_ = true;
            // Nothing is iterated yet: the nodes start where a default-goal
            // solve starts them, in binary64's stage.
            const slong default_working = static_cast<slong>(
                std::min<uint64_t>(options.max_work_bits,
                                   static_cast<uint64_t>(resources::precision_sum(cps_options{}.bits, 32))));
            for (auto &bits : iteration_bits_)
                bits = std::min(bits, default_working);
            if (policy_.trace) {
                TraceEvent event{"secular_default_goal_first", "attempt-refused", maximum_precision_, active(),
                                 &statistics, &search};
                event.engine = "secular";
                event.reason = "default-goal-first";
                policy_.trace(policy_.trace_data, event);
            }
            return CPS_STOP_NONE;
        }
    }
    if (policy_.secular_preliminary && !preliminary_done_) {
        preliminary(options, statistics);
        if (options.cancelled()) return CPS_STOP_CANCELLED;
        preliminary_done_ = true;
        if (!active()) return CPS_STOP_NONE;
    }
    if (growth_stage_ != GrowthStage::None)
        if (const auto completed = advance_growth(options, search, statistics))
            return *completed;
    // Equal finite proposals need separate interpolation nodes. This is a
    // heuristic restart, with all resulting root statements re-certified.
    separate_interpolation_nodes(proposals_, certified_.inactive, iteration_bits_);
    if (need_regeneration_) {
        switch (refresh_model(options, search, statistics, regeneration)) {
        case ModelRefresh::PacketReady:
            break;
        case ModelRefresh::StepComplete:
            return CPS_STOP_NONE;
        case ModelRefresh::Cancelled:
            return CPS_STOP_CANCELLED;
        }
    }
    research_packet_started_ = true;
    return advance_packet(options, search, statistics, interpolated_model, meter.packet);
}
} // namespace certpolysolve::detail
