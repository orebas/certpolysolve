// PiecewiseRootState::Impl's rounds (the class: piecewise_state.hpp): a round (start_round) streams the models at
// model_bits; each model (prepare_model) takes its count, carried points and fresh proposals, polished and certified
// (process_pending, process_candidate). finish_round observes and closes it; bounded deflated discovery seeds
// missing roots, then schedule_next_round preserves recovery obligations or grows model/arithmetic precision.
// Design: docs/piecewise-root-certification.md, docs/piecewise-ownership-invariants-20260920.md.
#include "methods/piecewise/piecewise_state.hpp"
#include "methods/piecewise/finishing/piecewise_finisher_source.hpp"
#include "methods/piecewise/piecewise_points.hpp"
#include "core/dovetail.hpp"
#include "core/telemetry.hpp"
#include "core/wide_int.hpp"
#include <bit>
#include <cmath>
#include <complex>
#include <numbers>


namespace certpolysolve::detail {

bool PiecewiseRootState::Impl::needed_candidate(size_t i) {
    if (guards.contains(round_points[i].value())) return false;
    if (config.owned_only_after_discovery && discovered()) {
        const auto id = owner_of(i);
        return id != PiecewiseOwnership::none && owner_candidate[id] == i;
    }
    if (ownership) {
        const auto id = ownership->containing(round_points[i].value());
        if (id && ownership->proof(*id).count == 1)
            return *id < owner_candidate.size() && owner_candidate[*id] == i;
    }
    return supplied || !seed_guards.contains(round_points[i].value(), i);
}
void PiecewiseRootState::Impl::seed_missing_roots(const cps_options &o) {
    const uint64_t owned = ownership->count();
    if (owned >= static_cast<uint64_t>(degree)) return;
    StageTimer timer(stats.deflation_ns);
    ++stats.deflation_rounds;
    const double budget = counted_work(stats.work) - round_start_work;
    const double start_work = counted_work(stats.work);
    if (!deflation_source)
        deflation_source = std::make_unique<PiecewiseFinisherSource>(std::vector<Coefficient>(
            input.coefficients.begin() + static_cast<ptrdiff_t>(zeros),
            input.coefficients.begin() + static_cast<ptrdiff_t>(length)), stats);
    using Complex = std::complex<double>;
    struct TwoPart { Complex high, low; };
    Float residual;
    const auto two_part = [&](const acb_t z) {
        TwoPart t;
        const auto split = [&](const arf_t x, double &high, double &low) {
            high = arf_get_d(x, ARF_RND_NEAR);
            arf_set_d(residual.value(), high);
            arf_sub(residual.value(), x, residual.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
            low = arf_get_d(residual.value(), ARF_RND_NEAR);
        };
        double rh, rl, ih, il;
        split(arb_midref(acb_realref(z)), rh, rl);
        split(arb_midref(acb_imagref(z)), ih, il);
        t.high = {rh, ih};
        t.low = {rl, il};
        return t;
    };
    const auto difference = [](const TwoPart &a, const TwoPart &b) { return (a.high - b.high) + (a.low - b.low); };
    std::vector<TwoPart> held;
    std::vector<double> weights;
    for (size_t id = 0; id < owner_candidate.size(); ++id) {
        if (!ownership->active(id)) continue;
        const auto &proof = ownership->proof(id);
        const size_t candidate = owner_candidate[id];
        ComplexBall point;
        if (candidate < carried.size() && ownership->contains(id, carried[candidate].value()))
            acb_get_mid(point.value(), carried[candidate].value());
        else
            acb_get_mid(point.value(), proof.box.value());
        held.push_back(two_part(point.value()));
        weights.push_back(static_cast<double>(proof.count));
    }
    // Fujiwara: every root has modulus at most 2 max(|c_(d-k) / c_d|^(1/k),
    // |c_0 / (2 c_d)|^(1/d)), here in log2.
    ComplexBall coefficient;
    Magnitude magnitude;
    const auto log2_magnitude = [&](slong k) {
        input.coefficients[zeros + static_cast<size_t>(k)].enclose(coefficient.value(), 64);
        acb_get_mag(magnitude.value(), coefficient.value());
        return mag_is_zero(magnitude.value()) ? -std::numeric_limits<double>::infinity()
                                              : mag_get_d_log2_approx(magnitude.value());
    };
    const double leading = log2_magnitude(degree);
    double log2_bound = -std::numeric_limits<double>::infinity();
    for (slong k = 1; k <= degree; ++k) {
        const double ratio = log2_magnitude(degree - k) - leading - (k == degree ? 1 : 0);
        log2_bound = std::max(log2_bound, ratio / static_cast<double>(k));
    }
    if (!std::isfinite(log2_bound)) return;
    const double radius = std::exp2(1 + log2_bound);
    // Free points on that circle at golden-angle steps, none on the real axis.
    const size_t missing = static_cast<size_t>(degree) - static_cast<size_t>(owned);
    std::vector<ComplexBall> points(missing);
    std::vector<TwoPart> parts(missing);
    std::vector<bool> converged(missing, false), stopped(missing, false), limited(missing, false);
    std::vector<double> nearest(missing, std::numeric_limits<double>::infinity());
    const double golden_angle = std::numbers::pi * (3 - std::sqrt(5.0));
    for (size_t i = 0; i < missing; ++i) {
        const auto start = std::polar(radius, golden_angle * (static_cast<double>(i) + 0.5));
        acb_set_d_d(points[i].value(), start.real(), start.imag());
        parts[i] = two_part(points[i].value());
    }
    ComplexBall value, derivative, newton, denominator, sum, correction;
    RealBall size, step;
    for (uint64_t sweep = 0; sweep < heuristics::initial_sweep_allowance; ++sweep) {
        bool active = false;
        for (size_t i = 0; i < missing; ++i) {
            if (converged[i] || stopped[i]) continue;
            if (o.cancelled() || counted_work(stats.work) - start_work > budget) {
                sweep = heuristics::initial_sweep_allowance;
                break;
            }
            active = true;
            if (!deflation_source->evaluate(value.value(), derivative.value(), points[i].value(), precision) ||
                acb_contains_zero(derivative.value())) {
                stopped[i] = true;
                continue;
            }
            // At the round's resolution with a step still above the model's:
            // a missing root lies about here, unresolved at this precision.
            if (acb_contains_zero(value.value())) { stopped[i] = limited[i] = true; continue; }
            acb_div(newton.value(), value.value(), derivative.value(), precision);
            Complex s = 0;
            nearest[i] = std::numeric_limits<double>::infinity();
            for (size_t j = 0; j < held.size(); ++j) {
                const Complex d = difference(parts[i], held[j]);
                nearest[i] = std::min(nearest[i], std::abs(d));
                s += weights[j] / d;
            }
            for (size_t k = 0; k < missing; ++k)
                if (k != i) s += 1.0 / difference(parts[i], parts[k]);
            stats.work.interactions += static_cast<double>(held.size() + missing - 1);
            if (!std::isfinite(s.real()) || !std::isfinite(s.imag())) { stopped[i] = true; continue; }
            acb_set_d_d(sum.value(), s.real(), s.imag());
            acb_mul(denominator.value(), newton.value(), sum.value(), precision);
            acb_sub_si(denominator.value(), denominator.value(), 1, precision);
            acb_neg(denominator.value(), denominator.value());
            acb_div(correction.value(), newton.value(), denominator.value(), precision);
            acb_sub(points[i].value(), points[i].value(), correction.value(), precision);
            acb_get_mid(points[i].value(), points[i].value());
            if (!acb_is_finite(points[i].value())) { stopped[i] = true; continue; }
            parts[i] = two_part(points[i].value());
            acb_abs(size.value(), points[i].value(), precision);
            arb_mul_2exp_si(size.value(), size.value(), -model_bits);
            acb_abs(step.value(), correction.value(), precision);
            if (arb_le(step.value(), size.value())) converged[i] = true;
        }
        ++stats.deflation_sweeps;
        if (!active) break;
    }
    // A converged point seeds its root unless an owner, a final guard or an
    // earlier seed of this call (within the model's resolution) has it. A
    // point stopped at the round's resolution only says a missing root is
    // about there: it seeds only where no other such seed, of this call or
    // an earlier one, is nearer than its nearest held root, so an
    // unresolved cluster gets one seed, not one per free point and rung.
    std::vector<ComplexBall> added;
    std::vector<TwoPart> earlier;
    for (const auto k : limited_seeds) earlier.push_back(two_part(carried[k].value()));
    RealBall gap, tolerance;
    for (const bool resolution_limited : {false, true}) {
        for (size_t i = 0; i < missing; ++i) {
            if (resolution_limited ? !limited[i] : !converged[i]) continue;
            auto &z = points[i];
            if (reflect && arb_is_negative(acb_imagref(z.value()))) acb_conj(z.value(), z.value());
            if (ownership->containing(z.value()) || guards.contains(z.value())) continue;
            bool repeated = false;
            if (!resolution_limited) {
                for (const auto &other : added) {
                    acb_sub(correction.value(), z.value(), other.value(), precision);
                    acb_abs(gap.value(), correction.value(), precision);
                    acb_abs(tolerance.value(), z.value(), precision);
                    arb_mul_2exp_si(tolerance.value(), tolerance.value(), -model_bits);
                    if (arb_le(gap.value(), tolerance.value())) { repeated = true; break; }
                }
            } else {
                const double reach = nearest[i];
                if (!std::isfinite(reach)) continue;
                const auto part = two_part(z.value());
                for (const auto &other : earlier)
                    if (std::abs(difference(part, other)) < reach) { repeated = true; break; }
                if (repeated) continue;
                earlier.push_back(part);
                limited_seeds.push_back(carried.size());
            }
            if (repeated) continue;
            added.push_back(z);
            carried.push_back(z);
            ++stats.deflation_seeds;
        }
    }
}
PiecewiseRootState::Impl::RoundStep PiecewiseRootState::Impl::finish_round(const cps_options &o) {
    observe_population(false);
    if (config.snapshot_rounds) snapshot(false, o, true);
    if (config.observe) {
        stats.discovery_count = ownership ? zeros + ownership->count() : 0;
        stats.final_accurate_count = covered;
        config.observe(stats, "round-complete");
    }
    models.reset();
    if (config.axis_only) axis->finish_round(o, result.roots);
    if (config.axis_only && (finished() || (config.axis_rounds && stats.rounds >= config.axis_rounds)))
        return RoundStep::Stop;
    if (force_reseed && !recovery_interrupted) {
        searched_model_bits = model_bits;
        searched_precision = precision;
    }
    force_reseed = false;
    const bool seeded = seed_round_discovery(o);
    schedule_next_round(o, seeded);
    return RoundStep::Continue;
}

bool PiecewiseRootState::Impl::seed_round_discovery(const cps_options &o) {
    // Discovery is what is left when the round owned no new root, or
    // when every owned root already meets the goal. Seeding, not
    // solving: the owned roots must outnumber the missing ones; with
    // fewer owned, the free points would be a global solve of most
    // roots, which the regime exit weighs against a global method.
    bool seeded = false;
    if (config.deflated_discovery && ownership && !supplied && !discovered() &&
        2 * ownership->count() > static_cast<uint64_t>(degree) &&
        (ownership->count() == previous_discovery || covered >= zeros + ownership->count()) &&
        (model_bits != deflated_model_bits || precision != deflated_precision)) {
        deflated_model_bits = model_bits;
        deflated_precision = precision;
        const uint64_t seeds_before = stats.deflation_seeds;
        seed_missing_roots(o);
        seeded = stats.deflation_seeds > seeds_before;
    }
    return seeded;
}

void PiecewiseRootState::Impl::schedule_next_round(const cps_options &o, bool seeded) {
    const bool improved = stats.candidate_improvements > previous_improvements;
    if (!improved && covered == previous_covered)
        ++stats.stationary_rounds;
    const bool candidate_round = stats.model_limited == previous_models &&
        stats.rounding_limited == previous_rounding && stats.candidate_limited > previous_candidates &&
        stats.accepted_moves > previous_moves;
    const bool recovery_requested = candidate_round && config.recover_discovery &&
        !supplied && covered == previous_covered;
    const bool searched = model_bits == searched_model_bits && precision == searched_precision;
    const bool deferred_stalled = proposal_deferred && !discovered() &&
        (!ownership || ownership->count() == previous_discovery) &&
        stats.model_limited == previous_models;
    // A restricted discovery obligation survives model-limited
    // rounds. Improve those models first; otherwise restore the
    // full chart before concluding that the carried pool is spent.
    const bool recovery_due = (recovery_requested || deferred_stalled) && !searched && !discovered();
    if (recovery_due && seeded && config.seeds_before_recovery) {
        // The deflated discovery has just seeded missing roots into the carried pool. A recovery
        // would start every chart's discovery again and polish every candidate it yields, nearly
        // all of them roots already owned. The seeds are tried first: the next round polishes
        // and certifies the pool at this rung. The recovery is owed, not dropped.
        // History: docs/code-history-20261021.md, piecewise_roots.cpp 13.
        recovery_owed = true;
    } else if (recovery_due || (recovery_owed && !searched && !discovered())) {
        // A stalled pool may be missing a root. Preserve the one
        // bounded discovery recovery before escalating precision,
        // even when candidate corrections have stopped improving.
        force_reseed = true;
        recovery_owed = false;
    } else if (candidate_round && (improved || covered > previous_covered)) {
        // Same-precision continuation needs improved corrections, not merely movement between
        // overlapping charts: scheduling, never a proof.
        // Continue useful refinement without repeating deterministic
        // discovery. Coordinate jitter alone does not enter here.
        if (recovery_requested && searched) ++stats.repeated_recoveries_avoided;
    } else {
        if (stats.rounding_limited - previous_rounding > stats.model_limited - previous_models)
            precision = clamped_rung(resources::precision_product(precision, 2), o);
        else {
            const slong refused_bits = model_bits;
            model_bits = resources::precision_product(model_bits, 2);
            // Per root, the doubled rung costs about as much as the
            // next one (measured on a flat input of high degree), so
            // a rung that cannot finish the most demanding
            // model-limited root is not worth building. History:
            // docs/code-history-20261021.md, piecewise_roots.cpp 14.
            if (config.predicted_refinement && discovered() && skip_doubled_rung(refused_bits, o)) {
                model_bits = resources::precision_product(model_bits, 2);
                ++stats.predicted_rungs_skipped;
            }
            precision = clamped_rung(std::max(precision,
                resources::precision_sum(
                    resources::precision_sum(resources::precision_product(model_bits, 2), o.bits), 64)), o);
        }
    }
}

bool PiecewiseRootState::Impl::skip_doubled_rung(slong refused_bits, const cps_options &o) const {
    uint64_t refused = 0, covered_by_doubling = 0;
    const auto doubling_bucket = static_cast<size_t>(std::bit_width(static_cast<uint64_t>(refused_bits) - 1));
    for (size_t k = 0; k < stats.certificate.round_shortfall_log2.size(); ++k) {
        refused += stats.certificate.round_shortfall_log2[k];
        if (k <= doubling_bucket) covered_by_doubling += stats.certificate.round_shortfall_log2[k];
    }
    if (!refused) return false;
    // Exact: covered <= (doubled / quadrupled) refused is covered x quadrupled <= doubled x refused, the two rungs'
    // arithmetic precisions being integers.
    const auto arithmetic = [&](uint64_t bits) { return 2 * bits + o.bits + 64; };
    const auto model_bits = static_cast<uint64_t>(refused_bits);
    const uint64_t doubled = arithmetic(2 * model_bits), quadrupled = arithmetic(4 * model_bits);
    return uint128{covered_by_doubling} * quadrupled <= uint128{doubled} * refused;
}
void PiecewiseRootState::Impl::start_round(const cps_options &o) {
    StageTimer timer(stats.setup_ns);
    ++stats.rounds;
    recovery_interrupted = false;
    if (config.axis_only) {
        if (!axis)
            axis.emplace(input, zeros, degree, stats);
        axis->start_round(model_bits, precision);
    }
    if (force_reseed) ++stats.recovery_rounds;
    report_precision();
    previous_models = stats.model_limited;
    previous_rounding = stats.rounding_limited;
    previous_candidates = stats.candidate_limited;
    previous_moves = stats.accepted_moves;
    previous_covered = covered;
    previous_improvements = stats.candidate_improvements;
    previous_discovery = ownership ? ownership->count() : 0;
    round_start_work = counted_work(stats.work);
    stats.certificate.round_model_shortfall_bits = 0;
    stats.certificate.round_shortfall_log2.fill(0);
    // Hints used in an earlier round can still explain a missing root.
    // Clear their recovery obligation only when ordinary discovery runs.
    if (force_reseed) proposal_deferred = false;
    if (progress_model_bits != model_bits || progress_precision != precision) {
        best_corrections.assign(carried.size(), std::numeric_limits<double>::infinity());
        progress_model_bits = model_bits;
        progress_precision = precision;
    }
    round_points = carried;
    pool_ready = supplied || !carried.empty();
    point_index = pool_ready ? std::make_unique<PiecewisePointIndex>(round_points) : nullptr;
    // Provisional guards are tied to the old coordinates. Rebuild them as
    // owned seeds move; final original-input certificates stay valid.
    seed_guards = GuardIndex{};
    live_pool = false;
    for (size_t i = 0; i < round_points.size() && !live_pool; ++i)
        live_pool = needed_candidate(i);
    ComplexPolynomial polynomial;
    ComplexBall coefficient;
    for (slong j = 0; j <= degree; ++j) {
        input.coefficients[zeros + static_cast<size_t>(j)].enclose(coefficient.value(), precision);
        acb_poly_set_coeff_acb(polynomial.value(), j, coefficient.value());
    }
    stats.work.enclosures += static_cast<double>(degree + 1) * operation_weight(precision);
    models = std::make_unique<PiecewiseModelCursor>(polynomial.value(), model_bits, precision, stats,
        [this](const RootModel &model) {
            // C1: only the charts centred on the axis (|Im c| well below the radius: sector 0 and S/2; a
            // neighbour's centre is about R / 2 off the axis).
            if (config.axis_only) {
                RealBall off, quarter;
                arb_abs(off.value(), acb_imagref(model.center.value()));
                arb_mul_2exp_si(quarter.value(), model.radius.value(), -2);
                if (!arb_lt(off.value(), quarter.value())) { ++stats.symmetry_skipped_models; return false; }
                return true; // every axis chart, every round: its count is this round's accounting
            }
            if (reflect) {
                RealBall top;
                arb_add(top.value(), acb_imagref(model.center.value()), model.radius.value(), 64);
                if (arb_is_negative(top.value())) { ++stats.symmetry_skipped_models; return false; }
            }
            if (!pool_ready || force_reseed) return true;
            for (auto i : point_index->query(model))
                if (config.owned_only_after_discovery && discovered() ? needed_candidate(i) :
                    !guards.contains(round_points[i].value())) return true;
            return false;
        }, config);
    (void)o;
}
void PiecewiseRootState::Impl::prepare_model(const RootModel &model, const cps_options &o) {
    const slong local_degree = acb_poly_degree(model.polynomial.value());
    if (local_degree < 1) return;
    if (config.axis_only) {
        // The axis pass owns copied charts and finishes at the round's end.
        axis->arrive(model, o);
        return;
    }
    std::optional<uint64_t> count;
    {
        StageTimer timer(stats.local_count_ns);
        count = count_root_model(model, precision);
        stats.work.count += static_cast<double>(local_degree + 1) * operation_weight(precision);
    }
    if (count) ++stats.local_counts;
    if (count && !*count) return;
    if (count && *count == 1) ++stats.one_root_models;
    else if (count) ++stats.other_counted_models;
    std::vector<size_t> selected;
    if (point_index) for (auto i : point_index->query(model))
        if (needed_candidate(i)) selected.push_back(i);
    bool reseed = !supplied && (!pool_ready || force_reseed || local_degree == 2 ||
        (count && *count > selected.size() && selected.empty()));
    if (reseed && ownership) {
        bool accounted = discovered();
        if (!accounted && count) {
            RealBall distance;
            Float lower;
            acb_get_abs_lbound_arf(lower.value(), model.center.value(), precision);
            arb_set_arf(distance.value(), lower.value());
            // A removed monomial must have no zero in this chart before a
            // model's Rouche count is an original-factor count.
            if (!model.lower_index || arb_gt(distance.value(), model.radius.value()))
                accounted = ownership->singletons() && ownership->contained_in(model, precision) == *count;
        }
        if (accounted) { reseed = false; ++stats.counted_discovery_retired; }
    }
    pending = std::make_unique<RootVector>(static_cast<slong>(selected.size()) + (reseed ? local_degree : 0));
    owners.assign(static_cast<size_t>(pending->size()), new_owner);
    pending_stable.assign(static_cast<size_t>(pending->size()), false);
    pending_count = 0;
    pending_cursor = 0;
    pending_polished = false;
    {
        StageTimer timer(stats.local_solve_ns);
        for (auto i : selected) {
            observe_selection(model, i);
            owners[static_cast<size_t>(pending_count)] = i;
            acb_sub((*pending)[pending_count], carried[i].value(), model.center.value(), precision);
            acb_div_arb((*pending)[pending_count], (*pending)[pending_count], model.radius.value(), precision);
            ++pending_count;
        }
        if (reseed) {
            RootVector fresh(local_degree);
            const auto old_core = stats.core_proposal_models;
            const slong n = proposals(fresh, model, o, result.statistics_v3, precision, stats, config,
                                      force_reseed);
            proposal_deferred = proposal_deferred || stats.core_proposal_models > old_core;
            for (slong j = 0; j < n; ++j)
                if (!held_copy(model, fresh[j])) acb_set((*pending)[pending_count++], fresh[j]);
            if (pool_ready) ++stats.reseeded_models;
        }
    }
}
void PiecewiseRootState::Impl::process_candidate(const RootModel &model, const cps_options &o, slong i) {
    candidate_jet.clear();
    const auto seed_before = stats.seed_certificate_ns, final_before = stats.final_certificate_ns;
    size_t owner = owners[static_cast<size_t>(i)];
    if (config.owned_only_after_discovery && discovered() && owner_of(owner) == PiecewiseOwnership::none)
        return;
    size_t candidate_id = owner;
    ComplexBall point;
    acb_mul_arb(point.value(), (*pending)[i], model.radius.value(), precision);
    acb_add(point.value(), point.value(), model.center.value(), precision);
    acb_get_mid(point.value(), point.value());
    if (!acb_is_finite(point.value())) return;
    if (guards.contains(point.value())) { ++stats.retired_points; return; }
    size_t owned = PiecewiseOwnership::none;
    if (ownership) {
        owned = owner_of(owner);
        if (owned != PiecewiseOwnership::none && !ownership->contains(owned, point.value())) {
            // Keep the old owned seed. The escaped point remains a new
            // unproved proposal and can still establish another proof.
            owner = new_owner;
            owned = PiecewiseOwnership::none;
            ++stats.ownership_escapes;
        }
        if (const auto other = ownership->containing(point.value());
            other && ownership->proof(*other).count == 1 &&
            (owner == new_owner || owner_candidate[*other] != owner)) {
            ++stats.owned_proposals_retired;
            return;
        }
    }
    std::optional<Root> reusable_guard;
    std::optional<Root> new_ownership;
    if (!supplied) {
        const size_t id = owner == new_owner ? carried.size() : owner;
        candidate_id = id;
        if (owned == PiecewiseOwnership::none && seed_guards.contains(point.value(), id)) {
            ++stats.seed_retired; return;
        }
        Root coarse;
        ChartFailure reason = ChartFailure::Candidate;
        bool guarded = false, proven = false;
        {
            StageTimer timer(stats.certification_ns);
            StageTimer observed(stats.seed_certificate_ns);
            if (owned != PiecewiseOwnership::none) {
                coarse = ownership->proof(owned);
                guarded = true;
                ++stats.ownership_reuses;
            } else {
                guarded = certify_guard(coarse, model, (*pending)[i], o, &reason) &&
                    (!zeros || !acb_contains_zero(coarse.box.value()));
                if (guarded && ownership) new_ownership = coarse;
            }
            proven = guarded;
            if (guarded && config.reuse_guards) reusable_guard = coarse;
            if (owned == PiecewiseOwnership::none)
                guarded = guarded && seed_guards.insert_point(std::move(coarse), point, id);
        }
        if (o.cancelled()) return;
        if (owner == new_owner) {
            if (!proven && (reason == ChartFailure::Domain ||
                (reason == ChartFailure::Candidate && live_pool))) {
                ++stats.seed_rejected;
                return;
            }
            // A finite in-chart Newton iterate need not yet support a
            // singleton guard, especially near a close pair. Keep such
            // iterates when seeding an empty pool. While usable carried
            // work remains, recovery must establish a coarse proof before
            // adding more candidate-limited iterates to that pool.
            // Retaining an iterate never establishes a root count.
            if (!proven && reason == ChartFailure::Candidate) ++stats.unproved_seeds_retained;
            carried.push_back(point);
        } else {
            if (!acb_equal(carried[id].value(), point.value())) ++stats.accepted_moves;
            carried[id] = point;
        }
        // If final certification is cancelled, resume this same owned
        // candidate rather than appending another copy of a fresh seed.
        owners[static_cast<size_t>(i)] = id;
        if (guarded) ++stats.seed_guards;
        if (new_ownership) remember_ownership(*new_ownership, id);
    }
    Root root;
    ChartFailure failure = ChartFailure::Candidate;
    bool certified;
    int raw_failure = static_cast<int>(failure);
    double correction_log2 = std::numeric_limits<double>::infinity();
    ChartCertificateDiagnostic diagnostic;
    const bool model_stopped = pending_stable[static_cast<size_t>(i)];
    const auto observe_candidate = [&](bool newly_certified, bool improved) {
        if (config.observe_candidate || config.observe_snapshot) {
            PiecewiseCandidateEvent event{candidate_id, static_cast<uint64_t>(model_bits),
                static_cast<uint64_t>(precision), candidate_id < round_points.size(), model_stopped,
                newly_certified, improved, certified ? -1 : static_cast<int>(failure),
                diagnostic.error_log2, diagnostic.margin_log2, correction_log2};
            event.component_terms_available = diagnostic.component_terms_available;
            event.component_log2 = diagnostic.component_log2;
            event.model = stats.models;
            event.local_degree = static_cast<uint64_t>(acb_poly_degree(model.polynomial.value()));
            event.owner = owned;
            event.seed_ns = stats.seed_certificate_ns - seed_before;
            event.final_ns = stats.final_certificate_ns - final_before;
            event.raw_failure = raw_failure;
            event.source = config.diagnostic_source_id;
            if (config.observe_snapshot && candidate_id < carried.size()) {
                candidate_history.resize(carried.size());
                auto &[attempts, last] = candidate_history[candidate_id];
                ++attempts;
                last = event;
            }
            if (config.observe_candidate) config.observe_candidate(event);
        }
    };
    {
        StageTimer timer(stats.certification_ns);
        StageTimer observed(stats.final_certificate_ns);
        ++stats.proposals;
        certified = certify_root_model(root, model, (*pending)[i], o, precision, &failure, &stats.certificate,
                                       &correction_log2,
                                       (config.observe_candidate || config.observe_snapshot) ? &diagnostic : nullptr,
                                       0, false, &stats.work.certificate, &candidate_jet);
        if (o.cancelled()) return;
        if (config.observe_certificate) {
            StageTimer timer(stats.diagnostic_ns);
            ComplexBall copy;acb_set(copy.value(),(*pending)[i]);
            config.observe_certificate(model,copy,o,precision,certified,certified?-1:static_cast<int>(failure),true);
        }
        if (certified && index.overlaps(root, result.roots)) {
            auto tighter = o;
            // Overlap alone cannot identify a root. Equal-count nested
            // certified rectangles do identify the same roots; otherwise
            // try smaller certificates, bounded by current arithmetic.
            auto same_roots = [&] {
                return index.overlaps(root, result.roots, [&](size_t existing) {
                    const auto &other = result.roots[existing];
                    return other.count == root.count &&
                        (contains_rectangle(other, root) || contains_rectangle(root, other));
                });
            };
            while (certified && index.overlaps(root, result.roots)) {
                if (same_roots()) { ++stats.overlap_identified; break; }
                if (!o.absolute_tolerance.zero() || precision <= 32 ||
                    tighter.bits + 8 > static_cast<uint64_t>(precision - 32)) break;
                tighter.bits += 8;
                ++stats.overlap_attempts;
                certified = certify_root_model(root, model, (*pending)[i], tighter, precision,
                                               &failure, &stats.certificate, &correction_log2, nullptr, 0, false,
                                               &stats.work.certificate, &candidate_jet);
                if (o.cancelled()) return;
            }
            if (certified && !index.overlaps(root, result.roots)) ++stats.overlap_resolved;
        }
    }
    raw_failure = certified ? -1 : static_cast<int>(failure);
    if (!certified) {
        bool improved = false;
        if (candidate_id != new_owner &&
            !std::isnan(correction_log2) && correction_log2 < std::numeric_limits<double>::infinity()) {
            best_corrections.resize(carried.size(), std::numeric_limits<double>::infinity());
            // A new binary magnitude level is robust to last-bit jitter
            // and chart changes; sub-bit coordinate movement alone is
            // insufficient reason to rebuild the same models forever.
            const double level = std::ceil(correction_log2);
            if (level < best_corrections[candidate_id]) {
                best_corrections[candidate_id] = level;
                // A fresh ID has no prior-round accuracy to improve.
                // Count progress only for a candidate that was already
                // available when this round began.
                if (candidate_id < round_points.size()) {
                    ++stats.candidate_improvements;
                    improved = true;
                }
            }
        }
        if (model_stopped) {
            ++stats.model_stopped_candidates;
        }
        if (failure == ChartFailure::ModelRemainder) ++stats.model_limited;
        else if (failure == ChartFailure::Rounding) ++stats.rounding_limited;
        else ++stats.candidate_limited;
        observe_candidate(false, improved);
        return;
    }
    ++stats.certificates;
    {
        StageTimer timer(stats.deduplication_ns);
        if (index.overlaps(root, result.roots)) { ++stats.duplicates; observe_candidate(false, false); return; }
        if (root.count > static_cast<uint64_t>(result.degree) - covered)
            throw std::logic_error("piecewise certificates exceed input degree");
        covered += root.count;
        index.insert(root, result.roots.size());
        result.roots.push_back(std::move(root));
    }
    observe_candidate(true, false);
    remember_ownership(result.roots.back(), candidate_id);
    mark_finished_owner(result.roots.back());
    {
        StageTimer timer(stats.certification_ns);
        Root guard;
        ++stats.guard_attempts;
        if (reusable_guard && guards.insert(std::move(*reusable_guard), result.roots.back())) {
            ++stats.guards;
            ++stats.reused_guards;
        } else if (certify_guard(guard, model, (*pending)[i], o) &&
            (!zeros || !acb_contains_zero(guard.box.value())) && guards.insert(std::move(guard), result.roots.back()))
            ++stats.guards;
        else guards.insert(result.roots.back(), result.roots.back());
    }
    publish_reflection();
}
void PiecewiseRootState::Impl::process_pending(const RootModel &model, const cps_options &o) {
    if (!pending_polished) {
        std::vector<uint64_t> owner_counts;
        if (config.observe_polish || config.observe_owned_refinement) {
            owner_counts.resize(static_cast<size_t>(pending_count));
            for (size_t j = 0; j < owner_counts.size(); ++j) {
                const auto id = owner_of(owners[j]);
                if (id != PiecewiseOwnership::none) {
                    const auto &proof=ownership->proof(id);
                    owner_counts[j] = proof.count;
                    if(config.observe_owned_refinement && proof.count==1) {
                        StageTimer diagnostic(stats.diagnostic_ns);
                        config.observe_owned_refinement(model,proof,carried[owners[j]],o,precision);
                    }
                }
            }
        }
        StageTimer timer(stats.local_solve_ns);
        StageTimer polishing(stats.polish_ns);
        polish(*pending, pending_count, model, o, precision, stats, &pending_stable, config, &owners, owner_counts);
        pending_polished = !o.cancelled();
    }
    while (pending_cursor < pending_count && !o.cancelled()) {
        const slong i = pending_cursor;
        process_candidate(model, o, i);
        // A cancelled proof leaves this item pending. Published proofs and
        // owned candidates remain usable when the same session resumes.
        if (o.cancelled()) break;
        ++pending_cursor;
        if (finished()) break;
    }
    if (pending_cursor == pending_count || finished()) pending.reset();
}

} // namespace certpolysolve::detail
