// Rouche's test on a chart's model (certify_root_model): the published rectangle about a candidate, from the goal's
// width (goal_half_width) halved until it passes or fails for a classified cause (the model's error, rounding, an
// inaccurate candidate); the model's jet at the candidate (of order certificate_jet_order or guard_jet_order),
// Taylor's integral remainder for the rest and the model's error (root_model_error); under the approximation goal a
// dominant order m > 1 counts a cluster. A model-remainder refusal records how much stronger the model must be
// (note_model_shortfall). count_root_model, Pellet's test at radius one (Rouche against the largest coefficient) with
// the model's error, only guides discovery.
//
// State: none; each call proves one rectangle or refuses it.
// Entry points: certify_root_model (piecewise_roots.cpp: the guards and final certificates of the piecewise state;
// certificate_attempt.cpp: certify_newton_point and the attempt's certificates, on linear models);
// count_root_model (piecewise_roots.cpp); goal_half_width (here and certificate_attempt.cpp).
// Not here: model construction (piecewise_models.cpp), model-error bounds (numeric/models/root_model.cpp), the
// piecewise state (piecewise_roots.cpp), the
// certificate attempt (certificate_attempt.cpp).
// Design notes: docs/piecewise-root-certification.md (the rectangle certificate), docs/mathematics.md ("The
// piecewise method").
// Moved verbatim from piecewise_roots.cpp in W8's S3 (docs/w8-plan-20261021.md).
#include "certify/chart_certificate.hpp"
#include "certify/chart_statistics.hpp"
#include "core/query.hpp"
#include "core/dovetail.hpp"
#include "core/heuristics.hpp"
#include "numeric/arithmetic/modulus_bound.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <optional>

namespace certpolysolve::detail {

std::optional<uint64_t> count_root_model(const RootModel &model, slong precision) {
    const auto degree = acb_poly_degree(model.polynomial.value());
    if (degree < 0)
        return {};
    // Rouche on the unit circle: the coefficient of largest size against the sum of the others' sizes and the
    // model's error. The others enter through upper bounds only, so FLINT's magnitudes do (an upper bound within
    // 2^-29 of the size, rounded up in every sum): a 64-bit bound of each coefficient's modulus and a ball sum
    // at the working precision would take most of the count's time to sharpen a comparison that the dominant
    // term wins or loses by far more. Its own lower bound stays at the working precision.
    // History: docs/code-history-20261021.md, piecewise_roots.cpp 1.
    const acb_srcptr coefficients = model.polynomial.value()->coeffs;
    Magnitude size, largest, others;
    slong dominant = 0;
    // One pass: `others` is the sum of every size seen but the largest so far.
    for (slong j = 0; j <= degree; ++j) {
        acb_get_mag(size.value(), coefficients + j);
        if (mag_cmp(size.value(), largest.value()) > 0) {
            mag_add(others.value(), others.value(), largest.value());
            mag_set(largest.value(), size.value());
            dominant = j;
        } else
            mag_add(others.value(), others.value(), size.value());
    }
    Float lower, bound;
    acb_get_abs_lbound_arf(lower.value(), coefficients + dominant, precision);
    arf_set_mag(bound.value(), others.value());
    RealBall lhs, rhs;
    arb_set_arf(lhs.value(), lower.value());
    arb_set_arf(rhs.value(), bound.value());
    arb_add(rhs.value(), rhs.value(), model.error.value(), precision);
    if (arb_gt(lhs.value(), rhs.value()))
        return static_cast<uint64_t>(dominant);
    return {};
}

// The half-width of the rectangle certify_root_model tries first around a
// centre: the goal's relative width with five guard bits, or an eighth of the
// absolute tolerance where that is larger.
void goal_half_width(arb_t radius, const acb_t center, const cps_options &options, slong precision) {
    upper(radius, center);
    arb_mul_2exp_si(radius, radius, -resources::precision_sum(options.bits, 5));
    if (!options.absolute_tolerance.zero()) {
        RealBall tolerance;
        options.absolute_tolerance.enclose(tolerance.value(), precision);
        arb_mul_2exp_si(tolerance.value(), tolerance.value(), -3);
        if (arb_gt(tolerance.value(), radius)) arb_set(radius, tolerance.value());
    }
}

// A model-remainder refusal says how much stronger the model must be: the
// refused error exceeds its allowance by this many bits.
void note_model_shortfall(ChartCertificateStatistics *statistics, const RealBall &error, const RealBall &allowance) {
    if (!statistics) return;
    Magnitude above, below;
    arb_get_mag(above.value(), error.value());
    arb_get_mag_lower(below.value(), allowance.value());
    if (mag_is_zero(below.value()) || mag_is_inf(above.value())) return;
    const double bits = mag_get_d_log2_approx(above.value()) - mag_get_d_log2_approx(below.value());
    if (!(bits >= 0) || !std::isfinite(bits)) return;
    const auto shortfall = static_cast<uint64_t>(std::ceil(bits)) + 1;
    statistics->round_model_shortfall_bits = std::max(statistics->round_model_shortfall_bits, shortfall);
    const auto bucket = std::min<size_t>(static_cast<size_t>(std::bit_width(shortfall - 1)),
                                         statistics->round_shortfall_log2.size() - 1);
    ++statistics->round_shortfall_log2[bucket];
}

bool certify_root_model(Root &root, const RootModel &model, const acb_t local_point, const cps_options &options,
                        slong precision, ChartFailure *failure, ChartCertificateStatistics *statistics,
                        double *correction_log2, ChartCertificateDiagnostic *diagnostic,
                        slong experimental_jet_order, bool defer_guard_jet, ChartCertificateWork *work,
                        CertificateJet *shared_jet) {
    if (experimental_jet_order < 0)
        throw std::invalid_argument("negative certificate jet order");
    if (correction_log2)
        *correction_log2 = std::numeric_limits<double>::infinity();
    auto fail = [&](ChartFailure why) {
        if (failure)
            *failure = why;
        return false;
    };
    if (options.cancelled() || !acb_is_finite(local_point) || !arb_is_positive(model.radius.value()))
        return fail(ChartFailure::Domain);
    RealBall local_distance, one;
    arb_one(one.value());
    upper(local_distance.value(), local_point);
    if (!arb_lt(local_distance.value(), one.value()))
        return fail(ChartFailure::Domain);
    ComplexBall center, local, delta;
    acb_mul_arb(center.value(), local_point, model.radius.value(), precision);
    acb_add(center.value(), center.value(), model.center.value(), precision);
    acb_get_mid(center.value(), center.value());
    if (!acb_is_finite(center.value()) || acb_is_zero(center.value()))
        return fail(ChartFailure::Domain);
    RealBall radius, magnitude;
    goal_half_width(radius.value(), center.value(), options, precision);
    acb_sub(local.value(), center.value(), model.center.value(), precision);
    acb_div_arb(local.value(), local.value(), model.radius.value(), precision);
    const slong degree = acb_poly_degree(model.polynomial.value());
    if (degree < 1)
        return false;
    // Small output boxes need only a cubic jet. Larger retirement guards use
    // more terms to retain cancellation near the root. Both retain a rigorous
    // derivative majorant for every remaining coefficient.
    const slong order = std::min<slong>(degree, experimental_jet_order ? experimental_jet_order :
                                               options.absolute_tolerance.zero() ? certificate_jet_order
                                                                                 : guard_jet_order);
    // Singleton guards can sometimes be refused from value/slope/model error
    // alone. Keep the complete jet for all multi-root predicates and detailed
    // diagnostics, whose component accounting requires those terms.
    const bool deferred = defer_guard_jet && options.goal != CPS_ROOT_APPROXIMATION &&
                          order > 1 && !statistics && !diagnostic;
    // shared_jet, when given, holds the longest jet computed so far at this point of this model at this precision
    // (the caller's promise: process_candidate certifies a coarse guard and then the final rectangle at one point).
    // Entry k of a jet depends on entries k and k - 1 alone, so a longer jet's head is the shorter jet bit for
    // bit and is taken from it; the ledger is charged as if it were computed.
    std::vector<ComplexBall> jet;
    const auto charge_jet = [&](slong length) {
        if (work)
            work->jets += static_cast<double>(length + 1) * static_cast<double>(degree + 1) * operation_weight(precision);
    };
    const auto evaluate_jet = [&](slong length) {
        charge_jet(length);
        if (shared_jet && shared_jet->jet.size() >= static_cast<size_t>(length + 1)) {
            jet.assign(shared_jet->jet.begin(), shared_jet->jet.begin() + (length + 1));
            if (statistics) ++statistics->jets_reused;
            return;
        }
        jet.clear();
        jet.resize(static_cast<size_t>(length + 1));
        for (slong j = degree; j >= 0; --j) {
            for (size_t k = static_cast<size_t>(length); k > 0; --k) {
                acb_mul(jet[k].value(), jet[k].value(), local.value(), precision);
                acb_add(jet[k].value(), jet[k].value(), jet[k - 1].value(), precision);
            }
            acb_mul(jet[0].value(), jet[0].value(), local.value(), precision);
            acb_add(jet[0].value(), jet[0].value(), model.polynomial.value()->coeffs + j, precision);
        }
        if (shared_jet && shared_jet->jet.size() < jet.size()) {
            shared_jet->jet = jet;
            shared_jet->extensible = false;
        }
    };
    // A singleton guard takes its jet by rungs (below): the same jet computed a column at a time. With h_k[j]
    // accumulator k after the coefficients of index j and above,
    //     h_k[j] = h_k[j + 1] local + h_(k-1)[j + 1],   h_0[j] = h_0[j + 1] local + c_j,   h_k[degree + 1] = 0,
    // which are the row-wise loop's operations in another order, so the same balls, and jet[k] = h_k[0]. The last
    // column is kept (shared_jet->top), and one more accumulator costs that accumulator's column alone.
    const bool laddered = deferred && shared_jet;
    const auto extend_jet = [&](slong length) {
        auto &kept = *shared_jet;
        const size_t rows = static_cast<size_t>(degree) + 2;
        if (!kept.extensible) kept.jet.clear();
        if (kept.top.size() < rows) {
            kept.top.resize(rows);
            kept.below.resize(rows);
        }
        while (kept.jet.size() < static_cast<size_t>(length + 1)) {
            const bool first = kept.jet.empty();
            std::swap(kept.top, kept.below);
            acb_zero(kept.top[static_cast<size_t>(degree) + 1].value());
            for (slong j = degree; j >= 0; --j) {
                const size_t row = static_cast<size_t>(j);
                acb_mul(kept.top[row].value(), kept.top[row + 1].value(), local.value(), precision);
                if (first)
                    acb_add(kept.top[row].value(), kept.top[row].value(), model.polynomial.value()->coeffs + j, precision);
                else
                    acb_add(kept.top[row].value(), kept.top[row].value(), kept.below[row + 1].value(), precision);
            }
            kept.jet.push_back(kept.top[0]);
            kept.extensible = true;
        }
        jet.assign(kept.jet.begin(), kept.jet.begin() + (length + 1));
    };
    bool full_jet_charged = false;
    if (laddered) {
        charge_jet(1);
        extend_jet(1);
    } else
        evaluate_jet(deferred ? 1 : order);
    if (correction_log2) {
        // Reuse the proof's Taylor jet: the correction in original coordinates
        // is radius * F(w)/F'(w). Only its magnitude level guides continuation;
        // no root count, containment or stopping certificate uses this estimate.
        Magnitude value_size, slope_size, scale_size;
        acb_get_mag(value_size.value(), jet[0].value());
        acb_get_mag_lower(slope_size.value(), jet[1].value());
        arb_get_mag(scale_size.value(), model.radius.value());
        if (!mag_is_zero(slope_size.value()))
            *correction_log2 = mag_get_d_log2_approx(value_size.value()) -
                mag_get_d_log2_approx(slope_size.value()) + mag_get_d_log2_approx(scale_size.value());
    }
    RealBall rounding, imaginary, constant;
    arb_get_rad_arb(rounding.value(), acb_realref(jet[0].value()));
    arb_get_rad_arb(imaginary.value(), acb_imagref(jet[0].value()));
    arb_add(rounding.value(), rounding.value(), imaginary.value(), precision);
    upper(constant.value(), jet[0].value());
    // Each halving is bounded by the available endpoint precision. Failure of
    // the constant/model floor requests better candidates or arithmetic instead
    // of shrinking indefinitely. The translated jet is reused throughout.
    for (slong shrink = 0; shrink <= precision && !options.cancelled(); ++shrink) {
        if (work)
            work->shrinks += operation_weight(precision);
        root.box = center;
        arb_add_error(acb_realref(root.box.value()), radius.value());
        arb_add_error(acb_imagref(root.box.value()), radius.value());
        root.count = root.multiplicity = 1;
        root.realness = CPS_REAL_UNKNOWN;
        root.publish(precision);
        if (!meets_goal(root, options, precision))
            return fail(ChartFailure::Rounding);
        acb_sub(delta.value(), root.box.value(), model.center.value(), precision);
        acb_div_arb(delta.value(), delta.value(), model.radius.value(), precision);
        upper(magnitude.value(), delta.value());
        bool domain = arb_lt(magnitude.value(), one.value()) &&
                      (!model.lower_index || !acb_contains_zero(root.box.value()));
        RealBall minimum, maximum, endpoint, distance;
        for (size_t k = 0; k < 4; ++k) {
            endpoint_ball(endpoint.value(), root.endpoints[k], precision);
            arb_sub(distance.value(), endpoint.value(),
                    k < 2 ? acb_realref(center.value()) : acb_imagref(center.value()), precision);
            if (k % 2 == 0) arb_neg(distance.value(), distance.value());
            if (!arb_is_positive(distance.value()))
                return fail(ChartFailure::Rounding);
            if (!k) minimum = distance;
            else arb_min(minimum.value(), minimum.value(), distance.value(), precision);
        }
        acb_sub(delta.value(), root.box.value(), center.value(), precision);
        upper(maximum.value(), delta.value());
        arb_div(minimum.value(), minimum.value(), model.radius.value(), precision);
        arb_div(maximum.value(), maximum.value(), model.radius.value(), precision);
        if (domain) {
            RealBall band_error, taylor_error;
            const auto model_error = root_model_error(model, root.box.value(), precision, &band_error, statistics ? &statistics->model_error : nullptr, &taylor_error);
            if (work) {
                work->bounds += operation_weight(precision);
                // The jet's tail below: one term per omitted coefficient.
                work->jets += static_cast<double>(degree - order) * operation_weight(precision);
            }
            if (jet.size() < static_cast<size_t>(order + 1)) {
                Float slope;
                RealBall linear, floor, allowance, power;
                acb_get_abs_lbound_arf(slope.value(), jet[1].value(), precision);
                arb_set_arf(linear.value(), slope.value());
                // Match the full predicate's operation order, including its
                // radius power, so borderline refusal classifications agree.
                arb_pow_ui(power.value(), minimum.value(), 1, precision);
                arb_mul(linear.value(), linear.value(), power.value(), precision);
                arb_add(floor.value(), constant.value(), model_error.value(), precision);
                // Every singleton Rouche RHS includes this floor. If it
                // already exceeds the linear LHS, the full jet cannot prove
                // this rectangle. Apply the same stopping classifications as
                // below; if none applies, retain the ordinary computation.
                if (arb_ge(floor.value(), linear.value())) {
                    arb_mul_2exp_si(allowance.value(), linear.value(), -2);
                    if (!arb_is_zero(model_error.value()) && !arb_lt(model_error.value(), allowance.value())) {
                        note_model_shortfall(statistics, model_error, allowance);
                        return fail(ChartFailure::ModelRemainder);
                    }
                    if (!arb_is_zero(rounding.value()) && !arb_lt(rounding.value(), allowance.value()))
                        return fail(ChartFailure::Rounding);
                    if (!arb_lt(constant.value(), allowance.value()))
                        return fail(ChartFailure::Candidate);
                }
            }
            if (diagnostic) {
                Magnitude bound;
                arb_get_mag(bound.value(), model_error.value());
                diagnostic->error_log2 = mag_get_d_log2_approx(bound.value());
            }
            // Rouche on the published rectangle with the first `used` orders of the jet in hand: true when it
            // proves a count. `linear` receives the singleton test's left side, which does not depend on the
            // order; `hopeless`, when asked for, whether that side fails even against the terms already exact
            // and the model's error, which no longer jet can change.
            RealBall linear;
            const auto rouche = [&](const slong used, bool *hopeless) -> bool {
                // Taylor integral remainder bounds every omitted jet coefficient
                // on the whole published rectangle, not just at its midpoint.
                RealBall radial, tail, term, power;
                upper(radial.value(), local.value());
                arb_add(radial.value(), radial.value(), maximum.value(), precision);
                // The binomials C(k, used + 1), k from the degree down: in a word up to word_binomial_degree,
                // in an fmpz above it (word_binomial).
                const bool in_word = degree <= word_binomial_degree;
                Integer binomial;
                uint128 word = 0;
                if (in_word) word = word_binomial(static_cast<ulong>(degree), static_cast<ulong>(used + 1));
                else fmpz_bin_uiui(binomial.value(), static_cast<ulong>(degree), static_cast<ulong>(used + 1));
                for (slong k = degree; k > used; --k) {
                    upper(term.value(), model.polynomial.value()->coeffs + k);
                    if (in_word) arb_mul_ui(term.value(), term.value(), static_cast<ulong>(word), precision);
                    else arb_mul_fmpz(term.value(), term.value(), binomial.value(), precision);
                    arb_mul(tail.value(), tail.value(), radial.value(), precision);
                    arb_add(tail.value(), tail.value(), term.value(), precision);
                    if (in_word) {
                        word = word * static_cast<ulong>(k - used - 1) / static_cast<ulong>(k);
                    } else {
                        fmpz_mul_si(binomial.value(), binomial.value(), k - used - 1);
                        fmpz_divexact_si(binomial.value(), binomial.value(), k);
                    }
                }
                arb_pow_ui(power.value(), maximum.value(), static_cast<ulong>(used + 1), precision);
                arb_mul(tail.value(), tail.value(), power.value(), precision);
                RealBall jet_tail;
                if ((statistics && model.error_components) || diagnostic) jet_tail = tail;
                arb_add(tail.value(), tail.value(), model_error.value(), precision);
                std::vector<RealBall> weighted(static_cast<size_t>(used + 1));
                arb_one(power.value());
                for (slong k = 0; k <= used; ++k) {
                    upper(term.value(), jet[static_cast<size_t>(k)].value());
                    arb_mul(weighted[static_cast<size_t>(k)].value(), term.value(), power.value(), precision);
                    arb_mul(power.value(), power.value(), maximum.value(), precision);
                }
                const slong last = options.goal == CPS_ROOT_APPROXIMATION ? used : 1;
                for (slong m = 1; m <= last; ++m) {
                    Float lower;
                    RealBall lhs, rhs = tail;
                    acb_get_abs_lbound_arf(lower.value(), jet[static_cast<size_t>(m)].value(), precision);
                    arb_set_arf(lhs.value(), lower.value());
                    arb_pow_ui(power.value(), minimum.value(), static_cast<ulong>(m), precision);
                    arb_mul(lhs.value(), lhs.value(), power.value(), precision);
                    for (slong k = 0; k <= used; ++k)
                        if (k != m)
                            arb_add(rhs.value(), rhs.value(), weighted[static_cast<size_t>(k)].value(), precision);
                    if (diagnostic && m == 1) {
                        Magnitude left, right;
                        arb_get_mag_lower(left.value(), lhs.value());
                        arb_get_mag(right.value(), rhs.value());
                        diagnostic->margin_log2 = mag_get_d_log2_approx(left.value()) - mag_get_d_log2_approx(right.value());
                        diagnostic->component_terms_available = model.error_components;
                        const std::array<const RealBall *, 5> terms{&band_error, &taylor_error,
                                                                   &jet_tail, &rounding, &constant};
                        for (size_t k = 0; k < terms.size(); ++k) {
                            Magnitude size;
                            arb_get_mag(size.value(), terms[k]->value());
                            diagnostic->component_log2[k] = mag_get_d_log2_approx(size.value());
                        }
                    }
                    if (m == 1 && hopeless) {
                        // What no longer jet can remove: the terms already exact and the model's error.
                        RealBall floor = model_error;
                        for (slong k = 0; k <= used; ++k)
                            if (k != m)
                                arb_add(floor.value(), floor.value(), weighted[static_cast<size_t>(k)].value(), precision);
                        *hopeless = !arb_gt(lhs.value(), floor.value());
                    }
                    if (arb_gt(lhs.value(), rhs.value())) {
                        root.count = static_cast<uint64_t>(m);
                        root.multiplicity = m == 1 ? 1 : 0;
                        root.certificate = "piecewise Rouche on published rectangle, count " + std::to_string(m);
                        return true;
                    }
                    if (m == 1 && statistics && model.error_components) {
                        ++statistics->error_samples;
                        RealBall allowance, without_band = jet_tail;
                        arb_mul_2exp_si(allowance.value(), lhs.value(), -2);
                        if (!arb_is_zero(band_error.value()) &&
                            !arb_lt(band_error.value(), allowance.value()))
                            ++statistics->omitted_band_limited;
                        if (!arb_is_zero(taylor_error.value()) &&
                            !arb_lt(taylor_error.value(), allowance.value()))
                            ++statistics->model_tail_limited;
                        if (!arb_is_zero(jet_tail.value()) && !arb_lt(jet_tail.value(), allowance.value()))
                            ++statistics->jet_tail_limited;
                        if (!arb_is_zero(rounding.value()) && !arb_lt(rounding.value(), allowance.value()))
                            ++statistics->arithmetic_limited;
                        arb_add(without_band.value(), without_band.value(), taylor_error.value(), precision);
                        for (slong k = 0; k <= used; ++k)
                            if (k != m) arb_add(without_band.value(), without_band.value(), weighted[static_cast<size_t>(k)].value(), precision);
                        // This counter is a hypothetical opportunity, not a proof:
                        // would the same rectangle pass if omitted bands vanished?
                        if (!arb_is_zero(band_error.value()) && arb_gt(lhs.value(), without_band.value()))
                            ++statistics->omitted_band_decisive;
                    }
                    if (m == 1) linear = lhs;
                }
                return false;
            };
            bool tested = false;
            if (jet.size() < static_cast<size_t>(order + 1)) {
                // The full jet is charged once a call, as the decisions that read the ledger were qualified.
                if (!full_jet_charged && laddered) charge_jet(order);
                full_jet_charged = true;
                if (!laddered)
                    evaluate_jet(order);
                else {
                    // The guard's ladder: the final certificate's order first (its jet is then the head of this
                    // one), then twice the entries each time, up to the guard's order. A rung that proves the
                    // rectangle ends it; so does a rung whose failure no longer jet can cure, and the
                    // classification below is the one the full order would reach (its test's right side is at
                    // least that rung's floor).
                    for (slong rung = std::max<slong>(std::min(order, certificate_jet_order),
                                                     static_cast<slong>(jet.size()) - 1);;
                         rung = std::min(order, 2 * rung + 1)) {
                        extend_jet(rung);
                        bool hopeless = false;
                        if (rouche(rung, rung < order ? &hopeless : nullptr)) return true;
                        if (hopeless || rung == order) break;
                    }
                    tested = true;
                }
            }
            if (!tested && rouche(order, nullptr)) return true;
            RealBall allowance;
            arb_mul_2exp_si(allowance.value(), linear.value(), -2);
            if (!arb_is_zero(model_error.value()) && !arb_lt(model_error.value(), allowance.value())) {
                note_model_shortfall(statistics, model_error, allowance);
                return fail(ChartFailure::ModelRemainder);
            }
            if (!arb_is_zero(rounding.value()) && !arb_lt(rounding.value(), allowance.value()))
                return fail(ChartFailure::Rounding);
            if (!arb_lt(constant.value(), allowance.value()))
                return fail(ChartFailure::Candidate);
        }
        arb_mul_2exp_si(radius.value(), radius.value(), -1);
        if (statistics) ++statistics->radius_shrinks;
    }
    if (statistics) ++statistics->width_limited;
    return fail(ChartFailure::Width);
}
} // namespace certpolysolve::detail
