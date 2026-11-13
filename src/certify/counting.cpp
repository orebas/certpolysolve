// Root counts in disks by Pellet's test with Graeffe root squaring (pellet_count), and the certificate of a counted
// square built from two such counts (certify_cluster_box); counting.hpp has DiskCount, CountContext and the contracts.
//
// pellet_count shifts p to the exact centre c and scales the radius r to one: q(w) = p(c + r w). If the modulus of one
// coefficient, of degree k, exceeds the sum of the others' in strict outward ball comparisons (dominant_coefficient),
// Rouche's theorem gives exactly k roots in |z - c| < r, counting multiplicity, and none on the circle, uniformly over
// the coefficient balls; a failed comparison is inconclusive, never a count of zero. Otherwise Graeffe steps square the
// roots, each staying on its side of the unit circle, and the test is repeated (count_shifted; up to the caller's
// max_graeffe steps). An adaptive context first tries mag_t bounds (magnitude_filter) and returns as soon as ball width obstructs
// a comparison, for the source to be re-enclosed at a higher precision. A direct test that width obstructed at a centre
// off both axes is retried once, under CountContext::real_shift, on rotated_real_shift: p(u (x + |c|)) with u = c/|c|
// has the same coefficient moduli and disk counts, and its shift by a real number does not re-box Arb's rectangular
// balls. Under a context, each Taylor shift and Graeffe step is priced before it runs (dovetail.hpp); one that the
// allowance (dovetailing) or the caller's cancellation does not cover suspends the test, and under a finite allowance a
// call with the same arguments resumes it. certify_cluster_box proves one count in the disks of radius r and 2r and
// checks that the exported square of half-side r lies strictly inside the outer disk, so the square holds that count.
// It claims the factor's multiplicity for one root only, except where c = 0 and the first count coefficients vanish
// exactly: x^count then divides p, and the region is the point 0 of multiplicity count (times the factor's). One root
// with p(c) = 0 is an exact point. Realness: nonreal off the axis, real for one root about a real centre of a real p.
//
// State: none here; a test's meter, allowance, suspension and last completed tests are in the caller's CountContext.
// Entry points (callers): pellet_count (subdivision.cpp, and certify_cluster_box); certify_cluster_box (clustering.cpp,
// subdivision.cpp, certification.cpp); rotated_real_shift (in src only pellet_count; the tests call it too).
// Not here: recover_clusters, declared in counting.hpp, is in clustering.cpp; Subdivision::step sets the allowance.
// Design notes: docs/mathematics.md (local counts), docs/selector-campaign-20260928.md (class G: the rotated real
// shift), docs/principled-audit-20260929.md (dovetailing as built).
#include "certify/counting.hpp"
#include "core/dovetail.hpp"
#include "core/service_meter.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <limits>

namespace certpolysolve::detail {
namespace {
slong magnitude_filter(const acb_poly_t p, slong expected) {
    mag_t lower, largest, term, sum;
    mag_init(lower);
    mag_init(largest);
    mag_init(term);
    mag_init(sum);
    slong candidate = expected >= 0 ? expected : 0;
    if (expected < 0)
        for (slong j = 0; j < acb_poly_length(p); ++j) {
            acb_get_mag_lower(lower, p->coeffs + j);
            if (mag_cmp(lower, largest) > 0) {
                candidate = j;
                mag_set(largest, lower);
            }
        }
    acb_get_mag_lower(lower, p->coeffs + candidate);
    for (slong j = 0; j < acb_poly_length(p); ++j)
        if (j != candidate) {
            acb_get_mag(term, p->coeffs + j);
            mag_add(sum, sum, term);
        }
    const bool proved = mag_cmp(lower, sum) > 0;
    mag_clear(lower);
    mag_clear(largest);
    mag_clear(term);
    mag_clear(sum);
    return proved ? candidate : -1;
}
bool comparison_needs_precision(const arb_t candidate, const arb_t others, slong precision) {
    if (arb_le(candidate, others))
        return false;
    // Near equality has a finite soft refusal; it need not be resolved.
    RealBall twice, thrice, others_twice, others_thrice;
    arb_mul_ui(twice.value(), candidate, 2, precision);
    arb_mul_ui(thrice.value(), candidate, 3, precision);
    arb_mul_ui(others_twice.value(), others, 2, precision);
    arb_mul_ui(others_thrice.value(), others, 3, precision);
    return !(arb_lt(twice.value(), others_thrice.value()) && arb_gt(thrice.value(), others_twice.value()));
}
DiskCount dominant_coefficient(const acb_poly_t polynomial, slong expected, slong precision, bool direct) {
    const slong length = acb_poly_length(polynomial);
    if (length == 0 || expected >= length)
        return {};
    RealVector magnitudes(length);
    Float lower, largest_lower;
    slong candidate = expected >= 0 ? expected : 0;
    for (slong j = 0; j < length; ++j) {
        acb_abs(magnitudes.data() + j, polynomial->coeffs + j, precision);
        if (!arb_is_finite(magnitudes.data() + j))
            return {-1, 0, true, direct, false, direct && expected <= 0};
        if (expected < 0) {
            arb_get_lbound_arf(lower.value(), magnitudes.data() + j, precision);
            if (arf_cmp(lower.value(), largest_lower.value()) > 0) {
                candidate = j;
                arf_set(largest_lower.value(), lower.value());
            }
        }
    }
    RealBall others;
    for (slong j = 0; j < length; ++j)
        if (j != candidate)
            arb_add(others.value(), others.value(), magnitudes.data() + j, precision);
    // Interval comparison is strict and outward. A failed comparison says
    // nothing about the count; no midpoint or residual is trusted here.
    if (arb_gt(magnitudes.data() + candidate, others.value()))
        return {candidate};
    const bool selected_needs_precision =
        comparison_needs_precision(magnitudes.data() + candidate, others.value(), precision);
    bool width_limited = direct && selected_needs_precision;
    if (expected < 0) {
        // Keep direct zero exclusion separate from the selected coefficient
        // and from the all-index request below. Uncertainty about a positive
        // count does not mean smaller cells are useless at a precision cap.
        // Reuse the magnitudes; no second Taylor shift or modulus pass.
        bool exclusion_limited = direct && candidate == 0 && selected_needs_precision;
        if (direct && candidate != 0) {
            RealBall zero_others;
            for (slong j = 1; j < length; ++j)
                arb_add(zero_others.value(), zero_others.value(), magnitudes.data() + j, precision);
            exclusion_limited = comparison_needs_precision(magnitudes.data(), zero_others.value(), precision);
        }
        width_limited = width_limited || exclusion_limited;
        // A decisive rejection of the selected lower bound says nothing
        // about another coefficient hidden by uncertainty. Refuse without
        // more precision only when every possible coefficient is below
        // 3/2 of the others: 5 max(u_i) < 3 sum(l_i) suffices, since
        // 2 u_i + 3 l_i <= 5 max(u_i). Clamp modulus lower bounds at zero.
        // The strict count above still decides every accepted certificate.
        Float upper, largest_upper, sum_lower;
        for (slong j = 0; j < length; ++j) {
            arb_get_lbound_arf(lower.value(), magnitudes.data() + j, precision);
            if (arf_sgn(lower.value()) > 0)
                arf_add(sum_lower.value(), sum_lower.value(), lower.value(), precision, ARF_RND_FLOOR);
            arb_get_ubound_arf(upper.value(), magnitudes.data() + j, precision);
            if (arf_cmp(upper.value(), largest_upper.value()) > 0)
                arf_set(largest_upper.value(), upper.value());
        }
        arf_mul_ui(largest_upper.value(), largest_upper.value(), 5, precision, ARF_RND_CEIL);
        arf_mul_ui(sum_lower.value(), sum_lower.value(), 3, precision, ARF_RND_FLOOR);
        DiskCount result{-1, 0, arf_cmp(largest_upper.value(), sum_lower.value()) >= 0, width_limited};
        result.exclusion_limited = exclusion_limited;
        return result;
    }
    DiskCount result{-1, 0, selected_needs_precision, width_limited};
    result.exclusion_limited = expected == 0 && width_limited;
    return result;
}
} // namespace

namespace {
// The Pellet predicate on shifted = p(center + y) (or a rotation of it, which
// leaves every coefficient modulus unchanged): q(w) = shifted(radius*w).
// A resumed test (CountContext::suspension) passes its scaled polynomial after
// `first_step` Graeffe steps, and the direct comparison's width flag.
DiskCount count_shifted(acb_poly_t shifted, const arf_t radius, slong precision, unsigned max_graeffe,
                        slong expected, CountContext *context, unsigned first_step = 0, bool scaled = false,
                        bool carried_width_limited = false, bool carried_exclusion_limited = false) {
    const slong degree = acb_poly_degree(shifted);
    if (!scaled) {
        RealBall power, scale;
        arb_one(power.value());
        arb_set_arf(scale.value(), radius);
        for (slong j = 0; j <= degree; ++j) {
            acb_mul_arb(shifted->coeffs + j, shifted->coeffs + j, power.value(), precision);
            arb_mul(power.value(), power.value(), scale.value(), precision);
        }
    }
    // Strict Pellet dominance on the unit circle gives a count by Rouche.
    // Graeffe squares the roots, preserving whether their modulus is below,
    // on, or above one, and preserving multiplicity.
    bool width_limited = carried_width_limited;
    bool exclusion_limited = carried_exclusion_limited;
    for (unsigned step = first_step; step <= max_graeffe; ++step) {
        if (context && context->adaptive) {
            const auto filtered = magnitude_filter(shifted, expected);
            if (filtered >= 0)
                return {filtered, step, false};
        }
        auto count = dominant_coefficient(shifted, expected, precision, step == 0);
        count.graeffe_steps = step;
        if (step == 0) {
            width_limited = count.width_limited;
            exclusion_limited = count.exclusion_limited;
        }
        count.width_limited = width_limited;
        count.exclusion_limited = exclusion_limited;
        if (count.count >= 0 || step == max_graeffe)
            return count;
        // Re-enclose the immutable source before further root squaring when
        // arithmetic uncertainty obstructs the comparison. This is scheduling
        // information only; an inconclusive predicate never excludes a cell.
        if (context && context->adaptive && count.needs_precision)
            return count;
        if (context) {
            const auto service_work = context->research_service_policy ?
                service_graeffe_work(shifted, precision) : ServicePolynomialWork{};
            const double cost =
                context->research_service_policy ? service_graeffe_price(service_work, precision) /
                    operation_weight(precision) :
                block_product_operations(shifted, precision);
            if (!context->covers(cost)) {
                // Stop before this step: the comparison at `step` is repeated
                // on resumption (it costs a pass over the coefficients).
                context->pending = cost;
                DiskCount stopped;
                stopped.suspended = true;
                stopped.graeffe_steps = step;
                stopped.width_limited = width_limited;
                stopped.exclusion_limited = exclusion_limited;
                return stopped;
            }
            ++context->statistics.graeffe_steps;
            context->metered += cost;
            service_observe(context->research_service_policy, ServiceOperation::Graeffe,
                            service_price(ServiceOperation::Graeffe, service_work.scalar * operation_weight(precision)), precision);
            service_observe(context->research_service_policy, ServiceOperation::GraeffeWords,
                            service_price(ServiceOperation::GraeffeWords,
                                service_work.transform * std::ceil(static_cast<double>(precision) / FLINT_BITS)), precision);
        }
        TraceScope service_trace(context ? context->research_service_policy : nullptr,
                                 "service_graeffe", static_cast<uint64_t>(precision), static_cast<uint64_t>(acb_poly_length(shifted)));
        acb_poly_graeffe_transform(shifted, shifted, precision);
    }
    return {};
}
} // namespace

void rotated_real_shift(acb_poly_t out, const acb_poly_t polynomial, const acb_t center, slong precision) {
    const slong degree = acb_poly_degree(polynomial);
    RealBall rho;
    acb_abs(rho.value(), center, precision);
    ComplexBall u, shift;
    acb_div_arb(u.value(), center, rho.value(), precision);
    // Powers by squaring and single products (logarithmic depth), so the
    // rotation factors stay within a few roundings of the unit circle.
    RootVector powers(degree + 1);
    _acb_vec_set_powers(powers.data(), u.value(), degree + 1, precision);
    acb_poly_set(out, polynomial);
    for (slong j = 0; j <= degree; ++j)
        acb_mul(out->coeffs + j, out->coeffs + j, powers[j], precision);
    acb_set_arb(shift.value(), rho.value());
    acb_poly_taylor_shift(out, out, shift.value(), precision);
}

namespace {
bool same_test(const CountContext::TestKey &key, const acb_poly_t polynomial, const acb_t center, const arf_t radius,
               slong precision, unsigned max_graeffe, slong expected) {
    return key.precision == precision && key.max_graeffe == max_graeffe && key.expected == expected &&
           arf_equal(arb_midref(key.radius.value()), radius) && acb_equal(key.center.value(), center) &&
           acb_poly_equal(key.source.value(), polynomial);
}
CountContext::TestKey test_key(const acb_poly_t polynomial, const acb_t center, const arf_t radius, slong precision,
                               unsigned max_graeffe, slong expected) {
    CountContext::TestKey key;
    acb_poly_set(key.source.value(), polynomial);
    acb_set(key.center.value(), center);
    arb_set_arf(key.radius.value(), radius);
    key.precision = precision;
    key.max_graeffe = max_graeffe;
    key.expected = expected;
    return key;
}
} // namespace

DiskCount pellet_count(const acb_poly_t polynomial, const acb_t center, const arf_t radius, slong precision,
                       unsigned max_graeffe, slong expected, CountContext *context) {
    const auto started = std::chrono::steady_clock::now();
    uint64_t shift_time = 0;
    struct Accounting {
        CountContext *context;
        std::chrono::steady_clock::time_point start;
        uint64_t &shift;
        ~Accounting() {
            if (context)
                context->statistics.counting_ns +=
                    static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                              std::chrono::steady_clock::now() - start)
                                              .count()) -
                    shift;
        }
    } accounting{context, started, shift_time};
    if (context)
        context->suspended = false;
    const slong degree = acb_poly_degree(polynomial);
    if (degree < 0 || expected < -1 || expected > degree || precision < 2 || !acb_is_finite(center) ||
        !acb_is_exact(center) || !arf_is_finite(radius) || arf_sgn(radius) <= 0 ||
        acb_contains_zero(polynomial->coeffs + degree))
        return {};
    // The caller owns the transform allowance, including subdivision's
    // degree-dependent bound. Precision/allowance stops remain inconclusive.
    // Under an allowance (dovetailing): a test completed earlier in this step
    // is not computed again, and a suspended one resumes where it stopped.
    const bool bounded = context && context->bounded();
    std::optional<CountContext::Suspension> resumed;
    if (bounded) {
        for (const auto &[key, result] : context->completed)
            if (same_test(key, polynomial, center, radius, precision, max_graeffe, expected))
                return result;
        if (context->suspension &&
            same_test(context->suspension->key, polynomial, center, radius, precision, max_graeffe, expected)) {
            resumed = std::move(context->suspension);
            context->suspension.reset();
        }
    }
    if (context && !resumed)
        ++context->statistics.disk_tests;
    const auto suspend = [&](int stage, bool ready, const DiskCount &stopped, acb_poly_t shifted,
                             const DiskCount &direct) {
        CountContext::Suspension state;
        state.key = resumed ? std::move(resumed->key) : test_key(polynomial, center, radius, precision, max_graeffe,
                                                                  expected);
        state.stage = stage;
        state.shifted_ready = ready;
        state.step = stopped.graeffe_steps;
        state.width_limited = stopped.width_limited;
        state.exclusion_limited = stopped.exclusion_limited;
        if (ready)
            acb_poly_swap(state.shifted.value(), shifted);
        state.direct = direct;
        context->suspension = std::move(state);
        context->suspended = true;
        DiskCount out;
        out.suspended = true;
        return out;
    };
    const auto finish = [&](const DiskCount &result) {
        if (bounded) {
            if (context->completed.size() == CountContext::completed_tests)
                context->completed.erase(context->completed.begin());
            context->completed.emplace_back(test_key(polynomial, center, radius, precision, max_graeffe, expected),
                                            result);
        }
        return result;
    };
    DiskCount direct;
    if (!resumed || resumed->stage == 0) {
        ComplexPolynomial transformed;
        const bool ready = resumed && resumed->shifted_ready;
        if (ready) {
            acb_poly_swap(transformed.value(), resumed->shifted.value());
        } else {
            const auto quote = context && context->research_service_policy
                ? service_taylor_quote(polynomial, center, precision) : ServiceTaylorQuote{};
            const double cost = !context ? 0 : context->research_service_policy
                ? quote.price() / operation_weight(precision) :
                context->research_identity_shift_meter && acb_is_zero(center)
                    ? static_cast<double>(acb_poly_length(polynomial)) :
                block_taylor_shift_operations(polynomial, precision);
            if (context) {
                if (!context->covers(cost)) {
                    context->pending = cost;
                    return suspend(0, false, DiskCount{}, transformed.value(), DiskCount{});
                }
            }
            {
                StageTimer timer(shift_time);
                TraceScope service_trace(context ? context->research_service_policy : nullptr,
                                         "service_taylor", static_cast<uint64_t>(precision), static_cast<uint64_t>(acb_poly_length(polynomial)));
                acb_poly_taylor_shift(transformed.value(), polynomial, center, precision);
            }
            if (context) {
                ++context->statistics.taylor_shifts;
                context->statistics.taylor_shift_ns += shift_time;
                context->metered += cost;
                if (context->research_service_policy)
                    quote.observe(context->research_service_policy, precision);
            }
        }
        direct = count_shifted(transformed.value(), radius, precision, max_graeffe, expected, context,
                               ready ? resumed->step : 0, ready, ready && resumed->width_limited,
                               ready && resumed->exclusion_limited);
        if (direct.suspended)
            return suspend(0, true, direct, transformed.value(), DiskCount{});
    } else {
        direct = resumed->direct;
    }
    // Arb's complex balls are rectangles: a product by a centre off both axes
    // re-boxes the error, up to sqrt(2) wider than a disk, and the Taylor
    // shift compounds that through the powers of the centre until the widths
    // can obstruct the dominance comparison (class G: near_multiple-193
    // rotated by (3+4i)/5). When the widths did obstruct it, retry once on
    // S(x) = p(u (x + |c|)) = p(c + u x), u = c / |c|: one rounding per
    // coefficient for the rotation, then a shift by a real number, which
    // scales real and imaginary parts without re-boxing. The rotation x -> u x
    // preserves the count in |x| < radius and every coefficient modulus the
    // test compares, and the balls enclose S for the exact u. The retry is
    // the costlier shift (full-precision |c| instead of a short dyadic), so
    // it is paid only where the direct test was arithmetic-limited.
    // "The widths did obstruct it": an adaptive test returns at the first
    // comparison that width obstructs, so its flag says so. A test that runs
    // every Graeffe step returns the last step's flag, which width always
    // obstructs at that depth; there the direct comparison's width is the
    // evidence.
    const bool obstructed = context && (context->adaptive ? direct.needs_precision : direct.width_limited);
    if (direct.count >= 0 || !obstructed || !context->real_shift ||
        arb_is_zero(acb_realref(center)) || arb_is_zero(acb_imagref(center)))
        return finish(direct);
    ComplexPolynomial rotated;
    const bool ready = resumed && resumed->stage == 1 && resumed->shifted_ready;
    if (ready) {
        acb_poly_swap(rotated.value(), resumed->shifted.value());
    } else {
        // The rotation's products (one a coefficient) and the real shift.
        ServiceTaylorQuote quote;
        if (context->research_service_policy) {
            ComplexBall real_center;
            acb_abs(acb_realref(real_center.value()), center, precision);
            quote = service_taylor_quote(polynomial, real_center.value(), precision, true);
        }
        const double rotation = service_price(ServiceOperation::Taylor, static_cast<double>(degree + 1));
        const double cost =
            context->research_service_policy ? rotation + quote.price() / operation_weight(precision) :
            static_cast<double>(degree + 1) + block_taylor_shift_operations(polynomial, precision);
        if (!context->covers(cost)) {
            context->pending = cost;
            return suspend(1, false, DiskCount{}, rotated.value(), direct);
        }
        const uint64_t before_rotation = shift_time;
        {
            StageTimer timer(shift_time);
            TraceScope service_trace(context->research_service_policy, "service_taylor", static_cast<uint64_t>(precision),
                                     static_cast<uint64_t>(acb_poly_length(polynomial)));
            rotated_real_shift(rotated.value(), polynomial, center, precision);
        }
        ++context->statistics.taylor_shifts;
        context->statistics.taylor_shift_ns += shift_time - before_rotation;
        if (context->research_service_policy) {
            // Repay exactly the admitted rotation plus shift, not just the shift.
            context->metered += cost;
            service_observe(context->research_service_policy, ServiceOperation::Taylor,
                            rotation * operation_weight(precision), precision);
            quote.observe(context->research_service_policy, precision);
        } else
            context->metered += block_taylor_shift_operations(polynomial, precision);
    }
    const auto retried = count_shifted(rotated.value(), radius, precision, max_graeffe, expected, context,
                                       ready ? resumed->step : 0, ready, ready && resumed->width_limited,
                                       ready && resumed->exclusion_limited);
    if (retried.suspended)
        return suspend(1, true, retried, rotated.value(), direct);
    return finish(retried.count >= 0 ? retried : direct);
}

bool certify_cluster_box(Root &out, const acb_poly_t polynomial, const acb_t center, const arf_t radius,
                         slong count, uint64_t factor_multiplicity, slong precision, unsigned max_graeffe,
                         CountContext *context, bool *needs_precision) {
    if (needs_precision)
        *needs_precision = false;
    if (count < 1 || factor_multiplicity == 0 ||
        static_cast<uint64_t>(count) > std::numeric_limits<uint64_t>::max() / factor_multiplicity)
        return false;
    const auto inner = pellet_count(polynomial, center, radius, precision, max_graeffe, count, context);
    if (inner.count != count) {
        if (needs_precision)
            *needs_precision = inner.needs_precision;
        return false;
    }
    Float outer_radius;
    arf_mul_2exp_si(outer_radius.value(), radius, 1);
    const auto outer =
        pellet_count(polynomial, center, outer_radius.value(), precision, max_graeffe, count, context);
    if (outer.count != count) {
        if (needs_precision)
            *needs_precision = outer.needs_precision;
        return false;
    }

    Root root;
    acb_set(root.box.value(), center);
    acb_add_error_arf(root.box.value(), radius);
    root.publish(precision);
    if (!acb_is_finite(root.box.value())) {
        if (needs_precision)
            *needs_precision = true;
        return false;
    }
    // The original disk lies in the square before outward export. Conversely,
    // check that even the rounded exported square lies strictly in the outer
    // disk. This matters when the center and radius have extreme scale gaps.
    ComplexBall offset;
    Float maximum_distance;
    acb_sub(offset.value(), root.box.value(), center, precision);
    acb_get_abs_ubound_arf(maximum_distance.value(), offset.value(), precision);
    if (!arf_is_finite(maximum_distance.value()) || arf_cmp(maximum_distance.value(), outer_radius.value()) >= 0) {
        // Both disk counts already succeeded. The exact radius-r square has
        // diagonal radius sqrt(2)r < 2r; only endpoint/norm rounding still
        // obstructs its exported certificate, so request fresh precision.
        if (needs_precision)
            *needs_precision = true;
        return false;
    }

    root.count = static_cast<uint64_t>(count) * factor_multiplicity;
    root.multiplicity = count == 1 ? factor_multiplicity : 0;
    root.certificate =
        inner.graeffe_steps || outer.graeffe_steps ? "pellet-graeffe-local-count" : "pellet-local-count";
    if (count == 1) {
        ComplexBall value;
        acb_poly_evaluate(value.value(), polynomial, center, precision);
        // An exact zero at the exact center, together with the strict one-root
        // count, proves a point root uniformly. Relative accuracy at zero
        // cannot be obtained merely by shrinking a positive-radius box.
        if (acb_is_zero(value.value())) {
            acb_set(root.box.value(), center);
            root.publish(precision);
            root.certificate = "pellet-exact-point";
        }
    }
    // At zero, exact vanishing coefficients prove divisibility by x^count
    // without squarefree preprocessing. Combined with the strict disk count,
    // this proves that the entire cluster is this one multiple point. Merely
    // shrinking its radius cannot establish multiplicity or relative accuracy.
    if (acb_is_zero(center) && count <= acb_poly_degree(polynomial)) {
        bool zero_factor = true;
        for (slong j = 0; j < count; ++j)
            zero_factor = zero_factor && acb_is_zero(polynomial->coeffs + j);
        if (zero_factor) {
            acb_zero(root.box.value());
            root.multiplicity = root.count;
            root.realness = CPS_REAL;
            root.publish(precision);
            root.certificate = "pellet-exact-zero-factor";
        }
    }
    if (!arb_contains_zero(acb_imagref(root.box.value())))
        root.realness = CPS_NONREAL;
    // A unique root in a conjugation-invariant disk of a real polynomial is
    // real. This remains true uniformly for a family of real coefficients.
    if (count == 1 && acb_poly_is_real(polynomial) && arb_is_zero(acb_imagref(center))) {
        arb_zero(acb_imagref(root.box.value()));
        root.realness = CPS_REAL;
        root.endpoints[2] = {};
        root.endpoints[3] = {};
    }
    out = std::move(root);
    return true;
}
} // namespace certpolysolve::detail
