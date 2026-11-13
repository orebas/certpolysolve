// Certified winding count with dyadic arc subdivision.
#include "certify/winding_count.hpp"
#include "core/heuristics.hpp"
#include "numeric/evaluate/monic_source.hpp"
#include <map>
#include <utility>
#include <vector>

namespace certpolysolve::detail {
namespace {
// The circle's parameter t in [0, 2) in ticks of 2^-winding_tick_log2: z(t) = centre + radius exp(pi i t). The
// smallest arc, a 2^program_winding_arc_log2 of the circle, spans a whole number of ticks.
constexpr int winding_tick_log2 = 27;
constexpr int64_t winding_ticks = int64_t{2} << winding_tick_log2;
static_assert(winding_tick_log2 + 1 + heuristics::program_winding_arc_log2 >= 1, "the smallest arc spans ticks");
static_assert(winding_ticks % heuristics::program_winding_arcs == 0, "the first arcs span whole ticks");

void circle_point(acb_t z, const acb_t centre, const arb_t radius, int64_t tick, slong precision) {
    RealBall t, sine, cosine;
    arb_set_si(t.value(), static_cast<slong>(tick));
    arb_mul_2exp_si(t.value(), t.value(), -winding_tick_log2);
    arb_sin_cos_pi(sine.value(), cosine.value(), t.value(), precision);
    acb_set_arb_arb(z, cosine.value(), sine.value());
    acb_mul_arb(z, z, radius, precision);
    acb_add(z, z, centre, precision);
}

// The larger of a ball's two radii.
void ball_radius(mag_t out, const acb_t z) { mag_max(out, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z))); }
} // namespace

WindingCount certified_root_count(MonicSource &source, const acb_t centre, const mag_t radius, slong precision,
                                  uint64_t max_arcs, const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics) {
    WindingCount out;
    const slong degree = source.degree();
    ComplexBall c;
    acb_get_mid(c.value(), centre);
    if (degree < 0 || !mag_is_finite(radius) || mag_is_zero(radius) || !acb_is_finite(c.value()))
        return out;
    RealBall rho;
    arf_set_mag(arb_midref(rho.value()), radius); // exact: the circle is at the magnitude's own value
    const auto evaluated = [&] {
        ++out.evaluations;
        ++statistics.evaluations;
        statistics.evaluation_terms += source.terms();
    };
    // The values at the arcs' ends by tick (the circle's two ends are one point), each checked once: a value whose
    // argument is not known to a few bits (a radius above a 2^-6 of its modulus, or 0 in the ball) ends the count,
    // since every arc that ends there would carry that doubt.
    std::map<int64_t, ComplexBall> values;
    ComplexBall z, unused;
    Magnitude spread, lower;
    const auto value_at = [&](int64_t tick) -> const ComplexBall * {
        tick %= winding_ticks;
        if (const auto found = values.find(tick); found != values.end())
            return &found->second;
        circle_point(z.value(), c.value(), rho.value(), tick, precision);
        evaluated();
        ComplexBall value;
        if (!source.evaluate(value.value(), unused.value(), z.value(), precision) || !acb_is_finite(value.value()))
            return nullptr;
        ball_radius(spread.value(), value.value());
        mag_mul_2exp_si(spread.value(), spread.value(), 6);
        acb_get_mag_lower(lower.value(), value.value());
        if (mag_cmp(spread.value(), lower.value()) >= 0) {
            out.imprecise = true;
            return nullptr;
        }
        return &values.emplace(tick, std::move(value)).first->second;
    };
    RealBall pi, total, change, half_length;
    arb_const_pi(pi.value(), precision);
    arb_zero(total.value());
    ComplexBall middle, end, ball, value_over, derivative_over, enclosure, difference, quotient;
    Magnitude reach, ends, part, lower_value, derivative_size;
    constexpr int64_t first = winding_ticks / heuristics::program_winding_arcs;
    constexpr int64_t smallest = winding_ticks >> -heuristics::program_winding_arc_log2;
    std::vector<std::pair<int64_t, int64_t>> arcs; // a stack: the arcs in increasing t, the first on top
    for (int64_t a = winding_ticks - first; a >= 0; a -= first)
        arcs.emplace_back(a, a + first);
    while (!arcs.empty()) {
        if (options.cancelled()) {
            out.cancelled = true;
            return out;
        }
        const auto [a, b] = arcs.back();
        arcs.pop_back();
        if (out.arcs >= max_arcs)
            return out;
        ++out.arcs;
        const int64_t m = a + (b - a) / 2;
        const ComplexBall *at_a = value_at(a), *at_b = value_at(b), *at_m = value_at(m);
        if (!at_a || !at_b || !at_m)
            return out;
        // B: about the arc's midpoint, half the arc's length (rho pi (t_b - t_a) / 2, at least the distance from
        // the midpoint to any point of the arc) plus the radii of the balls z(t_m), z(t_a) and z(t_b).
        arb_mul_si(half_length.value(), pi.value(), static_cast<slong>(b - a), precision);
        arb_mul(half_length.value(), half_length.value(), rho.value(), precision);
        arb_mul_2exp_si(half_length.value(), half_length.value(), -(winding_tick_log2 + 1));
        arb_get_mag(reach.value(), half_length.value());
        circle_point(middle.value(), c.value(), rho.value(), m, precision);
        ball_radius(part.value(), middle.value());
        mag_add(reach.value(), reach.value(), part.value());
        mag_zero(ends.value());
        for (const int64_t tick : {a, b}) {
            circle_point(end.value(), c.value(), rho.value(), tick, precision);
            ball_radius(part.value(), end.value());
            mag_max(ends.value(), ends.value(), part.value());
        }
        mag_add(reach.value(), reach.value(), ends.value());
        acb_get_mid(ball.value(), middle.value());
        arb_add_error_mag(acb_realref(ball.value()), reach.value());
        arb_add_error_mag(acb_imagref(ball.value()), reach.value());
        evaluated();
        if (!source.evaluate(value_over.value(), derivative_over.value(), ball.value(), precision))
            return out;
        // The image of the arc lies in the disk about p(z_m) of radius |p'(B)| h, h half the arc's length (each
        // point w of the arc is within h of the exact midpoint, and the segment between them lies in B): a disk
        // that misses 0 when |p(z_m)| > |p'(B)| h, which is tighter than the box product below. It also lies in
        // the mean-value box p(z_m) + p'(B) (B - z_m) and in p(B).
        acb_get_mag_lower(lower_value.value(), at_m->value());
        acb_get_mag(derivative_size.value(), derivative_over.value());
        arb_get_mag(part.value(), half_length.value());
        mag_mul(derivative_size.value(), derivative_size.value(), part.value());
        bool resolved = acb_is_finite(derivative_over.value()) && mag_cmp(derivative_size.value(), lower_value.value()) < 0;
        if (!resolved) {
            acb_sub(difference.value(), ball.value(), middle.value(), precision);
            acb_mul(enclosure.value(), derivative_over.value(), difference.value(), precision);
            acb_add(enclosure.value(), enclosure.value(), at_m->value(), precision);
            resolved = (acb_is_finite(enclosure.value()) && !acb_contains_zero(enclosure.value())) ||
                       (acb_is_finite(value_over.value()) && !acb_contains_zero(value_over.value()));
        }
        if (resolved) {
            acb_conj(quotient.value(), at_a->value());
            acb_mul(quotient.value(), quotient.value(), at_b->value(), precision);
            acb_arg(change.value(), quotient.value(), precision);
            // A change placed within 2^-3 radians either way; a wider ball is a change near a half-turn, which the
            // principal value cannot place: a shorter arc turns less.
            resolved = mag_cmp_2exp_si(arb_radref(change.value()), -3) < 0;
        }
        if (!resolved) {
            if (b - a <= smallest)
                return out;
            arcs.emplace_back(m, b);
            arcs.emplace_back(a, m);
            continue;
        }
        arb_add(total.value(), total.value(), change.value(), precision);
    }
    arb_div(total.value(), total.value(), pi.value(), precision);
    arb_mul_2exp_si(total.value(), total.value(), -1);
    Integer count;
    if (arb_get_unique_fmpz(count.value(), total.value()) && fmpz_sgn(count.value()) >= 0 &&
        fmpz_cmp_si(count.value(), degree) <= 0)
        out.roots = fmpz_get_ui(count.value());
    return out;
}

} // namespace certpolysolve::detail
