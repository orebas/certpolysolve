// The piecewise method's proposal engine on one chart's local model, its chart the unit disk: proposals finds candidate
// roots in binary64, and polish refines them by Newton's method in ball arithmetic.
//
// proposals: under interior_proposals, unless full_recovery, it runs on the model restricted to its proposal disk
// (restricted_root_model) and scales the points back; degrees one and two take closed forms in ball arithmetic.
// Otherwise, in binary64: chart_excluded (under chart_exclusion) skips a chart in whose kept disk Pellet's test at
// index 0 after Graeffe root-squaring steps finds no root; seed_count starts as many points as the Newton polygon
// counts roots within radius 2^seed_count_radius_log2, on a golden-angle spiral; a partial Ehrlich-Aberth iteration in
// total-step form moves each until its value is within rounding noise or its step is tiny; and the points inside
// chart_keep_radius are returned unpolished. polish takes up to heuristics::piecewise_polish_steps Newton steps a
// point, until the value is below a quarter of the model's error bound (under aligned_polish, root_model_error about
// the point), the derivative contains zero, the step has under 4 accurate bits or falls below 2^-(goal bits + 10), or
// the point leaves the chart; reuse_polish skips a point the model bound stopped earlier. guard_options sets a
// provisional ownership guard's query: the accuracy goal, absolute tolerance 2r/(n+1) for chart radius r and degree n.
//
// State: none; the statistics and the stable flags belong to the caller.
// Entry points: proposals, polish, guard_options, from piecewise_roots.cpp: PiecewiseRootState::Impl's prepare_model,
// process_pending and certify_guard, and the replays propose_piecewise_candidates and polish_piecewise_candidates.
// Not here: the certificates (certify_root_model, chart_certificate.cpp), the models and their error bounds
// (piecewise_models.cpp), and GuardIndex with the rectangle helpers (inline in piecewise_proposals.hpp).
// Design notes: docs/piecewise-root-certification.md (candidates, ownership and stopping),
// docs/piecewise-performance-findings-20260923.md (the interior proposals).
#include "methods/piecewise/piecewise_proposals.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include "core/dovetail.hpp"
#include "numeric/evaluate/blocked_polynomial_evaluator.hpp"
#include "core/telemetry.hpp"
#include <complex>
#include <numbers>

namespace certpolysolve::detail {
cps_options guard_options(const RootModel &model, const cps_options &options) {
    cps_options result = options;
    result.goal = CPS_ROOT_ACCURACY; // Provisional ownership guards contain one root.
    RealBall tolerance;
    arb_mul_2exp_si(tolerance.value(), model.radius.value(), 1);
    arb_div_si(tolerance.value(), tolerance.value(), acb_poly_degree(model.polynomial.value()) + 1, 64);
    Integer m, e;
    arf_get_fmpz_2exp(m.value(), e.value(), arb_midref(tolerance.value()));
    result.absolute_tolerance = Number(integer_string(m.value()) + "@" + integer_string(e.value()));
    return result;
}
namespace {
// A chart's binary64 proposal iteration keeps the points that end inside the
// disk |u| < chart_keep_radius (the chart is the unit disk; the margin admits
// a root on its edge).
inline constexpr double chart_keep_radius = 1.001;
// Points whose Horner recurrences the proposals' iteration runs side by side (proposals, below). A step of one
// chain is a product and two sums deep, about thirteen cycles, and has twelve operations; four chains keep the
// two multipliers and adders of the cores this runs on busy through those cycles, and two-wide vectors hold them.
inline constexpr size_t evaluation_lanes = 4;
// Arithmetic of the binary64 proposals (proposals, below): magnitudes compared squared, and Smith's division. A
// complex double's std::abs is a hypot call, and its division the C runtime's Annex G __divdc3 (infinities and NaNs
// handled in every operand), together a large part of the iteration's time. Proposals are never certificates: their
// rounding is free, their finiteness checked as before.
// History: docs/code-history-20261021.md, piecewise_proposals.cpp 1.
inline double magnitude_squared(std::complex<double> z) { return z.real() * z.real() + z.imag() * z.imag(); }
inline bool finite(std::complex<double> z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); }
// Smith's algorithm: the larger component of the divisor is never squared, so the quotient neither overflows nor
// underflows where it is representable. The divisor must be nonzero.
inline std::complex<double> quotient(std::complex<double> a, std::complex<double> b) {
    if (std::abs(b.imag()) < std::abs(b.real())) {
        const double r = b.imag() / b.real(), d = b.real() + b.imag() * r;
        return {(a.real() + a.imag() * r) / d, (a.imag() - a.real() * r) / d};
    }
    const double r = b.real() / b.imag(), d = b.imag() + b.real() * r;
    return {(a.real() * r + a.imag()) / d, (a.imag() * r - a.real()) / d};
}
// The points a chart's iteration starts: as many as its polynomial c has roots within 2^seed_count_radius_log2 of
// the centre, by the Newton polygon's count there (the index of the largest |c_k| R^k; the larger index on a tie).
// In the iteration's coordinates the kept disk is the unit disk, and a chart's own disk, the only place where its
// polynomial stands for the input, has radius 1 / proposal_radius, below 2 since the proposal radius exceeds a
// half: the count covers every root the model represents on any chart, and the kept disk's radius twice over.
// This count was measured against twice the dominant index plus sixteen, which starts many more points to keep a
// few more, and against the count in the chart's own disk, with none or one to spare, which loses roots, each a
// whole further round. A size is compared as its binary exponent plus k seed_count_radius_log2, then its
// fraction: exact, and no power can overflow. History: docs/code-history-20261021.md, piecewise_proposals.cpp 2.
inline constexpr int seed_count_radius_log2 = 1;
size_t seed_count(const std::vector<std::complex<double>> &c) {
    size_t counted = 0;
    slong best_exponent = std::numeric_limits<slong>::min();
    double best_fraction = 0;
    for (size_t k = 0; k < c.size(); ++k) {
        const double size = std::abs(c[k]);
        if (!(size > 0) || !std::isfinite(size)) continue;
        int exponent = 0;
        const double fraction = std::frexp(size, &exponent);
        const slong scaled = static_cast<slong>(exponent) + static_cast<slong>(k) * seed_count_radius_log2;
        if (scaled > best_exponent || (scaled == best_exponent && fraction >= best_fraction)) {
            best_exponent = scaled;
            best_fraction = fraction;
            counted = k;
        }
    }
    return counted;
}
// Whether the chart polynomial c (binary64 midpoints) has no root in the disk
// |u| <= chart_keep_radius, by Pellet's test at index 0 after root-squaring
// (Graeffe) steps: each step squares every root's modulus, so roots outside
// the disk move away doubly exponentially and the constant term comes to
// dominate. The test keeps a factor two for the rounding of the steps. A
// root a sixteenth of the radius outside takes about log2(16 ln(2n)) steps to
// be pushed out of reach, which bounds the steps; nearer ones are left to the
// iteration. Over half the charts of a large dense input hold no root, and
// the iteration spends much of its time on them for no point; these steps
// exclude most of them, a heuristic for proposals only (the certificates
// decide). Each step's multiply-adds are charged to the ledger as the
// iteration's terms are.
// History: docs/code-history-20261021.md, piecewise_proposals.cpp 3.
inline constexpr double chart_exclusion_reach = 16;
bool chart_excluded(std::vector<std::complex<double>> c, PiecewiseRootStatistics &observed) {
    const size_t n = c.size() - 1;
    if (n < 1)
        return false;
    double power = 1;
    for (auto &x : c) {
        x *= power;
        power *= chart_keep_radius;
    }
    const auto dominated = [&] {
        double others = 0;
        for (size_t j = 1; j <= n; ++j) others += std::abs(c[j]);
        return std::abs(c[0]) > 2 * others;
    };
    if (dominated())
        return true;
    const int steps = static_cast<int>(
        std::ceil(std::log2(chart_exclusion_reach * std::log(2 * static_cast<double>(n + 1)))));
    std::vector<std::complex<double>> squared(n + 1);
    for (int step = 0; step < steps; ++step) {
        // The coefficients of p(u) p(-u) in u^2, up to sign: c_k^2 + 2 sum_l
        // (-1)^l c_(k-l) c_(k+l).
        for (size_t k = 0; k <= n; ++k) {
            std::complex<double> sum = c[k] * c[k];
            const size_t reach = std::min(k, n - k);
            for (size_t l = 1; l <= reach; ++l)
                sum += (l % 2 ? -2.0 : 2.0) * c[k - l] * c[k + l];
            squared[k] = sum;
            observed.work.machine += static_cast<double>(reach + 1);
        }
        double top = 0;
        for (const auto &x : squared) top = std::max(top, std::abs(x));
        if (!(top > 0) || !std::isfinite(top))
            return false;
        for (size_t k = 0; k <= n; ++k) c[k] = squared[k] / top;
        if (dominated())
            return true;
    }
    return false;
}
// A chart's proposals: the points its binary64 iteration keeps (or the roots of a model of degree one or two),
// in chart coordinates and not polished. The caller polishes the ones it goes on with: a root lies in the kept
// disks of several charts, and every chart after the first proposes a copy that process_candidate retires;
// polished first, the copies would be most of the polished points.
// History: docs/code-history-20261021.md, piecewise_proposals.cpp 4.
} // namespace
slong proposals(RootVector &points, const RootModel &model, const cps_options &options,
                certpolysolve::detail::CandidateStatistics &stats, slong precision, PiecewiseRootStatistics &observed,
                const PiecewiseRootConfig &config, bool full_recovery) {
    if (config.interior_proposals && !full_recovery && arb_is_positive(model.proposal_radius.value())) {
        RootModel inner;
        const RealBall &fraction = model.proposal_radius;
        {
            StageTimer timer(observed.core_proposal_setup_ns);
            inner = restricted_root_model(model, fraction, precision);
            observed.work.restriction += 2 * static_cast<double>(acb_poly_length(model.polynomial.value())) *
                                         operation_weight(precision);
        }
        ++observed.core_proposal_models;
        // The rescaled polynomial is only a proposal representation: the
        // caller polishes and certifies in the full original chart. Bypass
        // further restriction in this recursive call.
        const slong count = proposals(points, inner, options, stats, precision, observed, config, true);
        {
            StageTimer timer(observed.core_proposal_setup_ns);
            for (slong j = 0; j < count; ++j)
                acb_mul_arb(points[j], points[j], fraction.value(), precision);
        }
        return count;
    }
    ExclusiveTimer machine_timer(observed.machine_ns, &observed.polish_ns);
    if (points.size() == 1) {
        acb_div(points[0], model.polynomial.value()->coeffs, model.polynomial.value()->coeffs + 1, precision);
        acb_neg(points[0], points[0]);
        acb_get_mid(points[0], points[0]);
        return 1;
    }
    if (points.size() == 2) {
        // Machine coefficients can erase a quadratic's tiny discriminant.
        // These are adaptive-precision proposals only; the original-model
        // rectangle certificate still decides whether either root is usable.
        const auto p = model.polynomial.value();
        ComplexBall discriminant, term, denominator;
        acb_sqr(discriminant.value(), p->coeffs + 1, precision);
        acb_mul(term.value(), p->coeffs, p->coeffs + 2, precision);
        acb_mul_2exp_si(term.value(), term.value(), 2);
        acb_sub(discriminant.value(), discriminant.value(), term.value(), precision);
        acb_sqrt(discriminant.value(), discriminant.value(), precision);
        acb_mul_2exp_si(denominator.value(), p->coeffs + 2, 1);
        acb_sub(points[0], discriminant.value(), p->coeffs + 1, precision);
        acb_neg(discriminant.value(), discriminant.value());
        acb_sub(points[1], discriminant.value(), p->coeffs + 1, precision);
        for (slong j = 0; j < 2; ++j) {
            acb_div(points[j], points[j], denominator.value(), precision);
            acb_get_mid(points[j], points[j]);
        }
        return 2;
    }
    // Partial Ehrlich--Aberth in the chart. A tropical coefficient index
    // estimates how many roots are relevant here. Temporary roundoff stopping
    // transfers proposals to ball polishing, never directly to a certificate.
    using Complex = std::complex<double>;
    std::vector<Complex> coefficients(static_cast<size_t>(points.size() + 1));
    double maximum = 0;
    for (slong j = 0; j <= points.size(); ++j) {
        const auto c = model.polynomial.value()->coeffs + j;
        coefficients[static_cast<size_t>(j)] = {arf_get_d(arb_midref(acb_realref(c)), ARF_RND_NEAR),
                                                arf_get_d(arb_midref(acb_imagref(c)), ARF_RND_NEAR)};
        maximum = std::max(maximum, std::abs(coefficients[static_cast<size_t>(j)]));
    }
    if (!(maximum > 0) || !std::isfinite(maximum))
        return 0;
    // Candidate-only truncation: the full ball polynomial still certifies and
    // polishes every point. No tail is omitted from the proof.
    while (coefficients.size() > 2 && std::abs(coefficients.back()) < maximum * 1e-19)
        coefficients.pop_back();
    if (config.chart_exclusion && chart_excluded(coefficients, observed)) {
        ++observed.excluded_charts;
        return 0;
    }
    std::vector<double> magnitudes;
    magnitudes.reserve(coefficients.size());
    for (const auto c : coefficients) magnitudes.push_back(std::abs(c));
    const size_t seeds = std::min(coefficients.size() - 1, std::max<size_t>(1, seed_count(coefficients)));
    observed.machine_seeds += seeds;
    std::vector<Complex> current, next, values, derivatives;
    std::vector<double> scales;
    std::vector<size_t> moving;
    std::vector<bool> settled;
    const double leading_magnitude = std::abs(coefficients.back());
    auto refine = [&] {
        next.resize(current.size());
        settled.assign(current.size(), false);
        for (unsigned iteration = 0; iteration < 96 && !options.cancelled(); ++iteration) {
            ++observed.machine_sweeps;
            bool active = false;
            // Value, slope and term scale at every point still moving, evaluation_lanes points at a time. A
            // point's Horner recurrence is one chain of dependent operations (a product and two sums deep a
            // term), which the processor cannot overlap with itself; chains of different points it can. Each
            // lane does the scalar recurrence's operations in its order, so every result is what it was. A
            // sweep reads the points of the sweep before and writes the next ones, so their order is free.
            moving.clear();
            for (size_t i = 0; i < current.size(); ++i) {
                next[i] = current[i];
                if (!settled[i]) moving.push_back(i);
            }
            active = !moving.empty();
            values.resize(moving.size());
            derivatives.resize(moving.size());
            scales.resize(moving.size());
            for (size_t first = 0; first < moving.size(); first += evaluation_lanes) {
                const size_t width = std::min(evaluation_lanes, moving.size() - first);
                double zr[evaluation_lanes], zi[evaluation_lanes], radius[evaluation_lanes], vr[evaluation_lanes],
                    vi[evaluation_lanes], dr[evaluation_lanes], di[evaluation_lanes], scale[evaluation_lanes];
                for (size_t lane = 0; lane < evaluation_lanes; ++lane) {
                    // A short last group repeats its last point in the spare lanes.
                    const Complex z = current[moving[first + std::min(lane, width - 1)]];
                    zr[lane] = z.real();
                    zi[lane] = z.imag();
                    radius[lane] = std::sqrt(magnitude_squared(z));
                    vr[lane] = coefficients.back().real();
                    vi[lane] = coefficients.back().imag();
                    dr[lane] = di[lane] = 0;
                    scale[lane] = leading_magnitude;
                }
                for (size_t j = coefficients.size() - 1; j > 0; --j) {
                    const double cr = coefficients[j - 1].real(), ci = coefficients[j - 1].imag();
                    const double size = magnitudes[j - 1];
                    for (size_t lane = 0; lane < evaluation_lanes; ++lane) {
                        // derivative = derivative z + value; value = value z + c; scale = scale |z| + |c|.
                        const double slope_real = dr[lane] * zr[lane] - di[lane] * zi[lane] + vr[lane];
                        const double slope_imaginary = dr[lane] * zi[lane] + di[lane] * zr[lane] + vi[lane];
                        const double value_real = vr[lane] * zr[lane] - vi[lane] * zi[lane] + cr;
                        const double value_imaginary = vr[lane] * zi[lane] + vi[lane] * zr[lane] + ci;
                        dr[lane] = slope_real;
                        di[lane] = slope_imaginary;
                        vr[lane] = value_real;
                        vi[lane] = value_imaginary;
                        scale[lane] = scale[lane] * radius[lane] + size;
                    }
                }
                for (size_t lane = 0; lane < width; ++lane) {
                    values[first + lane] = {vr[lane], vi[lane]};
                    derivatives[first + lane] = {dr[lane], di[lane]};
                    scales[first + lane] = scale[lane];
                }
            }
            for (size_t point = 0; point < moving.size(); ++point) {
                const size_t i = moving[point];
                const Complex z = current[i];
                const Complex value = values[point], derivative = derivatives[point];
                const double term_scale = scales[point];
                const double radius = std::sqrt(magnitude_squared(z));
                ++stats.evaluations;
                stats.evaluation_terms += coefficients.size();
                observed.work.machine += static_cast<double>(coefficients.size());
                if (derivative == Complex{})
                    continue;
                Complex q = quotient(value, derivative);
                if (!finite(q))
                    continue;
                const double noise = term_scale * (32 * static_cast<double>(coefficients.size()) + 16) *
                                     std::numeric_limits<double>::epsilon();
                const double small_step = 2e-15 * (1 + radius);
                if (magnitude_squared(value) < noise * noise || magnitude_squared(q) < small_step * small_step) {
                    settled[i] = true;
                    continue;
                }
                // The sum of 1 / (z - w) over the other points, each as the conjugate over the squared modulus:
                // one division a pair and no branch. A coincident point (squared modulus zero) is skipped, as
                // before; chart coordinates keep the squares in range.
                double interaction_real = 0, interaction_imaginary = 0;
                for (size_t j = 0; j < current.size(); ++j) {
                    const double x = z.real() - current[j].real(), y = z.imag() - current[j].imag();
                    const double size = x * x + y * y;
                    const double inverse = size > 0 ? 1 / size : 0;
                    interaction_real += x * inverse;
                    interaction_imaginary -= y * inverse;
                }
                const Complex interaction(interaction_real, interaction_imaginary);
                observed.work.interactions += static_cast<double>(current.size() - 1);
                const Complex denominator = 1.0 - q * interaction;
                if (denominator != Complex{})
                    q = quotient(q, denominator);
                if (const double size = magnitude_squared(q); size > .25)
                    q *= .5 / std::sqrt(size);
                if (finite(z - q))
                    next[i] = z - q;
                if (magnitude_squared(next[i]) > 9)
                    settled[i] = true;
            }
            current.swap(next);
            if (!active)
                break;
        }
        return std::all_of(settled.begin(), settled.end(), [](bool value) { return value; }) &&
            std::all_of(current.begin(), current.end(),
                        [](Complex z) { return magnitude_squared(z) < chart_keep_radius * chart_keep_radius; });
    };
    current.resize(seeds);
    for (size_t seed = 0; seed < seeds; ++seed) {
        const double radius = 1.1 * std::sqrt(static_cast<double>(seed) / static_cast<double>(seeds));
        current[seed] = std::polar(radius, static_cast<double>(seed) * std::numbers::pi * (3 - std::sqrt(5.0)));
    }
    refine();
    slong count = 0;
    for (size_t i = 0; i < current.size(); ++i) {
        const auto z = current[i];
        if (magnitude_squared(z) < chart_keep_radius * chart_keep_radius && count < points.size()) {
            if (!settled[i]) {
                ++observed.unsettled_points;
            }
            acb_set_d_d(points[count], z.real(), z.imag());
            ++count;
        }
    }
    observed.machine_points += static_cast<uint64_t>(count);
    return count;
}
void polish(RootVector &points, slong count, const RootModel &model, const cps_options &options,
            slong precision, PiecewiseRootStatistics &stats, std::vector<bool> *stable, const PiecewiseRootConfig &config,
            const std::vector<size_t> *owners, std::span<const uint64_t> owner_counts) {
    if (!owner_counts.empty() && owner_counts.size() != static_cast<size_t>(count))
        throw std::logic_error("polishing owner diagnostic size mismatch");
    if (config.observe_polish_batch && count>0) {
        StageTimer timer(stats.diagnostic_ns);
        config.observe_polish_batch(model,points,count,options,precision);
    }
    stats.polish_points += static_cast<uint64_t>(count);
    // Polish only chart-local proposals. This is local polynomial Newton work;
    // every accepted result is subsequently checked with the original-model
    // remainder. We neither solve nor evaluate the whole input per proposal.
    ComplexBall value, derivative, step;
    RealBall distance, one, value_bound, model_floor, global_floor, correction, tolerance;
    arb_one(one.value());
    arb_mul_2exp_si(model_floor.value(), model.error.value(), -2);
    global_floor = model_floor;
    arb_one(tolerance.value());
    arb_mul_2exp_si(tolerance.value(), tolerance.value(), -resources::precision_sum(options.bits, 10));
    std::unique_ptr<BlockedPolynomialEvaluator> evaluator;
    for (slong i = 0; i < count; ++i) {
        const size_t position = static_cast<size_t>(i);
        struct Observation {
            const PiecewiseRootConfig &config;
            const PiecewiseRootStatistics &stats;
            PiecewisePolishEvent event;
            uint64_t before;
            std::array<uint64_t, 4> stops_before;
            std::chrono::steady_clock::time_point start;
            ~Observation() {
                if (!config.observe_polish) return;
                event.evaluations = stats.polish_evaluations - before;
                event.nanoseconds = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - start).count());
                event.source = config.diagnostic_source_id;
                event.model_stops = stats.polish_model_stops - stops_before[0];
                event.rounding_stops = stats.polish_rounding_stops - stops_before[1];
                event.step_stops = stats.polish_step_stops - stops_before[2];
                event.reused = stats.polish_reused - stops_before[3];
                config.observe_polish(event);
            }
        } observation{config, stats,
            {owners ? static_cast<uint64_t>((*owners)[position]) : UINT64_MAX,
             stats.models, static_cast<uint64_t>(acb_poly_degree(model.polynomial.value())),
             stats.model_bits, static_cast<uint64_t>(precision), 0, 0,
             owner_counts.empty() ? 0 : owner_counts[position]}, stats.polish_evaluations,
             {stats.polish_model_stops, stats.polish_rounding_stops, stats.polish_step_stops, stats.polish_reused},
             config.observe_polish ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{}};
        if (config.reuse_polish && stable && (*stable)[position]) { ++stats.polish_reused; continue; }
        acb_abs(distance.value(), points[i], 64);
        if (!arb_is_finite(distance.value()) || !arb_lt(distance.value(), one.value())) {
            ++observation.event.domain_stops;
            continue;
        }
        RealBall previous_correction;
        bool correction_seen=false;
        uint64_t iteration=0;
        const uint64_t allowance = heuristics::piecewise_polish_steps;
        for (; iteration < allowance &&
                                     !options.cancelled(); ++iteration) {
            upper(distance.value(), points[i]);
            if (!arb_lt(distance.value(), one.value())) {
                ++observation.event.domain_stops;
                break;
            }
            ++stats.polish_evaluations;
            // A value and a derivative: two Horner terms per coefficient.
            stats.work.polish += 2 * static_cast<double>(acb_poly_length(model.polynomial.value())) *
                                 operation_weight(precision);
            {
            OptionalStageTimer detail(stats.polish_evaluation_ns, config.detailed_accounting);
            if (config.blocked_polish) {
                if (!evaluator) evaluator = std::make_unique<BlockedPolynomialEvaluator>(model.polynomial, precision);
                evaluator->evaluate(value.value(), derivative.value(), points[i]);
            } else acb_poly_evaluate2(value.value(), derivative.value(), model.polynomial.value(), points[i], precision);
            }
            upper(value_bound.value(), value.value());
            const bool skip_bound = config.reuse_polish_bound &&
                !arb_lt(value_bound.value(), global_floor.value());
            if (skip_bound) {
                // root_model_error returns the global bound or a provably
                // smaller local bound. If value < global cannot be established,
                // value < local cannot hold either. Refresh before any stop.
                // Thus this reuses the existing whole-chart bound without
                // changing a Newton step or the exact stopping decision.
                model_floor = global_floor;
                ++stats.polish_bounds_skipped;
            } else if (config.aligned_polish && model.error_envelope) {
                OptionalStageTimer detail(stats.polish_error_bound_ns, config.detailed_accounting);
                // Proposal stopping is a scheduling decision. Evaluate the
                // shared remainder on an outward candidate rectangle; final
                // certification recomputes it on the published rectangle.
                ComplexBall rectangle;
                acb_mul_arb(rectangle.value(), points[i], model.radius.value(), precision);
                acb_add(rectangle.value(), rectangle.value(), model.center.value(), precision);
                RealBall radius;
                upper(radius.value(), rectangle.value());
                arb_mul_2exp_si(radius.value(), radius.value(), -resources::precision_sum(options.bits, 5));
                if (!options.absolute_tolerance.zero()) {
                    RealBall absolute;
                    options.absolute_tolerance.enclose(absolute.value(), precision);
                    arb_mul_2exp_si(absolute.value(), absolute.value(), -3);
                    arb_max(radius.value(), radius.value(), absolute.value(), precision);
                }
                arb_add_error(acb_realref(rectangle.value()), radius.value());
                arb_add_error(acb_imagref(rectangle.value()), radius.value());
                const auto error = root_model_error(model, rectangle.value(), precision);
                arb_mul_2exp_si(model_floor.value(), error.value(), -2);
                ++stats.polish_local_bounds;
                stats.work.polish_bounds += operation_weight(precision);
            }
            if (arb_lt(value_bound.value(), model_floor.value())) {
                ++stats.polish_model_stops;
                // Repeating this same model/point/precision performs no update.
                // Reuse only this exact terminal observation, not small-step or
                // arithmetic-floor heuristics that may still change a point.
                if (stable) (*stable)[position] = true;
                break;
            }
            if (acb_contains_zero(derivative.value())) {
                ++observation.event.derivative_stops;
                break;
            }
            acb_div(step.value(), value.value(), derivative.value(), precision);
            if (!acb_is_finite(step.value()))
                break;
            // Keep the last finite Newton proposal even if it leaves this
            // chart or its correction has few reliable bits. New proposals
            // outside the chart are rejected below; owned seeds can move to a
            // neighboring chart. Clamping here preserves stale poor seeds.
            // This is candidate work, never a root-count or error guarantee.
            acb_sub(points[i], points[i], step.value(), precision);
            acb_get_mid(points[i], points[i]);
            upper(correction.value(), step.value());
            if (config.observe_polish) {
                if (correction_seen && !arb_lt(correction.value(),previous_correction.value()))
                    ++observation.event.nonimproving_steps;
                previous_correction=correction;
                correction_seen=true;
            }
            if (acb_rel_accuracy_bits(step.value()) < 4) { ++stats.polish_rounding_stops; break; }
            if (acb_is_zero(step.value()) || arb_lt(correction.value(), tolerance.value())) {
                ++stats.polish_step_stops; break;
            }
        }
        if(iteration==allowance && !options.cancelled())++observation.event.iteration_limit;
    }
}
} // namespace certpolysolve::detail
