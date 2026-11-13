// Proposal arithmetic shared by the methods: starting points for the root iterations and the multiprecision packets of
// sweeps that move them (candidate_packet, secular_candidate_packet). Every point leaves as an exact midpoint ball.
//
// initial_candidates (initialize) reads the Newton polygon, the upper convex hull of the points (j, h_j) with h_j the
// integer bit exponent of c_j: an edge from a to b puts b - a points on the circle of radius 2^((h_a - h_b)/(b - a)),
// at placement_angle's angles or, with binomial_seeds, at the arguments of the roots of c_a + c_b x^(b-a) turned by
// heuristics::binomial_seed_spacing_fraction of their spacing. scale_candidates enters a coefficient with no correct
// bit as its magnitude bound and reports whether every vertex is then resolved (determined).
//
// candidate_packet scales the variable and coefficients by powers of two and makes Ehrlich-Aberth sweeps in
// Gauss-Seidel order, with Newton's correction N = p/p' from a dense blocked Horner (MidpointPolynomialEvaluator) or
// from powers over a sparse support's gaps (evaluate), or, under CandidateEngine::Secular, secular_newton's N on
// residues regenerate_secular forms at the starts; classical_error_stopping stops a point at Horner's rounding floor
// (classical_evaluation_floor), and the research build's CandidateEngine::DurandKerner takes Weierstrass steps.
// secular_candidate_packet sweeps each point at its own precision, in Jacobi or Gauss-Seidel order, with N from a
// MonicSource or secular_newton on the model P = prod(x - b_j)(sum_j a_j/(x - b_j) - 1); under
// SecularModelPolicy::Continuation a point settles where the model's value is within its rounding error or the step is
// below the resolution (model_stop; NewtonPole where N is infinite). Either packet may take the interaction sums' far
// field from an InteractionTree (prepare_interactions); both add the near terms and take Aberth's correction
// N / (1 - N S) by aberth_step.hpp's add_interactions and aberth_delta.
//
// State: none at file scope. A caller's CandidateWorkspace::Impl keeps a packet's preparation (scaled coefficients,
// sparse support, dense evaluator, regenerated secular model) until the polynomial, precision or engine changes, and
// its stop flags for a packet resuming where it ended. No caller passes a SecularCandidateWorkspace.
// Entry points, by caller: secular_refinement.cpp calls scale_candidates, circle_candidates, candidate_packet
// (component_seeds, SecularRefinement::preliminary) and secular_candidate_packet (SecularRefinement::step); solver.cpp,
// secular_candidate_packet (program_advance) and initial_candidates; certificate_attempt.cpp (machine_newton_cover) and
// study_candidates.cpp, initial_candidates and candidate_packet; classical_isolation.cpp, candidate_packet and, for
// the research build's CandidateEngine::DoubleSeed, binary64_seeds (binary64 Aberth sweeps). Outside this file, only
// tests call regenerate_secular and midpoint_polynomial_evaluation.
// Not here: certification (the callers), binary64 packets (machine_candidates.cpp, study_candidates.cpp), stopping
// rules (model_stopping.hpp), the tree (interaction_tree.cpp), regenerate_residues (secular_refinement.cpp).
// Design notes: docs/mathematics.md ("Candidate experiments and certification", "Separate C1 experiments"),
// docs/solver-map.md ("Shared mathematical core"), docs/campaign-d5-model-stopping.md (the per-root model stop).
#include "numeric/candidates/candidates.hpp"
#include "core/query.hpp"
#include "platform/parallel.hpp"
#include <atomic>
#include "numeric/candidates/aberth_step.hpp"
#include "core/service_meter.hpp"
#include "core/heuristics.hpp"
#include "core/telemetry.hpp"
#include "numeric/evaluate/structured_source.hpp"
#include "numeric/candidates/interaction_tree.hpp"
#include "numeric/evaluate/midpoint_polynomial_evaluator.hpp"
#include "core/wide_int.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <flint/acf.h>
#include <limits>
#include <numbers>
#include <numeric>

namespace certpolysolve::detail {
namespace {
constexpr auto rounding = ARF_RND_NEAR;
struct Point {
    acf_t value;
    Point() { acf_init(value); }
    ~Point() { acf_clear(value); }
    Point(const Point &x) : Point() { acf_set(value, x.value); }
    Point(Point &&x) noexcept : Point() { acf_swap(value, x.value); }
    Point &operator=(Point x) noexcept {
        acf_swap(value, x.value);
        return *this;
    }
};
void zero(acf_t x) {
    arf_zero(acf_realref(x));
    arf_zero(acf_imagref(x));
}
void one(acf_t x) {
    arf_one(acf_realref(x));
    arf_zero(acf_imagref(x));
}
bool finite(const acf_t x) { return arf_is_finite(acf_realref(x)) && arf_is_finite(acf_imagref(x)); }
// Not "iszero": glibc's <math.h> declares a global template iszero(T) (T == 0),
// an exact match for a non-const acf_t argument, which would win every such
// call and compare the pointer with zero (always false): the zero tests on
// differences and denominators below would fall through to divisions by zero,
// caught only by the finite() checks after them.
// History: docs/code-history-20261021.md, candidates.cpp 1.
bool is_zero(const acf_t x) { return arf_is_zero(acf_realref(x)) && arf_is_zero(acf_imagref(x)); }
void midpoint(acf_t x, const acb_t y) {
    arf_set(acf_realref(x), arb_midref(acb_realref(y)));
    arf_set(acf_imagref(x), arb_midref(acb_imagref(y)));
}
void point_ball(acb_t x, const acf_t y) {
    arb_set_arf(acb_realref(x), acf_realref(y));
    arb_set_arf(acb_imagref(x), acf_imagref(y));
}
double magnitude_log2(const acf_t x) {
    if (is_zero(x))
        return -std::numeric_limits<double>::infinity();
    return static_cast<double>(
        std::max(arf_abs_bound_lt_2exp_si(acf_realref(x)), arf_abs_bound_lt_2exp_si(acf_imagref(x))));
}
void magnitude(mag_t out, const acf_t point) {
    Magnitude real, imag;
    arf_get_mag(real.value(), acf_realref(point));
    arf_get_mag(imag.value(), acf_imagref(point));
    mag_hypot(out, real.value(), imag.value());
}
// Cheap absolute-Horner cancellation estimate, used only to schedule candidate
// work. This is not an enclosure of a root or an acceptance predicate.
bool classical_evaluation_floor(const acf_t value, const acf_t point,
                                const std::vector<Magnitude> &coefficients, slong precision) {
    Magnitude radius, envelope, error, residual;
    magnitude(radius.value(), point);
    mag_set(envelope.value(), coefficients.back().value());
    for (size_t j = coefficients.size() - 1; j-- > 0;) {
        mag_mul(envelope.value(), envelope.value(), radius.value());
        mag_add(envelope.value(), envelope.value(), coefficients[j].value());
    }
    mag_mul_ui(error.value(), envelope.value(), static_cast<ulong>(coefficients.size()));
    mag_mul_2exp_si(error.value(), error.value(), 4 - precision);
    magnitude(residual.value(), value);
    return mag_cmp(residual.value(), error.value()) <= 0;
}
std::unique_ptr<InteractionTree> prepare_interactions(const std::vector<Point> &points, slong &scale) {
    double maximum = -std::numeric_limits<double>::infinity();
    for (const auto &p : points)
        maximum = std::max(maximum, magnitude_log2(p.value));
    // Range guard of an optional helper: the tree works on binary64 copies
    // scaled by 2^-scale. Past this magnitude no tree is built and the caller
    // sums the interactions directly. The literal only has to keep `scale` in
    // a slong; it is not a property of the input.
    if (!std::isfinite(maximum) || std::abs(maximum) >= 1.e12)
        return {};
    scale = static_cast<slong>(maximum);
    std::vector<std::complex<double>> scaled(points.size()), lows(points.size());
    // Each coordinate as two doubles, the second what the first leaves of it,
    // and what the two leave, bounded above (rounded away from zero): the
    // tree places its points to about 2^-106 of their size, and takes no far
    // field the copies cannot place (interaction_tree.cpp). With one double,
    // points 2^-52 of their size apart and closer would be neighbours, and a
    // cluster would pay exact interactions for them.
    // History: docs/code-history-20261021.md, candidates.cpp 2.
    std::vector<double> errors(points.size(), 0);
    Float coordinate, copy, residual, remainder;
    const auto rounded = [&](const arf_t part, double &low, double &error) {
        arf_mul_2exp_si(coordinate.value(), part, -scale);
        const double value = arf_get_d(coordinate.value(), ARF_RND_NEAR);
        arf_set_d(copy.value(), value);
        arf_sub(residual.value(), coordinate.value(), copy.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
        low = arf_get_d(residual.value(), ARF_RND_NEAR);
        arf_set_d(copy.value(), low);
        arf_sub(remainder.value(), residual.value(), copy.value(), std::numeric_limits<double>::digits, ARF_RND_UP);
        error += std::abs(arf_get_d(remainder.value(), ARF_RND_UP));
        return value;
    };
    for (size_t j = 0; j < points.size(); ++j) {
        double low_x = 0, low_y = 0;
        const double x = rounded(acf_realref(points[j].value), low_x, errors[j]);
        if (!std::isfinite(x) || !std::isfinite(low_x))
            return {};
        const double y = rounded(acf_imagref(points[j].value), low_y, errors[j]);
        if (!std::isfinite(y) || !std::isfinite(low_y))
            return {};
        scaled[j] = {x, y};
        lows[j] = {low_x, low_y};
    }
    // Underflowed coordinates coalesce in a nearby bucket. The tree retains
    // every index, and nearby interactions use the original arbitrary-
    // precision points. This geometry is candidate-only, never a certificate.
    return std::make_unique<InteractionTree>(scaled, errors, lows);
}
void power(acf_t out, const acf_t x, ulong n, slong prec) {
    Point base;
    acf_set(base.value, x);
    one(out);
    while (n) {
        if (n & 1)
            acf_mul(out, out, base.value, prec, rounding);
        n >>= 1;
        if (n)
            acf_mul(base.value, base.value, base.value, prec, rounding);
    }
}
void evaluate(acf_t value, acf_t derivative, const std::vector<Point> &coefficients,
              const std::vector<slong> &indices, const acf_t z, slong prec) {
    acf_set(value, coefficients.back().value);
    zero(derivative);
    if (indices.size() * 4 >= coefficients.size()) {
        for (size_t j = coefficients.size() - 1; j-- > 0;) {
            acf_mul(derivative, derivative, z, prec, rounding);
            acf_add(derivative, derivative, value, prec, rounding);
            acf_mul(value, value, z, prec, rounding);
            acf_add(value, value, coefficients[j].value, prec, rounding);
        }
    } else {
        Point zp, dzp, t;
        for (size_t j = 1; j < indices.size(); ++j) {
            const auto gap = static_cast<ulong>(indices[j - 1] - indices[j]);
            power(dzp.value, z, gap - 1, prec);
            acf_mul(zp.value, dzp.value, z, prec, rounding);
            arf_mul_ui(acf_realref(dzp.value), acf_realref(dzp.value), gap, prec, rounding);
            arf_mul_ui(acf_imagref(dzp.value), acf_imagref(dzp.value), gap, prec, rounding);
            acf_mul(derivative, derivative, zp.value, prec, rounding);
            acf_mul(t.value, value, dzp.value, prec, rounding);
            acf_add(derivative, derivative, t.value, prec, rounding);
            acf_mul(value, value, zp.value, prec, rounding);
            acf_add(value, value, coefficients[static_cast<size_t>(indices[j])].value, prec, rounding);
        }
        // The gaps end at the lowest nonzero coefficient, index m: value and
        // derivative are those of p(z)/z^m. Put the factor back, or the
        // Newton correction is that of a polynomial without p's zero root
        // (an inexact factor keeps a zero constant midpoint; exact input
        // strips its zero roots before).
        const auto m = static_cast<ulong>(indices.back());
        if (m) {
            power(t.value, z, m - 1, prec);
            acf_mul(dzp.value, derivative, z, prec, rounding);
            acf_set(zp.value, value);
            arf_mul_ui(acf_realref(zp.value), acf_realref(zp.value), m, prec, rounding);
            arf_mul_ui(acf_imagref(zp.value), acf_imagref(zp.value), m, prec, rounding);
            acf_add(dzp.value, dzp.value, zp.value, prec, rounding);
            acf_mul(derivative, dzp.value, t.value, prec, rounding);
            acf_mul(value, value, z, prec, rounding);
            acf_mul(value, value, t.value, prec, rounding);
        }
    }
}
// arg z of a midpoint, NaN for zero or non-finite. Both parts are scaled by a common power of two first, so a
// coefficient whose exponent lies outside binary64's range still has its argument.
double midpoint_argument(const acf_t z) {
    const arf_struct *re = acf_realref(z), *im = acf_imagref(z);
    if (!arf_is_finite(re) || !arf_is_finite(im) || (arf_is_zero(re) && arf_is_zero(im)))
        return std::numeric_limits<double>::quiet_NaN();
    const slong scale = arf_is_zero(re) ? arf_abs_bound_lt_2exp_si(im)
                        : arf_is_zero(im) ? arf_abs_bound_lt_2exp_si(re)
                                          : std::max(arf_abs_bound_lt_2exp_si(re), arf_abs_bound_lt_2exp_si(im));
    Float x, y;
    arf_mul_2exp_si(x.value(), re, -scale);
    arf_mul_2exp_si(y.value(), im, -scale);
    return std::atan2(arf_get_d(y.value(), ARF_RND_NEAR), arf_get_d(x.value(), ARF_RND_NEAR));
}
// vertices, when given, receives the Newton polygon's vertex indices. With binomial_seeds (a sparse polynomial,
// heuristics::binomial_seed_spacing_fraction), each edge's seeds lie at the roots of the edge's binomial: the
// polynomial is near that binomial on the edge's annulus, so each seed starts near a root, where the equally spaced
// start may put it half a spacing away. History: docs/code-history-20261021.md, candidates.cpp 3.
void initialize(std::vector<Point> &points, const std::vector<Point> &coefficients, slong prec,
                std::vector<slong> *vertices = nullptr, bool binomial_seeds = false) {
    std::vector<slong> hull;
    std::vector<double> heights;
    // Heights are the coefficients' integer bit exponents, carried as doubles:
    // exact up to 2^digits. Beyond that they are not the exponents any more,
    // and no seed scale can be derived from them.
    const double exact_height = std::ldexp(1.0, std::numeric_limits<double>::digits);
    bool exact_heights = true;
    for (slong j = 0; j < static_cast<slong>(coefficients.size()); ++j) {
        const double h = magnitude_log2(coefficients[static_cast<size_t>(j)].value);
        heights.push_back(h);
        if (!std::isfinite(h))
            continue;
        exact_heights = exact_heights && std::abs(h) <= exact_height;
        while (hull.size() >= 2) {
            const auto a = hull[hull.size() - 2], b = hull.back();
            if ((heights[static_cast<size_t>(b)] - heights[static_cast<size_t>(a)]) *
                    static_cast<double>(j - b) >
                (h - heights[static_cast<size_t>(b)]) * static_cast<double>(b - a))
                break;
            hull.pop_back();
        }
        hull.push_back(j);
    }
    if (vertices)
        *vertices = hull;
    // A coordinate shift can introduce zero roots even after the original
    // input's zero factors were removed. Account for that hull endpoint before
    // placing the remaining candidates on their coefficient-scale circles.
    size_t next = hull.empty() ? 0 : static_cast<size_t>(hull.front());
    for (size_t j = 0; j < next; ++j)
        zero(points[j].value);
    Float radius, direction;
    for (size_t k = 1; exact_heights && k < hull.size(); ++k) {
        const auto a = hull[k - 1], b = hull[k], count = b - a;
        // Radius 2^(drop/count), split exactly: heights are integer bit
        // exponents, so drop = quotient*count + remainder in integers. An exact
        // rescaling x -> 2^s x changes only the quotient, and the seeds (hence
        // every later floating-point decision) stay exact rescalings, where
        // rounding drop/count in floating point first would make q(x) and
        // q(8x) diverge.
        // History: docs/code-history-20261021.md, candidates.cpp 4.
        const auto drop = static_cast<int64_t>(heights[static_cast<size_t>(a)]) -
                          static_cast<int64_t>(heights[static_cast<size_t>(b)]);
        int64_t quotient = drop / count, remainder = drop % count;
        if (remainder < 0) {
            remainder += count;
            --quotient;
        }
        // The quotient is the radius's exponent as it stands. A clamp would
        // put the seeds of an input with very large roots, such as
        // x^3 + 2^(2^50+6) x - 1, far inside them, and the iteration would
        // then creep outward for ever.
        // History: docs/code-history-20261021.md, candidates.cpp 5.
        arf_set_d(radius.value(), std::exp2(static_cast<double>(remainder) / static_cast<double>(count)));
        arf_mul_2exp_si(radius.value(), radius.value(), static_cast<slong>(quotient));
        // The binomial c_a + c_b x^count has its roots at arg(-c_a / c_b) / count + 2 pi j / count. Where that
        // phase is not a number (a zero or non-finite midpoint), the equally spaced start.
        const double phase = binomial_seeds ? std::numbers::pi + midpoint_argument(coefficients[static_cast<size_t>(a)].value) -
                                                  midpoint_argument(coefficients[static_cast<size_t>(b)].value)
                                            : std::numeric_limits<double>::quiet_NaN();
        for (slong j = 0; j < count; ++j) {
            const double angle =
                std::isfinite(phase)
                    ? (phase + 2 * std::numbers::pi * (static_cast<double>(j) + heuristics::binomial_seed_spacing_fraction)) /
                          static_cast<double>(count)
                    : placement_angle(j, count, static_cast<slong>(k));
            arf_set_d(direction.value(), std::cos(angle));
            arf_mul(acf_realref(points[next].value), radius.value(), direction.value(), prec, rounding);
            arf_set_d(direction.value(), std::sin(angle));
            arf_mul(acf_imagref(points[next].value), radius.value(), direction.value(), prec, rounding);
            ++next;
        }
    }
    // Missing leading data is an unresolved numerical setup, and so are
    // heights that are not exact. These midpoint locations remain candidates;
    // they do not certify a zero root or its count.
    if (next != points.size() || !exact_heights)
        for (auto &p : points) {
            arf_nan(acf_realref(p.value));
            arf_nan(acf_imagref(p.value));
        }
}
// The Newton correction P/P' at z under the model. When it is infinite (see
// NewtonPole) the correction is unset, *pole says so, and the caller takes the
// Aberth limit; that needs model stopping, and without it the node is left
// where it is (the caller's invalid step), as it always was.
bool secular_newton(acf_t correction, const acf_t z, const std::vector<Point> &nodes,
                    const std::vector<Point> &residues, slong prec, ModelStop *stop = nullptr,
                    const std::vector<uint8_t> *unresolved = nullptr, NewtonPole *pole = nullptr) {
    Point difference, inverse, term, sum, derivative, poles, denominator;
    Float absolute_sum, absolute_term;
    if (stop) arf_one(absolute_sum.value());
    slong coincident = -1;
    for (size_t i = 0; i < nodes.size(); ++i)
        if (acf_equal(z, nodes[i].value)) {
            coincident = static_cast<slong>(i);
            break;
        }
    if (coincident >= 0) {
        // At b_k, P/P' = a_k / D, D = sum_{j!=k}(a_j+a_k)/(b_k-b_j)-1.
        // Evaluate this removable limit instead of dividing by a pole.
        const auto k = static_cast<size_t>(coincident);
        for (size_t j = 0; j < nodes.size(); ++j)
            if (j != k) {
                acf_sub(difference.value, z, nodes[j].value, prec, rounding);
                if (is_zero(difference.value))
                    return false;
                acf_add(term.value, residues[j].value, residues[k].value, prec, rounding);
                acf_approx_div(term.value, term.value, difference.value, prec, rounding);
                acf_add(sum.value, sum.value, term.value, prec, rounding);
            }
        arf_sub_ui(acf_realref(sum.value), acf_realref(sum.value), 1, prec, rounding);
        if (is_zero(sum.value)) {
            if (!stop || !pole)
                return false;
            // D's rounding error is needed only here: a second pass over the
            // same terms sums their magnitudes.
            for (size_t j = 0; j < nodes.size(); ++j)
                if (j != k) {
                    acf_sub(difference.value, z, nodes[j].value, prec, rounding);
                    acf_add(term.value, residues[j].value, residues[k].value, prec, rounding);
                    acf_approx_div(term.value, term.value, difference.value, prec, rounding);
                    arf_abs(absolute_term.value(), acf_realref(term.value));
                    arf_add(absolute_sum.value(), absolute_sum.value(), absolute_term.value(), 64, ARF_RND_CEIL);
                    arf_abs(absolute_term.value(), acf_imagref(term.value));
                    arf_add(absolute_sum.value(), absolute_sum.value(), absolute_term.value(), 64, ARF_RND_CEIL);
                }
            Point scale;
            arf_set(acf_realref(scale.value), absolute_sum.value());
            pole->at = true;
            pole->unresolved = unresolved && (*unresolved)[k];
            pole->margin = magnitude_log2(residues[k].value) - magnitude_log2(scale.value) -
                           evaluation_roundoff_log2(nodes.size(), prec);
            zero(correction);
            return true;
        }
        acf_approx_div(correction, residues[k].value, sum.value, prec, rounding);
        if (stop && finite(correction))
            *stop = unresolved && (*unresolved)[k] && !is_zero(correction)
                ? ModelStop::EvaluationFloor
                : model_stop(NAN, NAN, magnitude_log2(correction), magnitude_log2(z), nodes.size(), prec);
    } else {
        for (size_t j = 0; j < nodes.size(); ++j) {
            acf_sub(difference.value, z, nodes[j].value, prec, rounding);
            if (is_zero(difference.value))
                return false;
            acf_approx_inv(inverse.value, difference.value, prec, rounding);
            acf_add(poles.value, poles.value, inverse.value, prec, rounding);
            acf_mul(term.value, residues[j].value, inverse.value, prec, rounding);
            if (stop) {
                arf_abs(absolute_term.value(), acf_realref(term.value));
                arf_add(absolute_sum.value(), absolute_sum.value(), absolute_term.value(), 64, ARF_RND_CEIL);
                arf_abs(absolute_term.value(), acf_imagref(term.value));
                arf_add(absolute_sum.value(), absolute_sum.value(), absolute_term.value(), 64, ARF_RND_CEIL);
            }
            acf_add(sum.value, sum.value, term.value, prec, rounding);
            acf_mul(term.value, term.value, inverse.value, prec, rounding);
            acf_sub(derivative.value, derivative.value, term.value, prec, rounding);
        }
        arf_sub_ui(acf_realref(sum.value), acf_realref(sum.value), 1, prec, rounding);
        acf_mul(denominator.value, sum.value, poles.value, prec, rounding);
        acf_add(denominator.value, denominator.value, derivative.value, prec, rounding);
        if (is_zero(denominator.value))
            return false;
        acf_approx_div(correction, sum.value, denominator.value, prec, rounding);
        if (stop && finite(correction)) {
            Point scale;
            arf_set(acf_realref(scale.value), absolute_sum.value());
            *stop = model_stop(magnitude_log2(sum.value), magnitude_log2(scale.value),
                               magnitude_log2(correction), magnitude_log2(z), nodes.size(), prec);
        }
    }
    return finite(correction);
}
} // namespace

void midpoint_polynomial_evaluation(acb_t value, acb_t derivative, acb_srcptr coefficients, slong length,
                                    const acb_t z, slong prec) {
    std::vector<Point> midpoints(static_cast<size_t>(length));
    std::vector<slong> indices;
    for (slong j = length - 1; j >= 0; --j) {
        midpoint(midpoints[static_cast<size_t>(j)].value, coefficients + j);
        if (!is_zero(midpoints[static_cast<size_t>(j)].value))
            indices.push_back(j);
    }
    Point point, v, d;
    midpoint(point.value, z);
    if (!indices.empty())
        evaluate(v.value, d.value, midpoints, indices, point.value, prec);
    point_ball(value, v.value);
    point_ball(derivative, d.value);
}

bool initial_candidates(acb_ptr out, const acb_poly_t polynomial, slong prec, bool binomial_seeds) {
    const slong n = acb_poly_degree(polynomial);
    if (n < 1 || acb_contains_zero(polynomial->coeffs + n))
        return false;
    std::vector<Point> coefficients(static_cast<size_t>(n + 1)), points(static_cast<size_t>(n));
    for (slong i = 0; i <= n; ++i)
        midpoint(coefficients[static_cast<size_t>(i)].value, polynomial->coeffs + i);
    initialize(points, coefficients, prec, nullptr, binomial_seeds);
    for (slong i = 0; i < n; ++i) {
        if (!finite(points[static_cast<size_t>(i)].value))
            return false;
        point_ball(out + i, points[static_cast<size_t>(i)].value);
    }
    return true;
}

void circle_candidates(acb_ptr out, slong count, const acb_t center, const arb_t radius, slong prec) {
    RealBall direction;
    for (slong j = 0; j < count; ++j) {
        const double angle = placement_angle(j, count, 1);
        arb_set_d(direction.value(), std::cos(angle));
        arb_mul(acb_realref(out + j), direction.value(), radius, prec);
        arb_set_d(direction.value(), std::sin(angle));
        arb_mul(acb_imagref(out + j), direction.value(), radius, prec);
        acb_add(out + j, out + j, center, prec);
        acb_get_mid(out + j, out + j);
    }
}
bool scale_candidates(acb_ptr out, const acb_poly_t polynomial, slong prec, bool &determined,
                      UndeterminedReach *reach) {
    const slong n = acb_poly_degree(polynomial);
    if (n < 1 || acb_contains_zero(polynomial->coeffs + n))
        return false;
    std::vector<Point> coefficients(static_cast<size_t>(n + 1)), points(static_cast<size_t>(n));
    std::vector<uint8_t> unresolved(static_cast<size_t>(n + 1), 0);
    Magnitude bound;
    for (slong i = 0; i <= n; ++i) {
        const auto c = polynomial->coeffs + i;
        if (!acb_is_finite(c))
            return false;
        auto &coefficient = coefficients[static_cast<size_t>(i)];
        // No correct bit: the midpoint is rounding noise (Arb's block
        // products return it as exactly zero), and only the enclosure's size
        // is known. An exact zero is a known coefficient.
        if (!acb_is_zero(c) && acb_rel_accuracy_bits(c) <= 0) {
            unresolved[static_cast<size_t>(i)] = 1;
            acb_get_mag(bound.value(), c);
            arf_set_mag(acf_realref(coefficient.value), bound.value());
            arf_zero(acf_imagref(coefficient.value));
        } else
            midpoint(coefficient.value, c);
    }
    std::vector<slong> vertices;
    initialize(points, coefficients, prec, &vertices);
    determined = std::none_of(vertices.begin(), vertices.end(),
                              [&](slong v) { return unresolved[static_cast<size_t>(v)] != 0; });
    if (reach) {
        *reach = {};
        for (size_t e = 1; e < vertices.size(); ++e) {
            const auto a = static_cast<size_t>(vertices[e - 1]), b = static_cast<size_t>(vertices[e]);
            if (!unresolved[a] && !unresolved[b])
                continue;
            const double radius_log2 = (magnitude_log2(coefficients[a].value) - magnitude_log2(coefficients[b].value)) /
                                       static_cast<double>(b - a);
            if (radius_log2 > reach->radius_log2)
                *reach = {radius_log2, static_cast<slong>(b - a)};
        }
    }
    for (slong i = 0; i < n; ++i) {
        if (!finite(points[static_cast<size_t>(i)].value))
            return false;
        point_ball(out + i, points[static_cast<size_t>(i)].value);
    }
    return true;
}

struct SecularCandidateWorkspace::Impl {
    std::vector<Point> nodes, weights, points, next;
    std::vector<uint8_t> unresolved; // see residue_unresolved
    std::map<slong, std::pair<std::vector<Point>, std::vector<slong>>> midpoint_polynomials;
    std::vector<slong> precision;
    MonicSource *source = nullptr; // Borrowed only while the owning refinement lives.
    uint64_t revision = 0, ceiling = 0;
    bool ready = false;
};
SecularCandidateWorkspace::SecularCandidateWorkspace() : impl_(std::make_unique<Impl>()) {}
SecularCandidateWorkspace::~SecularCandidateWorkspace() = default;

CandidateProgress
secular_candidate_packet(acb_ptr out, acb_srcptr initial, std::span<const ComplexBall> node_balls,
                         std::span<const ComplexBall> residue_balls, std::span<const uint8_t> inactive,
                         std::span<const slong> root_precision, unsigned sweeps, const cps_options &options,
                         certpolysolve::detail::CandidateStatistics &statistics, MonicSource *direct_source, bool use_tree, bool jacobi,
                         SecularCandidateWorkspace *workspace, uint64_t revision,
                         const EnginePolicy *policy, ModelStoppingState *model) {
    StageTimer timer(statistics.candidate_ns);
    CandidateProgress progress;
    const size_t n = node_balls.size();
    if (!n || residue_balls.size() != n || root_precision.size() != n || inactive.size() != n)
        return progress;
    if (direct_source || !policy || policy->model_stopping == SecularModelPolicy::Legacy) model = nullptr;
    if (model) model->prepare(revision, 3, root_precision, inactive);
    SecularCandidateWorkspace::Impl temporary;
    auto &cache = workspace ? *workspace->impl_ : temporary;
    const bool representation_changed = !cache.ready || cache.revision != revision || cache.nodes.size() != n;
    const bool precision_changed = cache.precision.size() != root_precision.size() ||
        !std::equal(cache.precision.begin(), cache.precision.end(), root_precision.begin()) ||
        cache.ceiling != options.max_work_bits || cache.source != direct_source;
    {
        TraceScope trace(policy, "secular_preparation", 0, n, &statistics);
        if (representation_changed) {
            cache.nodes.resize(n);
            cache.weights.resize(n);
            cache.points.resize(n);
            cache.next.resize(n);
            cache.unresolved.resize(n);
            for (size_t i = 0; i < n; ++i) {
                midpoint(cache.nodes[i].value, node_balls[i].value());
                midpoint(cache.weights[i].value, residue_balls[i].value());
                cache.unresolved[i] = residue_unresolved(residue_balls[i].value()) ||
                                      (model && model->unsteered(revision, i));
            }
            cache.revision = revision;
        }
        if (representation_changed || precision_changed) cache.midpoint_polynomials.clear();
        if (workspace) {
            if (precision_changed) cache.precision.assign(root_precision.begin(), root_precision.end());
            cache.ceiling = options.max_work_bits;
            cache.source = direct_source;
            if (policy && policy->study_statistics) {
                auto &counts = *policy->study_statistics;
                if (representation_changed || precision_changed) ++counts.secular_preparations;
                else ++counts.secular_reuses;
            }
        }
        cache.ready = true;
    }
    auto &nodes = cache.nodes;
    auto &weights = cache.weights;
    auto &points = cache.points;
    auto &next = cache.next;
    for (size_t i = 0; i < n; ++i) {
        midpoint(points[i].value, initial + i);
    }
    // A point's scratch, and the counters of a sweep's points: one of each for the loop in order, one a block for
    // the split loop (W9), the counters added after the join.
    struct PointScratch {
        Point correction, difference, inverse, interaction, denominator, delta, value, derivative;
        ComplexBall point, value_ball, derivative_ball;
        std::vector<size_t> neighbours;
    };
    struct SweepCounts {
        uint64_t evaluations = 0, terms = 0, interactions = 0, model_evaluations = 0, evaluation_ns = 0,
                 interaction_ns = 0;
        double largest = -std::numeric_limits<double>::infinity();
        double service_evaluation = 0, service_interaction = 0;
        bool converged = true;
        void add(const SweepCounts &o) {
            evaluations += o.evaluations;
            terms += o.terms;
            interactions += o.interactions;
            service_evaluation += o.service_evaluation;
            service_interaction += o.service_interaction;
            model_evaluations += o.model_evaluations;
            evaluation_ns += o.evaluation_ns;
            interaction_ns += o.interaction_ns;
            largest = std::max(largest, o.largest);
            converged = converged && o.converged;
        }
    };
    auto *polynomial_source = dynamic_cast<PolynomialSource *>(direct_source);
    auto &midpoint_polynomials = cache.midpoint_polynomials;
    // Each node's last step relative to its point (log2), within this packet.
    std::vector<double> step_log(n, std::numeric_limits<double>::quiet_NaN());
    for (unsigned sweep = 0; sweep < sweeps; ++sweep) {
        std::unique_ptr<InteractionTree> tree;
        slong tree_scale = 0;
        if (use_tree && jacobi && n >= heuristics::interaction_tree_degree) {
            StageTimer interaction_timer(statistics.interaction_ns);
            tree = prepare_interactions(points, tree_scale);
        }
        // One point of the sweep, shared by the loop in order and the split one (W9, parallel.hpp): next[i] from
        // points[i], with the caller's scratch `w` and counters `c`.
        const auto step = [&](size_t i, PointScratch &w, SweepCounts &c) {
            auto &correction = w.correction, &difference = w.difference, &inverse = w.inverse,
                 &interaction = w.interaction, &denominator = w.denominator, &delta = w.delta, &value = w.value,
                 &derivative = w.derivative;
            auto &point = w.point, &value_ball = w.value_ball, &derivative_ball = w.derivative_ball;
            auto &neighbours = w.neighbours;
            acf_set(next[i].value, points[i].value);
            if (inactive[i] || (model && model->settled(i)))
                return;
            const slong prec = root_precision[i];
            // Count each active point at its own arithmetic precision, including
            // failed evaluations. The enclosing packet owns the final charge.
            struct PointMeter {
                SweepCounts &counts; const EnginePolicy *policy; slong precision;
                uint64_t terms, interactions;
                ~PointMeter() {
                    if (!policy || !policy->research_service_meter) return;
                    const double weight = operation_weight(precision);
                    counts.service_evaluation += static_cast<double>(counts.terms - terms) * weight;
                    counts.service_interaction += static_cast<double>(counts.interactions - interactions) * weight;
                }
            } point_meter{c, policy, prec, c.terms, c.interactions};
            bool valid;
            ModelStop reason = ModelStop::Active;
            NewtonPole pole;
            {
                StageTimer evaluation_timer(c.evaluation_ns);
                ++c.evaluations;
                c.terms += direct_source ? direct_source->terms() : n;
                if (direct_source) {
                    if (polynomial_source && polynomial_source->midpoint_evaluation() &&
                        n <= heuristics::small_complex_degree) {
                        auto [entry, inserted] = midpoint_polynomials.try_emplace(prec);
                        auto &[coefficients, indices] = entry->second;
                        if (inserted) {
                            TraceScope preparation(policy, "candidate_coefficient_preparation",
                                                   static_cast<uint64_t>(prec), n, &statistics);
                            const auto handle = polynomial_source->coefficient_enclosure(
                                std::min<slong>(resources::precision_sum(prec, 16), static_cast<slong>(options.max_work_bits)));
                            const auto p = handle->polynomial.value();
                            coefficients.resize(static_cast<size_t>(acb_poly_length(p)));
                            for (slong j = acb_poly_degree(p); j >= 0; --j) {
                                midpoint(coefficients[static_cast<size_t>(j)].value, p->coeffs + j);
                                if (!is_zero(coefficients[static_cast<size_t>(j)].value))
                                    indices.push_back(j);
                            }
                        }
                        valid = coefficients.size() == n + 1 && !indices.empty();
                        if (valid)
                            evaluate(value.value, derivative.value, coefficients, indices, points[i].value,
                                     prec);
                    } else {
                        point_ball(point.value(), points[i].value);
                        valid = direct_source->evaluate(
                            value_ball.value(), derivative_ball.value(), point.value(),
                            std::min<slong>(resources::precision_sum(prec, 16), static_cast<slong>(options.max_work_bits)));
                        midpoint(value.value, value_ball.value());
                        midpoint(derivative.value, derivative_ball.value());
                    }
                    valid = valid && !is_zero(derivative.value);
                    if (valid)
                        acf_approx_div(correction.value, value.value, derivative.value, prec, rounding);
                    valid = valid && finite(correction.value);
                } else {
                    if (model) ++c.model_evaluations;
                    valid = secular_newton(correction.value, points[i].value, nodes, weights, prec,
                                            model ? &reason : nullptr, &cache.unresolved, &pole);
                }
            }
            if (valid && model && reason != ModelStop::Active) {
                model->stop(i, reason);
                return;
            }
            zero(interaction.value);
            if (valid) {
                StageTimer interaction_timer(c.interaction_ns);
                bool far_ok = false;
                if (tree) {
                    std::complex<double> far;
                    double tail;
                    far_ok = tree->evaluate(i, far, neighbours, tail);
                    if (far_ok) {
                        arf_set_d(acf_realref(interaction.value), far.real());
                        arf_set_d(acf_imagref(interaction.value), far.imag());
                        arf_mul_2exp_si(acf_realref(interaction.value), acf_realref(interaction.value),
                                        -tree_scale);
                        arf_mul_2exp_si(acf_imagref(interaction.value), acf_imagref(interaction.value),
                                        -tree_scale);
                    }
                }
                if (!far_ok) {
                    zero(interaction.value);
                    neighbours.resize(n);
                    for (size_t j = 0; j < n; ++j)
                        neighbours[j] = j;
                }
                valid = add_interactions(interaction.value, points[i].value, points, neighbours, i, prec,
                                         difference.value, inverse.value, c.interactions);
            }
            if (valid && pole.at) {
                // Newton's correction is infinite here; the Aberth step is its
                // limit -1/A when that is determined (NewtonPole). Otherwise
                // this node's step is not: it settles at the evaluation floor.
                if (!pole.resolved(magnitude_log2(interaction.value))) {
                    model->stop(i, ModelStop::EvaluationFloor);
                    return;
                }
                valid = !is_zero(interaction.value);
                if (valid) {
                    acf_approx_inv(delta.value, interaction.value, prec, rounding);
                    acf_neg(delta.value, delta.value);
                }
            } else if (valid)
                valid = aberth_delta(delta.value, correction.value, interaction.value, denominator.value, prec);
            if (valid) {
                acf_sub(next[i].value, points[i].value, delta.value, prec, rounding);
                valid = finite(next[i].value);
                // The correction relative to the point. -1e15 stands in for
                // the height of a point that is exactly zero (its height is
                // -inf), so that any nonzero correction there counts as
                // large. A nonzero point below 2^-1e15 is measured against
                // the stand-in too, and looks converged early; its proposal
                // is still certified separately. No input is known to need
                // better (audit, "Other literals of the same kind").
                const double relative =
                    magnitude_log2(delta.value) - std::max(-1.e15, magnitude_log2(points[i].value));
                c.largest = std::max(c.largest, relative);
                c.converged = c.converged && relative < -static_cast<double>(prec - 12);
                // A converged cubic step needs no confirming sweep. Aberth's
                // steps shrink with order 3 near a simple root. When the two
                // last steps show that order (the ratio of their relative
                // sizes, never taken above 3), the next step is predicted below
                // the arithmetic's resolution, where model_stop would settle
                // the node without moving it. The node settles now, at the
                // position the confirming sweep would have left it.
                // History: docs/code-history-20261021.md, candidates.cpp 6.
                if (model && valid) {
                    const double resolution = 2.0 - static_cast<double>(prec);
                    const double before = step_log[i];
                    if (std::isfinite(before) && before < 0 && relative < before &&
                        std::min(relative / before, 3.0) * relative < resolution)
                        model->stop(i, ModelStop::SmallCorrection);
                    step_log[i] = relative;
                }
            }
            if (!valid) {
                acf_set(next[i].value, points[i].value);
                c.largest = std::numeric_limits<double>::infinity();
                c.converged = false;
            }
        };
        SweepCounts counts;
        bool cancelled = false;
        if (jacobi && !direct_source && n >= heuristics::secular_packet_split_points && parallel_width() > 1) {
            // Blocks of secular_packet_block points, fixed by n: each point reads the start of the sweep and writes
            // its own next[i], step_log[i] and model stop; the counters add and the largest correction is a
            // maximum, so they are merged after the join. Workers poll the token and the clock only. The stage
            // timers add each block's time, so with workers they count every thread's (more than the wall time);
            // nothing decides by them.
            const size_t block = heuristics::secular_packet_block, blocks = (n + block - 1) / block;
            std::vector<SweepCounts> parts(blocks);
            std::atomic<bool> interrupted{false};
            parallel_for(SplitSite::SecularPacket, blocks, [&](size_t b) {
                PointScratch w;
                for (size_t i = b * block; i < std::min(n, (b + 1) * block); ++i) {
                    if (options.interrupted()) {
                        interrupted.store(true, std::memory_order_relaxed);
                        return;
                    }
                    step(i, w, parts[b]);
                }
            });
            for (const auto &part : parts) counts.add(part);
            cancelled = interrupted.load(std::memory_order_relaxed) || options.cancelled();
        } else {
            PointScratch w;
            for (size_t i = 0; i < n && !cancelled; ++i) {
                cancelled = options.cancelled();
                if (cancelled) break;
                step(i, w, counts);
                if (!jacobi)
                    acf_set(points[i].value, next[i].value);
            }
        }
        progress.service_evaluation += counts.service_evaluation;
        progress.service_interaction += counts.service_interaction;
        statistics.evaluations += counts.evaluations;
        statistics.evaluation_terms += counts.terms;
        statistics.root_interactions += counts.interactions;
        statistics.evaluation_ns += counts.evaluation_ns;
        statistics.interaction_ns += counts.interaction_ns;
        if (model) model->evaluated(counts.model_evaluations);
        if (cancelled) {
            progress.cancelled = true;
            return progress;
        }
        const double largest = counts.largest;
        const bool converged = counts.converged;
        if (jacobi)
            points.swap(next); // Direct/tree A/B uses the same Jacobi ordering.
        ++statistics.candidate_sweeps;
        progress.correction_log2 = largest;
        progress.model_exhausted = model && model->exhausted(inactive);
        progress.model_immediate = model && model->immediate();
        progress.converged = model ? progress.model_exhausted : converged;
        if (progress.converged)
            break;
    }
    progress.finite = true;
    for (size_t i = 0; i < n; ++i) {
        progress.finite = progress.finite && finite(points[i].value);
        point_ball(out + i, points[i].value);
    }
    if (model) model->trace(*policy, statistics, inactive, "adaptive",
                            *std::max_element(root_precision.begin(), root_precision.end()));
    return progress;
}

#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
// The research build's DoubleSeed start (it uses this file's private arithmetic, so it stays here, fenced).
bool binary64_seeds(acb_ptr out, const acb_poly_t polynomial, const cps_options &options,
                    certpolysolve::detail::CandidateStatistics &statistics) {
    using D = std::complex<double>;
    const auto n = acb_poly_degree(polynomial);
    if (n < 2 || n > heuristics::binary64_seed_degree)
        return false;
    std::vector<Point> coefficients(static_cast<size_t>(n + 1)), initial(static_cast<size_t>(n));
    for (slong j = 0; j <= n; ++j)
        midpoint(coefficients[static_cast<size_t>(j)].value, polynomial->coeffs + j);
    const double average =
        (magnitude_log2(coefficients.front().value) - magnitude_log2(coefficients.back().value)) /
        static_cast<double>(n);
    if (!std::isfinite(average))
        return false;
    // x -> 2^scale x brings the roots' geometric mean to 1 for binary64. The
    // clamps are range guards: a scale they cut short leaves coefficients
    // that binary64 cannot hold, the conversion below fails, and the caller
    // takes its ordinary multiprecision start.
    const auto scale = static_cast<slong>(std::clamp(std::round(average), -1.e9, 1.e9));
    double maximum = -std::numeric_limits<double>::infinity();
    for (slong j = 0; j <= n; ++j)
        maximum = std::max(maximum, magnitude_log2(coefficients[static_cast<size_t>(j)].value) +
                                        static_cast<double>(scale) * static_cast<double>(j));
    if (!std::isfinite(maximum) || std::abs(maximum) > 1.e15)
        return false;
    const auto normalization = static_cast<slong>(maximum);
    const auto finite_d = [](D z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); };
    const auto convert = [&](const Point &point, D &z) {
        z = {arf_get_d(acf_realref(point.value), ARF_RND_NEAR),
             arf_get_d(acf_imagref(point.value), ARF_RND_NEAR)};
        return finite_d(z) && (z.real() != 0 || arf_is_zero(acf_realref(point.value))) &&
               (z.imag() != 0 || arf_is_zero(acf_imagref(point.value)));
    };
    std::vector<D> values(coefficients.size()), points(initial.size());
    for (slong j = 0; j <= n; ++j) {
        auto &c = coefficients[static_cast<size_t>(j)];
        slong product, exponent;
        if (multiply_overflows(scale, j, &product) || subtract_overflows(product, normalization, &exponent))
            return false;
        arf_mul_2exp_si(acf_realref(c.value), acf_realref(c.value), exponent);
        arf_mul_2exp_si(acf_imagref(c.value), acf_imagref(c.value), exponent);
        if (!convert(c, values[static_cast<size_t>(j)]))
            return false;
    }
    initialize(initial, coefficients, 53);
    for (size_t i = 0; i < points.size(); ++i)
        if (!convert(initial[i], points[i]))
            return false;
    for (unsigned sweep = 0; sweep < 32; ++sweep) {
        double largest = 0;
        for (size_t i = 0; i < points.size(); ++i) {
            if (options.cancelled())
                return false;
            const D z = points[i];
            D value = values.back(), derivative = 0, interaction = 0;
            {
                StageTimer timer(statistics.evaluation_ns);
                for (size_t j = values.size() - 1; j-- > 0;) {
                    derivative = derivative * z + value;
                    value = value * z + values[j];
                }
                ++statistics.evaluations;
                statistics.evaluation_terms += values.size();
            }
            if (!finite_d(value) || !finite_d(derivative) || derivative == D{})
                return false;
            {
                StageTimer timer(statistics.interaction_ns);
                for (size_t j = 0; j < points.size(); ++j)
                    if (j != i) {
                        if (z == points[j])
                            return false;
                        interaction += 1.0 / (z - points[j]);
                        ++statistics.root_interactions;
                    }
            }
            const D newton = value / derivative;
            const D delta = newton / (1.0 - newton * interaction);
            if (!finite_d(delta) || !finite_d(z - delta))
                return false;
            largest = std::max(largest, std::abs(delta) / std::max(std::abs(z), 1.e-300));
            points[i] = z - delta;
        }
        ++statistics.candidate_sweeps;
        if (largest < 0x1p-40)
            break;
    }
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = 0; j < i; ++j)
            if (points[i] == points[j])
                return false;
        acb_set_d_d(out + i, points[i].real(), points[i].imag());
        acb_mul_2exp_si(out + i, out + i, scale);
    }
    return true;
}
#endif

bool regenerate_secular(std::vector<ComplexBall> &residues, const acb_poly_t p,
                        std::span<const ComplexBall> nodes, slong prec, const cps_options &options,
                        certpolysolve::detail::CandidateStatistics &statistics) {
    StageTimer timer(statistics.regeneration_ns);
    if (nodes.size() != static_cast<size_t>(acb_poly_degree(p)) ||
        acb_contains_zero(p->coeffs + acb_poly_degree(p)))
        return false;
    std::vector<ComplexBall> proposed(nodes.size());
    ComplexBall product, difference, derivative;
    CoefficientEvaluator evaluator(p);
    std::vector<ComplexBall> products(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (options.cancelled())
            return false;
        size_t length = 1;
        acb_set(products[0].value(), p->coeffs + acb_poly_degree(p));
        for (size_t j = 0; j < nodes.size(); ++j)
            if (j != i) {
                acb_sub(difference.value(), nodes[i].value(), nodes[j].value(), prec);
                if (acb_contains_zero(difference.value()))
                    return false;
                acb_set(products[length++].value(), difference.value());
            }
        // A balanced product controls rectangular enclosure wrapping; a long
        // chain can lose every useful bit even when all nodes are distinct.
        while (length > 1) {
            size_t next = 0;
            for (size_t j = 0; j < length; j += 2) {
                if (j + 1 < length)
                    acb_mul(products[next].value(), products[j].value(), products[j + 1].value(), prec);
                else
                    acb_set(products[next].value(), products[j].value());
                ++next;
            }
            length = next;
        }
        acb_set(product.value(), products[0].value());
        {
            StageTimer evaluation_timer(statistics.evaluation_ns);
            if (!evaluator.evaluate(proposed[i].value(), derivative.value(), nodes[i].value(), prec))
                return false;
            ++statistics.evaluations;
            statistics.evaluation_terms += static_cast<uint64_t>(acb_poly_length(p));
        }
        acb_div(proposed[i].value(), proposed[i].value(), product.value(), prec);
        acb_neg(proposed[i].value(), proposed[i].value());
        if (!acb_is_finite(proposed[i].value()))
            return false;
    }
    residues = std::move(proposed);
    ++statistics.secular_regenerations;
    return true;
}

struct CandidateWorkspace::Impl {
    EnginePolicy policy;
    ResearchLocalPacketService local_packet_service = nullptr;
    void *local_packet_data = nullptr;
    std::vector<Point> coefficients,points,nodes,weights;
    std::vector<Magnitude> coefficient_magnitudes;
    std::vector<uint8_t> arithmetic_stopped;
    ComplexPolynomial original,scaled;
    std::unique_ptr<MidpointPolynomialEvaluator> midpoint_evaluator;
    std::vector<slong> indices;
    std::vector<size_t> neighbours;
    Point value,derivative,correction,difference,inverse,interaction,denominator,delta;
    slong scale=0,precision=0;
    CandidateEngine engine=CandidateEngine::Aberth;
    bool ready=false,secular_ready=false;
};
CandidateWorkspace::CandidateWorkspace(const EnginePolicy *policy):impl_(std::make_unique<Impl>()) {
    if (policy) impl_->policy=*policy;
}
CandidateWorkspace::~CandidateWorkspace()=default;
void CandidateWorkspace::research_set_local_packet_service(ResearchLocalPacketService service, void *data) {
    impl_->local_packet_service = service;
    impl_->local_packet_data = data;
}

CandidateProgress candidate_packet(acb_ptr out, const acb_poly_t polynomial, acb_srcptr initial, slong prec,
                                   unsigned sweeps, CandidateEngine engine, const cps_options &options,
                                   certpolysolve::detail::CandidateStatistics &statistics, std::span<const uint8_t> inactive,
                                   bool use_tree, CandidateWorkspace *workspace) {
    StageTimer timer(statistics.candidate_ns);
    const auto n = acb_poly_degree(polynomial);
    CandidateProgress progress;
    if (n < 1 || acb_contains_zero(polynomial->coeffs + n))
        return progress;
    CandidateWorkspace::Impl transient;
    auto &cache=workspace?*workspace->impl_:transient;
    auto &coefficients=cache.coefficients;
    auto &points=cache.points;
    auto &scaled=cache.scaled;
    auto &indices=cache.indices;
    const bool prepare=!cache.ready || cache.precision!=prec || cache.engine!=engine ||
        !acb_poly_equal(cache.original.value(),polynomial);
    if (prepare) {
    TraceScope trace(&cache.policy,"adaptive_preparation",static_cast<uint64_t>(prec),static_cast<uint64_t>(n),&statistics);
    cache.ready=false;cache.secular_ready=false;
    cache.precision=prec;cache.engine=engine;
    acb_poly_set(cache.original.value(),polynomial);
    coefficients.resize(static_cast<size_t>(n+1));points.resize(static_cast<size_t>(n));
    indices.clear();
    acb_poly_zero(scaled.value());
    Point first, last;
    midpoint(first.value, polynomial->coeffs);
    midpoint(last.value, polynomial->coeffs + n);
    const double average =
        (magnitude_log2(first.value) - magnitude_log2(last.value)) / static_cast<double>(n);
    // A normalisation for conditioning, exact in arf whatever its value: the
    // clamps here and on `normalization` below only keep the exponents in a
    // slong (the overflow checks follow); a cut scale leaves the polynomial
    // less well scaled, not wrong.
    cache.scale =
        std::isfinite(average) ? static_cast<slong>(std::clamp(std::round(average), -1.e9, 1.e9)) : 0;
    const slong scale=cache.scale;
    double maximum = -std::numeric_limits<double>::infinity();
    for (slong j = 0; j <= n; ++j) {
        midpoint(coefficients[static_cast<size_t>(j)].value, polynomial->coeffs + j);
        const auto h = magnitude_log2(coefficients[static_cast<size_t>(j)].value) +
                       static_cast<double>(scale) * static_cast<double>(j);
        if (std::isfinite(h))
            maximum = std::max(maximum, h);
    }
    const slong normalization =
        std::isfinite(maximum) ? static_cast<slong>(std::clamp(maximum, -1.e15, 1.e15)) : 0;

    ComplexBall coefficient;
    for (slong j = n; j >= 0; --j) {
        slong product, exponent;
        if (multiply_overflows(scale, j, &product) || subtract_overflows(product, normalization, &exponent))
            return progress;
        acb_mul_2exp_si(coefficient.value(), polynomial->coeffs + j, exponent);
        acb_poly_set_coeff_acb(scaled.value(), j, coefficient.value());
        midpoint(coefficients[static_cast<size_t>(j)].value, coefficient.value());
        if (!is_zero(coefficients[static_cast<size_t>(j)].value) || j == 0)
            indices.push_back(j);
    }
    if (cache.policy.classical_error_stopping) {
        cache.coefficient_magnitudes.resize(coefficients.size());
        for (size_t j = 0; j < coefficients.size(); ++j)
            acb_get_mag(cache.coefficient_magnitudes[j].value(), scaled.value()->coeffs + j);
    }
    cache.midpoint_evaluator.reset();
    if (indices.size() * 4 >= coefficients.size())
        cache.midpoint_evaluator = std::make_unique<MidpointPolynomialEvaluator>(scaled.value(), prec);
    cache.ready=true;
    }
    const slong scale=cache.scale;
    bool reuse_stopping = !prepare && initial && cache.arithmetic_stopped.size() == points.size();
    bool have_initial = initial != nullptr;
    Point previous;
    if (initial)
        for (slong i = 0; i < n; ++i) {
            if (reuse_stopping)
                acf_set(previous.value, points[static_cast<size_t>(i)].value);
            midpoint(points[static_cast<size_t>(i)].value, initial + i);
            arf_mul_2exp_si(acf_realref(points[static_cast<size_t>(i)].value),
                            acf_realref(points[static_cast<size_t>(i)].value), -scale);
            arf_mul_2exp_si(acf_imagref(points[static_cast<size_t>(i)].value),
                            acf_imagref(points[static_cast<size_t>(i)].value), -scale);
            have_initial = have_initial && finite(points[static_cast<size_t>(i)].value);
            if (reuse_stopping && !acf_equal(previous.value, points[static_cast<size_t>(i)].value))
                reuse_stopping = false;
        }
    if (!have_initial)
        initialize(points, coefficients, prec);
    for (const auto &point : points)
        if (!finite(point.value))
            return progress;
    // Separate coincident proposals locally; the perturbation is heuristic.
    // In each group of equal points the first (lowest index) stays and every
    // later one moves. The groups come from sorting the points (value, then
    // index), where comparing every pair costs a packet O(n^2) however few
    // points it moves (a repair round of the certificate attempt moves a few
    // hundred). A moved point that lands exactly on another proposal is left
    // there (a proposal-only matter).
    // History: docs/code-history-20261021.md, candidates.cpp 7.
    {
        std::vector<size_t> order(points.size());
        std::iota(order.begin(), order.end(), size_t{0});
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            if (const int real = arf_cmp(acf_realref(points[a].value), acf_realref(points[b].value)))
                return real < 0;
            if (const int imag = arf_cmp(acf_imagref(points[a].value), acf_imagref(points[b].value)))
                return imag < 0;
            return a < b;
        });
        std::vector<size_t> coincident;
        for (size_t k = 1; k < order.size(); ++k)
            if (acf_equal(points[order[k]].value, points[order[k - 1]].value))
                coincident.push_back(order[k]);
        for (const size_t i : coincident) {
            reuse_stopping = false;
            Float perturbation;
            arf_one(perturbation.value());
            // The perturbation is 2^(-prec/2) of the point's size. The
            // clamp keeps the shift in a slong; for a point beyond
            // 2^(+-1e15) the perturbation has the wrong scale (it is
            // still only a proposal; no input is known to reach it).
            const auto h = magnitude_log2(points[i].value);
            const auto shift = std::isfinite(h) ? static_cast<slong>(std::clamp(h, -1.e15, 1.e15)) : 0;
            arf_mul_2exp_si(perturbation.value(), perturbation.value(), shift - prec / 2);
            arf_mul_ui(perturbation.value(), perturbation.value(), static_cast<ulong>(i + 1), prec, rounding);
            arf_add(acf_imagref(points[i].value), acf_imagref(points[i].value), perturbation.value(), prec,
                    rounding);
        }
    }
    auto &nodes=cache.nodes;auto &weights=cache.weights;
    bool secular = engine == CandidateEngine::Secular;
    if (secular && !cache.secular_ready) {
        TraceScope trace(&cache.policy,"candidate_regeneration",static_cast<uint64_t>(prec),static_cast<uint64_t>(n),&statistics);
        std::vector<ComplexBall> node_balls(points.size()), residue_balls;
        for (size_t i = 0; i < points.size(); ++i)
            point_ball(node_balls[i].value(), points[i].value);
        secular = regenerate_secular(residue_balls, scaled.value(), node_balls, prec, options, statistics);
        if (!secular && cache.policy.study == StudyMethod::Secular) {
            // An isolated secular study cannot silently substitute polynomial
            // Aberth corrections when interpolation was inconclusive.
            for (size_t i=0;i<points.size();++i) {
                point_ball(out+static_cast<slong>(i),points[i].value);
                acb_mul_2exp_si(out+static_cast<slong>(i),out+static_cast<slong>(i),scale);
            }
            return progress;
        }
        if (secular) {
            cache.secular_ready=true;
            nodes = points;
            weights.resize(points.size());
            for (size_t i = 0; i < points.size(); ++i)
                midpoint(weights[i].value, residue_balls[i].value());
        }
    }
    auto &value=cache.value;auto &derivative=cache.derivative;auto &correction=cache.correction;
    auto &difference=cache.difference;auto &inverse=cache.inverse;auto &interaction=cache.interaction;
    auto &denominator=cache.denominator;auto &delta=cache.delta;
    const bool error_stopping = cache.policy.classical_error_stopping && engine == CandidateEngine::Aberth;
    // Inclusion may reorder proposals, and callers may translate or replace
    // them. Reuse temporary stopping only when every coordinate and its index
    // agree with the last packet, at the same polynomial and precision.
    auto &stopped = cache.arithmetic_stopped;
    if (error_stopping && !reuse_stopping)
        stopped.assign(points.size(), 0);
    size_t active = points.size();
    uint64_t guarded_evaluations = 0;
    for (unsigned sweep = 0; sweep < sweeps; ++sweep) {
        std::unique_ptr<InteractionTree> tree;
        slong tree_scale = 0;
        if (use_tree && !durand_kerner_engine(engine) &&
            points.size() >= static_cast<size_t>(heuristics::interaction_tree_degree)) {
            StageTimer interaction_timer(statistics.interaction_ns);
            tree = prepare_interactions(points, tree_scale);
        }
        auto &neighbours=cache.neighbours;
        double largest = -std::numeric_limits<double>::infinity();
        active = 0;
        for (size_t i = 0; i < points.size(); ++i) {
            if (options.cancelled()) {
                progress.cancelled = true;
                return progress;
            }
            if (have_initial && inactive.size() == points.size() && inactive[i])
                continue;
            if (error_stopping && stopped[i])
                continue;
            auto z = points[i].value;
            bool valid;
            {
                StageTimer evaluation_timer(statistics.evaluation_ns);
                ++statistics.evaluations;
                statistics.evaluation_terms += secular ? points.size()
                                               : indices.size() * 4 < coefficients.size()
                                                   ? indices.size()
                                                   : coefficients.size();
                if (secular)
                    valid = secular_newton(correction.value, z, nodes, weights, prec);
                else {
                    if (cache.midpoint_evaluator)
                        cache.midpoint_evaluator->evaluate(value.value, derivative.value, z);
                    else
                        evaluate(value.value, derivative.value, coefficients, indices, z, prec);
                    if (error_stopping) {
                        ++guarded_evaluations;
                        if (classical_evaluation_floor(value.value, z, cache.coefficient_magnitudes, prec)) {
                            stopped[i] = 1;
                            continue;
                        }
                    }
                    valid = durand_kerner_engine(engine) || !is_zero(derivative.value);
                    if (durand_kerner_engine(engine))
                        acf_set(correction.value, value.value);
                    else if (valid)
                        acf_approx_div(correction.value, value.value, derivative.value, prec, rounding);
                    valid = valid && finite(correction.value);
                }
            }
            ++active;
            if (!valid) {
                largest = std::numeric_limits<double>::infinity();
                continue;
            }
            zero(interaction.value);
            if (durand_kerner_engine(engine))
                acf_set(interaction.value, coefficients.back().value);
            {
                StageTimer interaction_timer(statistics.interaction_ns);
                bool far_ok = false;
                if (tree) {
                    std::complex<double> far;
                    double tail;
                    far_ok = tree->evaluate(i, far, neighbours, tail);
                    if (far_ok) {
                        arf_set_d(acf_realref(interaction.value), far.real());
                        arf_set_d(acf_imagref(interaction.value), far.imag());
                        arf_mul_2exp_si(acf_realref(interaction.value), acf_realref(interaction.value),
                                        -tree_scale);
                        arf_mul_2exp_si(acf_imagref(interaction.value), acf_imagref(interaction.value),
                                        -tree_scale);
                    }
                }
                if (!far_ok) {
                    neighbours.resize(points.size());
                    std::iota(neighbours.begin(), neighbours.end(), size_t(0));
                }
                // The far field uses this sweep's starting positions; nearby
                // interactions use the latest Gauss--Seidel positions. This
                // approximation changes proposals only. Inclusion is separate.
                if (durand_kerner_engine(engine)) {
                    // The research build's Durand-Kerner: the product c_n prod_{j != i} (z - z_j).
                    for (size_t j : neighbours)
                        if (j != i) {
                            acf_sub(difference.value, z, points[j].value, prec, rounding);
                            if (is_zero(difference.value)) {
                                valid = false;
                                break;
                            }
                            acf_mul(interaction.value, interaction.value, difference.value, prec, rounding);
                            ++statistics.root_interactions;
                        }
                } else
                    valid = add_interactions(interaction.value, z, points, neighbours, i, prec, difference.value,
                                             inverse.value, statistics.root_interactions);
            }
            if (!valid) {
                largest = std::numeric_limits<double>::infinity();
                continue;
            }
            bool stepped;
            if (durand_kerner_engine(engine)) {
                // Durand-Kerner's correction, Weierstrass's: the value over c_n prod_{j != i} (z - z_j).
                acf_set(denominator.value, interaction.value);
                stepped = !is_zero(denominator.value);
                if (stepped)
                    acf_approx_div(delta.value, correction.value, denominator.value, prec, rounding);
            } else
                stepped = aberth_delta(delta.value, correction.value, interaction.value, denominator.value, prec);
            if (!stepped) {
                largest = std::numeric_limits<double>::infinity();
                continue;
            }
            if (!finite(delta.value)) {
                largest = std::numeric_limits<double>::infinity();
                continue;
            }
            double relative = magnitude_log2(delta.value);
            const auto size = magnitude_log2(z);
            if (std::isfinite(size))
                relative -= size;
            largest = std::max(largest, relative);
            acf_sub(z, z, delta.value, prec, rounding);
        }
        ++statistics.candidate_sweeps;
        progress.correction_log2 = largest;
        if (cache.local_packet_service) {
            const ResearchLocalPacketProgress event{statistics.evaluation_terms, statistics.root_interactions,
                statistics.candidate_sweeps, static_cast<uint64_t>(prec), 0};
            if (cache.local_packet_service(cache.local_packet_data, event) == ResearchServiceAction::StopNumerical) {
                progress.cancelled = true;
                return progress;
            }
        }
        if (error_stopping && active == 0)
            break;
        if (!error_stopping && largest < -0.75 * static_cast<double>(prec)) {
            progress.converged = true;
            break;
        }
    }
    if (error_stopping) {
        progress.continuation = active ? CandidateContinuation::NextPacketAtRequestedPrecision
                                       : CandidateContinuation::PrecisionDecision;
        if (cache.policy.trace) {
            TraceEvent event{"classical_arithmetic_stopping", "packet", static_cast<uint64_t>(prec),
                             static_cast<uint64_t>(active), &statistics, nullptr};
            event.engine = "aberth";
            event.tier = "adaptive";
            event.reason = active ? "arithmetic-progress" : "arithmetic-floor";
            event.model_evaluations = guarded_evaluations;
            event.model_error_stops = static_cast<uint64_t>(std::count(stopped.begin(), stopped.end(), uint8_t(1)));
            cache.policy.trace(cache.policy.trace_data, event);
        }
    }
    progress.finite = true;
    for (slong i = 0; i < n; ++i) {
        point_ball(out + i, points[static_cast<size_t>(i)].value);
        acb_mul_2exp_si(out + i, out + i, scale);
        progress.finite = progress.finite && acb_is_finite(out + i);
    }
    return progress;
}
} // namespace certpolysolve::detail
