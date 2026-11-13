// Axis charts classify the real line by Bernstein bisection and refine proved
// one-root brackets against the original input. PiecewiseAxisPass owns the round;
// its outer driver streams models, prices work and decides the precision ladder.
// Design: docs/c1-axis-charts-note-20261015.md.
#include "methods/piecewise/axis/piecewise_axis.hpp"
#include "core/cover.hpp"
#include "methods/piecewise/piecewise_statistics.hpp"
#include "input/coefficients/polynomial_input.hpp"
#include "certify/real_sign_refinement.hpp"
#include "numeric/models/root_model.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <optional>
#include <utility>

namespace certpolysolve::detail {

struct PiecewiseAxisPass::Impl {
    const DensePolynomial &input;
    PiecewiseRootStatistics &stats;
    const size_t zeros;
    const slong degree;
    slong model_bits = 0, precision = 0;
    // This round's charts and whether its complete axis pass has answered.
    struct AxisChart {
        RootModel model;
    };
    std::vector<AxisChart> axis_charts;
    bool axis_complete = false;
    // The terms |a_j| r^j of the nonzero part at a radius, as balls (the precision only bounds them).
    std::vector<RealBall> axis_terms(const arb_t radius) const;
    // The term v dominating at radius r, |a_v| r^v > sum over j != v of |a_j| r^j, proved; -1 if none does.
    slong dominant_at(const arb_t radius) const;
    // The annulus r1 <= |z| <= r2 holds no root when one term dominates at both radii: over it each other term
    // is convex in log r, so their sum is largest at an end. Otherwise split at the geometric mean.
    bool annulus_root_free(const arb_t r1, const arb_t r2, unsigned depth) const;
    // ---------------------------------------------------------------------------------------- C1, one dimension
    // docs/c1-axis-charts-note-20261015.md, "Revision after the prototype". An axis chart's model q is real on the
    // axis (real input, real centre): f(w) = p(c + R w) / ((c + R w)^l 2^s) is real there and has p's zeros. On an
    // interval I of the chart's segment, with |f - q| <= eps(I) (root_model_error):
    //   |q| > eps on I                   -> no root of p in I;
    //   q(a), q(b) of strict opposite signs beyond eps, and |q'| > eps' on I with |f' - q'| <= eps'
    //                                     -> f is monotone on I: exactly one root of p in I, real.
    // eps' by Cauchy: |e'(w)| <= max_{|z - w| = s} |e(z)| / s, root_model_error on I fattened by s, over s.
    // `segment`: the part of each chart's axis diameter used (13/16 of it: more than the 2/3 its ring's axis stretch
    // takes, so neighbouring segments overlap); `fattening`: s. The model's error bound holds on the closed unit disk
    // only (root_model_error keeps its uniform bound outside it), so the fattened rectangle [m - h - s, m + h + s] x
    // [-s, s] must stay inside: (13/16 + 1/8)^2 + (1/8)^2 < 1. axis_interval checks it all the same.
    static constexpr double axis_segment = 13.0 / 16, axis_fattening = 1.0 / 8;
    // Bits beyond a root's goal width at which its bracket is refined, written and judged: one binary64 mantissa's
    // worth above the goal, so neither a rounded midpoint nor the published ball's radius reaches the goal's width.
    static constexpr slong axis_guard_bits = 64;
    // Probe counters of one pass (CERTPOLYSOLVE_PROBE_AXIS).
    struct AxisCounters {
        uint64_t values = 0, intervals = 0, parts = 0, deepest = 0, refine_steps = 0, dominance = 0;
        // Why a part was split: its value enclosure held 0 with |q(m)| below the error (noise) or above it (the
        // enclosure's width); its slope held 0 likewise.
        uint64_t value_noise = 0, value_wide = 0, slope_noise = 0, slope_wide = 0, ends_undecided = 0;
        uint64_t convex = 0, convex_free = 0, convex_roots = 0, convex_noise = 0; // convex parts; root-free; counted; noise
        uint64_t source_refinements = 0, source_evaluations = 0, box_retries = 0;
        double gaps = 0, classify = 0, refine = 0, conversion = 0, splits = 0;
    };
    mutable AxisCounters axis_counters;
    // Re q and its first three derivatives. For real w, Re q(w) is the polynomial of q's coefficients' real parts
    // exactly, and |f - Re q| <= |f - q| because f is real there; likewise for the derivatives.
    struct AxisReal {
        std::array<RealPolynomial, 4> d;
    };
    AxisReal axis_real(const AxisChart &chart) const;
    // z = c + R w for a real ball w, as a complex ball.
    void axis_point(ComplexBall &z, const AxisChart &chart, const arb_t w) const;
    // f at a point w (a ball of radius 0): Re q(w) with the model's local error.
    void axis_value(RealBall &out, const AxisChart &chart, const AxisReal &real, const arb_t w) const;
    // f's model on a piece [a, b] of the axis (chart units) in the Bernstein basis of degree K = deg q:
    // Re q(a + (b - a) u) = sum_k beta_k C(K, k) u^k (1 - u)^(K - k), u in [0, 1]. Re q lies between the least and
    // the largest beta_k there (the convex hull property), and its slope dq/dw between the least and the largest
    // K (beta_(k+1) - beta_k) / (b - a). Unlike a majorant of the monomial form, these bounds are sharp for the
    // charts' exponential-like models (rings away from |z| = 1), and splitting (de Casteljau) only averages.
    struct Bernstein {
        arb_ptr beta = nullptr;
        slong length = 0;
        explicit Bernstein(slong n) : beta(_arb_vec_init(n)), length(n) {}
        ~Bernstein() { if (beta) _arb_vec_clear(beta, length); }
        Bernstein(const Bernstein &) = delete;
        Bernstein &operator=(const Bernstein &) = delete;
        Bernstein(Bernstein &&other) noexcept : beta(other.beta), length(other.length) { other.beta = nullptr; other.length = 0; }
        Bernstein &operator=(Bernstein &&other) noexcept {
            std::swap(beta, other.beta);
            std::swap(length, other.length);
            return *this;
        }
    };
    // Working precision of the Bernstein form: the monomial-to-Bernstein transform may cancel up to the ratio of
    // q's majorant to q on the piece, at most 2^K for the binomial weights; K extra bits cover it.
    slong bernstein_precision(const AxisReal &real) const { return precision + std::max<slong>(arb_poly_degree(real.d[0].value()), 0); }
    Bernstein axis_bernstein(const AxisReal &real, const arb_t a, const arb_t b) const;
    // de Casteljau at u = 1/2: `node` is consumed into the right half; `left` receives the left half.
    void axis_split(Bernstein &node, Bernstein &left, slong prec) const;
    // The model's error on [a, b] (value) and, by Cauchy, on its slope (chart units): sup |f - q| over [a, b]
    // fattened by s, a rectangle inside the unit disk, over s. False when that rectangle leaves the disk.
    bool axis_errors(RealBall &value_error, RealBall &slope_error, const AxisChart &chart, const arb_t a, const arb_t b) const;
    static int sign_of(const arb_t x) { return real_sign_of(x); }
    // +1 (-1) when every entry exceeds `margin` (is below -margin), else 0.
    static int strict_sign(const arb_srcptr entries, slong count, const arb_t margin, slong prec);
    // A part where f is convex or concave: sigma f'' >= m2 > 0 on [a, b] (chart units, sigma = +-1), from the
    // Bernstein form's second differences less eps'' = 2 sup|f - q| / s^2 (Cauchy; `slope_error` is sup|f - q| / s).
    // For x0 in [a, b], sigma f(x) >= sigma f(x0) + sigma f'(x0)(x - x0) + m2 (x - x0)^2 / 2 >= sigma f(x0) -
    // f'(x0)^2 / (2 m2): no root where that is positive. Where sigma f(x0) < 0 instead, the set where sigma f < 0 is
    // an interval about x0 and the ends' signs count the roots: two (both ends of sigma f positive), one (the ends
    // differ) or none (both negative). This settles a pair beside the axis or a close real pair in one step, where
    // bisection needs parts about as narrow as the pair. x0 is the minimum of sigma q by Newton on q'.
    // 1 classified (roots appended), 0 not convex or not decided, -1 convex with sigma f(x0) inside the model's
    // error: no part about x0 is decided by this model, so the round stops here.
    int axis_convex(const AxisChart &chart, const AxisReal &real, const arb_t a, const arb_t b, const arb_srcptr beta,
                    slong K, const RealBall &slope_error, int sa, int sb,
                    std::vector<std::pair<RealBall, RealBall>> &roots, slong prec) const;
    // Classifies [a0, b0] (chart units) by Bernstein bisection; appends one-root intervals. False where some part
    // stays undecided: at the depth limit, or where f and f' are both within the model's error at a part's middle
    // (a double root, or a pair beside the axis this model cannot tell apart; more model bits might).
    bool axis_classify(const AxisChart &chart, const AxisReal &real, const arb_t a0, const arb_t b0,
                       std::vector<std::pair<RealBall, RealBall>> &roots, const cps_options &o) const;
    // Refines a one-root interval (chart units) by signs inside it, the middle first, then its quarter points, to
    // the goal's relative width; false where the model's accuracy stops it first (no decided sign at those points).
    bool axis_refine(const AxisChart &chart, const AxisReal &real, RealBall &a, RealBall &b, const cps_options &o) const;
    // The input's real polynomial p(z) / z^zeros, its exact coefficients enclosed at `bits`; rebuilt when the
    // precision a refinement needs changes.
    mutable RealPolynomial source_holder;
    mutable slong source_bits = 0;
    const arb_poly_struct *axis_source(slong bits) const;
    // The input's real polynomial for real_sign and real_refine (real_sign_refinement.hpp): the dense evaluation
    // of axis_source at the precision asked.
    struct AxisSourceEvaluator {
        const Impl &impl;
        void evaluate2(arb_t value, arb_t slope, const arb_t x, slong bits) const {
            arb_poly_evaluate2(value, slope, impl.axis_source(bits), x, bits);
        }
    };
    // From the round's precision to the query's cap, with the axis guard; the evaluations counted in the probe.
    RealSignBudget axis_sign_budget(const cps_options &o) const {
        return {precision, real_sign_cap(o), axis_guard_bits, axis_counters.source_evaluations};
    }
    // Refines a one-root interval [lo, hi] of the axis (z, exact ends of proved opposite signs) to the goal's
    // relative width with the input itself (real_refine). False where a sign stays undecided at the cap.
    bool source_refine(RealBall &lo, RealBall &hi, slong goal_bits, const cps_options &o) const;
    // -------------------------------------------------------------------------------------------- C1, the pass
    // The rings stream from the innermost out, each with its two axis charts (sector 0 on the positive axis, S/2 on
    // the negative). Each side's charts form a chain outwards from the gap about 0; a chart's piece is settled when
    // its outer neighbour arrives (the cut between them) or at the round's end (its segment's outer end). So the
    // pass runs while the rings are built, and a round that cannot finish stops at its first certain failure, not
    // after building every ring of the round. History: docs/code-history-20261021.md, piecewise_roots.cpp 9.
    static constexpr size_t axis_positive = 0, axis_negative = 1;
    // A published axis root: its bracket [low, high] (exact ends, p's signs proved opposite) inside its proved
    // one-root region [region_lo, region_hi] (exact, inner ends of the region's balls).
    struct AxisRoot { RealBall low, high, region_lo, region_hi; };

    struct AxisSide {
        std::vector<size_t> chain; // chart indices, innermost first
        RealBall start;            // the outermost chart's piece start, in its chart units (exact)
    };
    std::array<AxisSide, 2> axis_side;
    std::vector<AxisReal> axis_reals;
    std::vector<Root> axis_found;
    // One-root regions classified this round, refined and published only when the whole round succeeds (a round
    // that fails later would discard them).
    struct AxisPending { size_t chart; std::vector<std::pair<RealBall, RealBall>> roots; };
    std::vector<AxisPending> axis_pending;
    const char *axis_failure = nullptr;
    bool axis_inner_done = false;
    double axis_pass_seconds = 0;
    void axis_reset();
    void axis_fail(const char *why);
    size_t axis_side_of(size_t i) const;
    void axis_z(RealBall &z, size_t i, const arb_t w) const;
    void axis_w(RealBall &w, size_t i, const arb_t z) const;
    // p's sign at an exact point z by chart i's model, 0 when undecided. f = p z^-l 2^-scale, so at a negative z
    // p's sign is f's times (-1)^l. The model is evaluated over the ball w(z), so the sign holds at z itself.
    int axis_p_sign(size_t i, const arb_t z) const;
    // The gap about 0, between the two sides' innermost segments: no root in the disk of the larger of their inner
    // ends' radii, where the constant term dominates.
    void axis_inner();
    // Refines a root to the call's goal and publishes it: the goal's width first, then 32 bits more at a time
    // until the published box lies strictly inside the region (its representation, a ball whose radius keeps 30
    // bits, widens it by about 2^-30 of its width, which can reach the other root of a close pair) and meets
    // the goal by the library's own test (an absolute tolerance included). The bracket is refined with the
    // input itself unless the model already did (`refined`, at the first width only).
    // History: docs/code-history-20261021.md, piecewise_roots.cpp 10.
    const char *axis_finish_root(AxisRoot &state, bool refined, const cps_options &o, Root &out);
    // Refines and publishes the roots of chart i's piece (each with its proved one-root region) into axis_found.
    const char *axis_publish(size_t i, std::vector<std::pair<RealBall, RealBall>> &roots, const cps_options &o);
    // Settles the piece of chart a, its side's outermost so far: its outer end is the cut with b, its outer
    // neighbour, or with none its segment's outer end; then its classification and its roots. A cut lies in the
    // two segments' overlap at a point where both models decide p's sign. Each chart's piece is rounded outwards
    // at the cut (the ends of the balls w(z)), so the two pieces overlap in a sliver about z, which a's model must
    // prove root-free: no root is missed or counted twice. Across a gap the annulus must be root-free.
    void axis_settle(size_t side, size_t a, std::optional<size_t> b, const cps_options &o);
    // A chart arrives from the model stream (axis mode).
    void axis_arrive(const RootModel &model, const cps_options &o);
    // The round's end: the outermost pieces, the region beyond them, and the answer when nothing failed.
    void axis_solve(const cps_options &o, std::vector<Root> &destination);

    Impl(const DensePolynomial &p, size_t zero_order, slong nonzero_degree, PiecewiseRootStatistics &statistics)
        : input(p), stats(statistics), zeros(zero_order), degree(nonzero_degree) {}
};

PiecewiseAxisPass::PiecewiseAxisPass(const DensePolynomial &input, size_t zeros, slong degree,
                                   PiecewiseRootStatistics &stats)
    : impl_(std::make_unique<Impl>(input, zeros, degree, stats)) {}
PiecewiseAxisPass::~PiecewiseAxisPass() = default;
void PiecewiseAxisPass::start_round(slong model_bits, slong precision) {
    impl_->model_bits = model_bits;
    impl_->precision = precision;
    impl_->axis_charts.clear();
    impl_->axis_reset();
}
void PiecewiseAxisPass::arrive(const RootModel &model, const cps_options &options) {
    impl_->axis_arrive(model, options);
}
void PiecewiseAxisPass::finish_round(const cps_options &options, std::vector<Root> &destination) {
    impl_->axis_solve(options, destination);
}
bool PiecewiseAxisPass::failed() const { return impl_->axis_failure != nullptr; }
bool PiecewiseAxisPass::complete() const { return impl_->axis_complete; }

std::vector<RealBall> PiecewiseAxisPass::Impl::axis_terms(const arb_t radius) const {
    constexpr slong bound_bits = 128;
    std::vector<RealBall> terms(static_cast<size_t>(degree + 1));
    ComplexBall coefficient;
    RealBall power;
    arb_one(power.value());
    for (slong j = 0; j <= degree; ++j) {
        input.coefficients[zeros + static_cast<size_t>(j)].enclose(coefficient.value(), bound_bits);
        acb_abs(terms[static_cast<size_t>(j)].value(), coefficient.value(), bound_bits);
        arb_mul(terms[static_cast<size_t>(j)].value(), terms[static_cast<size_t>(j)].value(), power.value(), bound_bits);
        arb_mul(power.value(), power.value(), radius, bound_bits);
    }
    return terms;
}
slong PiecewiseAxisPass::Impl::dominant_at(const arb_t radius) const {
    ++axis_counters.dominance;
    constexpr slong bound_bits = 128;
    const auto terms = axis_terms(radius);
    slong v = 0;
    for (slong j = 1; j <= degree; ++j)
        if (arf_cmp(arb_midref(terms[static_cast<size_t>(j)].value()), arb_midref(terms[static_cast<size_t>(v)].value())) > 0)
            v = j;
    RealBall rest;
    for (slong j = 0; j <= degree; ++j)
        if (j != v) arb_add(rest.value(), rest.value(), terms[static_cast<size_t>(j)].value(), bound_bits);
    return arb_gt(terms[static_cast<size_t>(v)].value(), rest.value()) ? v : -1;
}
bool PiecewiseAxisPass::Impl::annulus_root_free(const arb_t r1, const arb_t r2, unsigned depth) const {
    const slong v1 = dominant_at(r1), v2 = dominant_at(r2);
    if (v1 >= 0 && v1 == v2) return true;
    if (!depth) return false;
    constexpr slong bound_bits = 128;
    RealBall middle;
    arb_mul(middle.value(), r1, r2, bound_bits);
    arb_sqrt(middle.value(), middle.value(), bound_bits);
    return annulus_root_free(r1, middle.value(), depth - 1) && annulus_root_free(middle.value(), r2, depth - 1);
}
PiecewiseAxisPass::Impl::AxisReal PiecewiseAxisPass::Impl::axis_real(const AxisChart &chart) const {
    AxisReal real;
    const acb_poly_struct *q = chart.model.polynomial.value();
    for (slong k = 0; k < acb_poly_length(q); ++k) arb_poly_set_coeff_arb(real.d[0].value(), k, acb_realref(q->coeffs + k));
    for (size_t i = 1; i < real.d.size(); ++i)
        arb_poly_derivative(real.d[i].value(), real.d[i - 1].value(), precision);
    return real;
}
void PiecewiseAxisPass::Impl::axis_point(ComplexBall &z, const AxisChart &chart, const arb_t w) const {
    acb_zero(z.value());
    arb_mul(acb_realref(z.value()), w, chart.model.radius.value(), precision);
    acb_add(z.value(), z.value(), chart.model.center.value(), precision);
}
void PiecewiseAxisPass::Impl::axis_value(RealBall &out, const AxisChart &chart, const AxisReal &real, const arb_t w) const {
    ++axis_counters.values;
    ComplexBall z;
    arb_poly_evaluate(out.value(), real.d[0].value(), w, precision);
    axis_point(z, chart, w);
    const RealBall error = root_model_error(chart.model, z.value(), precision, nullptr, nullptr, nullptr);
    arb_add_error(out.value(), error.value());
}
PiecewiseAxisPass::Impl::Bernstein PiecewiseAxisPass::Impl::axis_bernstein(const AxisReal &real, const arb_t a,
    const arb_t b) const {
    const slong K = std::max<slong>(arb_poly_degree(real.d[0].value()), 0), prec = bernstein_precision(real);
    Bernstein out(K + 1);
    {
        RealPolynomial shifted;
        arb_poly_taylor_shift(shifted.value(), real.d[0].value(), a, prec);
        RealBall width, power, binomial;
        arb_sub(width.value(), b, a, prec);
        arb_one(power.value());
        for (slong j = 0; j <= K; ++j) {
            // e_j = d_j / C(K, j), d_j = the shifted coefficient times (b - a)^j
            arb_poly_get_coeff_arb(out.beta + j, shifted.value(), j);
            arb_mul(out.beta + j, out.beta + j, power.value(), prec);
            arb_bin_uiui(binomial.value(), static_cast<ulong>(K), static_cast<ulong>(j), prec);
            arb_div(out.beta + j, out.beta + j, binomial.value(), prec);
            arb_mul(power.value(), power.value(), width.value(), prec);
        }
    }
    // beta_k = sum_(j <= k) C(k, j) e_j, by Pascal's rule in place
    for (slong i = 1; i <= K; ++i)
        for (slong j = K; j >= i; --j) arb_add(out.beta + j, out.beta + j, out.beta + j - 1, prec);
    return out;
}
void PiecewiseAxisPass::Impl::axis_split(Bernstein &node, Bernstein &left, slong prec) const {
    const slong n = node.length;
    arb_set(left.beta, node.beta);
    for (slong r = 1; r < n; ++r) {
        for (slong k = 0; k < n - r; ++k) {
            arb_add(node.beta + k, node.beta + k, node.beta + k + 1, prec);
            arb_mul_2exp_si(node.beta + k, node.beta + k, -1);
        }
        arb_set(left.beta + r, node.beta);
    }
    // node.beta[j] = beta^(K - j)_j: the right half's j-th coefficient.
}
bool PiecewiseAxisPass::Impl::axis_errors(RealBall &value_error, RealBall &slope_error, const AxisChart &chart,
    const arb_t a, const arb_t b) const {
    RealBall middle, half, interval;
    arb_add(middle.value(), a, b, precision);
    arb_mul_2exp_si(middle.value(), middle.value(), -1);
    arb_sub(half.value(), b, a, precision);
    arb_mul_2exp_si(half.value(), half.value(), -1);
    Magnitude h;
    arb_get_mag(h.value(), half.value());
    arb_set(interval.value(), middle.value());
    arb_add_error_mag(interval.value(), h.value());
    ComplexBall z;
    axis_point(z, chart, interval.value());
    value_error = root_model_error(chart.model, z.value(), precision, nullptr, nullptr, nullptr);
    ComplexBall local, fat;
    acb_zero(local.value());
    arb_set(acb_realref(local.value()), middle.value());
    Magnitude along, across;
    RealBall reach;
    arb_set_d(reach.value(), axis_fattening);
    arb_get_mag(across.value(), reach.value());
    arb_add(reach.value(), reach.value(), half.value(), precision);
    arb_get_mag(along.value(), reach.value());
    mag_add(arb_radref(acb_realref(local.value())), arb_radref(acb_realref(local.value())), along.value());
    mag_add(arb_radref(acb_imagref(local.value())), arb_radref(acb_imagref(local.value())), across.value());
    Float outward;
    acb_get_abs_ubound_arf(outward.value(), local.value(), precision);
    if (!arf_is_finite(outward.value()) || arf_cmp_ui(outward.value(), 1) > 0) return false;
    acb_mul_arb(fat.value(), local.value(), chart.model.radius.value(), precision);
    acb_add(fat.value(), fat.value(), chart.model.center.value(), precision);
    slope_error = root_model_error(chart.model, fat.value(), precision, nullptr, nullptr, nullptr);
    RealBall s;
    arb_set_d(s.value(), axis_fattening);
    arb_div(slope_error.value(), slope_error.value(), s.value(), precision);
    return true;
}
int PiecewiseAxisPass::Impl::strict_sign(const arb_srcptr entries, slong count, const arb_t margin, slong prec) {
    RealBall low;
    arb_neg(low.value(), margin);
    bool positive = true, negative = true;
    for (slong k = 0; k < count && (positive || negative); ++k) {
        if (positive && !arb_gt(entries + k, margin)) positive = false;
        if (negative && !arb_lt(entries + k, low.value())) negative = false;
    }
    (void)prec;
    return positive ? 1 : negative ? -1 : 0;
}
int PiecewiseAxisPass::Impl::axis_convex(const AxisChart &chart, const AxisReal &real, const arb_t a, const arb_t b,
    const arb_srcptr beta, slong K, const RealBall &slope_error, int sa, int sb,
    std::vector<std::pair<RealBall, RealBall>> &roots, slong prec) const {
    if (K < 2) return 0;
    RealBall width, scale, curvature_error, s, d2, m2, low_m2;
    arb_sub(width.value(), b, a, prec);
    arb_set_si(scale.value(), K * (K - 1));
    arb_div(scale.value(), scale.value(), width.value(), prec);
    arb_div(scale.value(), scale.value(), width.value(), prec);
    arb_set_d(s.value(), axis_fattening);
    arb_mul_2exp_si(curvature_error.value(), slope_error.value(), 1);
    arb_div(curvature_error.value(), curvature_error.value(), s.value(), prec);
    int sigma = 0;
    for (slong k = 0; k + 2 <= K; ++k) {
        arb_sub(d2.value(), beta + k + 2, beta + k + 1, prec);
        arb_sub(d2.value(), d2.value(), beta + k + 1, prec);
        arb_add(d2.value(), d2.value(), beta + k, prec);
        arb_mul(d2.value(), d2.value(), scale.value(), prec);
        const int here = arb_gt(d2.value(), curvature_error.value()) ? 1 :
                         (arb_neg(low_m2.value(), curvature_error.value()), arb_lt(d2.value(), low_m2.value())) ? -1 : 0;
        if (!here || (sigma && here != sigma)) return 0;
        sigma = here;
        // m2 = the least sigma d2 less eps''
        if (sigma < 0) arb_neg(d2.value(), d2.value());
        if (k == 0) arb_set(m2.value(), d2.value());
        else arb_min(m2.value(), m2.value(), d2.value(), prec);
    }
    arb_sub(m2.value(), m2.value(), curvature_error.value(), prec);
    if (!arb_is_positive(m2.value())) return 0;
    ++axis_counters.convex;
    // x0: Newton on q' from the middle, kept inside [a, b]; it stops when a step no longer moves x0.
    RealBall x, d1, dd, step, lo, hi;
    arb_add(x.value(), a, b, prec);
    arb_mul_2exp_si(x.value(), x.value(), -1);
    arb_set(lo.value(), a);
    arb_set(hi.value(), b);
    for (slong iteration = 0; iteration < prec; ++iteration) {
        arb_poly_evaluate2(d1.value(), dd.value(), real.d[1].value(), x.value(), prec);
        if (!arb_is_nonzero(dd.value())) break;
        arb_div(step.value(), d1.value(), dd.value(), prec);
        Float previous;
        arf_set(previous.value(), arb_midref(x.value()));
        arb_sub(x.value(), x.value(), step.value(), prec);
        mag_zero(arb_radref(x.value()));
        if (arf_cmp(arb_midref(x.value()), arb_midref(lo.value())) < 0) arf_set(arb_midref(x.value()), arb_midref(lo.value()));
        if (arf_cmp(arb_midref(x.value()), arb_midref(hi.value())) > 0) arf_set(arb_midref(x.value()), arb_midref(hi.value()));
        if (arf_equal(previous.value(), arb_midref(x.value()))) break;
    }
    // sigma f(x0) and |f'(x0)| with their errors
    RealBall value, slope, bound;
    arb_poly_evaluate2(value.value(), slope.value(), real.d[0].value(), x.value(), prec);
    ComplexBall z;
    axis_point(z, chart, x.value());
    const RealBall error = root_model_error(chart.model, z.value(), precision, nullptr, nullptr, nullptr);
    if (sigma < 0) arb_neg(value.value(), value.value());
    arb_abs(slope.value(), slope.value());
    arb_add(slope.value(), slope.value(), slope_error.value(), prec);
    arb_sqr(bound.value(), slope.value(), prec);
    arb_div(bound.value(), bound.value(), m2.value(), prec);
    arb_mul_2exp_si(bound.value(), bound.value(), -1);
    // lower bound of sigma f on [a, b]: sigma q(x0) - eps(x0) - (|q'(x0)| + eps')^2 / (2 m2)
    RealBall lower;
    arb_sub(lower.value(), value.value(), error.value(), prec);
    arb_sub(lower.value(), lower.value(), bound.value(), prec);
    if (arb_is_positive(lower.value())) { ++axis_counters.convex_free; return 1; }
    // sigma f(x0) < 0: count by the ends
    RealBall upper, floor;
    arb_add(upper.value(), value.value(), error.value(), prec);
    arb_sub(floor.value(), value.value(), error.value(), prec);
    if (!arb_is_negative(upper.value()) && !arb_is_positive(floor.value())) { ++axis_counters.convex_noise; return -1; }
    if (!arb_is_negative(upper.value()) || !sa || !sb) return 0;
    const int ea = sigma * sa, eb = sigma * sb; // signs of sigma f at the ends
    RealBall xa, xb;
    arb_set(xa.value(), a);
    arb_set(xb.value(), b);
    if (ea > 0) roots.emplace_back(xa, x);
    if (eb > 0) roots.emplace_back(x, xb);
    ++axis_counters.convex_roots;
    return 1;
}
bool PiecewiseAxisPass::Impl::axis_classify(const AxisChart &chart, const AxisReal &real, const arb_t a0,
    const arb_t b0, std::vector<std::pair<RealBall, RealBall>> &roots, const cps_options &o) const {
    // A part at depth d is 2^-d of the piece: below the round's arithmetic precision its ends and the model's
    // values no longer separate, so deeper parts wait for a round with more bits.
    const unsigned depth_limit = static_cast<unsigned>(precision);
    const slong prec = bernstein_precision(real), K = std::max<slong>(arb_poly_degree(real.d[0].value()), 0);
    struct Part { RealBall a, b; Bernstein form; unsigned depth; };
    std::vector<Part> stack;
    {
        const auto begin = std::chrono::steady_clock::now();
        stack.push_back({RealBall(), RealBall(), axis_bernstein(real, a0, b0), 0});
        axis_counters.conversion += std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    }
    arb_set(stack.back().a.value(), a0);
    arb_set(stack.back().b.value(), b0);
    RealBall value_error, slope_error, width, scale;
    Bernstein slopes(std::max<slong>(K, 1));
    while (!stack.empty()) {
        if (o.cancelled()) return false;
        Part part = std::move(stack.back());
        stack.pop_back();
        ++axis_counters.parts;
        axis_counters.deepest = std::max<uint64_t>(axis_counters.deepest, part.depth);
        if (!axis_errors(value_error, slope_error, chart, part.a.value(), part.b.value())) return false;
        const arb_srcptr beta = part.form.beta;
        if (strict_sign(beta, K + 1, value_error.value(), prec)) continue; // no root of p here
        // slopes K (beta_(k+1) - beta_k) / (b - a)
        arb_sub(width.value(), part.b.value(), part.a.value(), prec);
        arb_set_si(scale.value(), K);
        arb_div(scale.value(), scale.value(), width.value(), prec);
        for (slong k = 0; k < K; ++k) {
            arb_sub(slopes.beta + k, beta + k + 1, beta + k, prec);
            arb_mul(slopes.beta + k, slopes.beta + k, scale.value(), prec);
        }
        const auto end_sign = [&](arb_srcptr end) {
            RealBall v;
            arb_set(v.value(), end);
            arb_add_error(v.value(), value_error.value());
            return sign_of(v.value());
        };
        const int sa = end_sign(beta), sb = end_sign(beta + K);
        const bool monotone = K >= 1 && strict_sign(slopes.beta, K, slope_error.value(), prec);
        if (monotone && sa && sb) {
            // f monotone on [a, b]: one root where the ends' signs differ, none where they agree.
            if (sa != sb) roots.emplace_back(part.a, part.b);
            continue;
        }
        // Monotone with an end inside the model's error: at most one root, on a side this model cannot tell;
        // splitting keeps that end. More model bits decide it.
        if (monotone) { ++axis_counters.ends_undecided; return false; }
        if (const int convex = axis_convex(chart, real, part.a.value(), part.b.value(), beta, K, slope_error, sa, sb, roots, prec)) {
            if (convex < 0) return false;
            continue;
        }
        // A research probe, compiled only with CERTPOLYSOLVE_DEBUG_PROBES (the release library reads no environment).
#ifdef CERTPOLYSOLVE_DEBUG_PROBES
        static const bool deep_probe = std::getenv("CERTPOLYSOLVE_PROBE_AXIS_DEEP") != nullptr;
#else
        constexpr bool deep_probe = false;
#endif
        if (deep_probe && part.depth % 100 == 50) {
            double lo = 1e300, hi = -1e300, rad = 0;
            for (slong k = 0; k <= K; ++k) {
                const double v = arf_get_d(arb_midref(beta + k), ARF_RND_NEAR);
                lo = std::min(lo, v); hi = std::max(hi, v);
                rad = std::max(rad, mag_get_d(arb_radref(beta + k)));
            }
            double slo = 1e300, shi = -1e300, srad = 0;
            for (slong k = 0; k < K; ++k) {
                const double v = arf_get_d(arb_midref(slopes.beta + k), ARF_RND_NEAR);
                slo = std::min(slo, v); shi = std::max(shi, v);
                srad = std::max(srad, mag_get_d(arb_radref(slopes.beta + k)));
            }
            std::fprintf(stderr, "PROBE axisdeep depth=%u a=%.17g K=%ld beta[%.3e, %.3e] rad %.2e eps %.2e slopes[%.3e, %.3e] rad %.2e epsd %.2e ends %d %d\n",
                         part.depth, arf_get_d(arb_midref(part.a.value()), ARF_RND_NEAR), static_cast<long>(K), lo, hi, rad,
                         arf_get_d(arb_midref(value_error.value()), ARF_RND_NEAR), slo, shi, srad,
                         arf_get_d(arb_midref(slope_error.value()), ARF_RND_NEAR), sa, sb);
        }
        if (part.depth >= depth_limit) return false;
        Bernstein left(K + 1);
        {
            const auto begin = std::chrono::steady_clock::now();
            axis_split(part.form, left, prec);
            axis_counters.splits += std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        }
        // Noise at the middle: f and f' both within the model's error there.
        {
            RealBall magnitude, slope;
            arb_abs(magnitude.value(), left.beta + K);
            const bool value_noise = !arb_gt(magnitude.value(), value_error.value());
            if (K >= 1) {
                // s'(1/2) = K (beta^(K-1)_1 - beta^(K-1)_0) at u = 1/2: the left half's last two entries, doubled
                arb_sub(slope.value(), left.beta + K, left.beta + K - 1, prec);
                arb_mul_2exp_si(slope.value(), slope.value(), 1);
                arb_mul(slope.value(), slope.value(), scale.value(), prec);
                arb_abs(magnitude.value(), slope.value());
            }
            const bool slope_noise = K < 1 || !arb_gt(magnitude.value(), slope_error.value());
            if (value_noise) ++axis_counters.value_noise; else ++axis_counters.value_wide;
            if (slope_noise) ++axis_counters.slope_noise; else ++axis_counters.slope_wide;
            if (value_noise && slope_noise) return false;
        }
        RealBall middle;
        arb_add(middle.value(), part.a.value(), part.b.value(), precision);
        arb_mul_2exp_si(middle.value(), middle.value(), -1);
        Part low{part.a, middle, std::move(left), part.depth + 1};
        Part high{middle, part.b, std::move(part.form), part.depth + 1};
        stack.push_back(std::move(high));
        stack.push_back(std::move(low));
    }
    return true;
}
bool PiecewiseAxisPass::Impl::axis_refine(const AxisChart &chart, const AxisReal &real, RealBall &a, RealBall &b,
    const cps_options &o) const {
    RealBall va, vm, point, width, za, limit;
    axis_value(va, chart, real, a.value());
    const int sa = sign_of(va.value());
    if (!sa) return false;
    for (int step = 0; step < 4096 && !o.cancelled(); ++step) {
        ++axis_counters.refine_steps;
        // Done when R (b - a) <= |c + R a| 2^-(bits + 5).
        arb_sub(width.value(), b.value(), a.value(), precision);
        arb_mul(width.value(), width.value(), chart.model.radius.value(), precision);
        arb_mul(za.value(), a.value(), chart.model.radius.value(), precision);
        arb_add(za.value(), za.value(), acb_realref(chart.model.center.value()), precision);
        arb_abs(limit.value(), za.value());
        arb_mul_2exp_si(limit.value(), limit.value(), -static_cast<slong>(o.bits) - 5);
        if (arb_le(width.value(), limit.value())) return true;
        int sm = 0;
        for (const ulong quarters : {2ul, 1ul, 3ul}) {
            // point = a + (b - a) quarters / 4, exact
            arb_sub(point.value(), b.value(), a.value(), precision);
            arb_mul_ui(point.value(), point.value(), quarters, precision);
            arb_mul_2exp_si(point.value(), point.value(), -2);
            arb_add(point.value(), point.value(), a.value(), precision);
            axis_value(vm, chart, real, point.value());
            sm = sign_of(vm.value());
            if (sm) break;
        }
        if (!sm) return false;
        if (sm == sa) arb_set(a.value(), point.value());
        else arb_set(b.value(), point.value());
    }
    return false;
}
const arb_poly_struct *PiecewiseAxisPass::Impl::axis_source(slong bits) const {
    arb_poly_struct *source_poly = source_holder.value();
    if (source_bits != bits) {
        ComplexBall coefficient;
        arb_poly_zero(source_poly);
        for (slong j = 0; j <= degree; ++j) {
            input.coefficients[zeros + static_cast<size_t>(j)].enclose(coefficient.value(), bits);
            arb_poly_set_coeff_arb(source_poly, j, acb_realref(coefficient.value()));
        }
        source_bits = bits;
    }
    return source_poly;
}
bool PiecewiseAxisPass::Impl::source_refine(RealBall &lo, RealBall &hi, slong goal_bits, const cps_options &o) const {
    AxisSourceEvaluator p{*this};
    return real_refine(p, lo, hi, goal_bits, axis_sign_budget(o), o);
}
void PiecewiseAxisPass::Impl::axis_reset() {
    for (auto &side : axis_side) side.chain.clear();
    axis_reals.clear();
    axis_found.clear();
    axis_pending.clear();
    axis_failure = nullptr;
    axis_inner_done = false;
    axis_counters = AxisCounters();
    axis_pass_seconds = 0;
}
void PiecewiseAxisPass::Impl::axis_fail(const char *why) { if (!axis_failure) axis_failure = why; }
size_t PiecewiseAxisPass::Impl::axis_side_of(size_t i) const {
    return arb_is_negative(acb_realref(axis_charts[i].model.center.value())) ? axis_negative : axis_positive;
}
void PiecewiseAxisPass::Impl::axis_z(RealBall &z, size_t i, const arb_t w) const {
    arb_mul(z.value(), w, axis_charts[i].model.radius.value(), precision);
    arb_add(z.value(), z.value(), acb_realref(axis_charts[i].model.center.value()), precision);
}
void PiecewiseAxisPass::Impl::axis_w(RealBall &w, size_t i, const arb_t z) const {
    arb_sub(w.value(), z, acb_realref(axis_charts[i].model.center.value()), precision);
    arb_div(w.value(), w.value(), axis_charts[i].model.radius.value(), precision);
}
int PiecewiseAxisPass::Impl::axis_p_sign(size_t i, const arb_t z) const {
    RealBall w, v;
    axis_w(w, i, z);
    axis_value(v, axis_charts[i], axis_reals[i], w.value());
    const int sign = sign_of(v.value());
    return arb_is_negative(z) && (axis_charts[i].model.lower_index & 1) ? -sign : sign;
}
void PiecewiseAxisPass::Impl::axis_inner() {
    if (axis_inner_done || axis_side[axis_positive].chain.empty() || axis_side[axis_negative].chain.empty()) return;
    axis_inner_done = true;
    const auto begin = std::chrono::steady_clock::now();
    RealBall w, zp, zn, r;
    arb_set_d(w.value(), -axis_segment);
    axis_z(zp, axis_side[axis_positive].chain.front(), w.value());
    arb_neg(w.value(), w.value());
    axis_z(zn, axis_side[axis_negative].chain.front(), w.value());
    if (!arb_is_positive(zp.value()) || !arb_is_negative(zn.value())) {
        axis_fail("the innermost segments reach 0");
    } else {
        arb_neg(zn.value(), zn.value());
        arb_max(r.value(), zp.value(), zn.value(), precision);
        Float upper;
        arb_get_ubound_arf(upper.value(), r.value(), precision);
        arb_set_arf(r.value(), upper.value());
        if (dominant_at(r.value()) != 0) { ++stats.axis_gap_failures; axis_fail("inner disk"); }
    }
    axis_counters.gaps += std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
}
const char *PiecewiseAxisPass::Impl::axis_finish_root(AxisRoot &state, bool refined, const cps_options &o, Root &out) {
    constexpr slong goal_step = 32;
    const slong first_goal = static_cast<slong>(o.bits) + 5;
    for (slong goal = first_goal; ; goal += goal_step) {
        if (goal > static_cast<slong>(std::min<uint64_t>(o.max_work_bits, std::numeric_limits<slong>::max() / 2))) {
            ++stats.axis_unrefined;
            return "a root's box does not fit its region at the precision cap";
        }
        if (!refined || goal > first_goal) {
            ++axis_counters.source_refinements;
            if (!source_refine(state.low, state.high, goal, o)) {
                ++stats.axis_unrefined;
                return "a root not refined to the goal";
            }
        }
        // A box narrower than the working precision cannot be written or judged: work at the goal's width
        // and a guard (the refinement does the same).
        const slong work = std::max(precision, resources::precision_sum(goal, axis_guard_bits));
        Root root;
        acb_zero(root.box.value());
        arb_union(acb_realref(root.box.value()), state.low.value(), state.high.value(), work);
        root.count = root.multiplicity = 1;
        root.realness = CPS_REAL;
        root.certificate = "axis chart real model: a sign change inside a region of one root";
        root.publish(work);
        Float box_lo, box_hi;
        arb_get_lbound_arf(box_lo.value(), acb_realref(root.box.value()), work);
        arb_get_ubound_arf(box_hi.value(), acb_realref(root.box.value()), work);
        if (arf_cmp(box_lo.value(), arb_midref(state.region_lo.value())) > 0 &&
            arf_cmp(box_hi.value(), arb_midref(state.region_hi.value())) < 0 && meets_goal(root, o, work)) {
            out = std::move(root);
            return nullptr;
        }
        ++axis_counters.box_retries;
    }
}
const char *PiecewiseAxisPass::Impl::axis_publish(size_t i, std::vector<std::pair<RealBall, RealBall>> &roots,
    const cps_options &o) {
    const auto &chart = axis_charts[i];
    const auto &real = axis_reals[i];
    for (auto &[ra, rb] : roots) {
        const auto begin = std::chrono::steady_clock::now();
        const auto spent = [&] {
            axis_counters.refine += std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        };
        const auto to_z = [&](RealBall &z, const RealBall &w) {
            arb_mul(z.value(), w.value(), chart.model.radius.value(), precision);
            arb_add(z.value(), z.value(), acb_realref(chart.model.center.value()), precision);
        };
        // The proved one-root region in z, shrunk to the inner ends of its balls. A published box inside it
        // holds this root and no other: the box's own representation (a ball whose radius keeps 30 bits)
        // widens it by about 2^-30 of its width, which at the goal's width alone can reach the other root
        // of a close pair. History: docs/code-history-20261021.md, piecewise_roots.cpp 11.
        RealBall region_lo, region_hi;
        {
            RealBall z;
            Float end;
            to_z(z, ra);
            arb_get_ubound_arf(end.value(), z.value(), precision);
            arb_set_arf(region_lo.value(), end.value());
            to_z(z, rb);
            arb_get_lbound_arf(end.value(), z.value(), precision);
            arb_set_arf(region_hi.value(), end.value());
        }
        // The model refines while it decides signs (cheap: its degree); the input itself takes over there.
        const bool refined = axis_refine(chart, real, ra, rb, o);
        RealBall low, high;
        {
            RealBall z;
            Float end;
            to_z(z, ra);
            arb_get_lbound_arf(end.value(), z.value(), precision);
            arb_set_arf(low.value(), end.value());
            to_z(z, rb);
            arb_get_ubound_arf(end.value(), z.value(), precision);
            arb_set_arf(high.value(), end.value());
            arb_max(low.value(), low.value(), region_lo.value(), precision);
            arb_min(high.value(), high.value(), region_hi.value(), precision);
        }
        AxisRoot state{low, high, region_lo, region_hi};
        Root root;
        if (const char *why = axis_finish_root(state, refined, o, root)) {
            spent();
            return why;
        }
        axis_found.push_back(std::move(root));
        spent();
    }
    return nullptr;
}
void PiecewiseAxisPass::Impl::axis_settle(size_t side, size_t a, std::optional<size_t> b, const cps_options &o) {
    if (axis_failure) return;
    const auto begin = std::chrono::steady_clock::now();
    const auto spent = [&] { axis_pass_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count(); };
    const bool positive = side == axis_positive;
    // Outward is +w on the positive side and -w on the negative; a segment's outer end is w = +-sigma.
    RealBall outer, end_a, start_b;
    arb_set_d(outer.value(), positive ? axis_segment : -axis_segment);
    arb_set(end_a.value(), outer.value());
    if (b) {
        const auto gaps_begin = std::chrono::steady_clock::now();
        const auto gaps_spent = [&] { axis_counters.gaps += std::chrono::duration<double>(std::chrono::steady_clock::now() - gaps_begin).count(); };
        // Along this side, t = |z|: a's outer end, b's inner and outer ends.
        RealBall inner_w, za_out, zb_in, zb_out, ta_out, tb_in, tb_out;
        arb_neg(inner_w.value(), outer.value());
        axis_z(za_out, a, outer.value());
        axis_z(zb_in, *b, inner_w.value());
        axis_z(zb_out, *b, outer.value());
        arb_abs(ta_out.value(), za_out.value());
        arb_abs(tb_in.value(), zb_in.value());
        arb_abs(tb_out.value(), zb_out.value());
        if (arb_lt(tb_in.value(), ta_out.value())) {
            // An overlap [tb_in, min(ta_out, tb_out)]: its middle first, then points across it.
            RealBall upper;
            arb_min(upper.value(), ta_out.value(), tb_out.value(), precision);
            static constexpr int order[] = {8, 4, 12, 2, 6, 10, 14, 1, 3, 5, 7, 9, 11, 13, 15};
            bool found = false;
            RealBall z;
            for (const int k : order) {
                arb_sub(z.value(), upper.value(), tb_in.value(), precision);
                arb_mul_si(z.value(), z.value(), k, precision);
                arb_mul_2exp_si(z.value(), z.value(), -4);
                arb_add(z.value(), z.value(), tb_in.value(), precision);
                mag_zero(arb_radref(z.value()));
                if (!positive) arb_neg(z.value(), z.value());
                const int sa = axis_p_sign(a, z.value()), sb = axis_p_sign(*b, z.value());
                if (sa && sb && sa != sb) { ++stats.axis_sign_conflicts; gaps_spent(); spent(); axis_fail("sign conflict between neighbouring models"); return; }
                if (sa && sa == sb) { found = true; break; }
            }
            if (!found) { gaps_spent(); spent(); axis_fail("no decided boundary in an overlap"); return; }
            RealBall wa, wb;
            axis_w(wa, a, z.value());
            axis_w(wb, *b, z.value());
            Float end;
            // a's piece reaches past z, b's starts before it.
            if (positive) arb_get_ubound_arf(end.value(), wa.value(), precision);
            else arb_get_lbound_arf(end.value(), wa.value(), precision);
            arb_set_arf(end_a.value(), end.value());
            if (positive) arb_get_lbound_arf(end.value(), wb.value(), precision);
            else arb_get_ubound_arf(end.value(), wb.value(), precision);
            arb_set_arf(start_b.value(), end.value());
            // The sliver between b's start and a's end, as a ball in a's units: a's model must decide p's sign on it.
            RealBall z1, z2, hull, w_hull, v;
            axis_z(z1, a, end_a.value());
            axis_z(z2, *b, start_b.value());
            arb_union(hull.value(), z1.value(), z2.value(), precision);
            axis_w(w_hull, a, hull.value());
            axis_value(v, axis_charts[a], axis_reals[a], w_hull.value());
            if (!sign_of(v.value())) { gaps_spent(); spent(); axis_fail("a cut's sliver undecided"); return; }
        } else {
            // A gap [ta_out, tb_in] of radii: its annulus must be root-free.
            constexpr unsigned split_depth = 24;
            RealBall r1, r2;
            Float bound;
            arb_get_lbound_arf(bound.value(), ta_out.value(), precision);
            arb_set_arf(r1.value(), bound.value());
            arb_get_ubound_arf(bound.value(), tb_in.value(), precision);
            arb_set_arf(r2.value(), bound.value());
            if (!annulus_root_free(r1.value(), r2.value(), split_depth)) {
                ++stats.axis_gap_failures;
                gaps_spent();
                spent();
                axis_fail("a gap's annulus not root-free");
                return;
            }
            arb_set(start_b.value(), inner_w.value());
        }
        gaps_spent();
    }
    // a's piece as an increasing interval of its chart units.
    RealBall lo, hi;
    arb_set(lo.value(), positive ? axis_side[side].start.value() : end_a.value());
    arb_set(hi.value(), positive ? end_a.value() : axis_side[side].start.value());
    if (!arb_lt(lo.value(), hi.value())) {
        ++stats.axis_gap_failures;
        spent();
        axis_fail("pieces out of order (a segment overlapped twice or nested)");
        return;
    }
    std::vector<std::pair<RealBall, RealBall>> roots;
    bool classified;
    {
        const auto classify_begin = std::chrono::steady_clock::now();
        classified = axis_classify(axis_charts[a], axis_reals[a], lo.value(), hi.value(), roots, o);
        axis_counters.classify += std::chrono::duration<double>(std::chrono::steady_clock::now() - classify_begin).count();
    }
    if (!classified) {
        ++stats.axis_undecided;
        spent();
        axis_fail(o.cancelled() ? "cancelled" : "a piece undecided");
        return;
    }
    if (!roots.empty()) axis_pending.push_back({a, std::move(roots)});
    if (b) arb_set(axis_side[side].start.value(), start_b.value());
    spent();
}
void PiecewiseAxisPass::Impl::axis_arrive(const RootModel &model, const cps_options &o) {
    AxisChart chart;
    chart.model = model;
    axis_charts.push_back(std::move(chart));
    axis_reals.push_back(axis_real(axis_charts.back()));
    if (axis_failure) return;
    const size_t i = axis_charts.size() - 1;
    const size_t side = axis_side_of(i);
    auto &chain = axis_side[side].chain;
    if (chain.empty()) {
        // Its piece starts at its segment's inner end: w = -sigma on the positive side, +sigma on the negative.
        arb_set_d(axis_side[side].start.value(), side == axis_positive ? -axis_segment : axis_segment);
        chain.push_back(i);
        axis_inner();
        return;
    }
    const size_t previous = chain.back();
    chain.push_back(i);
    axis_settle(side, previous, i, o);
}
void PiecewiseAxisPass::Impl::axis_solve(const cps_options &o, std::vector<Root> &destination) {
    // A research probe, compiled only with CERTPOLYSOLVE_DEBUG_PROBES (the release library reads no environment).
#ifdef CERTPOLYSOLVE_DEBUG_PROBES
    static const bool probe = std::getenv("CERTPOLYSOLVE_PROBE_AXIS") != nullptr;
#else
    constexpr bool probe = false;
#endif
    if (axis_side[axis_positive].chain.empty() || axis_side[axis_negative].chain.empty()) axis_fail("a side without axis charts");
    if (!axis_inner_done) axis_fail("inner disk");
    for (const size_t side : {axis_positive, axis_negative})
        if (!axis_failure) axis_settle(side, axis_side[side].chain.back(), std::nullopt, o);
    if (!axis_failure) {
        // Beyond the outermost segments: the leading term dominates at the smaller of their outer ends' radii.
        RealBall w, zp, zn, r;
        arb_set_d(w.value(), axis_segment);
        axis_z(zp, axis_side[axis_positive].chain.back(), w.value());
        arb_neg(w.value(), w.value());
        axis_z(zn, axis_side[axis_negative].chain.back(), w.value());
        arb_neg(zn.value(), zn.value());
        arb_min(r.value(), zp.value(), zn.value(), precision);
        Float lower;
        arb_get_lbound_arf(lower.value(), r.value(), precision);
        arb_set_arf(r.value(), lower.value());
        if (dominant_at(r.value()) != degree) { ++stats.axis_gap_failures; axis_fail("outer region"); }
    }
    if (!axis_failure) {
        const auto begin = std::chrono::steady_clock::now();
        for (auto &pending : axis_pending)
            if (const char *why = axis_publish(pending.chart, pending.roots, o)) { axis_fail(why); break; }
        axis_pass_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    }
    const size_t published = axis_failure ? 0 : axis_found.size();
    if (!axis_failure) {
        for (auto &root : axis_found)
            destination.push_back(std::move(root));
        axis_found.clear();
    }
    axis_complete = axis_failure == nullptr;
    if (probe)
        std::fprintf(stderr, "PROBE axis1d round=%lu model_bits=%ld precision=%ld charts=%zu roots=%zu construction_total=%.3f seconds=%.3f "
                     "(gaps %.3f classify %.3f [conversion %.3f splits %.3f] refine %.3f) values=%lu intervals=%lu parts=%lu deepest=%lu "
                     "refine_steps=%lu source refinements %lu evaluations %lu box retries %lu dominance=%lu split: value noise %lu wide %lu, ends undecided %lu, slope noise %lu wide %lu; convex %lu free %lu counted %lu: %s\n",
                     static_cast<unsigned long>(stats.rounds), static_cast<long>(model_bits), static_cast<long>(precision),
                     axis_charts.size(), published, static_cast<double>(stats.construction_ns) / 1e9, axis_pass_seconds,
                     axis_counters.gaps, axis_counters.classify, axis_counters.conversion, axis_counters.splits, axis_counters.refine,
                     static_cast<unsigned long>(axis_counters.values), static_cast<unsigned long>(axis_counters.intervals),
                     static_cast<unsigned long>(axis_counters.parts), static_cast<unsigned long>(axis_counters.deepest),
                     static_cast<unsigned long>(axis_counters.refine_steps),
                     static_cast<unsigned long>(axis_counters.source_refinements), static_cast<unsigned long>(axis_counters.source_evaluations),
                     static_cast<unsigned long>(axis_counters.box_retries),
                     static_cast<unsigned long>(axis_counters.dominance),
                     static_cast<unsigned long>(axis_counters.value_noise), static_cast<unsigned long>(axis_counters.value_wide),
                     static_cast<unsigned long>(axis_counters.ends_undecided), static_cast<unsigned long>(axis_counters.slope_noise),
                     static_cast<unsigned long>(axis_counters.slope_wide), static_cast<unsigned long>(axis_counters.convex),
                     static_cast<unsigned long>(axis_counters.convex_free), static_cast<unsigned long>(axis_counters.convex_roots),
                     axis_failure ? axis_failure : "complete");
}

} // namespace certpolysolve::detail
