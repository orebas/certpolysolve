// The all-real seed pass (real_line_seeds): binary64 starting points on the real line for an integer polynomial whose
// roots are all real, from one Chebyshev interpolant of f = p 2^-phi on the interval of its roots. The seeds only start
// an iteration whose certified route decides every answer; a pass that stops gives none, only a RealLineSeedStop.
//
// A root at zero is divided out and seeded by 0, and an input that violates Newton's inequalities
// (newton_inequality_failure, exact) has a non-real root and is refused. The stages are numbered as in real_line_seeds.
// (1) Laguerre's method from outside a Fujiwara-type bound (Pass::root_bound, Pass::laguerre) reaches the extreme
// roots, which on a real-rooted polynomial it never passes, and its discriminant tests Laguerre's inequality
// (n - 1) p'^2 >= n p p''. (2) The range of log2 |p| at range_samples points sets the node count N and the precisions,
// or ends the pass (WideRange). (3) One change of basis gives the Chebyshev coefficients
// (Pass::chebyshev_from_monomials), its cancellation bounded by the same computation on absolute values, after a
// storage prediction (Memory), and (4) one DCT-I in ball arithmetic their values v_j at the nodes cos(pi j / N)
// (ball_cosine_transform), with a floor on their certain sign changes (FewCells). (5) The envelope phi of log2 |v_j| is
// a running median (running_median), smoothed by a Gaussian in the Chebyshev angle and shifted to lie on or above every
// local maximum. (6) The values f_j = v_j 2^-phi_j are rounded to binary64, (4) redone at a higher precision until
// every radius is below 2^-value_bits of 2^phi_j. (7) Their Chebyshev coefficients follow (chebyshev_coefficients), N
// doubled until the tail is small, and (8) the interpolant's roots are its sign changes on a finer grid
// (chebyshev_values), refined by the Illinois method on Clenshaw's recurrence (chebyshev_evaluate), and near-touching
// minima; more than n are noise (NoisyInterpolant). (9) Among the cells between nodes, one without a seed whose ends
// have certain opposite signs gets the secant's zero, a seed alone in such a cell is kept, and any other is tested by
// p's certain signs on its two sides (Pass::certain_sign); at most n are kept, verified ones preferred.
//
// State: none. Each call builds its own Pass (p, p', p''), balls and Binary64Transform (transform_of), and keeps none.
// Entry points: real_line_seeds, called by solver.cpp (seed_from_real_line) and real_line_routing.cpp
// (every_root_seeded), and real_line_seed_stop_name, by solver.cpp; both also by tests/test_real_line_seeds.cpp and
// benchmarks/real-roots/real_line_seeds_driver.cpp, and the test alone calls newton_inequality_failure,
// Binary64Transform, chebyshev_coefficients, chebyshev_values and chebyshev_evaluate.
// Not here: the limits from a call's options and Auto's check before routing (real_line_seed_limits, every_root_seeded:
// real_line_routing.cpp); the flag for a factor's start (cps_options::real_line_seeds, set in real_only_route.cpp and
// session.cpp); the degree floor and the cell share (heuristics.hpp); the same-bits rule (real_line_seeds.hpp).
// Design notes: docs/real-line-seed-pass-design-20261020.md, docs/real-roots-real-line-seeds-note-20261018.md (why an
// interpolant of p itself fails) and docs/allroots-seed-routing-note-20261020.md (every_root_seeded).
#include "numeric/candidates/real_line_seeds.hpp"
#include <flint/acb_dft.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/arb_poly.h>
#include <flint/fmpq.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <optional>
#include <stdexcept>

// The stages follow prototype v4. History: docs/code-history-20261021.md, real_line_seeds.cpp 1.

namespace certpolysolve::detail {
namespace {

using Clock = std::chrono::steady_clock;
uint64_t elapsed_ns(Clock::time_point start) {
    return static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count());
}

// A ball value is accepted at this many correct bits: above binary64's 53, so that a binary64 iterate, a sign or a
// log2 read from it is right to its last bit (v3's rule).
constexpr slong correct_bits = 60;
// The first precision of an evaluation; it doubles until correct_bits are met. Below about four times correct_bits
// the cancellation of Horner's scheme at these degrees always forces a second try (v3). A sequence of evaluations
// (Laguerre's steps toward one end, the per-seed tests) starts each one at the precision the previous one needed:
// approaching a root, the cancellation only grows.
constexpr slong first_evaluation_bits = 256;
// Laguerre's iteration toward an end stops in two phases: when a step is below 2^-20 of the root bound (both ends,
// which gives the interval's width W), then when it is below 2^-35 of W. The interval is widened by 2^-30 of W
// (below); in the slowest linear convergence measured (the clustered family) an end that stopped below 2^-35 W is
// nearer its root than the widening. From outside, Laguerre's iterates on a polynomial with only real roots never pass
// the extreme root, so the interval contains every root wherever the iteration stops. History:
// docs/code-history-20261021.md, real_line_seeds.cpp 2.
constexpr int laguerre_rough_bits = 20;
constexpr int laguerre_end_bits = 35;
// Steps a side. From outside, Laguerre's method converges monotonically for a polynomial with only real roots;
// into a cluster at the end of the root range (the arcsine families) it converges linearly. 500 is about 2.5 times
// the most measured on the prototype's inputs; past it the pass gives up (no seeds), it does not fail. History:
// docs/code-history-20261021.md, real_line_seeds.cpp 3.
constexpr int laguerre_max_steps = 500;
// The interval is the two extreme roots widened by 2^-30 of their distance, so that they lie strictly inside it and
// the interpolant does not depend on how its ends round. Wider costs nodes: 1% wider and Chebyshev's own polynomial
// reaches 2^216 at the ends (docs/real-roots-real-line-seeds-note-20261018.md). History:
// docs/code-history-20261021.md, real_line_seeds.cpp 4.
constexpr double interval_widening = 0x1p-30;
// The pass ends before the change of basis when log2 |p| ranges over more than 4 bits a degree across the interval:
// the interpolant would need more than 1.3 (1 + 0.7 4) = 5 nodes a root (the node count below), and the roots spread
// over scales the nodes cannot follow. Measured on the screening (M6): the inputs whose seeds paid lie below 4 bits a
// degree (the Laguerre and Hermite families, whose weights e^(-x/2) and e^(-x^2/2) give their range), and those above
// it failed the cell check after paying the change of basis and the transform. History:
// docs/code-history-20261021.md, real_line_seeds.cpp 5.
constexpr double max_range_bits_per_degree = 4;
// Ball samples of log2 |p| along the interval before the change of basis: they set the precisions and the node
// count (the range is measured again at every node afterwards).
constexpr int range_samples = 33;
// Precision margin on top of the range and the cancellation bound (v4).
constexpr slong guard_bits = 192;
// A node's ball radius must be below 2^-54 of its envelope: one bit beyond binary64, so the scaled value f_j is
// exact to its last bit.
constexpr slong value_bits = 54;
// The series of 2^-phi has about ln 2 times the range of phi (in bits) terms (the Bessel coefficients of
// e^(a cos theta) are significant up to k about a), so the interpolant needs about n + 0.7 range nodes; N is the
// power of two at or above 1.3 times that (v4, measured against v3's observed N). History:
// docs/code-history-20261021.md, real_line_seeds.cpp 6.
constexpr double bessel_factor = 0.7;
constexpr double nodes_margin = 1.3;
// N is doubled until the interpolant's tail is met, at most until 16 times the degree.
constexpr slong max_nodes_per_root = 16;
// The tail: the last tenth of the Chebyshev coefficients below 2^-45 of the largest (v3).
constexpr int tail_bits = 45;
constexpr double tail_start = 0.9;
// The envelope: a running median of log2 |p| over a half-window of 2 mean node spacings between roots (its median
// is the envelope minus 1 bit whatever the local density), smoothed by a Gaussian of 8 such spacings in the
// Chebyshev angle (the measured balance: a narrower one can double N, a wider one leave dense regions below it),
// cut where the multiplier is below 2^-60. A zero value counts as a dip 64 bits below the smallest. History:
// docs/code-history-20261021.md, real_line_seeds.cpp 7.
constexpr double median_spacings = 2;
constexpr double sigma_spacings = 8;
constexpr int truncate_bits = 60;
constexpr double zero_value_dip_bits = 64;
// phi is evaluated at the nodes at 128 bits: a binary64 phi of size 10^4 carries 2^-40 of non-smooth noise, which
// would show in the interpolant's tail (v3's trap).
constexpr slong scale_bits = 128;
// The roots of the interpolant: sign changes on a grid 16 times finer than the nodes (fine enough for jacobi's
// close roots), refined by 8 Illinois steps, which are superlinear on a cell 1/16 of a node spacing. History:
// docs/code-history-20261021.md, real_line_seeds.cpp 8.
constexpr slong grid_per_node = 16;
constexpr int illinois_steps = 8;
// The interpolant is evaluated with its coefficients above 2^-53 of the largest: the rest is noise.
constexpr int series_bits = 53;
// A sign change counts only where one side is above 2^-46 of the largest coefficient: below it the series is
// noise (where roots are denser than the nodes, the nodes' certified signs give the seeds instead).
constexpr int noise_bits = 46;
// Near-touching minima: a local minimum of |f| without a sign change, below 2^-6 of the local amplitude (blocks of
// 4 grid cells per node, with their neighbours) where that amplitude is above the noise by 8 bits, is refined on
// 256 points over its two cells; no sign change there and below 2^-20 of the amplitude is a touching pair (one seed).
constexpr slong amplitude_block_cells = 4 * grid_per_node;
constexpr int dip_bits = 6;
constexpr int amplitude_margin_bits = 8;
constexpr slong refine_points = 256;
constexpr int touch_bits = 20;
// A seed's own bracket for its test: s (1 +- 2^-10), the oracle measurement's criterion. History:
// docs/code-history-20261021.md, real_line_seeds.cpp 9.
constexpr int seed_half_width_bits = 10;
// The memory model of the pass's peak, in units of one vector of balls: L (P/8 + 64) bytes for the change of basis
// (L coefficients of P bits, about 64 bytes of struct and radius each) and 2N (P/8 + 64) for the ball transform.
// 40 and 12 bound every peak measured on the driver's inputs, of the change of basis (FLINT's block multiplication
// holds both operands and the double-width product as integers, besides the two levels of blocks and the powers of
// w) and of the transform (input, output, FLINT's roots of unity and its work vectors, each of two components).
// History: docs/code-history-20261021.md, real_line_seeds.cpp 10.
constexpr uint64_t ball_overhead_bytes = 64;
constexpr uint64_t conversion_vectors = 40;
constexpr uint64_t transform_vectors = 12;

// A binary64 value of a ball's midpoint, nearest.
double mid_double(arb_srcptr x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }

// log2 of the magnitude of a ball's midpoint, -inf at zero: the exponent exactly, the mantissa's log2 from arb.
double log2_abs_mid(arb_srcptr x) {
    if (arf_is_zero(arb_midref(x))) return -std::numeric_limits<double>::infinity();
    Float mantissa;
    Integer exponent;
    arf_frexp(mantissa.value(), exponent.value(), arb_midref(x));
    RealBall t, ln2;
    arb_set_arf(t.value(), mantissa.value());
    arb_abs(t.value(), t.value());
    arb_log(t.value(), t.value(), 64);
    arb_const_log2(ln2.value(), 64);
    arb_div(t.value(), t.value(), ln2.value(), 64);
    return fmpz_get_d(exponent.value()) + mid_double(t.value());
}

// log2 of a ball's radius, -inf for an exact ball, computed as log2_abs_mid is (FLINT's own approximate log2 of a
// radius may go through libm, whose bits can depend on the CPU, and two decisions read this value).
double log2_radius(arb_srcptr x) {
    if (mag_is_zero(arb_radref(x))) return -std::numeric_limits<double>::infinity();
    RealBall r;
    arf_set_mag(arb_midref(r.value()), arb_radref(x));
    return log2_abs_mid(r.value());
}

double log2_double(double x) {
    RealBall t;
    arb_set_d(t.value(), x);
    return log2_abs_mid(t.value());
}

double arb_exp_double(double x) {
    RealBall t;
    arb_set_d(t.value(), x);
    arb_exp(t.value(), t.value(), 64);
    return mid_double(t.value());
}

double arb_pi_double() {
    RealBall t;
    arb_const_pi(t.value(), 64);
    return mid_double(t.value());
}

bool is_power_of_two(slong n) { return n > 0 && (n & (n - 1)) == 0; }
slong next_power_of_two(double x) {
    slong n = 1;
    while (static_cast<double>(n) < x) n *= 2;
    return n;
}

struct Pass {
    const fmpz_poly_struct *p;
    const RealLineSeedLimits &limits;
    RealLineSeedStatistics &stats;
    slong n;
    IntegerPolynomial d1, d2; // p', p''

    Pass(const fmpz_poly_struct *poly, const RealLineSeedLimits &l, RealLineSeedStatistics &s)
        : p(poly), limits(l), stats(s), n(fmpz_poly_degree(poly)) {
        fmpz_poly_derivative(d1.value(), p);
        fmpz_poly_derivative(d2.value(), d1.value());
    }

    bool cancelled() const { return limits.cancelled && limits.cancelled(); }

    // f(x) as a ball with correct_bits relative accuracy (or an exact zero), else false: the precision would pass
    // the limit.
    bool evaluate(arb_t value, const fmpz_poly_struct *f, double x, slong &start) {
        RealBall point;
        arb_set_d(point.value(), x);
        for (slong prec = start; prec <= limits.max_precision; prec *= 2) {
            start = prec;
            arb_fmpz_poly_evaluate_arb(value, f, point.value(), prec);
            ++stats.ball_evaluations;
            stats.largest_evaluation_precision = std::max(stats.largest_evaluation_precision, prec);
            if (arb_is_zero(value) || arb_rel_accuracy_bits(value) >= correct_bits) return true;
        }
        return false;
    }

    // The sign of p(x), certain, or 0 when it cannot be decided within the limit (or p(x) = 0 exactly).
    int certain_sign(double x, slong &start) {
        RealBall point, value;
        arb_set_d(point.value(), x);
        for (slong prec = start; prec <= limits.max_precision; prec *= 2) {
            start = prec;
            arb_fmpz_poly_evaluate_arb(value.value(), p, point.value(), prec);
            ++stats.ball_evaluations;
            if (arb_is_positive(value.value())) return 1;
            if (arb_is_negative(value.value())) return -1;
            if (arb_is_zero(value.value())) return 0;
        }
        return 0;
    }

    // Laguerre's method from x toward the nearest root. For a polynomial with only real roots, started outside them,
    // it converges monotonically to the extreme root on that side. Its discriminant is (n - 1) (n H - G^2), with
    // G = p'/p and H = G^2 - p''/p, and n H - G^2 = ((n - 1) p'^2 - n p p'') / p^2 >= 0 at every real point for such
    // a polynomial (Cauchy-Schwarz on sum 1/(x - x_i)): a certainly negative discriminant proves a non-real root.
    // Steps until one is below 2^-bits of `scale`; `start` holds the three evaluations' warm precisions.
    RealLineSeedStop laguerre(double &x, double scale, int bits, slong (&start)[3]) {
        RealBall pv, qv, rv, g, h, t, disc, den;
        for (int step = 0; step < laguerre_max_steps; ++step) {
            if (cancelled()) return RealLineSeedStop::Cancelled;
            ++stats.laguerre_steps;
            if (!evaluate(pv.value(), p, x, start[0]) || !evaluate(qv.value(), d1.value(), x, start[1]) ||
                !evaluate(rv.value(), d2.value(), x, start[2]))
                return RealLineSeedStop::Precision;
            if (arb_is_zero(pv.value())) return RealLineSeedStop::Seeds; // x is a root exactly
            constexpr slong step_bits = 2 * correct_bits + 8; // the step's own arithmetic: its inputs carry 60 bits
            arb_div(g.value(), qv.value(), pv.value(), step_bits);
            arb_div(t.value(), rv.value(), pv.value(), step_bits);
            arb_sqr(h.value(), g.value(), step_bits);
            arb_sub(h.value(), h.value(), t.value(), step_bits);       // H = G^2 - p''/p
            arb_mul_si(disc.value(), h.value(), n, step_bits);
            arb_sqr(t.value(), g.value(), step_bits);
            arb_sub(disc.value(), disc.value(), t.value(), step_bits); // n H - G^2
            if (arb_is_negative(disc.value())) return RealLineSeedStop::LaguerreInequality;
            arb_mul_si(disc.value(), disc.value(), n - 1, step_bits);
            if (arb_is_positive(disc.value()))
                arb_sqrt(disc.value(), disc.value(), step_bits);
            else
                arb_zero(disc.value());
            if (arb_is_positive(g.value()))
                arb_add(den.value(), g.value(), disc.value(), step_bits);
            else
                arb_sub(den.value(), g.value(), disc.value(), step_bits);
            arb_set_si(t.value(), n);
            arb_div(t.value(), t.value(), den.value(), step_bits);
            const double delta = mid_double(t.value());
            if (!std::isfinite(delta)) return RealLineSeedStop::Numerical;
            x -= delta;
            if (!std::isfinite(x)) return RealLineSeedStop::Numerical;
            if (std::fabs(delta) <= std::ldexp(scale, -bits)) return RealLineSeedStop::Seeds;
        }
        return RealLineSeedStop::LaguerreStalled;
    }

    // A Fujiwara-type bound on the roots' moduli, rounded up to a power of two: 2 max_k |a_{n-k}/a_n|^(1/k) <=
    // 2 max_k 2^((bits(a_{n-k}) - bits(a_n) + 1) / k).
    double root_bound() const {
        const slong lead = static_cast<slong>(fmpz_bits(p->coeffs + n));
        slong exponent = std::numeric_limits<slong>::min();
        for (slong k = 1; k <= n; ++k) {
            const fmpz *c = p->coeffs + n - k;
            if (fmpz_is_zero(c)) continue;
            const slong numerator = static_cast<slong>(fmpz_bits(c)) - lead + 1;
            // ceil(numerator / k) for either sign
            const slong e = numerator >= 0 ? (numerator + k - 1) / k : -((-numerator) / k);
            exponent = std::max(exponent, e);
        }
        return std::ldexp(2.0, static_cast<int>(exponent));
    }

    // The Chebyshev coefficients c_0..c_n of p on [centre - half, centre + half] at `prec` bits, from
    //     z^n p(centre + half (z + 1/z) / 2) = sum_k a_k z^(n-k) w(z)^k,   w(z) = (half/2) z^2 + centre z + half/2,
    // whose coefficient of z^(n+k) is c_k / 2 (c_0 for k = 0): a divide-and-conquer sum of polynomial products,
    // H(i, 2m) = H(i, m) z^m + w^m H(i + m, m) with H(i, 1) = a_i, so that H(0, L) = z^(L-1-n) (that sum). With
    // `absolute`, the same sum of |a_k| and |w|'s coefficients: no cancellation, so its largest coefficient bounds
    // the terms the true sum adds. Empty when cancelled between levels.
    std::vector<RealBall> chebyshev_from_monomials(double centre, double half, slong prec, bool absolute) {
        const slong length = next_power_of_two(static_cast<double>(n + 1));
        std::vector<RealPolynomial> blocks(static_cast<size_t>(length));
        RealBall c;
        for (slong i = 0; i <= n; ++i) {
            arb_set_fmpz(c.value(), p->coeffs + i);
            if (absolute) arb_abs(c.value(), c.value());
            arb_poly_set_coeff_arb(blocks[static_cast<size_t>(i)].value(), 0, c.value());
        }
        RealPolynomial power; // w^m
        arb_set_d(c.value(), half / 2);
        arb_poly_set_coeff_arb(power.value(), 0, c.value());
        arb_poly_set_coeff_arb(power.value(), 2, c.value());
        arb_set_d(c.value(), absolute ? std::fabs(centre) : centre);
        arb_poly_set_coeff_arb(power.value(), 1, c.value());
        RealPolynomial product, shifted;
        for (slong m = 1; m < length; m *= 2) {
            if (cancelled()) return {};
            std::vector<RealPolynomial> next(static_cast<size_t>(length / (2 * m)));
            for (slong i = 0, b = 0; i < length; i += 2 * m, ++b) {
                arb_poly_mul(product.value(), power.value(), blocks[static_cast<size_t>(i / m + 1)].value(), prec);
                arb_poly_shift_left(shifted.value(), blocks[static_cast<size_t>(i / m)].value(), m);
                arb_poly_add(next[static_cast<size_t>(b)].value(), shifted.value(), product.value(), prec);
            }
            blocks = std::move(next);
            if (2 * m < length) arb_poly_mul(power.value(), power.value(), power.value(), prec);
        }
        std::vector<RealBall> out(static_cast<size_t>(n + 1));
        const auto *top = blocks[0].value();
        arb_poly_get_coeff_arb(out[0].value(), top, length - 1);
        for (slong k = 1; k <= n; ++k) {
            arb_poly_get_coeff_arb(out[static_cast<size_t>(k)].value(), top, length - 1 + k);
            arb_mul_2exp_si(out[static_cast<size_t>(k)].value(), out[static_cast<size_t>(k)].value(), 1);
        }
        return out;
    }
};

// v_j = sum_{k=0}^{K} c_k cos(pi j k / N), j = 0..N (K < N), in ball arithmetic at `prec`: the DFT of length 2N of the
// even extension e (e_k = e_{2N-k} = c_k) is X_j = c_0 + 2 sum_{k>=1} c_k cos(pi j k / N), so v_j = (X_j + c_0) / 2.
// Each c_k is rounded to prec first (its radius kept).
std::vector<RealBall> ball_cosine_transform(const std::vector<RealBall> &c, slong nodes, slong prec) {
    const slong length = 2 * nodes;
    RootVector e(length), x(length);
    const slong K = static_cast<slong>(c.size()) - 1;
    if (K >= nodes) throw std::logic_error("ball_cosine_transform: more coefficients than nodes");
    for (slong k = 0; k <= K; ++k) {
        arb_set_round(acb_realref(e.data() + k), c[static_cast<size_t>(k)].value(), prec);
        if (k > 0) arb_set(acb_realref(e.data() + length - k), acb_realref(e.data() + k));
    }
    acb_dft(x.data(), e.data(), length, prec);
    std::vector<RealBall> v(static_cast<size_t>(nodes + 1));
    RealBall c0;
    arb_set_round(c0.value(), c[0].value(), prec);
    for (slong j = 0; j <= nodes; ++j) {
        arb_add(v[static_cast<size_t>(j)].value(), acb_realref(x.data() + j), c0.value(), prec);
        arb_mul_2exp_si(v[static_cast<size_t>(j)].value(), v[static_cast<size_t>(j)].value(), -1);
    }
    return v;
}

// The running median of y over windows of 2 w + 1 values, the ends padded with the end values. Sorting exact
// binary64 values: the same bits everywhere.
std::vector<double> running_median(const std::vector<double> &y, slong w) {
    const slong count = static_cast<slong>(y.size());
    std::vector<double> out(y.size()), window(static_cast<size_t>(2 * w + 1));
    for (slong i = 0; i < count; ++i) {
        for (slong d = -w; d <= w; ++d) {
            const slong j = std::clamp<slong>(i + d, 0, count - 1);
            window[static_cast<size_t>(d + w)] = y[static_cast<size_t>(j)];
        }
        std::nth_element(window.begin(), window.begin() + w, window.end());
        out[static_cast<size_t>(i)] = window[static_cast<size_t>(w)];
    }
    return out;
}

} // namespace

const char *real_line_seed_stop_name(RealLineSeedStop stop) {
    switch (stop) {
    case RealLineSeedStop::Seeds: return "seeds";
    case RealLineSeedStop::Degree: return "degree";
    case RealLineSeedStop::NewtonInequality: return "Newton's inequality";
    case RealLineSeedStop::LaguerreInequality: return "Laguerre's inequality";
    case RealLineSeedStop::LaguerreStalled: return "Laguerre's iteration did not converge";
    case RealLineSeedStop::FewCells: return "too few sign changes";
    case RealLineSeedStop::NoisyInterpolant: return "a noisy interpolant";
    case RealLineSeedStop::WideRange: return "too wide a range";
    case RealLineSeedStop::Memory: return "memory";
    case RealLineSeedStop::Precision: return "precision";
    case RealLineSeedStop::Cancelled: return "cancelled";
    case RealLineSeedStop::Numerical: return "numerical";
    }
    return "unknown";
}

namespace {
enum class NewtonEstimate { Holds, Fails, Close };

// Newton's inequality at k, a_k^2 k (n - k) >= a_(k-1) a_(k+1) (k + 1) (n - k + 1), with a_(k-1) a_(k+1) > 0, from
// each coefficient's leading bits: |a| = d 2^e with d in [1/2, 1) to within one ulp (fmpz_get_d_2exp). With
// E = 2 e_k - e_(k-1) - e_(k+1) the two sides compare as d_k^2 k (n - k) 2^E against d_(k-1) d_(k+1) (k + 1)
// (n - k + 1). Beyond newton_estimate_far_bits the exponent decides alone (each other factor lies between 1/4 and
// n^2). Otherwise both sides are formed by IEEE products and an exact power-of-two scaling (no libm transcendental,
// as the pass's rule asks), whose rounding stays far below newton_estimate_margin; a relative difference within it
// is left to the exact products (Close), so the decision is the exact comparison's on every CPU.
constexpr double newton_estimate_margin = 0x1p-40;
constexpr slong newton_estimate_far_bits = 96;

NewtonEstimate newton_estimate(const fmpz *before, const fmpz *at, const fmpz *after, slong k, slong n) {
    slong e_before, e_at, e_after;
    const double d_before = std::fabs(fmpz_get_d_2exp(&e_before, before));
    const double d_at = std::fabs(fmpz_get_d_2exp(&e_at, at));
    const double d_after = std::fabs(fmpz_get_d_2exp(&e_after, after));
    const slong exponent = 2 * e_at - e_before - e_after;
    if (exponent > newton_estimate_far_bits) return NewtonEstimate::Holds;
    if (exponent < -newton_estimate_far_bits) return NewtonEstimate::Fails;
    const double left = std::ldexp(d_at * d_at * (static_cast<double>(k) * static_cast<double>(n - k)),
                                   static_cast<int>(exponent));
    const double right = d_before * d_after * (static_cast<double>(k + 1) * static_cast<double>(n - k + 1));
    if (left > right * (1 + newton_estimate_margin)) return NewtonEstimate::Holds;
    if (left < right * (1 - newton_estimate_margin)) return NewtonEstimate::Fails;
    return NewtonEstimate::Close;
}
} // namespace

slong newton_inequality_failure(const fmpz_poly_struct *p) {
    const slong n = fmpz_poly_degree(p);
    Integer left, right;
    for (slong k = 1; k < n; ++k) {
        const fmpz *before = p->coeffs + k - 1, *at = p->coeffs + k, *after = p->coeffs + k + 1;
        if (fmpz_is_zero(at)) {
            // A zero coefficient of a real-rooted polynomial with p(0) != 0 sits between coefficients of opposite signs.
            if (fmpz_sgn(before) * fmpz_sgn(after) >= 0) return k;
            continue;
        }
        // a_(k-1) a_(k+1) <= 0: the right side is at most 0 and the left side positive, so the inequality holds.
        if (fmpz_sgn(before) * fmpz_sgn(after) <= 0) continue;
        // Both sides positive: the leading bits decide, and the exact products only where the two are close (the
        // products of large coefficients cost more than the rest of a refused pass; dev/docs/slowdowns-20261029.md).
        const auto estimate = newton_estimate(before, at, after, k, n);
        if (estimate == NewtonEstimate::Holds) continue;
        if (estimate == NewtonEstimate::Fails) return k;
        fmpz_mul(left.value(), at, at);
        fmpz_mul_si(left.value(), left.value(), k * (n - k));
        fmpz_mul(right.value(), before, after);
        fmpz_mul_si(right.value(), right.value(), (k + 1) * (n - k + 1));
        if (fmpz_cmp(left.value(), right.value()) < 0) return k;
    }
    return -1;
}

Binary64Transform::Binary64Transform(slong max_length) : max_length_(max_length) {
    if (max_length < 2 || !is_power_of_two(max_length))
        throw std::invalid_argument("Binary64Transform: the length must be a power of two, at least 2");
    const slong half = max_length / 2;
    cos_.resize(static_cast<size_t>(half));
    sin_.resize(static_cast<size_t>(half));
    // cos and sin of 2 pi k / L for k < L/2, from arb: computed for k <= L/4 and mirrored (cos(pi - a) = -cos a,
    // sin(pi - a) = sin a), which halves the arb calls and gives the same bits as computing each.
    RealBall s, c;
    const slong quarter = max_length / 4;
    for (slong k = 0; k <= std::min(quarter, half - 1); ++k) {
        Rational angle; // 2 pi k / L as pi times this
        fmpq_set_si(angle.value(), 2 * k, static_cast<ulong>(max_length));
        arb_sin_cos_pi_fmpq(s.value(), c.value(), angle.value(), 64);
        cos_[static_cast<size_t>(k)] = mid_double(c.value());
        sin_[static_cast<size_t>(k)] = mid_double(s.value());
    }
    for (slong k = quarter + 1; k < half; ++k) {
        cos_[static_cast<size_t>(k)] = -cos_[static_cast<size_t>(half - k)];
        sin_[static_cast<size_t>(k)] = sin_[static_cast<size_t>(half - k)];
    }
}

void Binary64Transform::forward(std::vector<double> &re, std::vector<double> &im) const {
    const slong length = static_cast<slong>(re.size());
    if (!is_power_of_two(length) || length > max_length_ || im.size() != re.size())
        throw std::invalid_argument("Binary64Transform::forward: bad length");
    // Bit reversal, then iterative radix-2 butterflies (decimation in time), each butterfly in a fixed order.
    for (slong i = 1, j = 0; i < length; ++i) {
        slong bit = length >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            std::swap(re[static_cast<size_t>(i)], re[static_cast<size_t>(j)]);
            std::swap(im[static_cast<size_t>(i)], im[static_cast<size_t>(j)]);
        }
    }
    for (slong len = 2; len <= length; len *= 2) {
        const slong half = len / 2, stride = max_length_ / len;
        for (slong start = 0; start < length; start += len) {
            for (slong k = 0; k < half; ++k) {
                const double wr = cos_[static_cast<size_t>(k * stride)];
                const double wi = -sin_[static_cast<size_t>(k * stride)];
                const size_t a = static_cast<size_t>(start + k), b = static_cast<size_t>(start + k + half);
                const double vr = re[b] * wr - im[b] * wi;
                const double vi = re[b] * wi + im[b] * wr;
                re[b] = re[a] - vr;
                im[b] = im[a] - vi;
                re[a] = re[a] + vr;
                im[a] = im[a] + vi;
            }
        }
    }
}

std::vector<double> Binary64Transform::dct1(const std::vector<double> &x) const {
    const slong m = static_cast<slong>(x.size()) - 1;
    if (!is_power_of_two(m) || 2 * m > max_length_) throw std::invalid_argument("Binary64Transform::dct1: bad length");
    std::vector<double> re(static_cast<size_t>(2 * m)), im(static_cast<size_t>(2 * m), 0.0);
    for (slong k = 0; k <= m; ++k) re[static_cast<size_t>(k)] = x[static_cast<size_t>(k)];
    for (slong k = 1; k < m; ++k) re[static_cast<size_t>(2 * m - k)] = x[static_cast<size_t>(k)];
    forward(re, im);
    re.resize(static_cast<size_t>(m + 1));
    return re;
}

double Binary64Transform::cos_pi(slong k, slong m) const {
    // cos(pi k / m) = cos(2 pi (k stride) / L) with stride = L / (2 m); for k > m/2... k up to m is k stride up to L/2,
    // one past the table: cos(pi) = -1.
    const slong stride = max_length_ / (2 * m);
    const slong index = k * stride;
    if (index == max_length_ / 2) return -1.0;
    return cos_[static_cast<size_t>(index)];
}

std::vector<double> chebyshev_coefficients(const Binary64Transform &t, const std::vector<double> &values) {
    const slong m = static_cast<slong>(values.size()) - 1;
    std::vector<double> c = t.dct1(values);
    for (auto &v : c) v = v / static_cast<double>(m);
    c[0] = c[0] / 2;
    c[static_cast<size_t>(m)] = c[static_cast<size_t>(m)] / 2;
    return c;
}

std::vector<double> chebyshev_values(const Binary64Transform &t, const std::vector<double> &coefficients, slong m) {
    std::vector<double> x(static_cast<size_t>(m + 1), 0.0);
    const slong count = std::min<slong>(static_cast<slong>(coefficients.size()), m + 1);
    for (slong k = 0; k < count; ++k) x[static_cast<size_t>(k)] = coefficients[static_cast<size_t>(k)];
    for (slong k = 1; k < m; ++k) x[static_cast<size_t>(k)] = x[static_cast<size_t>(k)] / 2;
    return t.dct1(x);
}

double chebyshev_evaluate(const std::vector<double> &c, slong length, double x) {
    // Clenshaw: b_k = c_k + 2 x b_{k+1} - b_{k+2}; the sum is c_0 + x b_1 - b_2.
    double b1 = 0, b2 = 0;
    const double twice = 2 * x;
    for (slong k = length - 1; k >= 1; --k) {
        const double b0 = c[static_cast<size_t>(k)] + twice * b1 - b2;
        b2 = b1;
        b1 = b0;
    }
    return c[0] + x * b1 - b2;
}

RealLineSeeds real_line_seeds(const fmpz_poly_struct *p, const RealLineSeedLimits &limits) {
    // A root at zero (a factor seeded in centred coordinates whose centre is a root, say) is known exactly: the pass
    // runs on p / x^v and the point 0 seeds it (once: equal start points are no use to the iteration).
    if (fmpz_poly_degree(p) >= 1 && fmpz_is_zero(p->coeffs)) {
        slong v = 0;
        while (v < fmpz_poly_length(p) && fmpz_is_zero(p->coeffs + v)) ++v;
        IntegerPolynomial divided;
        fmpz_poly_shift_right(divided.value(), p, v);
        RealLineSeeds result = real_line_seeds(divided.value(), limits);
        result.statistics.degree = fmpz_poly_degree(p);
        if (result.statistics.stop == RealLineSeedStop::Seeds) {
            result.seeds.insert(std::upper_bound(result.seeds.begin(), result.seeds.end(), 0.0), 0.0);
            result.seeds.erase(std::unique(result.seeds.begin(), result.seeds.end()), result.seeds.end());
        }
        return result;
    }
    RealLineSeeds result;
    auto &stats = result.statistics;
    const auto started = Clock::now();
    const slong n = fmpz_poly_degree(p);
    stats.degree = n;
    auto stop = [&](RealLineSeedStop why) {
        stats.stop = why;
        stats.total_ns = elapsed_ns(started);
        result.seeds.clear();
        return result;
    };
    if (n < 2 || fmpz_is_zero(p->coeffs)) return stop(RealLineSeedStop::Degree);
    stats.newton_failure = newton_inequality_failure(p);
    if (stats.newton_failure >= 0) return stop(RealLineSeedStop::NewtonInequality);
    Pass pass(p, limits, stats);
    if (pass.cancelled()) return stop(RealLineSeedStop::Cancelled);

    // 1. The interval: the extreme roots by Laguerre's method from outside the root bound.
    auto t0 = Clock::now();
    const double bound = pass.root_bound();
    double hi = bound, lo = -bound;
    slong warm_hi[3] = {first_evaluation_bits, first_evaluation_bits, first_evaluation_bits};
    slong warm_lo[3] = {first_evaluation_bits, first_evaluation_bits, first_evaluation_bits};
    if (auto why = pass.laguerre(hi, bound, laguerre_rough_bits, warm_hi); why != RealLineSeedStop::Seeds) return stop(why);
    if (auto why = pass.laguerre(lo, bound, laguerre_rough_bits, warm_lo); why != RealLineSeedStop::Seeds) return stop(why);
    if (!(hi > lo)) return stop(RealLineSeedStop::Numerical);
    const double rough_width = hi - lo;
    if (auto why = pass.laguerre(hi, rough_width, laguerre_end_bits, warm_hi); why != RealLineSeedStop::Seeds) return stop(why);
    if (auto why = pass.laguerre(lo, rough_width, laguerre_end_bits, warm_lo); why != RealLineSeedStop::Seeds) return stop(why);
    stats.lower_root = lo;
    stats.upper_root = hi;
    stats.interval_ns = elapsed_ns(t0);
    const double width = hi - lo;
    if (!(width > 0)) return stop(RealLineSeedStop::Numerical);
    const double a = lo - width * interval_widening, b = hi + width * interval_widening;
    const double centre = (a + b) / 2, half = (b - a) / 2;

    // 2. The range of log2 |p| over the interval, from ball samples at the Chebyshev points of the second kind.
    t0 = Clock::now();
    const slong sample_intervals = range_samples - 1;
    // One binary64 transform, rebuilt only when a stage needs a longer one (its twiddles serve every shorter length
    // and every cosine of the nodes and the grid).
    std::optional<Binary64Transform> transform_storage;
    auto transform_of = [&](slong length) -> const Binary64Transform & {
        if (!transform_storage || transform_storage->max_length() < length) transform_storage.emplace(length);
        return *transform_storage;
    };
    double sample_max = -std::numeric_limits<double>::infinity(), sample_min = std::numeric_limits<double>::infinity();
    {
        RealBall value;
        slong warm = first_evaluation_bits;
        for (slong j = 0; j <= sample_intervals; ++j) {
            if (pass.cancelled()) return stop(RealLineSeedStop::Cancelled);
            const double x = centre + half * transform_of(2 * sample_intervals).cos_pi(j, sample_intervals);
            if (!pass.evaluate(value.value(), p, x, warm)) return stop(RealLineSeedStop::Precision);
            const double l = log2_abs_mid(value.value());
            if (!std::isfinite(l)) continue;
            sample_max = std::max(sample_max, l);
            sample_min = std::min(sample_min, l);
        }
    }
    if (!std::isfinite(sample_max)) return stop(RealLineSeedStop::Numerical);
    const double range = sample_max - sample_min;
    stats.range_bits = range;
    stats.range_ns = elapsed_ns(t0);
    if (range > max_range_bits_per_degree * static_cast<double>(n)) return stop(RealLineSeedStop::WideRange);
    slong nodes = next_power_of_two(nodes_margin * (static_cast<double>(n) + bessel_factor * range));
    slong transform_bits = static_cast<slong>(std::ceil(range)) + guard_bits;

    // 3. The Chebyshev coefficients by the change of basis, at the precision its cancellation bound asks.
    t0 = Clock::now();
    std::vector<RealBall> chebyshev;
    {
        const auto absolute = pass.chebyshev_from_monomials(centre, half, 64, true);
        if (absolute.empty()) return stop(RealLineSeedStop::Cancelled);
        double term_bound = -std::numeric_limits<double>::infinity();
        for (const auto &c : absolute) term_bound = std::max(term_bound, log2_abs_mid(c.value()));
        const double loss_bound = std::max(0.0, term_bound - sample_max);
        stats.cancellation_bits = loss_bound;
        slong conversion_bits = static_cast<slong>(std::ceil(loss_bound + range)) + value_bits + guard_bits;
        const slong length = next_power_of_two(static_cast<double>(n + 1));
        const uint64_t conversion_bytes =
            conversion_vectors * static_cast<uint64_t>(length) * (static_cast<uint64_t>(conversion_bits) / 8 + ball_overhead_bytes);
        const uint64_t transform_bytes =
            transform_vectors * static_cast<uint64_t>(2 * nodes) * (static_cast<uint64_t>(transform_bits) / 8 + ball_overhead_bytes);
        stats.predicted_bytes = std::max(conversion_bytes, transform_bytes);
        if (stats.predicted_bytes > limits.max_bytes) return stop(RealLineSeedStop::Memory);
        for (;;) {
            if (conversion_bits > limits.max_precision) return stop(RealLineSeedStop::Precision);
            stats.conversion_precision = conversion_bits;
            chebyshev = pass.chebyshev_from_monomials(centre, half, conversion_bits, false);
            if (chebyshev.empty()) return stop(RealLineSeedStop::Cancelled);
            double largest = -std::numeric_limits<double>::infinity(), widest = -std::numeric_limits<double>::infinity();
            for (const auto &c : chebyshev) {
                largest = std::max(largest, log2_abs_mid(c.value()));
                widest = std::max(widest, log2_radius(c.value()));
            }
            // Enough when the coefficients' radii sit the range plus value_bits plus 64 bits below the largest.
            if (largest - widest >= range + static_cast<double>(value_bits) + 64) break;
            const double lost = static_cast<double>(conversion_bits) - (largest - widest);
            conversion_bits = static_cast<slong>(std::ceil(lost + range)) + value_bits + guard_bits;
        }
    }
    stats.conversion_ns = elapsed_ns(t0);

    // 4. The values at the N + 1 nodes cos(pi j / N) by the ball transform; their certified signs are counted at once.
    std::vector<RealBall> values;
    auto transform_values = [&]() -> bool {
        const auto t1 = Clock::now();
        values = ball_cosine_transform(chebyshev, nodes, transform_bits);
        stats.transform_ns += elapsed_ns(t1);
        return !pass.cancelled();
    };
    if (transform_bits > limits.max_precision) return stop(RealLineSeedStop::Precision);
    if (!transform_values()) return stop(RealLineSeedStop::Cancelled);
    auto sign_changes = [&]() {
        slong changes = 0, previous = 0;
        for (const auto &v : values) {
            const slong s = arb_is_positive(v.value()) ? 1 : (arb_is_negative(v.value()) ? -1 : 0);
            if (s != 0 && previous != 0 && s != previous) ++changes;
            if (s != 0) previous = s;
        }
        return changes;
    };
    if (limits.min_cell_fraction > 0 &&
        static_cast<double>(sign_changes()) < limits.min_cell_fraction * static_cast<double>(n))
        return stop(RealLineSeedStop::FewCells);

    // 5. The envelope phi of log2 |p| at the nodes (binary64), as Chebyshev coefficients.
    t0 = Clock::now();
    std::vector<double> envelope;
    {
        const slong count = nodes + 1;
        std::vector<double> logs(static_cast<size_t>(count));
        double smallest = std::numeric_limits<double>::infinity();
        for (slong j = 0; j < count; ++j) {
            logs[static_cast<size_t>(j)] = log2_abs_mid(values[static_cast<size_t>(j)].value());
            if (std::isfinite(logs[static_cast<size_t>(j)])) smallest = std::min(smallest, logs[static_cast<size_t>(j)]);
        }
        if (!std::isfinite(smallest)) return stop(RealLineSeedStop::Numerical);
        std::vector<double> y(logs);
        for (auto &v : y)
            if (!std::isfinite(v)) v = smallest - zero_value_dip_bits;
        const double spacing = static_cast<double>(nodes) / static_cast<double>(n);
        const slong w = std::max<slong>(2, static_cast<slong>(std::floor(median_spacings * spacing + 0.5)));
        const double sigma = sigma_spacings * spacing * arb_pi_double() / static_cast<double>(nodes);
        const auto median = running_median(y, w);
        const auto &tr = transform_of(2 * nodes);
        auto coefficients = chebyshev_coefficients(tr, median);
        slong kept = 0;
        for (slong k = 0; k < count; ++k) {
            const double multiplier = arb_exp_double(-0.5 * (sigma * static_cast<double>(k)) * (sigma * static_cast<double>(k)));
            if (multiplier > std::ldexp(1.0, -truncate_bits)) kept = k;
            coefficients[static_cast<size_t>(k)] = coefficients[static_cast<size_t>(k)] * multiplier;
            if (pass.cancelled()) return stop(RealLineSeedStop::Cancelled);
        }
        coefficients.resize(static_cast<size_t>(kept + 1));
        const auto phi = chebyshev_values(tr, coefficients, nodes);
        // Raised to the largest residual of the local maxima (the ends count when above their one neighbour).
        double raise = -std::numeric_limits<double>::infinity();
        for (slong j = 0; j < count; ++j) {
            const double l = logs[static_cast<size_t>(j)];
            if (!std::isfinite(l)) continue;
            const bool left = j == 0 || l >= logs[static_cast<size_t>(j - 1)];
            const bool right = j == count - 1 || l >= logs[static_cast<size_t>(j + 1)];
            if (left && right) raise = std::max(raise, l - phi[static_cast<size_t>(j)]);
        }
        if (!std::isfinite(raise)) return stop(RealLineSeedStop::Numerical);
        coefficients[0] = coefficients[0] + raise;
        stats.envelope_terms = kept + 1;
        stats.envelope_raise_bits = raise;
        envelope = std::move(coefficients);
    }
    stats.envelope_ns = elapsed_ns(t0);

    // 6. f_j = v_j 2^-phi_j in binary64, phi_j from the envelope's coefficients in balls at scale_bits (exact inputs);
    // the node balls' radii checked against the envelope, the transform redone wider if they are not small enough.
    std::vector<double> f;
    auto scale = [&]() -> RealLineSeedStop {
        const auto t1 = Clock::now();
        std::vector<RealBall> envelope_balls(envelope.size());
        for (size_t k = 0; k < envelope.size(); ++k) arb_set_d(envelope_balls[k].value(), envelope[k]);
        if (static_cast<slong>(envelope.size()) > nodes) envelope_balls.resize(static_cast<size_t>(nodes));
        const auto phi = ball_cosine_transform(envelope_balls, nodes, scale_bits);
        f.assign(static_cast<size_t>(nodes + 1), 0.0);
        double worst = -std::numeric_limits<double>::infinity();
        RealBall ln2, t;
        Integer whole;
        arb_const_log2(ln2.value(), scale_bits);
        for (slong j = 0; j <= nodes; ++j) {
            if ((j & 1023) == 0 && pass.cancelled()) return RealLineSeedStop::Cancelled;
            const auto *e = phi[static_cast<size_t>(j)].value();
            const auto *v = values[static_cast<size_t>(j)].value();
            // 2^-phi = 2^-k 2^(k - phi) with k = floor(phi): the exponential only of a number in (-1, 0].
            arf_get_fmpz(whole.value(), arb_midref(e), ARF_RND_FLOOR);
            arb_set_fmpz(t.value(), whole.value());
            arb_sub(t.value(), t.value(), e, scale_bits);
            arb_mul(t.value(), t.value(), ln2.value(), scale_bits);
            arb_exp(t.value(), t.value(), scale_bits);
            fmpz_neg(whole.value(), whole.value());
            arb_mul_2exp_fmpz(t.value(), t.value(), whole.value());
            RealBall scaled;
            arb_set_arf(scaled.value(), arb_midref(v));
            arb_mul(scaled.value(), scaled.value(), t.value(), scale_bits);
            f[static_cast<size_t>(j)] = mid_double(scaled.value());
            worst = std::max(worst, log2_radius(v) - mid_double(e));
        }
        stats.scale_ns += elapsed_ns(t1);
        if (worst <= -static_cast<double>(value_bits)) return RealLineSeedStop::Seeds;
        transform_bits += static_cast<slong>(std::ceil(worst)) + value_bits + 64;
        if (transform_bits > limits.max_precision) return RealLineSeedStop::Precision;
        if (!transform_values()) return RealLineSeedStop::Cancelled;
        return RealLineSeedStop::Numerical; // signals "redo"
    };
    for (;;) {
        const auto why = scale();
        if (why == RealLineSeedStop::Seeds) break;
        if (why != RealLineSeedStop::Numerical) return stop(why);
    }
    stats.transform_precision = transform_bits;

    // 7. The interpolant's coefficients; N doubled (a new ball transform, the change of basis reused) until the tail
    // is met.
    std::vector<double> c;
    for (;;) {
        c = chebyshev_coefficients(transform_of(2 * nodes), f);
        double largest = 0, tail = 0;
        for (slong k = 0; k <= nodes; ++k) {
            const double m = std::fabs(c[static_cast<size_t>(k)]);
            largest = std::max(largest, m);
            if (static_cast<double>(k) >= tail_start * static_cast<double>(nodes)) tail = std::max(tail, m);
        }
        if (!(largest > 0)) return stop(RealLineSeedStop::Numerical);
        stats.tail_bits = log2_double(tail / largest);
        if (tail <= std::ldexp(largest, -tail_bits) || nodes >= max_nodes_per_root * n) break;
        const uint64_t doubled_bytes = transform_vectors * static_cast<uint64_t>(4 * nodes) *
                                       (static_cast<uint64_t>(transform_bits) / 8 + ball_overhead_bytes);
        stats.predicted_bytes = std::max(stats.predicted_bytes, doubled_bytes);
        if (doubled_bytes > limits.max_bytes) return stop(RealLineSeedStop::Memory);
        nodes *= 2;
        if (!transform_values()) return stop(RealLineSeedStop::Cancelled);
        for (;;) {
            const auto why = scale();
            if (why == RealLineSeedStop::Seeds) break;
            if (why != RealLineSeedStop::Numerical) return stop(why);
        }
    }
    stats.nodes = nodes;

    // 8. The interpolant's roots: sign changes above the noise on a grid 16 times finer, refined by Illinois steps;
    // near-touching minima refined locally.
    t0 = Clock::now();
    std::vector<double> seeds_t;
    std::vector<int> kind; // 0 interpolant, 1 touching pair
    {
        double fmax = 0;
        for (const double v : c) fmax = std::max(fmax, std::fabs(v));
        slong significant = 0;
        for (slong k = 0; k <= nodes; ++k)
            if (std::fabs(c[static_cast<size_t>(k)]) > std::ldexp(fmax, -series_bits)) significant = k;
        const slong terms = significant + 1;
        const double noise = std::ldexp(fmax, -noise_bits);
        const slong m = grid_per_node * nodes;
        const auto &tr = transform_of(2 * m);
        auto grid = chebyshev_values(tr, c, m); // at cos(pi j / m), j = 0..m: decreasing t
        std::reverse(grid.begin(), grid.end());
        std::vector<double> t(static_cast<size_t>(m + 1));
        for (slong j = 0; j <= m; ++j) t[static_cast<size_t>(j)] = tr.cos_pi(m - j, m);
        auto illinois = [&](double left, double right) {
            double fl = chebyshev_evaluate(c, terms, left), fr = chebyshev_evaluate(c, terms, right);
            int side = 0;
            for (int s = 0; s < illinois_steps; ++s) {
                double mid = (left * fr - right * fl) / (fr - fl);
                if (!(std::isfinite(mid) && mid > left && mid < right)) mid = (left + right) / 2;
                const double fm = chebyshev_evaluate(c, terms, mid);
                const bool same = (fm > 0) == (fl > 0);
                if (same) {
                    if (side == 1) fr = fr / 2;
                    left = mid;
                    fl = fm;
                    side = 1;
                } else {
                    if (side == -1) fl = fl / 2;
                    right = mid;
                    fr = fm;
                    side = -1;
                }
            }
            return std::fabs(fl) < std::fabs(fr) ? left : right;
        };
        for (slong j = 0; j < m; ++j) {
            if (pass.cancelled()) return stop(RealLineSeedStop::Cancelled);
            const double g0 = grid[static_cast<size_t>(j)], g1 = grid[static_cast<size_t>(j + 1)];
            if (((g0 > 0 && g1 < 0) || (g0 < 0 && g1 > 0)) && std::max(std::fabs(g0), std::fabs(g1)) > noise) {
                seeds_t.push_back(illinois(t[static_cast<size_t>(j)], t[static_cast<size_t>(j + 1)]));
                kind.push_back(0);
                ++stats.interpolant_seeds;
                // p has at most n real roots, and f = p 2^-phi changes sign only at them: an interpolant that changes
                // sign more often is noise (roots far denser than the nodes, or values the scale did not resolve).
                if (stats.interpolant_seeds > n) return stop(RealLineSeedStop::NoisyInterpolant);
            }
        }
        // The local amplitude: block maxima of |grid| over blocks of amplitude_block_cells, each with its neighbours.
        const slong blocks = (m + 1 + amplitude_block_cells - 1) / amplitude_block_cells;
        std::vector<double> block(static_cast<size_t>(blocks), 0.0);
        for (slong j = 0; j <= m; ++j)
            block[static_cast<size_t>(j / amplitude_block_cells)] =
                std::max(block[static_cast<size_t>(j / amplitude_block_cells)], std::fabs(grid[static_cast<size_t>(j)]));
        auto local = [&](slong j) {
            const slong b0 = j / amplitude_block_cells;
            double l = block[static_cast<size_t>(b0)];
            if (b0 > 0) l = std::max(l, block[static_cast<size_t>(b0 - 1)]);
            if (b0 + 1 < blocks) l = std::max(l, block[static_cast<size_t>(b0 + 1)]);
            return l;
        };
        auto sgn = [](double v) { return v > 0 ? 1 : (v < 0 ? -1 : 0); };
        for (slong j = 1; j < m; ++j) {
            const double a0 = std::fabs(grid[static_cast<size_t>(j - 1)]), a1 = std::fabs(grid[static_cast<size_t>(j)]),
                         a2 = std::fabs(grid[static_cast<size_t>(j + 1)]);
            const int s0 = sgn(grid[static_cast<size_t>(j - 1)]), s1 = sgn(grid[static_cast<size_t>(j)]),
                      s2 = sgn(grid[static_cast<size_t>(j + 1)]);
            const double amplitude = local(j);
            if (!(a1 <= a0 && a1 <= a2 && s0 == s2 && s1 == s0 && a1 < std::ldexp(amplitude, -dip_bits) &&
                  amplitude > std::ldexp(fmax, -(noise_bits - amplitude_margin_bits))))
                continue;
            // refine over its two cells
            const double left = t[static_cast<size_t>(j - 1)], right = t[static_cast<size_t>(j + 1)];
            std::vector<double> ts(static_cast<size_t>(refine_points + 1)), fs(static_cast<size_t>(refine_points + 1));
            for (slong r = 0; r <= refine_points; ++r) {
                ts[static_cast<size_t>(r)] = left + (right - left) * (static_cast<double>(r) / static_cast<double>(refine_points));
                fs[static_cast<size_t>(r)] = chebyshev_evaluate(c, terms, ts[static_cast<size_t>(r)]);
            }
            bool changed = false;
            for (slong r = 0; r < refine_points; ++r) {
                const double f0 = fs[static_cast<size_t>(r)], f1 = fs[static_cast<size_t>(r + 1)];
                if (((f0 > 0 && f1 < 0) || (f0 < 0 && f1 > 0)) && std::max(std::fabs(f0), std::fabs(f1)) > noise) {
                    seeds_t.push_back(illinois(ts[static_cast<size_t>(r)], ts[static_cast<size_t>(r + 1)]));
                    kind.push_back(0);
                    ++stats.interpolant_seeds;
                    if (stats.interpolant_seeds > n) return stop(RealLineSeedStop::NoisyInterpolant);
                    changed = true;
                }
            }
            if (!changed) {
                slong best = 0;
                for (slong r = 1; r <= refine_points; ++r)
                    if (std::fabs(fs[static_cast<size_t>(r)]) < std::fabs(fs[static_cast<size_t>(best)])) best = r;
                if (std::fabs(fs[static_cast<size_t>(best)]) < std::ldexp(amplitude, -touch_bits)) {
                    seeds_t.push_back(ts[static_cast<size_t>(best)]);
                    kind.push_back(1);
                    ++stats.touching_seeds;
                }
            }
        }
    }
    stats.roots_ns = elapsed_ns(t0);

    // 9. The cells between the nodes (points with a certain sign, increasing x).
    t0 = Clock::now();
    std::vector<double> seeds(seeds_t.size());
    for (size_t i = 0; i < seeds_t.size(); ++i) seeds[i] = centre + half * seeds_t[i];
    std::vector<double> points(static_cast<size_t>(nodes + 1));
    std::vector<int> signs(static_cast<size_t>(nodes + 1));
    std::vector<const RealBall *> point_values(static_cast<size_t>(nodes + 1));
    const auto &cells_transform = transform_of(2 * nodes);
    for (slong j = 0; j <= nodes; ++j) {
        // node j is cos(pi j / N): increasing x is j = N..0
        const slong i = nodes - j;
        points[static_cast<size_t>(i)] = centre + half * cells_transform.cos_pi(j, nodes);
        const auto *v = values[static_cast<size_t>(j)].value();
        signs[static_cast<size_t>(i)] = arb_is_positive(v) ? 1 : (arb_is_negative(v) ? -1 : 0);
        point_values[static_cast<size_t>(i)] = &values[static_cast<size_t>(j)];
    }
    const slong cells = nodes;
    std::vector<char> change(static_cast<size_t>(cells), 0);
    for (slong i = 0; i < cells; ++i) {
        change[static_cast<size_t>(i)] = signs[static_cast<size_t>(i)] * signs[static_cast<size_t>(i + 1)] < 0;
        stats.cells_with_sign_change += change[static_cast<size_t>(i)];
    }
    auto cell_of = [&](double x) {
        const auto it = std::upper_bound(points.begin(), points.end(), x);
        return std::clamp<slong>(static_cast<slong>(it - points.begin()) - 1, 0, cells - 1);
    };
    std::vector<slong> counts(static_cast<size_t>(cells), 0);
    for (const double s : seeds) ++counts[static_cast<size_t>(cell_of(s))];
    // Node seeds: a sign-change cell the interpolant left empty gets the secant through its two values.
    {
        RealBall xi, xj, vi, vj, t;
        for (slong i = 0; i < cells; ++i) {
            if (!change[static_cast<size_t>(i)] || counts[static_cast<size_t>(i)] != 0) continue;
            arb_set_d(xi.value(), points[static_cast<size_t>(i)]);
            arb_set_d(xj.value(), points[static_cast<size_t>(i + 1)]);
            arb_set_arf(vi.value(), arb_midref(point_values[static_cast<size_t>(i)]->value()));
            arb_set_arf(vj.value(), arb_midref(point_values[static_cast<size_t>(i + 1)]->value()));
            // x_i - v_i (x_j - x_i) / (v_j - v_i)
            arb_sub(t.value(), xj.value(), xi.value(), scale_bits);
            arb_mul(t.value(), t.value(), vi.value(), scale_bits);
            arb_sub(vj.value(), vj.value(), vi.value(), scale_bits);
            arb_div(t.value(), t.value(), vj.value(), scale_bits);
            arb_sub(t.value(), xi.value(), t.value(), scale_bits);
            seeds.push_back(mid_double(t.value()));
            kind.push_back(2);
            ++stats.node_seeds;
        }
    }
    {
        std::vector<size_t> order(seeds.size());
        for (size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::stable_sort(order.begin(), order.end(), [&](size_t x, size_t y) { return seeds[x] < seeds[y]; });
        std::vector<double> s2(seeds.size());
        std::vector<int> k2(seeds.size());
        for (size_t i = 0; i < order.size(); ++i) {
            s2[i] = seeds[order[i]];
            k2[i] = kind[order[i]];
        }
        seeds.swap(s2);
        kind.swap(k2);
    }
    std::fill(counts.begin(), counts.end(), 0);
    std::vector<slong> cell(seeds.size());
    for (size_t i = 0; i < seeds.size(); ++i) {
        cell[i] = cell_of(seeds[i]);
        ++counts[static_cast<size_t>(cell[i])];
    }
    std::vector<double> certified, unverified;
    slong test_warm = first_evaluation_bits;
    for (size_t i = 0; i < seeds.size(); ++i) {
        if (pass.cancelled()) return stop(RealLineSeedStop::Cancelled);
        const slong ci = cell[i];
        const double s = seeds[i];
        if (change[static_cast<size_t>(ci)] && counts[static_cast<size_t>(ci)] == 1) {
            certified.push_back(s);
            ++stats.alone;
            continue;
        }
        if (!limits.seed_tests) {
            unverified.push_back(s);
            ++stats.unverified;
            continue;
        }
        const double full = std::ldexp(std::fabs(s), -seed_half_width_bits);
        double delta = full;
        if (i > 0) delta = std::min(delta, (s - seeds[i - 1]) / 2);
        if (i + 1 < seeds.size()) delta = std::min(delta, (seeds[i + 1] - s) / 2);
        ++stats.tested;
        const int left = pass.certain_sign(s - delta, test_warm), right = pass.certain_sign(s + delta, test_warm);
        if (left * right < 0) {
            certified.push_back(s);
            ++stats.certified_by_test;
        } else if (delta < full) {
            unverified.push_back(s); // a close pair's seed: its bracket was cut by a neighbour
            ++stats.unverified;
        } else {
            ++stats.dropped;
        }
    }
    // At most n seeds: the certified ones first (distinct roots), then unverified ones while there is room.
    std::vector<double> kept = certified;
    for (const double s : unverified) {
        if (static_cast<slong>(kept.size()) >= n) break;
        kept.push_back(s);
    }
    if (static_cast<slong>(kept.size()) > n) kept.resize(static_cast<size_t>(n));
    std::sort(kept.begin(), kept.end());
    // Equal seeds are no use to the iteration (its points must be distinct): one of each.
    kept.erase(std::unique(kept.begin(), kept.end()), kept.end());
    stats.cells_ns = elapsed_ns(t0);
    stats.stop = RealLineSeedStop::Seeds;
    stats.total_ns = elapsed_ns(started);
    result.seeds = std::move(kept);
    return result;
}

} // namespace certpolysolve::detail
