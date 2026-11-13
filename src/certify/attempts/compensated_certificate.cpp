// The certificate attempt's first tier (compensated_certificate.hpp has the method; the proof is in
// dev/docs/small-input-fast-path-20261026.md, "The proof").
//
// The rounding model. Every operation below is an IEEE binary64 operation rounded to nearest (the library is built
// with -ffp-contract=off: a product is fused only where std::fma is written). In the normal range each result is the
// exact result times (1 + delta), |delta| <= u = 2^-53. Below it a sum or difference is exact, and a product, a
// fused product or a quotient is off by at most the smallest subnormal, 2^-1074, in absolute terms. Overflow cannot
// pass silently: every quantity a test compares is checked finite, and a comparison with NaN is false.
//
// What the code relies on:
// - TwoSum (machine::DoubleDouble::sum) is exact: a + b = s + e. TwoProd, p = ab and e = fma(a, b, -p), is exact in
//   the normal range and off by at most 2^-1074 below it.
// - Complex arithmetic by the textbook formulas: a product's error is at most sqrt(2) gamma_2 |a||b| (Higham,
//   Accuracy and Stability of Numerical Algorithms, 2nd ed., Lemma 3.5), plus 2^-1074 for each of its four real
//   products below the normal range; a sum's componentwise relative errors give at most u |a + b|.
// - gamma_k = k u / (1 - k u), and gamma_a + gamma_b + gamma_a gamma_b <= gamma_(a + b).
// - Upper bounds on nonnegative quantities are computed in binary64 and inflated by the number of roundings on their
//   longest path (inflated); lower bounds on positive quantities deflated likewise.
//
// State: the prepared coefficients (CompensatedCertificate). Entry points: certify, evaluate,
// second_derivative_bound (certificate_attempt.cpp; tests/test_compensated_certificate.cpp).
#include "certify/attempts/compensated_certificate.hpp"
#include "numeric/arithmetic/arf_binary64.hpp"
#include "certify/chart_certificate.hpp"
#include "core/heuristics.hpp"
#include "numeric/arithmetic/machine_conversion.hpp"
#include "numeric/arithmetic/machine_arithmetic.hpp"
#include "core/resource_limits.hpp"
#include "core/raii.hpp"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <utility>

namespace certpolysolve::detail {
namespace {
constexpr double unit_roundoff = 0x1p-53;
constexpr double smallest_subnormal = std::numeric_limits<double>::denorm_min(); // 2^-1074
// The largest degree the rounding counts below are proved for: gamma_bound needs k u <= 2^-30 for k up to 4 n + 5,
// and inflated needs m (m + 1) u <= 1 for m up to 4 n + 10; 2^20 leaves both with room.
constexpr int64_t largest_degree = int64_t(1) << 20;

// An upper bound on a nonnegative quantity that `rounded` operations, each within a relative u, computed as `value`
// (absolute errors below the normal range are added separately): the quantity is at most value (1 - u)^-rounded <=
// value (1 + (rounded + 1) u) for rounded (rounded + 1) u <= 1, and the product below, rounded once more, stays above
// that. One smallest subnormal more covers that product's own rounding below the normal range.
double inflated(double value, uint64_t rounded) {
    return value * (1 + static_cast<double>(rounded + 2) * 0x1p-52) + smallest_subnormal;
}
// A lower bound on a positive quantity that `rounded` such operations computed as `value`: at least value (1 +
// u)^-rounded >= value (1 - rounded u), and the product below, rounded once more, stays below that (and one smallest
// subnormal less below the normal range). The result may be 0 or negative: callers then decline.
double deflated(double value, uint64_t rounded) {
    return value * (1 - static_cast<double>(rounded + 2) * 0x1p-52) - smallest_subnormal;
}
// gamma_k from above: k u is exact, gamma_k <= k u (1 + 2 k u) <= k u (1 + 2^-29) for k u <= 2^-30, and the rounded
// product by (1 + 2^-28) stays above that.
double gamma_bound(uint64_t k) { return (static_cast<double>(k) * unit_roundoff) * (1 + 0x1p-28); }
// One operation on nonnegative numbers, bounded from above or below in every range: in the normal range the result is
// within a relative u, which inflated and deflated cover; below it a sum is exact and a product is off by at most
// 2^-1075, which their added or removed subnormal covers. The bounds after the evaluation's loop are built from these,
// so that no absolute error there needs counting.
double add_up(double a, double b) { return inflated(a + b, 1); }
double mul_up(double a, double b) { return inflated(a * b, 1); }
double sub_down(double a, double b) { return deflated(a - b, 1); }
double mul_down(double a, double b) { return deflated(a * b, 1); }
// An upper bound on y^n for y >= 0, by binary powering with every product rounded outward; infinite where it leaves
// binary64's range.
double power_up(double y, uint64_t n) {
    double result = 1, base = y;
    for (; n; n >>= 1) {
        if (n & 1)
            result = mul_up(result, base);
        if (n > 1)
            base = mul_up(base, base);
    }
    return result;
}

// sqrt(a^2 + b^2) on parts scaled so that the larger lies in [1, 2): the squares' sum then lies in [1, 8] and no
// square falls below the normal range where it matters (a part 2^-511 below the other adds at most 2^-1074 to a sum
// of at least 1). Four roundings: two squares, the sum, the square root. `exponent` is the scale.
double scaled_modulus(double a, double b, int &exponent) {
    const double largest = std::max(std::abs(a), std::abs(b));
    if (largest == 0) {
        exponent = 0;
        return 0;
    }
    exponent = std::ilogb(largest);
    const double as = std::ldexp(a, -exponent), bs = std::ldexp(b, -exponent);
    return std::sqrt(as * as + bs * bs);
}
double modulus_upper(double a, double b) {
    int exponent = 0;
    const double scaled = scaled_modulus(a, b, exponent);
    // The scaled parts' own rounding (the smaller part, where it falls below the normal range) is far inside the
    // margin of one more counted rounding.
    return std::ldexp(inflated(scaled, 5), exponent) + smallest_subnormal;
}
double modulus_lower(double a, double b) {
    int exponent = 0;
    const double scaled = scaled_modulus(a, b, exponent);
    return std::ldexp(deflated(scaled, 5), exponent) - smallest_subnormal;
}
// An arf as a double rounded up, never 0 for a nonzero value (arf_get_d may give 0 below binary64's range).
double upper_double(const arf_t value) {
    if (arf_is_zero(value))
        return 0;
    Float magnitude;
    arf_abs(magnitude.value(), value);
    const double d = arf_to_binary64(magnitude.value(), ARF_RND_UP);
    return d > 0 ? d : smallest_subnormal;
}
double upper_double(const mag_t value) {
    Float f;
    arf_set_mag(f.value(), value);
    return upper_double(f.value());
}
} // namespace

CompensatedCertificate::CompensatedCertificate(const acb_poly_t polynomial) {
    const slong n = acb_poly_degree(polynomial);
    if (n < 1 || n > largest_degree)
        return;
    // The scale: every component's modulus below 2^scale_, so that the scaled midpoints lie below 1.
    bool first = true;
    for (slong k = 0; k <= n; ++k) {
        const acb_srcptr c = polynomial->coeffs + k;
        if (!acb_is_finite(c))
            return;
        for (const arb_srcptr part : {acb_realref(c), acb_imagref(c)}) {
            Float bound;
            arb_get_abs_ubound_arf(bound.value(), part, 64);
            if (arf_is_zero(bound.value()))
                continue;
            const slong e = arf_abs_bound_lt_2exp_si(bound.value());
            if (e == -ARF_PREC_EXACT || e == ARF_PREC_EXACT)
                return; // an exponent beyond what slong carries
            scale_ = first ? e : std::max<int64_t>(scale_, e);
            first = false;
        }
    }
    if (first)
        return; // every coefficient is zero
    std::vector<Term> terms(static_cast<size_t>(n) + 1);
    // Whether some coefficient carries an error before the rounding allowance: a ball's radius, or what the two words
    // leave out (error below is that sum inflated, never 0 since the inflation adds a subnormal).
    bool carried = false;
    for (slong k = 0; k <= n; ++k) {
        const acb_srcptr c = polynomial->coeffs + k;
        Term &t = terms[static_cast<size_t>(k)];
        double errors[2] = {0, 0};
        for (int part = 0; part < 2; ++part) {
            const arb_srcptr component = part ? acb_imagref(c) : acb_realref(c);
            Float scaled, residual, word;
            arf_mul_2exp_si(scaled.value(), arb_midref(component), -scale_);
            Magnitude error;
            mag_mul_2exp_si(error.value(), arb_radref(component), -scale_);
            // The two words by rounding to nearest, the exact residual into the error (below binary64's range a
            // word is 0 or subnormal and the residual holds the rest).
            const double hi = arf_to_binary64(scaled.value(), ARF_RND_NEAR);
            arf_set_d(word.value(), hi);
            arf_sub(residual.value(), scaled.value(), word.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
            const double lo = arf_to_binary64(residual.value(), ARF_RND_NEAR);
            arf_set_d(word.value(), lo);
            arf_sub(residual.value(), residual.value(), word.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
            if (!std::isfinite(hi) || !std::isfinite(lo))
                return;
            Magnitude rest;
            arf_get_mag(rest.value(), residual.value());
            mag_add(error.value(), error.value(), rest.value());
            errors[part] = upper_double(error.value());
            (part ? t.im_hi : t.re_hi) = hi;
            (part ? t.im_lo : t.re_lo) = lo;
        }
        // |delta| <= |delta_re| + |delta_im|, one rounded sum.
        carried = carried || errors[0] > 0 || errors[1] > 0;
        t.error = inflated(errors[0] + errors[1], 1);
        if (!std::isfinite(t.error))
            return;
    }
    // Every polynomial in the balls has degree n: |a_n| from below stays positive.
    const Term &lead = terms.back();
    const double lead_lower = modulus_lower(lead.re_hi, lead.im_hi);
    const double lead_rest =
        inflated(((std::abs(lead.re_lo) + std::abs(lead.im_lo)) + lead.error), 2);
    if (!(lead_lower > lead_rest))
        return;
    // The majorants of p' and p'': k A_k and k (k - 1) A_k with A_k >= |a_k| (|z| <= |Re z| + |Im z|, a sum of five
    // terms).
    std::vector<double> slope(static_cast<size_t>(n) + 1, 0), second(static_cast<size_t>(n) + 1, 0);
    for (slong k = 1; k <= n; ++k) {
        const Term &t = terms[static_cast<size_t>(k)];
        const double size = inflated(
            (((std::abs(t.re_hi) + std::abs(t.im_hi)) + std::abs(t.re_lo)) + std::abs(t.im_lo)) + t.error, 4);
        slope[static_cast<size_t>(k)] = mul_up(static_cast<double>(k), size);
        if (k >= 2)
            second[static_cast<size_t>(k)] = mul_up(static_cast<double>(k) * static_cast<double>(k - 1), size);
    }
    errors_ = carried;
    terms_ = std::move(terms);
    first_ = std::move(slope);
    second_ = std::move(second);
}

CompensatedEvaluation CompensatedCertificate::evaluate(double xr, double xi) const {
    CompensatedEvaluation out;
    out.x_re = xr;
    out.x_im = xi;
    if (!available() || !std::isfinite(xr) || !std::isfinite(xi) || (xr == 0 && xi == 0))
        return out;
    const uint64_t n = terms_.size() - 1;
    const double modulus = modulus_upper(xr, xi);
    if (!std::isfinite(modulus) || !(modulus > 0))
        return out;
    // max(1, |x|)^k <= amplification for every k up to the degree: the most an absolute error made below binary64's
    // normal range can grow on its way to the result. Nothing else is admitted or refused here: an overflow anywhere
    // below leaves a value or a bound non-finite, and the tier declines.
    const double amplification = modulus >= 1 ? power_up(modulus, n) : 1;
    if (!std::isfinite(amplification))
        return out;
    out.modulus = modulus;
    out.amplification = amplification;
    const Term &top = terms_.back();
    // The compensated Horner scheme. With S_k the exact partial sums for the coefficients hi + lo (S_n = a_n,
    // S_k = S_(k+1) x + a_k) and s_k the rounded ones, every step is split exactly:
    //     s_(k+1) x + a_k = s_k + t_k,  t_k = sigma + h1 - h2 + h5 + lo_re + i (sigma' + h3 + h4 + h6 + lo_im)
    // (Cameron and Graillat 2024, Theorems 2.4 and 2.5), up to 2^-1074 for each of the four TwoProd below the normal
    // range. Hence S_k - s_k = sum_(j >= k) t_j x^(j - k), and p(x) = s_0 + E(x), E(x) = sum_j t_j x^j, which the
    // correction c evaluates by Horner's scheme on the rounded t_j.
    double sr = top.re_hi, si = top.im_hi;
    double cr = top.re_lo, ci = top.im_lo; // t_n = the leading coefficient's low word, exactly
    double dr = 0, di = 0;                 // the derivative recurrence on the rounded s
    // The bounds' accumulators, nonnegative, rounded to nearest (their roundings are counted below):
    //   terms        = sum_j T_j X^j,             T_j >= |t_j| (the moduli of t_j's terms, by |Re| + |Im|)
    //   terms_slope  = sum_j j T_j X^(j - 1)
    //   values       = sum_(k < n) (|Re s_(k+1)| + |Im s_(k+1)|) X^k
    //   errors       = sum_k error_k X^k
    //   errors_slope = sum_k k error_k X^(k - 1)
    double terms = std::abs(top.re_lo) + std::abs(top.im_lo), terms_slope = 0, values = 0;
    double errors = top.error, errors_slope = 0;
    const double X = modulus;
    for (uint64_t k = n; k-- > 0;) {
        const Term &a = terms_[k];
        // d = d x + s_(k+1)
        {
            const double pr = dr * xr - di * xi, pi = dr * xi + di * xr;
            dr = pr + sr;
            di = pi + si;
        }
        values = values * X + (std::abs(sr) + std::abs(si));
        // s x as w + x' + y + z exactly (TwoProdCmplx), then w + hi = s' + sigma (TwoSumCmplx).
        const double g1 = sr * xr, h1 = std::fma(sr, xr, -g1);
        const double g2 = si * xi, h2 = std::fma(si, xi, -g2);
        const double g3 = sr * xi, h3 = std::fma(sr, xi, -g3);
        const double g4 = si * xr, h4 = std::fma(si, xr, -g4);
        const auto w_re = machine::DoubleDouble::sum(g1, -g2);
        const auto w_im = machine::DoubleDouble::sum(g3, g4);
        const auto s_re = machine::DoubleDouble::sum(w_re.hi, a.re_hi);
        const auto s_im = machine::DoubleDouble::sum(w_im.hi, a.im_hi);
        const double tr = (((h1 - h2) + w_re.lo) + s_re.lo) + a.re_lo;
        const double ti = (((h3 + h4) + w_im.lo) + s_im.lo) + a.im_lo;
        const double T = (((((((((std::abs(h1) + std::abs(h2)) + std::abs(w_re.lo)) + std::abs(s_re.lo)) +
                               std::abs(a.re_lo)) + std::abs(h3)) + std::abs(h4)) + std::abs(w_im.lo)) +
                           std::abs(s_im.lo)) + std::abs(a.im_lo));
        // c = c x + t
        {
            const double pr = cr * xr - ci * xi, pi = cr * xi + ci * xr;
            cr = pr + tr;
            ci = pi + ti;
        }
        terms_slope = terms_slope * X + terms;
        terms = terms * X + T;
        errors_slope = errors_slope * X + errors;
        errors = errors * X + a.error;
        sr = s_re.hi;
        si = s_im.hi;
    }
    out.value_re = sr;
    out.value_im = si;
    out.correction_re = cr;
    out.correction_im = ci;
    out.slope_re = dr;
    out.slope_im = di;
    // The bounds, for every g in the balls (P the polynomial of the words hi + lo, |g - P| coefficientwise <= error):
    //   |g(x) - s_0 - c_0| <= sum_k error_k |x|^k + |E(x) - c_0|,
    //   |E(x) - c_0| <= gamma_4 sum_j T_j |x|^j (the rounded t_j, sums of five terms a part)
    //                   + gamma_(4n+1) sum_j |t^_j| |x|^j (Horner's scheme in complex arithmetic: n steps of a product
    //                   within (1 + u)^3 and a sum within (1 + u)) <= gamma_(4n+5) sum_j T_j |x|^j;
    //   |g'(x) - d_0| <= sum_k k error_k |x|^(k-1) + sum_(k<n) |S_(k+1) - s_(k+1)| |x|^k + gamma_(4n+1) sum_(k<n)
    //                   |s_(k+1)| |x|^k, and sum_(k<n) |S_(k+1) - s_(k+1)| |x|^k <= sum_j j T_j |x|^(j-1).
    // Roundings: a path through an accumulator has at most 2n of its own (a product and a sum a step) and at most 2n
    // of the accumulator it adds (terms_slope adds terms), plus T_j's nine sums: at most 4n + 10.
    const uint64_t depth = 4 * n + 10;
    // Absolute errors below the normal range: a step makes at most 17 of them (four TwoProd, the correction's and the
    // derivative's four real products each, five accumulator products), each at most 2^-1074, and a step's error
    // reaches the result multiplied by at most |x|^k <= amplification, or k |x|^(k - 1) <= n amplification in the
    // slope's sums. 17 n (n + 1) of them bound every case (17 n (n + 1) 2^-1074 is exact: an integer below 2^45 times
    // a power of two).
    const double allowance =
        mul_up(17.0 * static_cast<double>(n) * static_cast<double>(n + 1) * smallest_subnormal, amplification);
    out.value_error = add_up(
        add_up(mul_up(gamma_bound(4 * n + 5), inflated(terms, depth)), inflated(errors, depth)), allowance);
    out.slope_error = add_up(add_up(add_up(inflated(terms_slope, depth),
                                           mul_up(gamma_bound(4 * n + 1), inflated(values, depth))),
                                    inflated(errors_slope, depth)),
                             allowance);
    out.ok = std::isfinite(out.value_error) && std::isfinite(out.slope_error) && std::isfinite(sr) &&
             std::isfinite(si) && std::isfinite(cr) && std::isfinite(ci) && std::isfinite(dr) && std::isfinite(di);
    return out;
}

double CompensatedCertificate::second_derivative_bound(double modulus, double reach) const {
    const uint64_t n = terms_.size() - 1;
    if (n < 2)
        return 0;
    const double Y = add_up(modulus, reach); // >= |z| on the disk
    if (!std::isfinite(Y))
        return std::numeric_limits<double>::infinity();
    const double amplification = Y >= 1 ? power_up(Y, n) : 1; // >= max(1, Y)^k for every k up to the degree
    if (!std::isfinite(amplification))
        return std::numeric_limits<double>::infinity();
    double acc = second_[n];
    for (uint64_t k = n; k-- > 2;)
        acc = acc * Y + second_[k];
    // 2 (n - 2) roundings; n - 2 products below the normal range at most, each off by 2^-1075 and reaching the sum
    // times at most max(1, Y)^n (n 2^-1074 is exact).
    const double allowance = mul_up(static_cast<double>(n) * smallest_subnormal, amplification);
    return add_up(inflated(acc, 2 * n), allowance);
}

bool CompensatedCertificate::certify(Root &root, const acb_t point, const cps_options &options, slong precision,
                                     CompensatedEvaluation *seen) const {
    if (seen)
        *seen = {};
    // A goal at or beyond the two-word precision is never within the tier's reach.
    if (!available() || options.bits >= static_cast<uint64_t>(machine::double_double_bits))
        return false;
    const double xr = arf_to_binary64(arb_midref(acb_realref(point)), ARF_RND_NEAR);
    const double xi = arf_to_binary64(arb_midref(acb_imagref(point)), ARF_RND_NEAR);
    CompensatedRetry retry;
    if (!may_certify(xr, xi, options))
        return false;
    if (prove(root, xr, xi, options, precision, &retry, seen))
        return true;
    for (unsigned k = 0; k < heuristics::compensated_newton_retries && retry.available; ++k) {
        const CompensatedRetry from = retry;
        retry = {};
        if (may_certify(from.re, from.im, options) && prove(root, from.re, from.im, options, precision, &retry))
            return true;
    }
    return false;
}

// A necessary condition of prove's test, from what the evaluation's value bound must contain. With X = |x| from above
// (evaluate's own), S = sum_k error_k X^k is computed by evaluate's operations in evaluate's order, so its value bound,
// which adds S inflated, is at least S. prove's test on the circle of radius 2r needs value_error < |D| (2r - step
// error) <= 2 r U, with U >= |D|. Hence a point where S >= 2 R U, R >= r, cannot pass, and the evaluation is not made.
//   U: Horner's scheme in complex arithmetic (a product within sqrt(2) gamma_2, a sum within u: Higham, Lemma 3.5)
//      keeps every computed partial sum s_k within (1 + gamma_(4 (n - k))) sum_(j >= k) |a_j| X^(j - k) in modulus,
//      and the derivative recurrence on them adds as much again: |D| <= (1 + gamma_(8n + 1)) sum_k k A_k X^(k - 1),
//      A_k >= |a_k| (first_), plus the products' absolute errors below the normal range, which evaluate's slope
//      allowance bounds.
//   R: prove's r is at most |z1| 2^-(bits + 5), |z1| bounded at 64 bits (goal_half_width, without an absolute
//      tolerance), z1 = x - q. Where |q| <= X, |z1| <= 2X. Where |q| > X, the test on the circle of radius r needs
//      K (|q| + r)^2 / 2 < |D| r, K >= K0 = second_derivative_bound(X, 0) (the test's K is the same majorant over a
//      larger disk), and r <= 2 c |q| (c = 2^-(bits + 5) and the bound's rounding): then |q| < 4 c U / K0 and
//      r < 8 c^2 U / K0.
bool CompensatedCertificate::may_certify(double xr, double xi, const cps_options &options) const {
    if (!errors_ || !options.absolute_tolerance.zero() || !available())
        return true;
    const double X = modulus_upper(xr, xi);
    if (!std::isfinite(X) || !(X > 0))
        return true; // evaluate declines here itself
    const uint64_t n = terms_.size() - 1;
    double S = terms_.back().error;
    for (uint64_t k = n; k-- > 0;)
        S = S * X + terms_[k].error;
    const double amplification = X >= 1 ? power_up(X, n) : 1;
    if (!std::isfinite(amplification))
        return true;
    double slope = first_[n];
    for (uint64_t k = n; k-- > 1;)
        slope = slope * X + first_[k];
    // 2 (n - 1) roundings in the majorant's Horner scheme; the absolute errors as evaluate's slope allowance.
    const double allowance =
        mul_up(17.0 * static_cast<double>(n) * static_cast<double>(n + 1) * smallest_subnormal, amplification);
    const double U = add_up(mul_up(1 + gamma_bound(8 * n + 1), inflated(slope, 2 * n)), allowance);
    // c >= 2^-(bits + 5) (1 + 2^-63), goal_half_width's factor with its 64-bit bound on |z1|.
    const auto guard = static_cast<uint64_t>(resources::precision_sum(options.bits, 5));
    const double c = std::ldexp(1 + 0x1p-52, -static_cast<int>(std::min<uint64_t>(guard, std::numeric_limits<int>::max())));
    const double near = mul_up(2 * X, c);
    const double K0 = second_derivative_bound(X, 0);
    // The ratio first: 8 c^2 U itself can fall below binary64's range, where rounding up would make it a subnormal.
    const double far =
        K0 > 0 ? mul_up(mul_up(8 * c, c), inflated(U / K0, 1)) : std::numeric_limits<double>::infinity();
    const double R = std::max(near, far);
    // The value's bound holds S and, added after it, the same allowance for absolute errors below the normal range.
    return !(std::max(S, allowance) >= mul_up(2 * R, U));
}

bool CompensatedCertificate::prove(Root &root, double xr, double xi, const cps_options &options, slong precision,
                                   CompensatedRetry *retry, CompensatedEvaluation *seen) const {
    if (!available() || options.bits >= static_cast<uint64_t>(machine::double_double_bits))
        return false;
    const auto e = evaluate(xr, xi);
    if (seen)
        *seen = e;
    if (!e.ok)
        return false;
    // The Newton step q ~ V / D, V = s_0 + c_0 rounded, D = d_0: the textbook quotient on parts scaled by 2^-f so that
    // D's larger part lies in [1, 2) (|D|^2 then lies in [1, 8]). |V / D - q| <= (sqrt(2) gamma_4 (1 + u) + u) |V/D|
    // (Higham, Lemma 3.5, and the rounding of V) <= 7 u |q| (1 + 8 u), and below the normal range the scaled V's
    // parts and the quotient's products add at most a few 2^-1074 (eight are allowed).
    const double vr = e.value_re + e.correction_re, vi = e.value_im + e.correction_im;
    const double largest = std::max(std::abs(e.slope_re), std::abs(e.slope_im));
    if (!(largest > 0) || !std::isfinite(largest) || !std::isfinite(vr) || !std::isfinite(vi))
        return false;
    const int f = std::ilogb(largest);
    const double sr = std::ldexp(e.slope_re, -f), si = std::ldexp(e.slope_im, -f);
    const double wr = std::ldexp(vr, -f), wi = std::ldexp(vi, -f);
    const double den = sr * sr + si * si;
    const double qr = (wr * sr + wi * si) / den, qi = (wi * sr - wr * si) / den;
    if (!std::isfinite(qr) || !std::isfinite(qi))
        return false;
    const double step = modulus_upper(qr, qi);
    const double step_error = add_up(mul_up(8 * unit_roundoff, step), 8 * smallest_subnormal);
    // |D| from below.
    const double slope_lower = modulus_lower(e.slope_re, e.slope_im);
    if (!(slope_lower > 0))
        return false;
    // The centre z1 = x - q, exactly as two words a part.
    const auto z_re = machine::DoubleDouble::sum(xr, -qr);
    const auto z_im = machine::DoubleDouble::sum(xi, -qi);
    ComplexBall center;
    {
        Float word;
        arf_set_d(arb_midref(acb_realref(center.value())), z_re.hi);
        arf_set_d(word.value(), z_re.lo);
        arf_add(arb_midref(acb_realref(center.value())), arb_midref(acb_realref(center.value())), word.value(),
                ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_set_d(arb_midref(acb_imagref(center.value())), z_im.hi);
        arf_set_d(word.value(), z_im.lo);
        arf_add(arb_midref(acb_imagref(center.value())), arb_midref(acb_imagref(center.value())), word.value(),
                ARF_PREC_EXACT, ARF_RND_DOWN);
    }
    if (acb_is_zero(center.value()))
        return false;
    // The goal's half-width about z1, as a double from below (a smaller rectangle meets the goal as well). A half-width
    // below binary64's normal range is left to the ball certificate: above it, 5r/4 below is rounded by a relative u.
    RealBall goal;
    goal_half_width(goal.value(), center.value(), options, precision);
    Float goal_lower;
    arb_get_lbound_arf(goal_lower.value(), goal.value(), precision);
    const double r = arf_to_binary64(goal_lower.value(), ARF_RND_DOWN);
    if (!(r >= DBL_MIN) || !std::isfinite(r))
        return false;
    // The majorant of |g''| over the disk about x that holds both circles: |z - x| <= |q| + 2r.
    const double reach_out = add_up(step, 2 * r);
    const double second = second_derivative_bound(e.modulus, reach_out);
    if (!std::isfinite(second))
        return false;
    // Rouché's theorem on |z - z1| = rho against the line l(z) = V + D (z - x), whose one zero is x - V/D, within
    // step_error of z1: on the circle |z - x| <= |q| + rho = reach, and
    //   |g(z) - l(z)| <= value_error + slope_error reach + second reach^2 / 2 < |D| (rho - step_error) <= |l(z)|
    // makes g's zeros in the open disk as many as l's: one.
    const auto rouche = [&](double rho) {
        if (!(rho > step_error))
            return false;
        const double reach = add_up(step, rho);
        const double room = sub_down(rho, step_error);
        if (!(room > 0))
            return false;
        const double lhs = add_up(add_up(e.value_error, mul_up(e.slope_error, reach)),
                                  mul_up(mul_up(0.5, second), mul_up(reach, reach)));
        const double rhs = mul_down(slope_lower, room);
        return std::isfinite(lhs) && std::isfinite(rhs) && lhs < rhs;
    };
    // The goal's circles first. Where they fail and the second-derivative term is why (it grows as the radius
    // squared: at a coarse goal r passes |D| / K however close the point is), circles of a radius rho below r, from
    // the same evaluation: K bounds |g''| on the disk the goal's circles used, which holds the smaller ones, and the
    // rectangle of half-width rho meets the goal as well (heuristics::compensated_circle_curvature_share). rho is the
    // power of two halfway, in exponent, between `low`, below which the inner circle's room |D| (rho - step error)
    // cannot hold the errors and the step's own terms (the rho^2 term left out), and `high`, above which the
    // second-derivative term takes more than the share of the outer circle's |D| 2 rho. Both are only a choice:
    // rouche decides.
    double half_width = r;
    bool passed = rouche(r) && rouche(2 * r);
    const double high =
        passed || !(second > 0) ? 0 : heuristics::compensated_circle_curvature_share * slope_lower / second;
    if (std::isfinite(high) && high < r && high >= 2 * DBL_MIN) {
        double rho = std::ldexp(1.0, std::ilogb(high));
        const double denominator = slope_lower - e.slope_error - second * step;
        const double low =
            (e.value_error + e.slope_error * step + 0.5 * second * step * step + slope_lower * step_error) /
            denominator;
        if (denominator > 0 && std::isfinite(low) && low >= DBL_MIN && low < high) {
            const int lowest = std::ilogb(low) + 1, highest = std::ilogb(high);
            if (lowest <= highest) {
                const int sum = lowest + highest; // halved toward minus infinity
                rho = std::ldexp(1.0, sum >= 0 ? sum / 2 : -((1 - sum) / 2));
            }
        }
        if (rouche(rho) && rouche(2 * rho)) {
            half_width = rho;
            passed = true;
        }
    }
    if (!passed) {
        // Kantorovich's h = K |q| / |D| <= 1/2, with K the majorant of |g''| the test used (over the disk of radius
        // |q| + 2r about x, where Kantorovich's own is 2|q|: a choice of where to evaluate next, which proves or
        // refuses by itself): a root lies within 2|q| of x and Newton's iteration from x converges to it, so the point
        // nearest the step's end is closer to it than x. A point of a cluster has h of at least about 1.
        if (retry && (z_re.hi != xr || z_im.hi != xi) && mul_up(second, step) <= 0.5 * slope_lower) {
            retry->re = z_re.hi;
            retry->im = z_im.hi;
            retry->available = true;
        }
        return false;
    }
    // With h the half-width the circles passed (r, or rho): the rectangle z1 +- h (1 + i) holds the disk of radius h
    // and lies in the open disk of radius 2h; published, it must still lie there: inside the square of half-width
    // 5h/4, whose corners are 5 sqrt(2)/4 h < 2h from z1.
    Root proved;
    proved.box = center;
    Float half;
    arf_set_d(half.value(), half_width);
    arb_add_error_arf(acb_realref(proved.box.value()), half.value());
    arb_add_error_arf(acb_imagref(proved.box.value()), half.value());
    proved.count = proved.multiplicity = 1;
    proved.realness = CPS_REAL_UNKNOWN;
    proved.publish(precision);
    // The rounded 5h/4 and the radius's upward rounding to a magnitude (2^-29 relative) keep its corners below
    // 1.77 h. The box's own radius is h rounded up the same way: it holds the disk of radius h, and its corners lie
    // within sqrt(2) h (1 + 2^-29) < 2h.
    ComplexBall container = center;
    arf_set_d(half.value(), 1.25 * half_width);
    arb_add_error_arf(acb_realref(container.value()), half.value());
    arb_add_error_arf(acb_imagref(container.value()), half.value());
    if (!box_contains(container.value(), proved.box.value()) || !meets_goal(proved, options, precision))
        return false;
    proved.certificate =
        "original-coefficient binary64 compensated evaluation, Rouche on two circles about the Newton step";
    root = std::move(proved);
    return true;
}

bool CompensatedCertificate::refine(acb_ptr points, int64_t n, const std::vector<int64_t> &moving, int64_t goal_bits,
                                    uint64_t &evaluated, unsigned max_steps) const {
    using DD = machine::DoubleDouble;
    using Z = machine::Complex<DD>;
    if (!available() || moving.empty() || n != degree()) return false;
    std::vector<Z> z(static_cast<size_t>(n));
    for (int64_t i = 0; i < n; ++i)
        if (!from_ball(z[static_cast<size_t>(i)], points + i)) return false;
    std::vector<Z> a(terms_.size());
    for (size_t k = 0; k < terms_.size(); ++k)
        a[k] = Z(DD::sum(terms_[k].re_hi, terms_[k].re_lo), DD::sum(terms_[k].im_hi, terms_[k].im_lo));
    // The nudge: z (1 + i t), t its own small angle, so no two points keep a conjugate symmetry.
    for (size_t m = 0; m < moving.size(); ++m) {
        Z &w = z[static_cast<size_t>(moving[m])];
        const DD t(heuristics::double_double_nudge * static_cast<double>(m + 1));
        w = w + Z(-(t * w.im), t * w.re);
    }
    // A converged point: its correction below 2^-(goal + margin) of its modulus.
    const double tolerance = std::ldexp(1.0, -static_cast<int>(goal_bits + heuristics::double_double_margin_bits));
    // The reach test (heuristics::double_double_reach_margin_bits): a point stops where its cluster, as the
    // double-double points show it, is beyond what the stop can be met at.
    const double reach_log2 = static_cast<double>(machine::double_double_bits) -
                              static_cast<double>(goal_bits + heuristics::double_double_margin_bits) +
                              heuristics::double_double_reach_margin_bits;
    // threshold[k - 1] = (2^-reach_log2)^(2/k): a point with k neighbours within s, s^k < 2^-reach_log2, is out of
    // reach; compared as squares, s^2 < threshold[k - 1] (relative to the point's squared modulus).
    const size_t group = std::min<size_t>(heuristics::double_double_reach_group - 1, static_cast<size_t>(n - 1));
    std::vector<double> threshold(group);
    for (size_t k = 1; k <= group; ++k) threshold[k - 1] = std::exp2(-2 * reach_log2 / static_cast<double>(k));
    std::vector<double> squares;
    const auto out_of_reach = [&](int64_t i) {
        const Z &w = z[static_cast<size_t>(i)];
        const double re = static_cast<double>(w.re), im = static_cast<double>(w.im);
        const double size = re * re + im * im;
        if (!(size > 0) || !std::isfinite(size) || group == 0) return false;
        squares.clear();
        for (int64_t j = 0; j < n; ++j) {
            if (j == i) continue;
            const Z gap = w - z[static_cast<size_t>(j)];
            const double gr = static_cast<double>(gap.re), gi = static_cast<double>(gap.im);
            squares.push_back(gr * gr + gi * gi);
        }
        std::partial_sort(squares.begin(), squares.begin() + static_cast<std::ptrdiff_t>(group), squares.end());
        for (size_t k = 1; k <= group; ++k) // k neighbours: a cluster of k + 1 points within sqrt(squares[k - 1])
            if (squares[k - 1] < threshold[k - 1] * size) return true;
        return false;
    };
    const Z one(DD(1.0));
    bool moved = false;
    std::vector<Z> next(moving.size());
    std::vector<uint8_t> ok(moving.size()), done(moving.size(), 0);
    for (unsigned step = 0; step < max_steps; ++step) {
        if (step >= heuristics::double_double_reach_first_step &&
            (step - heuristics::double_double_reach_first_step) % heuristics::double_double_reach_every == 0) {
            bool every = true;
            for (size_t m = 0; m < moving.size(); ++m) {
                if (!done[m] && out_of_reach(moving[m])) done[m] = 1;
                every = every && done[m];
            }
            if (every) break;
        }
        bool converged = true;
        for (size_t m = 0; m < moving.size(); ++m) {
            ok[m] = 0;
            // A point whose last correction was below the tolerance stays where it is (its neighbours' steps are
            // then taken against it, as against the points held).
            if (done[m]) continue;
            ++evaluated;
            const Z w = z[static_cast<size_t>(moving[m])];
            Z p = a.back(), d{};
            for (size_t k = a.size() - 1; k-- > 0;) {
                d = d * w + p;
                p = p * w + a[k];
            }
            if (!machine::finite(p) || !machine::finite(d) || (d.re.hi == 0 && d.im.hi == 0)) continue;
            const Z newton = p / d;
            Z sum{};
            for (int64_t j = 0; j < n; ++j) {
                if (j == moving[m]) continue;
                const Z gap = w - z[static_cast<size_t>(j)];
                if (gap.re.hi == 0 && gap.im.hi == 0) continue;
                sum += one / gap;
            }
            const Z denominator = one - newton * sum;
            if (denominator.re.hi == 0 && denominator.im.hi == 0) continue;
            const Z correction = newton / denominator;
            next[m] = w - correction;
            if (!machine::finite(next[m])) continue;
            ok[m] = 1;
            const double size = std::hypot(static_cast<double>(w.re), static_cast<double>(w.im));
            if (std::hypot(static_cast<double>(correction.re), static_cast<double>(correction.im)) > tolerance * size)
                converged = false;
            else
                done[m] = 1;
        }
        for (size_t m = 0; m < moving.size(); ++m)
            if (ok[m]) {
                z[static_cast<size_t>(moving[m])] = next[m];
                moved = true;
            }
        if (converged) break;
    }
    if (!moved) return false;
    for (const int64_t i : moving) {
        to_arf(arb_midref(acb_realref(points + i)), z[static_cast<size_t>(i)].re);
        to_arf(arb_midref(acb_imagref(points + i)), z[static_cast<size_t>(i)].im);
        mag_zero(arb_radref(acb_realref(points + i)));
        mag_zero(arb_radref(acb_imagref(points + i)));
    }
    return true;
}
} // namespace certpolysolve::detail
