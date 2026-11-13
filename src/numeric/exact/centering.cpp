// Exact translation of an exact polynomial, p(x) to a nonzero constant multiple of p(x + c) (center_integer,
// center_gaussian, ExactTaylorCache), and the cost heuristics, never certificates, that the centring rules read.
//
// The centre: round_center rounds each part to a short dyadic, or to zero when its ball contains zero;
// center_materialization_fits and center_owned_bits bound its size. The height rule keeps a translate whose monic
// height (the largest coefficient's bit length less the leading one's) falls by at least bitcount(n):
// center_height_allows_gain is its necessary condition, translated_height_lower_bound bounds the translate's height
// from below by one ball evaluation of p(c), its constant term up to a factor, and gaussian_center_improves_height
// decides it for Gaussian translates on exact squared moduli. The scale-invariant rules read a HeightProfile of
// log2 |coefficient| from bit lengths (height_profile): height_at is the monic height at scale 2^t, root_scale the
// log2 geometric mean of the nonzero roots' moduli, worst_height_excess the translate's largest excess over all
// scales (at the Newton polygons' breakpoints and limits), and unseparated_origin_roots weighs, by Pellet's test at
// the origin (certified_group_scale), the roots near the origin that a shift by c carries away.
// The shift: center_integer charges a bound on every coefficient of B^n p(x + A/B) to the budget, then builds it by
// FLINT's Taylor shift and divides out the content; center_gaussian shifts by the real part, then by the imaginary
// part between exact quarter turns (gaussian_integer_shift).
//
// State: none, but each ExactTaylorCache owns its source and the translate at its last centre, or that it failed.
// Entry points: all of centering.hpp, from factor_state.cpp (prepare_secular_center, refine_center, attach_origin, and
// owned_exact_taylor, whose cache structured_source.cpp's PolynomialSource reads) and solver.cpp (factor_step).
// Not here: the rules that accept or refuse a centre (those two files) and their EnginePolicy switches (engine_policy.hpp).
// Design notes: docs/mathematics.md, docs/controller-coordinates-findings-20260924.md (the Gaussian shift, the height
// rule), docs/selector-campaign-20260928.md (the scale-invariant rules), docs/final-report-20261011/README.md (item 6).
#include "numeric/exact/centering.hpp"
#include <algorithm>
#include <optional>
#include <cmath>
#include <vector>
#include <climits>
#include <flint/arb_fmpz_poly.h>
namespace certpolysolve::detail {
namespace {
// The word convention is the existing exact-operation estimate's one word
// plus coefficient bits / FLINT_BITS. These are scheduling proxies, not an
// instruction count or a bound on FLINT's arithmetic. Measuring both sides
// of an executed shift keeps denominator/center growth and cancellation in
// view without repricing exact integers at the enclosing ball precision.
double integer_polynomial_words(const fmpz_poly_t polynomial) {
    double words = 0;
    for (slong j = 0; j < fmpz_poly_length(polynomial); ++j)
        words += 1 + static_cast<double>(fmpz_bits(polynomial->coeffs + j) / FLINT_BITS);
    return words;
}
double gaussian_polynomial_words(const GaussianPolynomial &polynomial) {
    double words = 0;
    for (const auto &coefficient : polynomial)
        for (const auto part : {coefficient.real.value(), coefficient.imag.value()})
            words += 2 + static_cast<double>(fmpz_bits(fmpq_numref(part)) / FLINT_BITS) +
                         static_cast<double>(fmpz_bits(fmpq_denref(part)) / FLINT_BITS);
    return words;
}
double enclosed_polynomial_words(const acb_poly_t polynomial) {
    double words = 0;
    for (slong j = 0; j < acb_poly_length(polynomial); ++j)
        for (const auto part : {acb_realref(polynomial->coeffs + j), acb_imagref(polynomial->coeffs + j)})
            // One word for the magnitude radius, plus the midpoint's words.
            words += 2 + static_cast<double>(arf_bits(arb_midref(part)) / FLINT_BITS);
    return words;
}
void research_integer_shift(fmpz_poly_t polynomial, const fmpz_t shift, double *work) {
    const slong degree = fmpz_poly_degree(polynomial);
    const bool charged = work && degree > 0 && !fmpz_is_zero(shift);
    const double before = charged ? integer_polynomial_words(polynomial) : 0;
    fmpz_poly_taylor_shift(polynomial, polynomial, shift);
    if (charged)
        *work += static_cast<double>(degree) * std::max(before, integer_polynomial_words(polynomial));
}
// P(t+A+iB) via real FLINT shifts. For the imaginary shift use
// Q(u)=P(iu), shift Q(u+B), then substitute u=-it. Quarter turns are exact
// signed swaps, so this delegates all large polynomial algebra to FLINT.
void gaussian_integer_shift(std::vector<Integer> &re, std::vector<Integer> &im,
                            const fmpz_t real_shift, const fmpz_t imaginary_shift
                            , double *research_work
                            ) {
    IntegerPolynomial real, imaginary;
    const slong length=static_cast<slong>(re.size());
    for (slong j=0;j<length;++j) {
        fmpz_poly_set_coeff_fmpz(real.value(),j,re[static_cast<size_t>(j)].value());
        fmpz_poly_set_coeff_fmpz(imaginary.value(),j,im[static_cast<size_t>(j)].value());
    }
    if (!fmpz_is_zero(real_shift)) {
        research_integer_shift(real.value(), real_shift, research_work);
        research_integer_shift(imaginary.value(), real_shift, research_work);
    }
    const auto rotate = [&](bool inverse) {
        Integer a,b;
        for (slong j=0;j<length;++j) {
            fmpz_poly_get_coeff_fmpz(a.value(),real.value(),j);
            fmpz_poly_get_coeff_fmpz(b.value(),imaginary.value(),j);
            const slong turn=inverse ? (4-j%4)%4 : j%4;
            if (turn==1) { fmpz_swap(a.value(),b.value());fmpz_neg(a.value(),a.value()); }
            else if (turn==2) { fmpz_neg(a.value(),a.value());fmpz_neg(b.value(),b.value()); }
            else if (turn==3) { fmpz_swap(a.value(),b.value());fmpz_neg(b.value(),b.value()); }
            fmpz_poly_set_coeff_fmpz(real.value(),j,a.value());
            fmpz_poly_set_coeff_fmpz(imaginary.value(),j,b.value());
        }
    };
    if (!fmpz_is_zero(imaginary_shift)) {
        rotate(false);
        research_integer_shift(real.value(), imaginary_shift, research_work);
        research_integer_shift(imaginary.value(), imaginary_shift, research_work);
        rotate(true);
    }
    for (slong j=0;j<length;++j) {
        fmpz_poly_get_coeff_fmpz(re[static_cast<size_t>(j)].value(),real.value(),j);
        fmpz_poly_get_coeff_fmpz(im[static_cast<size_t>(j)].value(),imaginary.value(),j);
    }
}

void gaussian_norms(const GaussianPolynomial &p, Rational &maximum, Rational &leading) {
    Rational value, imaginary;
    for (const auto &c : p) {
        fmpq_mul(value.value(), c.real.value(), c.real.value());
        fmpq_mul(imaginary.value(), c.imag.value(), c.imag.value());
        fmpq_add(value.value(), value.value(), imaginary.value());
        if (fmpq_cmp(value.value(), maximum.value()) > 0) fmpq_set(maximum.value(), value.value());
    }
    fmpq_set(leading.value(), value.value());
}
}
bool center_height_allows_gain(const fmpz_poly_t p) {
    // Monic coefficient height cannot be negative. If it is already smaller
    // than the existing required gain, no translated polynomial can pass the
    // admission rule. Avoid constructing a provably rejected exact shift.
    const slong degree=fmpz_poly_degree(p);
    return degree>0 && FLINT_ABS(fmpz_poly_max_bits(p))-static_cast<slong>(fmpz_bits(fmpz_poly_lead(p))) >=
        static_cast<slong>(FLINT_BIT_COUNT(static_cast<ulong>(degree)));
}
std::optional<slong> translated_height_lower_bound(const fmpz_poly_t p, const fmpq_t center) {
    const slong n = fmpz_poly_degree(p);
    if (n < 1) return std::nullopt;
    // Only the binary exponent of p(center) is used, so two words of working precision are ample; a
    // ball that still contains zero gives no bound, and the caller then builds the translate as before.
    const slong working = 2 * FLINT_BITS;
    RealBall point, value;
    arb_set_fmpq(point.value(), center, working);
    arb_fmpz_poly_evaluate_arb(value.value(), p, point.value(), working);
    Magnitude lower;
    arb_get_mag_lower(lower.value(), value.value());
    if (mag_is_zero(lower.value()) || mag_is_inf(lower.value()) || !fmpz_fits_si(MAG_EXPREF(lower.value())))
        return std::nullopt;
    // A magnitude m 2^(e - MAG_BITS) has m >= 2^(MAG_BITS - 1), so log2|p(center)| >= e - 1.
    const slong log2_value = fmpz_get_si(MAG_EXPREF(lower.value())) - 1;
    // The translate is B^n p(x + A/B) / g for center = A/B and its content g. Its constant term is
    // N/g and its leading one L/g with N = B^n p(center), L = B^n lead(p). A bit length is
    // floor(log2) + 1, and bits(N/g) - bits(L/g) >= bits(N) - bits(L) - 1 for any common divisor g.
    // With bits(N) >= n log2 B + log2_value and bits(L) <= n log2 B + bits(lead) + 1 (equality less
    // one when B is a power of two), the translate's height is at least the value returned.
    const slong lead_bits = static_cast<slong>(fmpz_bits(fmpz_poly_lead(p)));
    const slong rounding = 2; // one bit for the common divisor, one for B^n when B is not a power of two
    if (log2_value < WORD_MIN / 2 || log2_value > WORD_MAX / 2) return std::nullopt;
    return log2_value - lead_bits - rounding;
}
namespace {
slong lowest_present(const HeightProfile &p) {
    return static_cast<slong>(std::find(p.present.begin(), p.present.end(), true) - p.present.begin());
}
// Upper hull of the points (k, bits_k), in increasing k.
std::vector<std::pair<slong, slong>> upper_hull(const HeightProfile &p) {
    std::vector<std::pair<slong, slong>> hull;
    for (slong k = 0; k < static_cast<slong>(p.bits.size()); ++k) {
        if (!p.present[static_cast<size_t>(k)]) continue;
        const std::pair<slong, slong> q{k, p.bits[static_cast<size_t>(k)]};
        while (hull.size() >= 2) {
            const auto &a = hull[hull.size() - 2], &b = hull.back();
            if (static_cast<double>(b.second - a.second) * static_cast<double>(q.first - b.first) >
                static_cast<double>(q.second - b.second) * static_cast<double>(b.first - a.first))
                break;
            hull.pop_back();
        }
        hull.push_back(q);
    }
    return hull;
}
// Scales where the dominant term of max_k (bits_k + t k) changes: consecutive
// hull terms k0 < k1 tie at t = (bits0 - bits1) / (k1 - k0).
void hull_breakpoints(const HeightProfile &p, std::vector<double> &out) {
    const auto hull = upper_hull(p);
    for (size_t i = 0; i + 1 < hull.size(); ++i)
        out.push_back(static_cast<double>(hull[i].second - hull[i + 1].second) /
                      static_cast<double>(hull[i + 1].first - hull[i].first));
}
slong part_bits(const fmpq_t q) {
    // log2|q| within a bit: the bit lengths of numerator and denominator.
    return static_cast<slong>(fmpz_bits(fmpq_numref(q))) - static_cast<slong>(fmpz_bits(fmpq_denref(q)));
}
} // namespace
HeightProfile height_profile(const fmpz_poly_t p) {
    HeightProfile out;
    const slong length = fmpz_poly_length(p);
    out.bits.assign(static_cast<size_t>(std::max<slong>(length, 0)), 0);
    out.present.assign(out.bits.size(), false);
    for (slong k = 0; k < length; ++k) {
        out.present[static_cast<size_t>(k)] = !fmpz_is_zero(p->coeffs + k);
        out.bits[static_cast<size_t>(k)] = static_cast<slong>(fmpz_bits(p->coeffs + k));
    }
    return out;
}
HeightProfile height_profile(const GaussianPolynomial &p) {
    // |re + i im| lies between the larger part and sqrt(2) times it, so the
    // larger part's estimate is within 1.5 bits of log2 of the modulus.
    HeightProfile out;
    out.bits.assign(p.size(), 0);
    out.present.assign(p.size(), false);
    for (size_t k = 0; k < p.size(); ++k) {
        const bool re = !fmpq_is_zero(p[k].real.value()), im = !fmpq_is_zero(p[k].imag.value());
        out.present[k] = re || im;
        if (re && im) out.bits[k] = std::max(part_bits(p[k].real.value()), part_bits(p[k].imag.value()));
        else if (re || im) out.bits[k] = part_bits(re ? p[k].real.value() : p[k].imag.value());
    }
    return out;
}
double height_at(const HeightProfile &p, double t) {
    const slong n = static_cast<slong>(p.bits.size()) - 1;
    if (n <= 0 || !p.present.back()) return 0;
    double maximum = -HUGE_VAL;
    for (slong k = 0; k <= n; ++k)
        if (p.present[static_cast<size_t>(k)])
            maximum = std::max(maximum, static_cast<double>(p.bits[static_cast<size_t>(k)]) + t * static_cast<double>(k));
    return maximum - (static_cast<double>(p.bits.back()) + t * static_cast<double>(n));
}
double root_scale(const HeightProfile &p) {
    const slong n = static_cast<slong>(p.bits.size()) - 1, low = lowest_present(p);
    if (n <= 0 || !p.present.back() || low >= n) return 0;
    return static_cast<double>(p.bits[static_cast<size_t>(low)] - p.bits.back()) / static_cast<double>(n - low);
}
double worst_height_excess(const HeightProfile &after, const HeightProfile &before) {
    if (after.bits.size() != before.bits.size() || after.bits.size() < 2 || !after.present.back() ||
        !before.present.back())
        return HUGE_VAL;
    const slong low_after = lowest_present(after), low_before = lowest_present(before);
    // Below every breakpoint both heights are lines of slope (lowest order - n).
    if (low_after < low_before) return HUGE_VAL;
    std::vector<double> scales;
    hull_breakpoints(after, scales);
    hull_breakpoints(before, scales);
    double worst = 0; // the t -> +inf limit, where both heights are zero
    for (const double t : scales) worst = std::max(worst, height_at(after, t) - height_at(before, t));
    if (low_after == low_before) {
        const size_t low = static_cast<size_t>(low_after);
        worst = std::max(worst, static_cast<double>((after.bits[low] - after.bits.back()) -
                                                    (before.bits[low] - before.bits.back())));
    }
    return worst;
}
namespace {
// Pellet at the origin for `group` roots, radius 2^t: |a_K| 2^(tK) above the sum of the other terms. The
// profile places log2|a_k| within [bits_k - 1, bits_k + 1.5] (height_profile), and the other terms are at
// most n, so each must stay below the K-th term by 2.5 + log2 n bits. Returns the smallest such t, if any.
std::optional<double> certified_group_scale(const HeightProfile &p, slong group) {
    const slong n = static_cast<slong>(p.bits.size()) - 1;
    const double slack = 2.5 + std::log2(static_cast<double>(n));
    const double own = static_cast<double>(p.bits[static_cast<size_t>(group)]);
    double low = -HUGE_VAL, high = HUGE_VAL;
    for (slong k = 0; k <= n; ++k) {
        if (k == group || !p.present[static_cast<size_t>(k)])
            continue;
        const double other = static_cast<double>(p.bits[static_cast<size_t>(k)]);
        if (k < group)
            low = std::max(low, (other - own + slack) / static_cast<double>(group - k));
        else
            high = std::min(high, (own - other - slack) / static_cast<double>(k - group));
    }
    if (!(low < high))
        return std::nullopt;
    return low;
}
} // namespace
slong unseparated_origin_roots(const HeightProfile &p, double center_log2, double reach_bits) {
    const slong n = static_cast<slong>(p.bits.size()) - 1;
    if (n < 2 || !p.present.back() || !std::isfinite(center_log2))
        return 0;
    const auto hull = upper_hull(p);
    const auto scale = [&](size_t v) {
        return static_cast<double>(hull[v - 1].second - hull[v].second) /
               static_cast<double>(hull[v].first - hull[v - 1].first);
    };
    // The hull's slopes are the roots' scales, increasing along it: hull[v].first roots lie at or below
    // scale(v).
    slong carried = 0;
    for (size_t v = 1; v < hull.size() && scale(v) < center_log2 - reach_bits; ++v)
        carried = hull[v].first;
    if (carried < 2)
        return 0;
    for (size_t v = 1; v + 1 < hull.size() && scale(v) < center_log2 - 1; ++v) {
        if (hull[v].first < carried)
            continue;
        const auto low = certified_group_scale(p, hull[v].first);
        if (low && *low < center_log2 - 1)
            return 0;
    }
    return carried;
}
bool center_height_allows_gain(const GaussianPolynomial &p) {
    if (p.size()<2) return false;
    Rational maximum,leading;
    gaussian_norms(p,maximum,leading);
    if (fmpq_is_zero(leading.value())) return false;
    fmpq_mul_2exp(leading.value(),leading.value(),2*FLINT_BIT_COUNT(static_cast<ulong>(p.size()-1)));
    return fmpq_cmp(maximum.value(),leading.value())>=0;
}
bool gaussian_center_improves_height(const GaussianPolynomial &before, const GaussianPolynomial &after) {
    // Compare squared moduli exactly, so complex scaling of the polynomial
    // cannot alter admission. This is a cost heuristic, never a certificate.
    if (before.size() < 2 || before.size() != after.size()) return false;
    Rational old_max, old_leading, new_max, new_leading, left, right;
    gaussian_norms(before, old_max, old_leading);
    gaussian_norms(after, new_max, new_leading);
    if (fmpq_is_zero(old_leading.value()) || fmpq_is_zero(new_leading.value())) return false;
    fmpq_mul(left.value(), old_max.value(), new_leading.value());
    fmpq_mul(right.value(), new_max.value(), old_leading.value());
    fmpq_mul_2exp(right.value(), right.value(), 2 * FLINT_BIT_COUNT(static_cast<ulong>(before.size()-1)));
    return fmpq_cmp(left.value(), right.value()) >= 0;
}

uint64_t center_owned_bits(const fmpz_poly_t polynomial) {
    uint64_t owned = sizeof(IntegerPolynomial) * CHAR_BIT;
    for (slong i = 0; i < fmpz_poly_length(polynomial); ++i)
        owned = resources::add(owned, resources::add(sizeof(fmpz) * CHAR_BIT, fmpz_bits(polynomial->coeffs + i)));
    return owned;
}
uint64_t center_owned_bits(const GaussianPolynomial &polynomial) {
    uint64_t owned = sizeof(GaussianPolynomial) * CHAR_BIT;
    for (const auto &c : polynomial)
        for (const auto q : {c.real.value(), c.imag.value()})
            owned = resources::add(owned, resources::add(sizeof(Rational) * CHAR_BIT,
                resources::add(fmpz_bits(fmpq_numref(q)), fmpz_bits(fmpq_denref(q)))));
    return owned;
}
bool center_materialization_fits(const acb_t target, uint64_t limit) {
    if (!acb_is_finite(target)) return false;
    for (const auto coordinate : {acb_realref(target), acb_imagref(target)}) {
        const auto midpoint = arb_midref(coordinate);
        if (arf_is_zero(midpoint)) continue;
        Integer exponent;
        arf_abs_bound_lt_2exp_fmpz(exponent.value(), midpoint);
        fmpz_abs(exponent.value(), exponent.value());
        if (!fmpz_abs_fits_ui(exponent.value())) return false;
        const auto magnitude = static_cast<uint64_t>(fmpz_get_ui(exponent.value()));
        const auto bits = resources::add(magnitude, resources::nonnegative(arf_bits(midpoint)));
        if (bits == UINT64_MAX || bits > limit) return false;
        limit -= bits;
    }
    return true;
}
bool ExactTaylorCache::enclose(acb_poly_t out, const acb_t target, slong precision) {
    if (!acb_is_finite(target) || !acb_is_exact(target))
        return false;
    if (!attempted_ || !acb_equal(target, center_.value())) {
        acb_set(center_.value(), target);
        attempted_ = true;
        available_ = false;
        const uint64_t n = std::max(static_cast<uint64_t>(fmpz_poly_length(integer_.value())),
                                    static_cast<uint64_t>(gaussian_.size()));
        if (!operations_ || resources::multiply(n, n) > operations_) return false;
        const uint64_t owned = integer_.value()->length ? center_owned_bits(integer_.value()) : center_owned_bits(gaussian_);
        if (!center_materialization_fits(target, allow_expansion_ ? budget_ : std::min(budget_, owned))) return false;
        GaussianRational center;
        arf_get_fmpq(center.real.value(), arb_midref(acb_realref(target)));
        arf_get_fmpq(center.imag.value(), arb_midref(acb_imagref(target)));
        use_integer_ = fmpz_poly_degree(integer_.value()) >= 1 && arb_is_zero(acb_imagref(target));
        if (use_integer_)
            available_ = center_integer(shifted_integer_.value(), integer_.value(), center.real.value(), budget_
                                        , nullptr, research_service_meter_ ? &research_work_ : nullptr
                                        );
        else {
            if (gaussian_.empty()) {
                gaussian_.resize(static_cast<size_t>(fmpz_poly_length(integer_.value())));
                for (size_t j = 0; j < gaussian_.size(); ++j)
                    fmpq_set_fmpz(gaussian_[j].real.value(), integer_.value()->coeffs + j);
            }
            available_ = center_gaussian(shifted_gaussian_, gaussian_, center, budget_
                                         , nullptr, research_service_meter_ ? &research_work_ : nullptr
                                         );
        }
    }
    if (available_) {
        if (use_integer_)
            acb_poly_set_fmpz_poly(out, shifted_integer_.value(), precision);
        else
            enclose_gaussian(out, shifted_gaussian_, precision);
        if (research_service_meter_)
            // The exact translation can be reused at another precision;
            // only this actual conversion/copy is paid on such a hit.
            research_work_ += (use_integer_ ? integer_polynomial_words(shifted_integer_.value())
                                           : gaussian_polynomial_words(shifted_gaussian_)) +
                              enclosed_polynomial_words(out);
    }
    return available_;
}
void round_center(acb_t center, slong bits) {
    if (!acb_is_finite(center))
        return;
    for (arb_ptr part : {acb_realref(center), acb_imagref(center)}) {
        if (arb_contains_zero(part))
            arb_zero(part);
        else
            arf_set_round(arb_midref(part), arb_midref(part), bits, ARF_RND_NEAR);
    }
}
bool center_integer(fmpz_poly_t out, const fmpz_poly_t p, const fmpq_t center, uint64_t budget,
                    PreparationDiagnostic *diagnostic
                    , double *research_work
                    ) {
    if (diagnostic) *diagnostic = {};
    const auto refusal = [&] {
        if (diagnostic) diagnostic->reject("centering-growth-budget", budget);
        return false;
    };
    const slong n = fmpz_poly_degree(p);
    if (n < 1)
        return false;
    const auto numerator_bits = static_cast<uint64_t>(fmpz_bits(fmpq_numref(center)));
    const auto denominator_bits = static_cast<uint64_t>(fmpz_bits(fmpq_denref(center)));
    uint64_t remaining = budget;
    std::vector<slong> support;
    for (slong k = 0; k <= n; ++k)
        if (!fmpz_is_zero(p->coeffs + k))
            support.push_back(k);
    const uint64_t summation_bits =
        static_cast<uint64_t>(FLINT_BIT_COUNT(static_cast<ulong>(support.size())));
    // Bound EACH cleared-denominator coefficient:
    // sum(k>=j) a_k binom(k,j) A^(k-j) B^(n-k+j).
    // binom(k,j)<=2^k and the term count give a rigorous integer bit bound.
    // This retains sparse/graded growth instead of charging the largest
    // possible coefficient to every position. The pre-substitution shift has
    // no larger bounds since B>=1. Division by content only reduces storage.
    for (slong j = 0; j <= n; ++j) {
        uint64_t largest = 0;
        for (const auto k : support)
            if (k >= j) {
                uint64_t term = static_cast<uint64_t>(fmpz_bits(p->coeffs + k));
                for (const auto [exponent, bits] :
                     {std::pair<uint64_t, uint64_t>{static_cast<uint64_t>(k), 1},
                      {static_cast<uint64_t>(k - j), numerator_bits},
                      {static_cast<uint64_t>(n - k + j), denominator_bits}}) {
                    if (term > remaining || (bits && exponent > (remaining - term) / bits))
                        return refusal();
                    term += exponent * bits;
                }
                largest = std::max(largest, term);
            }
        if (summation_bits > remaining || largest > remaining - summation_bits)
            return refusal();
        remaining -= largest + summation_bits;
    }
    IntegerPolynomial transformed;
    Integer power, coefficient, content;
    fmpz_one(power.value());
    // R(t)=B^n p(t/B), then R(t+A), then substitute t=B*x.
    for (slong j = n; j >= 0; --j) {
        fmpz_mul(coefficient.value(), p->coeffs + j, power.value());
        fmpz_poly_set_coeff_fmpz(transformed.value(), j, coefficient.value());
        fmpz_mul(power.value(), power.value(), fmpq_denref(center));
    }
    research_integer_shift(transformed.value(), fmpq_numref(center), research_work);
    fmpz_one(power.value());
    for (slong j = 0; j <= n; ++j) {
        fmpz_mul(transformed.value()->coeffs + j, transformed.value()->coeffs + j, power.value());
        fmpz_mul(power.value(), power.value(), fmpq_denref(center));
    }
    fmpz_poly_content(content.value(), transformed.value());
    fmpz_poly_scalar_divexact_fmpz(out, transformed.value(), content.value());
    if (diagnostic) {
        diagnostic->reason = "accepted";
        diagnostic->accepted = true;
        diagnostic->required_bits_lower_bound = budget - remaining;
        for (slong j = 0; j <= n; ++j) {
            diagnostic->input_bits += static_cast<uint64_t>(fmpz_bits(p->coeffs + j));
            diagnostic->output_bits += static_cast<uint64_t>(fmpz_bits(out->coeffs + j));
        }
    }
    return true;
}
bool center_gaussian(GaussianPolynomial &out, const GaussianPolynomial &p, const GaussianRational &center,
                     uint64_t budget, PreparationDiagnostic *diagnostic
                     , double *research_work
                     ) {
    if (p.size() < 2 || p.size() > static_cast<size_t>(WORD_MAX))
        return false;
    const size_t n = p.size() - 1;
    Integer common, denominator, ar, ai, power, temp, content;
    fmpz_one(common.value());
    uint64_t largest = 0;
    for (const auto &z : p)
        for (const auto q : {z.real.value(), z.imag.value()}) {
            fmpz_lcm(common.value(), common.value(), fmpq_denref(q));
            largest = std::max(largest, static_cast<uint64_t>(fmpz_bits(fmpq_numref(q))));
            if (static_cast<uint64_t>(fmpz_bits(common.value())) > budget / 2 / (n + 1))
                return false;
        }
    fmpz_lcm(denominator.value(), fmpq_denref(center.real.value()), fmpq_denref(center.imag.value()));
    fmpz_divexact(temp.value(), denominator.value(), fmpq_denref(center.real.value()));
    fmpz_mul(ar.value(), temp.value(), fmpq_numref(center.real.value()));
    fmpz_divexact(temp.value(), denominator.value(), fmpq_denref(center.imag.value()));
    fmpz_mul(ai.value(), temp.value(), fmpq_numref(center.imag.value()));
    const auto shift_bits = static_cast<uint64_t>(
        std::max({fmpz_bits(ar.value()), fmpz_bits(ai.value()), fmpz_bits(denominator.value())}));
    const auto per_coefficient = budget / 2 / (n + 1);
    const auto denominator_bits = static_cast<uint64_t>(fmpz_bits(common.value()));
    if (largest > per_coefficient || denominator_bits > per_coefficient - largest)
        return false;
    // The existing conservative bound, largest + common-denominator bits +
    // 2*n*(shift_bits+3), also covers the two real shifts: denominator powers
    // have total degree at most n, and binomial factors and coefficient sums
    // are bounded by the remaining linear-in-n allowance. Quarter turns add
    // no bits. This bounds coefficients, not total allocator workspace/RSS.
    const auto shift_allowance = (per_coefficient - largest - denominator_bits) / 2 / n;
    if (shift_allowance < 3 || shift_bits > shift_allowance - 3)
        return false;
    std::vector<Integer> re(n + 1), im(n + 1);
    fmpz_one(power.value());
    for (size_t k = n + 1; k-- > 0;) {
        fmpz_divexact(temp.value(), common.value(), fmpq_denref(p[k].real.value()));
        fmpz_mul(re[k].value(), temp.value(), fmpq_numref(p[k].real.value()));
        fmpz_mul(re[k].value(), re[k].value(), power.value());
        fmpz_divexact(temp.value(), common.value(), fmpq_denref(p[k].imag.value()));
        fmpz_mul(im[k].value(), temp.value(), fmpq_numref(p[k].imag.value()));
        fmpz_mul(im[k].value(), im[k].value(), power.value());
        fmpz_mul(power.value(), power.value(), denominator.value());
    }
    gaussian_integer_shift(re,im,ar.value(),ai.value()
                           , research_work
                           );
    fmpz_one(power.value());
    for (size_t j = 0; j <= n; ++j) {
        fmpz_mul(re[j].value(), re[j].value(), power.value());
        fmpz_mul(im[j].value(), im[j].value(), power.value());
        fmpz_gcd(content.value(), content.value(), re[j].value());
        fmpz_gcd(content.value(), content.value(), im[j].value());
        fmpz_mul(power.value(), power.value(), denominator.value());
    }
    GaussianPolynomial proposed(n + 1);
    for (size_t j = 0; j <= n; ++j) {
        fmpz_divexact(re[j].value(), re[j].value(), content.value());
        fmpz_divexact(im[j].value(), im[j].value(), content.value());
        fmpq_set_fmpz(proposed[j].real.value(), re[j].value());
        fmpq_set_fmpz(proposed[j].imag.value(), im[j].value());
    }
    if (diagnostic)
        // Numerator bits before and after, as center_integer reports them.
        for (size_t j = 0; j <= n; ++j) {
            diagnostic->input_bits += static_cast<uint64_t>(fmpz_bits(fmpq_numref(p[j].real.value()))) +
                                      static_cast<uint64_t>(fmpz_bits(fmpq_numref(p[j].imag.value())));
            diagnostic->output_bits += static_cast<uint64_t>(fmpz_bits(re[j].value())) +
                                       static_cast<uint64_t>(fmpz_bits(im[j].value()));
        }
    out = std::move(proposed);
    return true;
}
} // namespace certpolysolve::detail
