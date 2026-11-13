// The algebra of certified covers, declared with most of its contracts in cover.hpp: exporting a region
// (Root::publish), disjointness, merging, realness by conjugation, completing a cover at a stop, and the goal tests.
//
// Root::publish rounds the box outward to exact dyadic endpoints and replaces it by a containing ball, so that later
// tests judge the exported rectangle. disjoint, overlapping_roots and reflection_met test pairs of regions
// (visit_real_neighbours sweeps the real axis). prove_real_by_reflection applies conjugation
// symmetry: in a complete disjoint cover for real coefficients (the caller's precondition), a region of one distinct
// root whose mirror image meets no other region holds that root's conjugate, so the root is real (keep_real_ends, for
// a cover already published and judged, keeps the real ends instead of publishing the region again). merge_overlapping
// unites overlapping regions until none overlap, adding counts, which is sound for coprime factors.
// complete_cover_at_stop puts the roots no region holds into one square about the origin, which absorbs the regions
// inside it: its half-side halves from the power of two above Cauchy's bound (64 rungs at most) while no region
// straddles it and certify_cluster_box proves its count, the degree less the regions outside; without such a square the
// cover is Cauchy's square alone (cauchy_cluster: every root has |z| < 1 + max |a_k/a_n|). The goal tests are
// sufficient: with E the sum of the coordinate radii and L the larger lower bound of |Re z| and |Im z|, E <= 2^-bits L
// bounds the centre's relative error at every point of the box; accuracy_met also needs one distinct root with its
// multiplicity, approximation_met also tries the squared Euclidean form, and meets_goal picks the goal's test.
// cover_structure, unresolved_radius_log2, cover_ahead and prefer_cover rank certified covers and prove nothing.
//
// State: none.
// Entry points (callers): Root::publish, meets_goal and disjoint (fifteen, thirteen and eleven files);
// merge_overlapping and prove_real_by_reflection (solver.cpp, solver_program.cpp, prepared_piecewise.cpp);
// complete_cover_at_stop (prepared_piecewise.cpp, piecewise_roots.cpp); overlapping_roots (prepared_piecewise.cpp,
// certificate_attempt.cpp, secular_refinement.cpp); prefer_cover (solver.cpp, solver_program.cpp, solver_native.cpp,
// secular_refinement.cpp); cover_structure, unresolved_radius_log2, cover_ahead and cauchy_cluster (solver.cpp);
// accuracy_met (real_refinement.cpp, sturm.cpp); approximation_met (clustering.cpp).
// Not here: disk counts (counting.cpp), power covers (power_composition.cpp), secular's classify_cover_realness.
// Design notes: docs/mathematics.md (coverage, Cauchy's fallback, export), docs/work-queue-20260930.md (item 3.2).
#include "core/cover.hpp"
#include "certify/cover_sweep.hpp"
#include "numeric/arithmetic/arf_binary64.hpp"
#include "core/query.hpp"
#include "certify/counting.hpp"
#include "core/wide_int.hpp"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <numeric>
#include <tuple>

namespace certpolysolve::detail {
namespace {
void publish_component(arb_t component, Endpoint &lower, Endpoint &upper, slong precision) {
    Float lo, hi;
    // Finite precision here prevents a huge exponent gap from expanding an
    // otherwise compact ball into an enormous integer during serialization.
    arb_get_interval_arf(lo.value(), hi.value(), component, precision);
    lower = Endpoint(lo.value());
    upper = Endpoint(hi.value());
    arb_set_interval_arf(component, lo.value(), hi.value(), precision);
}
void enclosure_error(arb_t error, const acb_t box, slong precision) {
    RealBall extra;
    arb_zero(error);
    arf_set_mag(arb_midref(error), arb_radref(acb_realref(box)));
    arf_set_mag(arb_midref(extra.value()), arb_radref(acb_imagref(box)));
    arb_add(error, error, extra.value(), precision);
}
} // namespace

namespace {
// The decimal digits of an unsigned 128-bit integer, as fmpz_get_str writes them, into a caller's buffer backwards
// from its end: chunks of the largest power of ten a 64-bit word holds (10^19), the lowest first, each but the
// highest padded to 19 digits; within a word two digits a division (the table of "00" to "99"). The text is built
// in place and copied into its string once: with a std::string a piece and the sign prepended, a degree-32 solve's
// 128 ends cost about 955 instructions each, most of it the insertion and a second allocation
// (dev/docs/easy-panel-coverage-20261027.md).
constexpr int decimal_digits_per_word = 19;
constexpr uint64_t decimal_word = 10'000'000'000'000'000'000ull; // 10^decimal_digits_per_word
constexpr int max_decimal_digits_128 = 39;                        // 2^128 - 1 has 39 digits
constexpr char two_digits[] = "00010203040506070809101112131415161718192021222324252627282930313233343536373839"
                              "40414243444546474849505152535455565758596061626364656667686970717273747576777879"
                              "8081828384858687888990919293949596979899";
// The digits of `value`, at least `width` of them (zero-padded), ending at `end`; returns the first.
char *word_digits(uint64_t value, char *end, int width) {
    char *first = end;
    while (value >= 100) {
        const auto pair = static_cast<size_t>(value % 100);
        value /= 100;
        *--first = two_digits[2 * pair + 1];
        *--first = two_digits[2 * pair];
    }
    if (value >= 10) {
        *--first = two_digits[2 * value + 1];
        *--first = two_digits[2 * value];
    } else {
        *--first = static_cast<char>('0' + value);
    }
    while (end - first < width) *--first = '0';
    return first;
}
char *decimal_digits(uint128 value, char *end) {
    char *first = end;
    while (value >= decimal_word) {
        const auto chunk = static_cast<uint64_t>(value % decimal_word);
        value /= decimal_word;
        first = word_digits(chunk, first, decimal_digits_per_word);
    }
    return word_digits(static_cast<uint64_t>(value), first, 1);
}
} // namespace

// arf_get_fmpz_2exp's mantissa and exponent, in decimal, without an fmpz where they fit in words: a mantissa of one or
// two limbs (an end of up to 128 bits, as every end published at up to 112 bits of precision) and a word-size
// exponent. The odd mantissa is the limbs shifted right by their trailing zeros, the exponent the arf's less the
// limbs' bits and plus that shift (FLINT 3.6, arf/get.c). An fmpz above 62 bits takes an mpz from FLINT's pool, whose
// first use in a process initializes a block of 4,064 (about 1.1 million instructions: with ends as arf, a dense
// degree-8 solve uses no other mpz, so a one-shot solve skipped nothing else; dev/docs/easy-panel-rerun-20261027.md).
// Larger mantissas and exponents take FLINT's path. test_endpoint_text compares both with it.
EndpointText Endpoint::text() const {
    constexpr slong word_limbs = 2; // the mantissa limbs a uint128 holds
    nn_srcptr limbs;
    slong count;
    if (arf_is_special(value_)) return {}; // zero (a published end is finite): "0" and "0", as arf_get_fmpz_2exp
    ARF_GET_MPN_READONLY(limbs, count, value_);
    if (count <= word_limbs && !COEFF_IS_MPZ(*ARF_EXPREF(value_))) {
        // The lowest limb is not zero (an arf has no zero limbs at the bottom), so the shift is a bit position, 0 to
        // 63; flint_ctz gives an int under GCC and Clang, a 64-bit ulong under MSVC.
        const int shift = static_cast<int>(flint_ctz(limbs[0]));
        uint128 mantissa = limbs[0];
        if (count == 2) mantissa |= static_cast<uint128>(limbs[1]) << FLINT_BITS;
        mantissa >>= shift;
        // |the arf's exponent| <= COEFF_MAX < 2^62 and the subtrahend is at most 128: no overflow
        const slong exponent = *ARF_EXPREF(value_) - (count * FLINT_BITS - shift);
        char text[1 + max_decimal_digits_128]; // a sign and the digits
        char *const end = text + sizeof text;
        char *first = decimal_digits(mantissa, end);
        if (ARF_SGNBIT(value_)) *--first = '-';
        return {std::string(first, end), std::to_string(exponent)};
    }
    Integer mantissa, exponent;
    get(mantissa.value(), exponent.value());
    return {integer_string(mantissa.value()), integer_string(exponent.value())};
}
std::string Endpoint::mantissa_text() const { return text().mantissa; }
std::string Endpoint::exponent_text() const { return text().exponent; }

void publish_text(cps_result &result) {
    result.exported.resize(result.roots.size());
    for (size_t i = 0; i < result.roots.size(); ++i)
        for (size_t k = 0; k < result.exported[i].size(); ++k)
            result.exported[i][k] = result.roots[i].endpoints[k].text();
}

cps_result *published_result(cps_result result) {
    auto out = std::make_unique<cps_result>(std::move(result));
    publish_text(*out);
    return out.release();
}

void Root::publish(slong precision) {
    if (!acb_is_finite(box.value()))
        throw std::logic_error("cannot publish a nonfinite root enclosure");
    publish_component(acb_realref(box.value()), endpoints[0], endpoints[1], resources::precision_sum(precision, 16));
    publish_component(acb_imagref(box.value()), endpoints[2], endpoints[3], resources::precision_sum(precision, 16));
}

namespace {
// A box's parts as binary64 intervals rounded outward: the midpoint between the binary64 number below it
// (arf_get_d's floor, beyond binary64's range too, as arf_to_binary64 gives it) and the next one up, the radius from
// above (mag_get_d), one step outward after each rounded sum (next_down and next_up, std::nextafter on the bits). Two
// boxes whose binary64 intervals lie apart in one part lie apart exactly, so
// acb_overlaps and arb_overlaps would say so; the exact tests run only on the pairs these bounds do not separate.
// The real part is bounded for every region (the sweep's order); the imaginary part only for a pair the sweep visits.
struct OutwardBox {
    double re_lo, re_hi;
    mutable double im_lo = 0, im_hi = 0;
    mutable bool imaginary = false;
};
void outward_part(double &lo, double &hi, const arb_struct *x) {
    const double radius = mag_get_d(arb_radref(x));
    const double below = arf_to_binary64(arb_midref(x), ARF_RND_FLOOR);
    lo = next_down(below - radius);
    hi = next_up(next_up(below) + radius);
}
void imaginary_bounds(const OutwardBox &b, const Root &root) {
    if (!b.imaginary) {
        outward_part(b.im_lo, b.im_hi, acb_imagref(root.box.value()));
        b.imaginary = true;
    }
}
// visit(i, j, a, b) for every pair of regions whose real intervals may meet, each pair once, until a visit returns
// false: a sweep over the real axis, the regions sorted by their lower real bounds, each visited with those before it
// whose upper bound it has not passed. The pairs left out have disjoint real intervals, so their boxes, and a box and
// the reflection of another in the real axis, cannot meet. Every visitor's answer is independent of the order of the
// visits (disjoint, overlapping_roots and reflection_met only collect). Until 25 October 2026 fewer than 128 regions
// compared every pair, the most expensive part of a solve at degree 64 (dev/docs/easy-inputs-profile-20261025.md).
// Since 26 October the sweep runs on binary64 bounds (a and b, the pair's OutwardBox) wherever every region's bounds
// are finite there, and passes them on for a first test; otherwise on 64-bit bounds, as before (a and b null).
template <class Visitor>
bool visit_real_neighbours(const std::vector<Root> &roots, Visitor visit) {
    if (roots.size() < 2) return true;  // no pair to visit
    std::vector<OutwardBox> outward(roots.size());
    bool finite = true;
    for (size_t i = 0; i < roots.size() && finite; ++i) {
        auto &b = outward[i];
        outward_part(b.re_lo, b.re_hi, acb_realref(roots[i].box.value()));
        finite = std::isfinite(b.re_lo) && std::isfinite(b.re_hi);
    }
    std::vector<size_t> order(roots.size());
    std::iota(order.begin(), order.end(), 0);
    std::vector<size_t> active;
    if (finite) {
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return outward[a].re_lo < outward[b].re_lo; });
        for (auto i : order) {
            active.erase(std::remove_if(active.begin(), active.end(),
                                        [&](size_t j) { return outward[j].re_hi < outward[i].re_lo; }),
                         active.end());
            for (auto j : active)
                if (!visit(i, j, &outward[i], &outward[j])) return false;
            active.push_back(i);
        }
        return true;
    }
    struct Bounds {
        Float lo, hi;
    };
    std::vector<Bounds> bounds(roots.size());
    for (size_t i = 0; i < roots.size(); ++i) {
        // Rounded OUTWARD bounds only prune impossible overlaps. Keeping
        // 64-bit bounds avoids expanding compressed exponent differences.
        arb_get_interval_arf(bounds[i].lo.value(), bounds[i].hi.value(), acb_realref(roots[i].box.value()), 64);
    }
    std::sort(order.begin(), order.end(),
              [&](size_t a, size_t b) { return arf_cmp(bounds[a].lo.value(), bounds[b].lo.value()) < 0; });
    for (auto i : order) {
        active.erase(
            std::remove_if(active.begin(), active.end(),
                           [&](size_t j) { return arf_cmp(bounds[j].hi.value(), bounds[i].lo.value()) < 0; }),
            active.end());
        for (auto j : active)
            if (!visit(i, j, nullptr, nullptr)) return false;
        active.push_back(i);
    }
    return true;
}
template <class Visitor>
bool visit_overlaps(const std::vector<Root> &roots, Visitor visit) {
    return visit_real_neighbours(roots, [&](size_t i, size_t j, const OutwardBox *a, const OutwardBox *b) {
        if (a) {
            if (a->re_hi < b->re_lo || b->re_hi < a->re_lo) return true;
            imaginary_bounds(*a, roots[i]);
            imaginary_bounds(*b, roots[j]);
            if (a->im_hi < b->im_lo || b->im_hi < a->im_lo) return true;
        }
        return !acb_overlaps(roots[i].box.value(), roots[j].box.value()) || visit(i, j);
    });
}
} // namespace

bool disjoint(const std::vector<Root> &roots, slong precision) {
    if (std::any_of(roots.begin(), roots.end(), [](const Root &root) {
            return !acb_is_finite(root.box.value());
        })) return false;
    if (auto sweep = CoverSweep::prepare(roots, precision)) return sweep->disjoint();
    return visit_overlaps(roots, [](size_t, size_t) { return false; });
}

std::vector<bool> reflection_met(const std::vector<Root> &roots, slong precision) {
    if (auto sweep = CoverSweep::prepare(roots, precision); sweep && sweep->disjoint())
        return sweep->reflection_met();
    std::vector<bool> met(roots.size(), false);
    RealBall reflected;
    visit_real_neighbours(roots, [&](size_t i, size_t j, const OutwardBox *a, const OutwardBox *b) {
        // The reflection of box i meets box j when their real intervals meet
        // and -Im(box i) meets Im(box j). Negating both sides of the second
        // gives the same statement with i and j exchanged, so one test
        // answers for both regions. The binary64 bounds decide first where
        // they show the real intervals, or -Im(box i) and Im(box j), apart.
        if (met[i] && met[j])
            return true;
        if (a) {
            if (a->re_hi < b->re_lo || b->re_hi < a->re_lo) return true;
            imaginary_bounds(*a, roots[i]);
            imaginary_bounds(*b, roots[j]);
            if (-a->im_lo < b->im_lo || b->im_hi < -a->im_hi) return true;
        }
        if (!arb_overlaps(acb_realref(roots[i].box.value()), acb_realref(roots[j].box.value())))
            return true;
        arb_neg(reflected.value(), acb_imagref(roots[i].box.value()));
        if (arb_overlaps(reflected.value(), acb_imagref(roots[j].box.value())))
            met[i] = met[j] = true;
        return true;
    });
    return met;
}

void prove_real_by_reflection(std::vector<Root> &cover, slong precision, bool keep_real_ends) {
    // A complete disjoint cover gives a representation-independent realness
    // proof. For a unique distinct root, its conjugate must lie in this same
    // box if no other box meets the reflection.
    //
    // Which reflections are met is computed once, for the cover as it stands,
    // by a sweep (reflection_met) where every pair of regions was compared.
    // Replacing a region proved real by its real segment, below, does not
    // change the answer for a later one. A region made real holds a root z
    // and, the cover being complete, its conjugate (no other region met the
    // reflection), so its imaginary interval contains zero. For a later
    // region i: if Im(box i) contains zero, the segment and the old box both
    // meet the reflection of box i exactly when the real intervals meet; if
    // it does not, the conjugates of box i's roots are not real and lie in
    // regions that are not replaced, which meet the reflection either way.
    const auto met = reflection_met(cover);
    for (size_t i = 0; i < cover.size(); ++i) {
        auto &root = cover[i];
        if (root.realness != CPS_REAL_UNKNOWN ||
            !(root.count == 1 || (root.multiplicity && root.count == root.multiplicity)))
            continue;
        if (!met[i]) {
            arb_zero(acb_imagref(root.box.value()));
            root.realness = CPS_REAL;
            root.certificate += "+conjugation";
            // A published region's box is the ball about its exported ends, so its real ends hold the root, real,
            // and they stay as published: publishing again would round them outward once more.
            if (keep_real_ends)
                root.endpoints[2] = root.endpoints[3] = Endpoint{};
            else
                root.publish(precision);
        }
    }
}

bool merge_overlapping(std::vector<Root> &cover, slong precision) {
    // The usual cover has no overlap, and then the pairwise loop below changes
    // nothing: the sweep finds that out without comparing every pair.
    if (disjoint(cover))
        return true;
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < cover.size(); ++i)
            for (size_t j = i + 1; j < cover.size();) {
                if (!acb_overlaps(cover[i].box.value(), cover[j].box.value())) {
                    ++j;
                    continue;
                }
                acb_union(cover[i].box.value(), cover[i].box.value(), cover[j].box.value(),
                          resources::precision_sum(precision, 16));
                cover[i].count += cover[j].count;
                cover[i].multiplicity = 0;
                cover[i].realness = cover[i].realness == CPS_REAL && cover[j].realness == CPS_REAL
                                        ? CPS_REAL
                                        : CPS_REAL_UNKNOWN;
                cover[i].certificate = "complete-factor-cover-union";
                cover[i].publish(precision);
                cover.erase(cover.begin() + static_cast<std::ptrdiff_t>(j));
                changed = true;
            }
    }
    return false;
}
bool complete_cover_at_stop(std::vector<Root> &cover, const acb_poly_t polynomial, slong precision,
                            const cps_options *cancellation) {
    const slong degree = acb_poly_degree(polynomial);
    uint64_t certified = 0;
    for (const auto &root : cover) certified += root.count;
    if (degree < 1 || certified >= static_cast<uint64_t>(degree)) return degree >= 0;
    Root global;
    if (!cauchy_cluster(global, polynomial, precision)) return false;
    // 2^e for the smallest e with Cauchy's bound below it: the square of that
    // half-side holds every root.
    Float bound, radius;
    arb_get_ubound_arf(bound.value(), acb_realref(global.box.value()), precision);
    const slong top = arf_abs_bound_lt_2exp_si(bound.value());
    ComplexBall center;
    std::vector<Root> outside, best;
    bool found = false;
    // The disk tests' count context only carries the caller's cancellation: no
    // allowance, so a test that is not cancelled runs exactly as without one, and
    // its counters stay unreported, as the stop cover's tests always were.
    certpolysolve::detail::SearchStatistics unreported{};
    CountContext counting(unreported);
    counting.cancellation = cancellation;
    // The ladder stops at the first half-side that no longer qualifies: a
    // smaller square cannot hold the uncertified roots a larger one did not.
    for (slong e = top; e > top - 64; --e) {
        if (cancellation && cancellation->cancelled()) break;
        arf_one(radius.value());
        arf_mul_2exp_si(radius.value(), radius.value(), e);
        // The exported region is the square about the centre of half-side
        // `radius` (certify_cluster_box); a certified region is outside when it
        // does not meet that square, inside when the square contains it.
        ComplexBall square;
        acb_add_error_arf(square.value(), radius.value());
        outside.clear();
        uint64_t outside_count = 0;
        bool straddles = false;
        for (const auto &root : cover) {
            if (!acb_overlaps(square.value(), root.box.value())) {
                outside.push_back(root);
                outside_count += root.count;
            } else if (!acb_contains(square.value(), root.box.value())) {
                straddles = true;
                break;
            }
        }
        if (straddles) break;
        const slong expected = degree - static_cast<slong>(outside_count);
        if (expected <= 0) break;
        Root disk;
        if (!certify_cluster_box(disk, polynomial, center.value(), radius.value(), expected, 1, precision,
                                 default_graeffe_steps, &counting))
            break;
        disk.certificate = "stop: uncertified roots in a counted square";
        best = outside;
        best.push_back(std::move(disk));
        found = true;
    }
    if (found) {
        cover = std::move(best);
        return true;
    }
    global.certificate = "stop: cauchy-global-count";
    cover = {std::move(global)};
    return true;
}
std::vector<bool> overlapping_roots(const std::vector<Root> &roots) {
    const bool nonfinite = std::any_of(roots.begin(), roots.end(), [](const Root &root) {
        return !acb_is_finite(root.box.value());
    });
    std::vector<bool> conflict(roots.size(), nonfinite);
    if (!nonfinite) visit_overlaps(roots, [&](size_t i, size_t j) {
        conflict[i] = conflict[j] = true;
        return true;
    });
    return conflict;
}

// |z - midpoint|_2 <= r_real + r_imag, while |z|_2 >= the larger of the two coordinate absolute lower bounds. Thus
// this is a sufficient TRUE-relative error test, not merely a midpoint-relative one.
bool relative_accuracy_ball(const acb_t box, uint64_t bits, slong precision) {
    RealBall error, allowance;
    enclosure_error(error.value(), box, precision);
    if (arb_is_zero(error.value()))
        return true;
    Float real_lower, imag_lower;
    arb_get_abs_lbound_arf(real_lower.value(), acb_realref(box), precision);
    arb_get_abs_lbound_arf(imag_lower.value(), acb_imagref(box), precision);
    arb_set_arf(allowance.value(),
                arf_cmp(real_lower.value(), imag_lower.value()) >= 0 ? real_lower.value() : imag_lower.value());
    arb_mul_2exp_si(allowance.value(), allowance.value(), -static_cast<slong>(bits));
    return arb_le(error.value(), allowance.value());
}

// The same test in binary64, decided before the ball form (the design pass of 26 October 2026: the goal checks were
// a twelfth of a degree-8 solve, five a root, nearly all far inside the bound). It holds only where the ball form
// holds at any precision of 53 bits or more, so every answer is the same:
// - E_hi >= the ball sum's upper end. Each radius is a mag (a 30-bit mantissa), exact in binary64 inside its range
//   and replaced from above by 2^-1000 below it; their sum, rounded to nearest, is inflated by 2^-48, past both its
//   own rounding and the ball sum's at 53 bits or more (each within a relative 2^-52).
// - L_lo <= the ball path's lower bound of each |part|, |mid| - rad rounded down at 53 bits or more (within a
//   relative 2^-52 of the exact one): here |mid| rounded toward zero, less the radius, rounded to nearest, then
//   deflated by 2^-49.
// - 2^-bits L_lo is formed only where it stays in binary64's normal range, where the scaling is exact.
// Midpoints or radii beyond 2^1000 in magnitude, and precisions below 53 bits, are left to the ball form (false).
bool relative_accuracy_binary64(const acb_t box, uint64_t bits, slong precision) {
    constexpr slong range = 1000;
    if (precision < 53 || bits > static_cast<uint64_t>(range))
        return false;
    double radius[2], lower[2];
    for (int part = 0; part < 2; ++part) {
        const arb_struct *x = part ? acb_imagref(box) : acb_realref(box);
        const mag_struct *r = arb_radref(x);
        if (mag_is_zero(r))
            radius[part] = 0;
        else if (mag_cmp_2exp_si(r, -range) < 0)
            radius[part] = std::ldexp(1.0, -static_cast<int>(range));
        else if (mag_cmp_2exp_si(r, range) >= 0)
            return false;
        else
            radius[part] = mag_get_d(r);
        const arf_struct *m = arb_midref(x);
        if (arf_is_zero(m) || arf_cmpabs_2exp_si(m, -range) < 0) {
            lower[part] = 0;
            continue;
        }
        if (arf_cmpabs_2exp_si(m, range) >= 0)
            return false;
        const double magnitude = std::abs(arf_to_binary64(m, ARF_RND_DOWN));
        const double difference = magnitude - radius[part];
        lower[part] = difference > 0 ? difference * (1 - 0x1p-49) : 0;
    }
    const double error = (radius[0] + radius[1]) * (1 + 0x1p-48);
    const double allowance = std::ldexp(std::max(lower[0], lower[1]), -static_cast<int>(bits));
    return allowance >= DBL_MIN && error <= allowance;
}

namespace {
// |m_outer - m_inner| + r_inner <= r_outer with the left side rounded up (ARF_RND_UP rounds away from zero): a pass
// means the exact inequality holds, which is arb_contains's question for finite balls.
bool part_contained(const arb_t outer, const arb_t inner) {
    constexpr slong bound_bits = 64;
    Float gap, radius;
    arf_sub(gap.value(), arb_midref(outer), arb_midref(inner), bound_bits, ARF_RND_UP);
    arf_abs(gap.value(), gap.value());
    arf_set_mag(radius.value(), arb_radref(inner));
    arf_add(gap.value(), gap.value(), radius.value(), bound_bits, ARF_RND_UP);
    arf_set_mag(radius.value(), arb_radref(outer));
    return arf_cmp(gap.value(), radius.value()) <= 0;
}
} // namespace

bool box_contains(const acb_t outer, const acb_t inner) {
    if (acb_is_finite(outer) && acb_is_finite(inner) && part_contained(acb_realref(outer), acb_realref(inner)) &&
        part_contained(acb_imagref(outer), acb_imagref(inner)))
        return true;
    return acb_contains(outer, inner) != 0;
}

namespace {
bool coordinate_accuracy_met(const Root &root, uint64_t bits, const Number &absolute_tolerance,
                             slong precision) {
    if (!acb_is_finite(root.box.value()))
        return false;
    if (relative_accuracy_binary64(root.box.value(), bits, precision) ||
        relative_accuracy_ball(root.box.value(), bits, precision))
        return true;
    if (!absolute_tolerance.zero()) {
        RealBall error, allowance;
        enclosure_error(error.value(), root.box.value(), precision);
        absolute_tolerance.enclose(allowance.value(), precision);
        if (arb_le(error.value(), allowance.value()))
            return true;
    }
    return false;
}
} // namespace

bool accuracy_met(const Root &root, uint64_t bits, const Number &absolute_tolerance, slong precision) {
    // Preserve the qualified individual-root stopping rule. A counted cluster
    // is not one distinct root, even when its diameter is sufficiently small.
    return root.multiplicity == root.count &&
           coordinate_accuracy_met(root, bits, absolute_tolerance, precision);
}

bool approximation_met(const Root &root, uint64_t bits, const Number &absolute_tolerance, slong precision) {
    if (coordinate_accuracy_met(root, bits, absolute_tolerance, precision))
        return true;
    if (!acb_is_finite(root.box.value()))
        return false;

    // Squared Euclidean geometry avoids a square root and the extra slack of
    // the inexpensive coordinate bound. For the whole rectangle,
    //   |z - center|^2 <= r_real^2 + r_imag^2,
    //   |z|^2 >= lower_abs(real)^2 + lower_abs(imag)^2.
    // Interval comparisons require the error's upper bound to be no larger
    // than the allowance's lower bound; arithmetic rounding stays explicit.
    RealBall error_squared, allowance_squared, component;
    arb_get_rad_arb(error_squared.value(), acb_realref(root.box.value()));
    arb_sqr(error_squared.value(), error_squared.value(), precision);
    arb_get_rad_arb(component.value(), acb_imagref(root.box.value()));
    arb_sqr(component.value(), component.value(), precision);
    arb_add(error_squared.value(), error_squared.value(), component.value(), precision);

    Float lower;
    arb_get_abs_lbound_arf(lower.value(), acb_realref(root.box.value()), precision);
    arb_set_arf(allowance_squared.value(), lower.value());
    arb_sqr(allowance_squared.value(), allowance_squared.value(), precision);
    arb_get_abs_lbound_arf(lower.value(), acb_imagref(root.box.value()), precision);
    arb_set_arf(component.value(), lower.value());
    arb_sqr(component.value(), component.value(), precision);
    arb_add(allowance_squared.value(), allowance_squared.value(), component.value(), precision);
    // Two shifts avoid overflowing a signed precision type by forming 2*bits.
    arb_mul_2exp_si(allowance_squared.value(), allowance_squared.value(), -static_cast<slong>(bits));
    arb_mul_2exp_si(allowance_squared.value(), allowance_squared.value(), -static_cast<slong>(bits));
    if (arb_le(error_squared.value(), allowance_squared.value()))
        return true;
    if (!absolute_tolerance.zero()) {
        absolute_tolerance.enclose(allowance_squared.value(), precision);
        arb_sqr(allowance_squared.value(), allowance_squared.value(), precision);
        if (arb_le(error_squared.value(), allowance_squared.value()))
            return true;
    }
    return false;
}

bool cluster_radius_met(const Root &root, const Number &radius, slong precision) {
    if (radius.sign() <= 0 || !acb_is_finite(root.box.value()))
        return false;
    RealBall error, allowance;
    enclosure_error(error.value(), root.box.value(), precision);
    radius.enclose(allowance.value(), precision);
    return arb_le(error.value(), allowance.value());
}

bool meets_goal(const Root &root, const cps_options &options, slong precision) {
    switch (engine_goal(options)) {
    case CPS_ROOT_ACCURACY:
        return accuracy_met(root, options.bits, options.absolute_tolerance, precision);
    case CPS_ROOT_APPROXIMATION:
        // Real-domain counts currently come from individual real-root proofs.
        return options.real_only
                   ? accuracy_met(root, options.bits, options.absolute_tolerance, precision)
                   : approximation_met(root, options.bits, options.absolute_tolerance, precision);
    case CPS_CLUSTER_RADIUS:
        return cluster_radius_met(root, options.cluster_radius, precision);
    case CPS_ROOT_ISOLATION:
        // One distinct root with its multiplicity established, and no width requirement.
        return root.multiplicity == root.count && acb_is_finite(root.box.value());
    case CPS_GOAL_DEFAULT: // engine_goal resolves it to CPS_ROOT_ACCURACY
        break;
    }
    return false;
}

std::tuple<uint64_t, uint64_t, int64_t> cover_structure(const std::vector<Root> &roots) {
    uint64_t unresolved = 0, largest = 0;
    for (const auto &root : roots)
        if (root.count != root.multiplicity) {
            unresolved += root.count;
            largest = std::max(largest, root.count);
        }
    return std::tuple{unresolved, largest, -static_cast<int64_t>(roots.size())};
}

double unresolved_radius_log2(const std::vector<Root> &roots) {
    double largest = -std::numeric_limits<double>::infinity();
    for (const auto &root : roots)
        if (root.count != root.multiplicity) {
            Magnitude radius;
            mag_add(radius.value(), arb_radref(acb_realref(root.box.value())),
                    arb_radref(acb_imagref(root.box.value())));
            if (!mag_is_zero(radius.value()))
                largest = std::max(largest, mag_get_d_log2_approx(radius.value()));
        }
    return largest;
}

bool cover_ahead(const std::vector<Root> &a, const std::vector<Root> &b) {
    const auto a_structure = cover_structure(a), b_structure = cover_structure(b);
    return a_structure != b_structure ? a_structure < b_structure
                                      : unresolved_radius_log2(a) < unresolved_radius_log2(b);
}

bool prefer_cover(const std::vector<Root> &candidate, const std::vector<Root> &previous, slong precision) {
    // Quality ordering for already certified complete covers. Prefer resolving
    // multiplicities, then smaller unresolved clusters, then more regions
    // (cover_structure, smaller is better). With equal structure, never
    // replace a cover by one with a worse largest error radius. This is not
    // used to establish mathematical correctness.
    const auto candidate_structure = cover_structure(candidate), previous_structure = cover_structure(previous);
    if (candidate_structure != previous_structure)
        return candidate_structure < previous_structure;
    auto maximum_radius = [&](arf_t maximum, const std::vector<Root> &roots) {
        RealBall error;
        Float upper;
        arf_zero(maximum);
        for (const auto &root : roots) {
            enclosure_error(error.value(), root.box.value(), precision);
            arb_get_ubound_arf(upper.value(), error.value(), precision);
            if (arf_cmp(upper.value(), maximum) > 0)
                arf_set(maximum, upper.value());
        }
    };
    Float candidate_radius, previous_radius;
    maximum_radius(candidate_radius.value(), candidate);
    maximum_radius(previous_radius.value(), previous);
    return arf_cmp(candidate_radius.value(), previous_radius.value()) <= 0;
}

bool cauchy_cluster(Root &root, const acb_poly_t polynomial, slong precision) {
    const auto degree = acb_poly_degree(polynomial);
    if (degree < 1 || acb_contains_zero(polynomial->coeffs + degree))
        return false;
    RealBall absolute, radius;
    ComplexBall ratio;
    Float upper, maximum;
    arf_zero(maximum.value());
    for (slong k = 0; k < degree; ++k) {
        acb_div(ratio.value(), polynomial->coeffs + k, polynomial->coeffs + degree, precision);
        acb_abs(absolute.value(), ratio.value(), precision);
        if (!arb_is_finite(absolute.value()))
            return false;
        arb_get_ubound_arf(upper.value(), absolute.value(), precision);
        if (arf_cmp(upper.value(), maximum.value()) > 0)
            arf_set(maximum.value(), upper.value());
    }
    // Cauchy's bound: every root has |z| < 1 + max_{k<n}|a_k/a_n|.
    // Outward rounding and the containing square preserve the exact count n,
    // uniformly for every polynomial in the coefficient balls. This is a
    // conservative GLOBAL cluster, not Ccluster-style local subdivision.
    arb_set_arf(radius.value(), maximum.value());
    arb_add_ui(radius.value(), radius.value(), 1, precision);
    acb_zero(root.box.value());
    arb_add_error(acb_realref(root.box.value()), radius.value());
    arb_add_error(acb_imagref(root.box.value()), radius.value());
    root.count = static_cast<uint64_t>(degree);
    root.multiplicity = degree == 1 ? 1u : 0u;
    root.certificate = "cauchy-global-count";
    root.publish(precision);
    return acb_is_finite(root.box.value());
}
void sort_roots(std::vector<Root> &roots) {
    std::sort(roots.begin(), roots.end(), [](const Root &a, const Root &b) {
        const auto c = arf_cmp(arb_midref(acb_realref(a.box.value())), arb_midref(acb_realref(b.box.value())));
        return c < 0 || (c == 0 && arf_cmp(arb_midref(acb_imagref(a.box.value())),
                                           arb_midref(acb_imagref(b.box.value()))) < 0);
    });
}
} // namespace certpolysolve::detail
