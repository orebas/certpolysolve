// The evaluation at a point with short parts (limb_horner_evaluator.hpp) against exact arithmetic. Its value
// becomes the numerator of a secular residue that certificates rest on, so the property tested is enclosure:
// the exact value of the midpoint polynomial, and of polynomials at the corners of the coefficients' balls, lie
// in the returned ball. The second property is the header's bound: the radius is at most (n + 1) 2^-precision of
// the sum of the terms' sizes plus the coefficients' own radii carried through, and zero where nothing is
// dropped. The cases are the ones of the derivation: every sign of the point's parts, real and imaginary
// points, points far from one in size, parts far apart (a multiplier at a limb offset, of two or three limbs),
// parts of two limbs, coefficients whose exponents make the window jump by several limbs, cancellation at and
// beside roots, and the refusals.
#include "numeric/evaluate/limb_horner_evaluator.hpp"
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
namespace d = certpolysolve::detail;
static void require(bool value, const std::string &why) {
    if (!value) throw std::runtime_error(why);
}
namespace {
struct Random {
    uint64_t state;
    uint64_t next() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        return state;
    }
    double unit() { return static_cast<double>(next() >> 11) * 0x1p-53; } // [0, 1)
    double signed_unit() { return 2 * unit() - 1; }
    double signed_half() { return (next() & 1 ? -1 : 1) * (0.5 + unit() / 2); } // in size within [1/2, 1)
    long below(long bound) { return static_cast<long>(next() % static_cast<uint64_t>(bound)); }
};
// A number of `bits` random bits (the lowest set) in size within [2^(e - 1), 2^e), e uniform within `spread`
// either way of `centre`; zero one time in sixteen where `zeros`.
void random_float(arf_t out, Random &random, slong bits, long spread, long centre = 0, bool zeros = true) {
    arf_zero(out);
    if (zeros && random.below(16) == 0) return;
    d::Float limb;
    for (slong have = 0; have < bits; have += 64) {
        arf_set_ui(limb.value(), random.next() | 1);
        arf_mul_2exp_si(limb.value(), limb.value(), have);
        arf_add(out, out, limb.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
    }
    arf_set_round(out, out, bits, ARF_RND_DOWN);
    arf_mul_2exp_si(out, out, centre + (spread ? random.below(2 * spread + 1) - spread : 0) - arf_abs_bound_lt_2exp_si(out));
    if (random.next() & 1) arf_neg(out, out);
}
// Ball coefficients: midpoints of `bits` bits, radii 2^-radius_bits of each part's size (none for zero).
d::ComplexPolynomial random_polynomial(Random &random, slong degree, slong bits, long spread, bool real, slong radius_bits) {
    d::ComplexPolynomial p;
    d::ComplexBall c;
    d::Magnitude radius;
    for (slong k = 0; k <= degree; ++k) {
        for (arb_ptr part : {acb_realref(c.value()), acb_imagref(c.value())}) {
            arb_zero(part);
            if (real && part == acb_imagref(c.value())) continue;
            random_float(arb_midref(part), random, bits, spread);
            if (k == degree && arf_is_zero(arb_midref(part))) arf_one(arb_midref(part));
            if (radius_bits) {
                arf_get_mag(radius.value(), arb_midref(part));
                mag_mul_2exp_si(radius.value(), radius.value(), -radius_bits);
                mag_set(arb_radref(part), radius.value());
            }
        }
        acb_poly_set_coeff_acb(p.value(), k, c.value());
    }
    return p;
}
// The exact value at an exact point of the polynomial whose coefficients are the balls' midpoints moved by
// their radii with the signs `corner` gives (no move for corner zero).
void exact_value(acb_t out, const d::ComplexPolynomial &p, const acb_t point, uint64_t corner) {
    const slong length = acb_poly_length(p.value());
    d::ComplexPolynomial exact;
    d::ComplexBall c;
    d::Float radius;
    Random signs{corner | 1};
    for (slong k = 0; k < length; ++k) {
        acb_get_mid(c.value(), p.value()->coeffs + k);
        if (corner) {
            arf_set_mag(radius.value(), arb_radref(acb_realref(p.value()->coeffs + k)));
            if (signs.next() & 1) arf_neg(radius.value(), radius.value());
            arf_add(arb_midref(acb_realref(c.value())), arb_midref(acb_realref(c.value())), radius.value(),
                    ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_set_mag(radius.value(), arb_radref(acb_imagref(p.value()->coeffs + k)));
            if (signs.next() & 1) arf_neg(radius.value(), radius.value());
            arf_add(arb_midref(acb_imagref(c.value())), arb_midref(acb_imagref(c.value())), radius.value(),
                    ARF_PREC_EXACT, ARF_RND_DOWN);
        }
        acb_poly_set_coeff_acb(exact.value(), k, c.value());
    }
    _acb_poly_evaluate_horner(out, exact.value()->coeffs, acb_poly_length(exact.value()), point, ARF_PREC_EXACT);
    require(acb_is_exact(out), "the reference value is not exact");
}
// Enclosure and the radius bound at one point. Whether the value came out exact.
bool check(const d::ComplexPolynomial &p, const acb_t point, slong precision, const std::string &what) {
    d::LimbHornerEvaluator evaluator(p, precision);
    require(evaluator.usable(), what + ": the coefficients were refused");
    d::ComplexBall value, exact;
    require(d::LimbHornerEvaluator::serves(point), what + ": the point was refused");
    require(evaluator.evaluate_value(value.value(), point), what + ": no value");
    // Enclosure: the midpoint polynomial and four corners of the coefficients' balls.
    for (uint64_t corner : {0u, 2u, 4u, 6u, 8u}) {
        exact_value(exact.value(), p, point, corner);
        require(acb_contains(value.value(), exact.value()),
                what + ": the value does not enclose the exact one (corner " + std::to_string(corner) + ")");
    }
    // The bound: (n + 1) 2^-precision sum |a_k| |x|^k and twice the radii's own sum.
    const slong length = acb_poly_length(p.value());
    d::Magnitude modulus, power, size, sizes, radii, term, allowed, radius;
    acb_get_mag(modulus.value(), point);
    mag_one(power.value());
    for (slong k = 0; k < length; ++k) {
        const acb_struct *a = p.value()->coeffs + k;
        arf_get_mag(size.value(), arb_midref(acb_realref(a)));
        arf_get_mag(term.value(), arb_midref(acb_imagref(a)));
        mag_add(size.value(), size.value(), term.value());
        mag_addmul(sizes.value(), size.value(), power.value());
        mag_add(term.value(), arb_radref(acb_realref(a)), arb_radref(acb_imagref(a)));
        mag_addmul(radii.value(), term.value(), power.value());
        mag_mul(power.value(), power.value(), modulus.value());
    }
    mag_add(sizes.value(), sizes.value(), radii.value());
    mag_mul_ui(allowed.value(), sizes.value(), static_cast<ulong>(length));
    mag_mul_2exp_si(allowed.value(), allowed.value(), -precision);
    mag_mul_2exp_si(radii.value(), radii.value(), 1);
    mag_add(allowed.value(), allowed.value(), radii.value());
    mag_max(radius.value(), arb_radref(acb_realref(value.value())), arb_radref(acb_imagref(value.value())));
    require(mag_cmp(radius.value(), allowed.value()) <= 0, what + ": the radius is above the derived bound");
    return acb_is_exact(value.value());
}
bool check(const d::ComplexPolynomial &p, double xr, double xi, slong precision, const std::string &what) {
    d::ComplexBall point;
    acb_set_d_d(point.value(), xr, xi);
    return check(p, point.value(), precision, what);
}
void random_cases() {
    Random random{0x9e3779b97f4a7c15ULL};
    d::ComplexBall point;
    for (slong precision : {64, 117, 234, 468, 936, 1872})
        for (slong degree : {0, 1, 2, 17, 96})
            for (long spread : {0L, 40L, 3000L})
                for (int kind = 0; kind < 8; ++kind) {
                    const bool real = kind & 1;
                    const slong bits = kind & 2 ? precision : 53;
                    const slong radius_bits = kind & 4 ? precision : 0;
                    const auto p = random_polynomial(random, degree, bits, spread, real, radius_bits);
                    const std::string what = "precision " + std::to_string(precision) + ", degree " +
                                             std::to_string(degree) + ", spread " + std::to_string(spread) +
                                             ", kind " + std::to_string(kind);
                    // Every sign of the two parts, in the unit square.
                    for (int signs = 0; signs < 4; ++signs)
                        check(p, (signs & 1 ? -1 : 1) * (0.25 + random.unit()), (signs & 2 ? -1 : 1) * (0.25 + random.unit()),
                              precision, what + ", signs " + std::to_string(signs));
                    // Real and imaginary points, either sign; the origin.
                    check(p, random.signed_unit(), 0, precision, what + ", real point");
                    check(p, 0, random.signed_unit(), precision, what + ", imaginary point");
                    check(p, 0, 0, precision, what + ", the origin");
                    // Points far from one in size: the window slides up or down at every step.
                    check(p, std::ldexp(random.signed_half(), 300), std::ldexp(random.signed_half(), 297), precision,
                          what + ", large point");
                    check(p, std::ldexp(random.signed_half(), -300), std::ldexp(random.signed_half(), -303), precision,
                          what + ", small point");
                    // Parts apart: a multiplier that straddles two limbs, one moved onto a limb boundary, ones
                    // at a limb offset; the real part the larger, then the imaginary.
                    for (int apart : {11, 12, 30, 60, 64, 75, 100, 190, 900})
                        for (int larger = 0; larger < 2; ++larger) {
                            const double big = random.signed_half(), small = std::ldexp(random.signed_half(), -apart);
                            check(p, larger ? small : big, larger ? big : small, precision,
                                  what + ", parts " + std::to_string(apart) + " bits apart, larger " + std::to_string(larger));
                        }
                    // Parts of two limbs (the nodes of a last regeneration), together and apart.
                    for (int apart : {0, 40, 130}) {
                        random_float(arb_midref(acb_realref(point.value())), random, 100 + random.below(29), 0, 0, false);
                        random_float(arb_midref(acb_imagref(point.value())), random, 100 + random.below(29), 0, -apart, false);
                        check(p, point.value(), precision, what + ", two-limb parts " + std::to_string(apart) + " bits apart");
                    }
                }
}
// Cancellation: a polynomial with known binary64 roots, its coefficients rounded to the precision, at the roots
// and beside them. The exact polynomial lies in the coefficients' balls, so its value (zero at a root) must be
// enclosed, and beside a root the enclosure and the bound hold as anywhere.
void roots() {
    Random random{0x2545f4914f6cdd1dULL};
    for (slong precision : {117, 468, 1872})
        for (slong count : {3, 24}) {
            std::vector<d::ComplexBall> points(static_cast<size_t>(count));
            d::RootVector roots(count);
            for (slong j = 0; j < count; ++j) {
                // Every third root beside the real axis: a multiplier at a limb offset.
                acb_set_d_d(points[static_cast<size_t>(j)].value(), random.signed_unit(),
                            j % 3 ? random.signed_unit() : std::ldexp(random.signed_unit(), -70));
                acb_set(roots[j], points[static_cast<size_t>(j)].value());
            }
            d::ComplexPolynomial exact, rounded;
            acb_poly_product_roots(exact.value(), roots.data(), count, ARF_PREC_EXACT);
            acb_poly_set_round(rounded.value(), exact.value(), precision);
            d::LimbHornerEvaluator evaluator(rounded, precision);
            require(evaluator.usable(), "rounded coefficients were refused");
            d::ComplexBall value, reference, beside;
            for (slong j = 0; j < count; ++j) {
                const std::string what = "precision " + std::to_string(precision) + ", " + std::to_string(count) +
                                         " roots, root " + std::to_string(j);
                require(evaluator.evaluate_value(value.value(), points[static_cast<size_t>(j)].value()), what + ": no value");
                require(acb_contains_zero(value.value()), what + ": the value at a root does not enclose zero");
                // Beside the root by 2^-30 of the real part.
                const double xr = arf_get_d(arb_midref(acb_realref(points[static_cast<size_t>(j)].value())), ARF_RND_DOWN);
                const double xi = arf_get_d(arb_midref(acb_imagref(points[static_cast<size_t>(j)].value())), ARF_RND_DOWN);
                acb_set_d_d(beside.value(), xr + std::ldexp(std::fabs(xr), -30), xi);
                require(evaluator.evaluate_value(value.value(), beside.value()), what + ": no value beside the root");
                _acb_poly_evaluate_horner(reference.value(), exact.value()->coeffs, acb_poly_length(exact.value()),
                                          beside.value(), ARF_PREC_EXACT);
                require(acb_contains(value.value(), reference.value()),
                        what + ": the value beside a root does not enclose the exact one");
                check(rounded, beside.value(), precision, what + ", beside");
            }
        }
}
void exact_cases() {
    // (x - 3)(x + 5)(x - 7) = x^3 - 5 x^2 - 29 x + 105 with integer coefficients: every step is exact.
    d::ComplexPolynomial p;
    for (auto [k, a] : {std::pair<slong, long>{0, 105}, {1, -29}, {2, -5}, {3, 1}}) acb_poly_set_coeff_si(p.value(), k, a);
    d::LimbHornerEvaluator evaluator(p, 117);
    d::ComplexBall point, value, expected;
    for (auto [x, y] : {std::pair<double, long>{3, 0}, {-5, 0}, {7, 0}, {2, 35}, {-0.5, 118}}) {
        acb_set_d(point.value(), x);
        require(evaluator.evaluate_value(value.value(), point.value()), "no value of the integer cubic");
        if (x == -0.5) { // 105 + 14.5 - 1.25 - 0.125
            acb_set_d(expected.value(), 118.125);
        } else
            acb_set_si(expected.value(), y);
        require(acb_is_exact(value.value()) && acb_equal(value.value(), expected.value()),
                "the integer cubic is not exact at " + std::to_string(x));
    }
    // At a complex point the same: (2 + i)^3 - 5 (2 + i)^2 - 29 (2 + i) + 105 = 2 + 11 i - 15 - 20 i - 58 - 29 i + 105.
    acb_set_d_d(point.value(), 2, 1);
    require(evaluator.evaluate_value(value.value(), point.value()), "no value of the integer cubic at 2 + i");
    acb_set_d_d(expected.value(), 34, -38);
    require(acb_is_exact(value.value()) && acb_equal(value.value(), expected.value()),
            "the integer cubic is not exact at 2 + i");
    // And with the multipliers at limb offsets, straddling limbs, and of two limbs, where the window holds
    // every bit of the exact value (three powers of parts up to 200 bits apart, at 936 bits).
    d::ComplexBall wide;
    arf_set_ui_2exp_si(arb_midref(acb_imagref(wide.value())), 5, -103);
    arf_add_ui(arb_midref(acb_realref(wide.value())), arb_midref(acb_imagref(wide.value())), 3, ARF_PREC_EXACT, ARF_RND_DOWN);
    for (auto [xr, xi, name] : {std::tuple<double, double, const char *>{3, 0x1p-70, "a limb offset"},
                                {0x1p-70, -3, "a limb offset in the real part"},
                                {-(1 + 0x1p-52), (1 + 0x1p-52) * 0x1p-28, "a straddling multiplier"},
                                {1 + 0x1p-52, -0x1p-200, "a multiplier moved onto a limb boundary"}}) {
        require(check(p, xr, xi, 936, std::string("the integer cubic at ") + name),
                std::string("the integer cubic is not exact at ") + name);
    }
    require(check(p, wide.value(), 936, "the integer cubic at a two-limb part"),
            "the integer cubic is not exact at a two-limb part");
    // A real polynomial at a real point has a real value.
    Random random{7};
    const auto real = random_polynomial(random, 40, 200, 30, true, 200);
    d::LimbHornerEvaluator real_evaluator(real, 234);
    acb_set_d(point.value(), -1.375);
    require(real_evaluator.evaluate_value(value.value(), point.value()) && arb_is_zero(acb_imagref(value.value())),
            "a real polynomial at a real point has an imaginary part");
}
void refusals() {
    Random random{11};
    const auto p = random_polynomial(random, 8, 53, 0, false, 0);
    d::LimbHornerEvaluator evaluator(p, 117);
    d::ComplexBall point, value;
    acb_set_d_d(value.value(), 5, 6);
    const auto refused = [&](const std::string &what) {
        d::ComplexBall five_six;
        acb_set_d_d(five_six.value(), 5, 6);
        require(!d::LimbHornerEvaluator::serves(point.value()) && !evaluator.evaluate_value(value.value(), point.value()) &&
                    acb_equal(value.value(), five_six.value()),
                what + " was served, or the output was written");
    };
    // An inexact point.
    acb_set_d_d(point.value(), 0.5, 0.25);
    arb_add_error_2exp_si(acb_realref(point.value()), -80);
    refused("an inexact point");
    // A part of more than two limbs: 129 bits.
    acb_set_d_d(point.value(), 0.5, 0.25);
    d::Float tiny;
    arf_one(tiny.value());
    arf_mul_2exp_si(tiny.value(), tiny.value(), -130);
    arf_add(arb_midref(acb_imagref(point.value())), arb_midref(acb_imagref(point.value())), tiny.value(), ARF_PREC_EXACT,
            ARF_RND_DOWN);
    refused("a part of 129 bits");
    // A part of exactly two limbs is served.
    acb_set_d_d(point.value(), 0.5, 0.25);
    arf_mul_2exp_si(tiny.value(), tiny.value(), 1);
    arf_add(arb_midref(acb_imagref(point.value())), arb_midref(acb_imagref(point.value())), tiny.value(), ARF_PREC_EXACT,
            ARF_RND_DOWN);
    require(d::LimbHornerEvaluator::serves(point.value()), "a part of 128 bits was refused");
    // A part outside the admitted range.
    acb_set_d_d(point.value(), std::ldexp(1.0, d::LimbHornerEvaluator::exponent_range + 2), 1);
    refused("a part outside the range");
    // Coefficients: a midpoint longer than the window less a limb; a ball that is not finite.
    const auto longer = random_polynomial(random, 4, 64 * d::LimbHornerEvaluator::limbs_for(64), 0, false, 0);
    require(!d::LimbHornerEvaluator(longer, 64).usable(), "a midpoint as long as the window was taken");
    auto infinite = p;
    arb_pos_inf(acb_realref(infinite.value()->coeffs + 3));
    require(!d::LimbHornerEvaluator(infinite, 117).usable(), "an infinite coefficient was taken");
    require(d::LimbHornerEvaluator::limbs_for(117) == 3 && d::LimbHornerEvaluator::limbs_for(936) == 16 &&
                d::LimbHornerEvaluator::limbs_for(122) == 3 && d::LimbHornerEvaluator::limbs_for(123) == 4,
            "limbs_for is not the least window of the derivation");
}
// The window's limbs are counted in int: limbs_for declines before adding to a precision that could pass slong
// or counting limbs that pass int, and the constructor sizes nothing for a declined precision.
// A window one limb past the largest would be 8 GiB: a bad_alloc fails the test.
void boundaries() {
    using E = d::LimbHornerEvaluator;
    constexpr slong largest = std::numeric_limits<slong>::max();
    require(E::limbs_for(0) == 0 && E::limbs_for(-1) == 0 && E::limbs_for(std::numeric_limits<slong>::min()) == 0,
            "a precision below one bit took a window");
    require(E::limbs_for(E::maximum_precision) == E::maximum_limbs, "the largest precision does not take the largest window");
    require(E::limbs_for(E::maximum_precision - 64) == E::maximum_limbs - 1, "the rule a limb below the largest window");
    require(E::limbs_for(E::maximum_precision + 1) == 0 && E::limbs_for(largest) == 0,
            "a precision past the largest window was taken");
    Random random{13};
    const auto p = random_polynomial(random, 4, 53, 0, false, 0);
    d::ComplexBall point, value, five_six;
    acb_set_d_d(point.value(), 0.5, 0.25);
    acb_set_d_d(five_six.value(), 5, 6);
    for (slong precision : {slong{0}, E::maximum_precision + 1, largest}) {
        d::LimbHornerEvaluator evaluator(p, precision);
        acb_set_d_d(value.value(), 5, 6);
        require(!evaluator.usable() && !evaluator.evaluate_value(value.value(), point.value()) &&
                    acb_equal(value.value(), five_six.value()),
                "a declined precision was served, or the output was written");
    }
}
} // namespace
int main() {
    try {
        boundaries();
        refusals();
        exact_cases();
        roots();
        random_cases();
    } catch (const std::exception &error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
    std::cout << "limb Horner evaluator: ok\n";
    return 0;
}
