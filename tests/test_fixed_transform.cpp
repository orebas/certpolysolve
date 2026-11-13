// The fixed-point transform (fixed_transform.hpp) against FLINT's ball transform at a far higher precision.
// Its balls become model coefficients that certificates rest on, so the property tested is enclosure: every
// output rectangle must contain the exact transform of every point of the input balls, with the bound the header
// derives and no more. The inputs are random with mixed sizes and the extreme cases of the derivation: components
// at the top of the range (no overflow), a single frequency (one output takes the whole growth), zeros.
#include "numeric/transforms/fixed_transform.hpp"
#include "core/resource_limits.hpp"
#include <flint/acb_dft.h>
#include <algorithm>
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
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return state >> 11;
    }
    // A dyadic number with `bits` random mantissa bits and the given exponent of two, of either sign.
    void dyadic(arf_t x, slong bits, slong exponent) {
        arf_zero(x);
        d::Float limb;
        for (slong done = 0; done < bits; done += 32) {
            arf_mul_2exp_si(x, x, 32);
            arf_set_ui(limb.value(), static_cast<ulong>(next() & 0xffffffffu));
            arf_add(x, x, limb.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
        }
        arf_mul_2exp_si(x, x, exponent - bits);
        if (next() & 1) arf_neg(x, x);
    }
};
// The exact transform of the midpoints (FLINT's ball transform on exact inputs at `reference` bits: its balls are
// 2^-reference relative, far below anything compared here) and rho, the sum of the inputs' radii.
struct Reference {
    d::RootVector values;
    d::RealBall radii;
    Reference(const d::RootVector &input, int exponent, slong reference) : values(input.size()) {
        d::RootVector midpoints(input.size());
        d::RealBall radius;
        for (slong j = 0; j < input.size(); ++j) {
            acb_get_mid(midpoints[j], input[j]);
            arb_get_rad_arb(radius.value(), acb_realref(input[j]));
            arb_add(radii.value(), radii.value(), radius.value(), reference);
            arb_get_rad_arb(radius.value(), acb_imagref(input[j]));
            arb_add(radii.value(), radii.value(), radius.value(), reference);
        }
        acb_dft_rad2_t plan;
        acb_dft_rad2_init(plan, exponent, reference);
        acb_dft_rad2_precomp(values.data(), midpoints.data(), plan, reference);
        acb_dft_rad2_clear(plan);
    }
};
// Every output rectangle holds the disk of radius rho about the exact transform of the midpoints, and is no wider
// than the derived bound allows: the arithmetic part for inputs below 2^e, rho, and the rounding of the midpoint
// to `precision` bits.
void check(const d::RootVector &input, int exponent, int limbs, slong precision, const std::string &what) {
    const slong reference = 64 * static_cast<slong>(limbs) + 256;
    d::FixedTransform transform(exponent);
    d::FixedTransformStatistics stats;
    const auto out = transform.apply_with(input, precision, limbs, stats);
    require(out.has_value() && stats.accepted == 1 && stats.limbs == static_cast<uint64_t>(limbs),
            what + ": the transform declined");
    const Reference exact(input, exponent, reference);
    slong size = -ARF_PREC_EXACT;
    for (slong j = 0; j < input.size(); ++j)
        for (const arb_struct *x : {acb_realref(input[j]), acb_imagref(input[j])})
            if (!arf_is_zero(arb_midref(x))) size = std::max(size, arf_abs_bound_lt_2exp_si(arb_midref(x)));
    if (size == -ARF_PREC_EXACT) size = 0;
    d::Magnitude arithmetic;
    d::FixedTransform::arithmetic_bound(arithmetic.value(), size, exponent, limbs);
    d::RealBall allowed, distance, radius, slack;
    d::Float bound;
    arf_set_mag(bound.value(), arithmetic.value());
    arb_set_arf(allowed.value(), bound.value());
    arb_add(allowed.value(), allowed.value(), exact.radii.value(), reference);
    for (slong k = 0; k < input.size(); ++k)
        for (int part = 0; part < 2; ++part) {
            const arb_struct *got = part ? acb_imagref(out->values[k]) : acb_realref(out->values[k]);
            const arb_struct *want = part ? acb_imagref(exact.values[k]) : acb_realref(exact.values[k]);
            require(arb_is_finite(got), what + ": an output is not finite");
            // |midpoint - exact| + rho <= radius: the disk of the true outputs lies inside the rectangle.
            arb_set_arf(distance.value(), arb_midref(got));
            arb_sub(distance.value(), distance.value(), want, reference);
            arb_abs(distance.value(), distance.value());
            arb_add(distance.value(), distance.value(), exact.radii.value(), reference);
            arb_get_rad_arb(radius.value(), got);
            require(arb_le(distance.value(), radius.value()),
                    what + ": output " + std::to_string(k) + " does not enclose the exact transform");
            // The radius is the derived bound and the midpoint's rounding to the working precision (at most
            // 2^-precision of it); a thousandth more for the upward rounding of each magnitude addition.
            arb_set_arf(slack.value(), arb_midref(got));
            arb_abs(slack.value(), slack.value());
            arb_mul_2exp_si(slack.value(), slack.value(), 1 - precision);
            arb_add(slack.value(), slack.value(), allowed.value(), reference);
            arb_mul_2exp_si(distance.value(), allowed.value(), -10);
            arb_add(slack.value(), slack.value(), distance.value(), reference);
            require(arb_le(radius.value(), slack.value()),
                    what + ": output " + std::to_string(k) + " is wider than the derived bound");
        }
}
// Random inputs: mantissas of `bits` bits, sizes spread over `spread` binary orders below 2^top, a share of
// zeros, and radii of 2^-radius relative where radius is positive.
void fill(d::RootVector &input, Random &random, slong bits, slong top, slong spread, slong radius) {
    for (slong j = 0; j < input.size(); ++j)
        for (arb_struct *x : {acb_realref(input[j]), acb_imagref(input[j])}) {
            arb_zero(x);
            if (random.next() % 8 == 0) continue;
            const slong exponent = top - static_cast<slong>(random.next() % static_cast<uint64_t>(spread + 1));
            random.dyadic(arb_midref(x), bits, exponent);
            if (radius > 0 && random.next() % 2) arb_add_error_2exp_si(x, exponent - radius);
        }
}
void random_inputs() {
    Random random{20261008};
    for (int exponent = 1; exponent <= 11; ++exponent)
        for (int limbs = d::FixedTransform::minimum_limbs; limbs <= d::FixedTransform::maximum_limbs; ++limbs) {
            if (exponent > 8 && limbs > 5) continue; // the long words are the same code; keep the suite quick
            d::RootVector input(slong{1} << exponent);
            const slong precision = 64 * static_cast<slong>(limbs) - 40;
            const std::string name = "length 2^" + std::to_string(exponent) + ", " + std::to_string(limbs) + " limbs";
            // The same order of size everywhere; sizes spread far beyond the word; a huge and a tiny scale.
            fill(input, random, 150, 0, 2, 0);
            check(input, exponent, limbs, precision, name + ", flat");
            fill(input, random, 200, 37, 400, 0);
            check(input, exponent, limbs, precision, name + ", spread");
            fill(input, random, 64, 100000, 90, 170);
            check(input, exponent, limbs, precision, name + ", large with radii");
            fill(input, random, 300, -200000, 30, 250);
            check(input, exponent, limbs, precision, name + ", small with radii");
        }
}
// The extremes of the derivation. Every component just below 2^e with one sign: output zero is N times the
// input, the largest sum the headroom allows for. A single frequency with every component near the top: one
// output takes the modulus N sqrt(2) 2^e, through every level's rotation.
void extreme_inputs() {
    for (int exponent : {1, 2, 5, 9, 12})
        for (int limbs : {2, 3, 4}) {
            const slong N = slong{1} << exponent, precision = 64 * static_cast<slong>(limbs) - 40;
            const std::string name = "length 2^" + std::to_string(exponent) + ", " + std::to_string(limbs) + " limbs";
            d::RootVector input(N);
            for (slong j = 0; j < N; ++j) {
                // 1 - 2^-400: below one, and every bit of the word set after truncation.
                arb_one(acb_realref(input[j]));
                arb_set(acb_imagref(input[j]), acb_realref(input[j]));
                d::Float tiny;
                arf_one(tiny.value());
                arf_mul_2exp_si(tiny.value(), tiny.value(), -400);
                arf_sub(arb_midref(acb_realref(input[j])), arb_midref(acb_realref(input[j])), tiny.value(),
                        ARF_PREC_EXACT, ARF_RND_DOWN);
                arf_set(arb_midref(acb_imagref(input[j])), arb_midref(acb_realref(input[j])));
            }
            check(input, exponent, limbs, precision, name + ", every component at the top");
            for (slong j = 0; j < N; ++j)
                if (j & 1) acb_neg(input[j], input[j]);
            check(input, exponent, limbs, precision, name + ", alternating signs");
            // x_j = (1 + i)(1 - 2^-30) exp(2 pi i j f / N): output f is N (1 + i)(1 - 2^-30) up to rounding.
            d::RootVector roots(N);
            _acb_vec_unit_roots(roots.data(), N, N, 64 * static_cast<slong>(limbs) + 64);
            const slong frequency = N > 2 ? 3 % N : 1;
            for (slong j = 0; j < N; ++j) {
                d::ComplexBall value;
                acb_set_d_d(value.value(), 1 - 0x1p-30, 1 - 0x1p-30);
                acb_mul(value.value(), value.value(), roots[(j * frequency) % N], 64 * static_cast<slong>(limbs) + 64);
                acb_mul_2exp_si(value.value(), value.value(), -1); // each component of the product below one
                acb_get_mid(input[j], value.value());
            }
            check(input, exponent, limbs, precision, name + ", one frequency");
        }
}
void word_size_rule() {
    using T = d::FixedTransform;
    // 64 K >= precision + 2 L + headroom_bits + ceil(log2(error_units_per_point)) = precision + 2 L + 6.
    require(T::limbs_for(149, 15) == 3, "149 bits at length 2^15 take three limbs (185 bits)");
    require(T::limbs_for(181, 15) == 4, "181 bits at length 2^15 take four limbs (217 bits)");
    require(T::limbs_for(156, 15) == 3 && T::limbs_for(157, 15) == 4, "the rule's edge at three limbs");
    require(T::limbs_for(53, 4) == T::minimum_limbs, "a small request takes the smallest word");
    require(T::limbs_for(64 * T::maximum_limbs, 15) == 0, "a request beyond the largest word is declined");
    // With that word the arithmetic part of the bound is at most 2^(e - precision).
    for (slong precision : {64, 128, 149, 181, 298, 400})
        for (int exponent : {1, 7, 15, 20}) {
            const int limbs = T::limbs_for(precision, exponent);
            if (!limbs) continue;
            d::Magnitude bound;
            T::arithmetic_bound(bound.value(), 0, exponent, limbs);
            require(mag_cmp_2exp_si(bound.value(), -precision) <= 0, "the word chosen misses the working precision");
            if (limbs > T::minimum_limbs) {
                T::arithmetic_bound(bound.value(), 0, exponent, limbs - 1);
                require(mag_cmp_2exp_si(bound.value(), -precision) > 0, "a smaller word would have done");
            }
        }
}
// The constructor's range is checked before anything is shifted by the exponent or sized from it, and the
// word-size rule declines before adding to a precision that could pass slong. The exponents
// rejected include one whose table of positions could not be allocated (63) and ones whose shift would be
// undefined (-1, 64): a bad_alloc, or anything but the invalid argument, fails the test. The largest exponent
// accepted, 31, takes an 8 GiB table and is tested through the helper alone.
void boundaries() {
    using T = d::FixedTransform;
    for (int exponent : {-1, 0, 32, 63, 64}) {
        bool rejected = false;
        try {
            T transform(exponent);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "exponent " + std::to_string(exponent) + " was not rejected as an invalid argument");
    }
    // The helper at the exponents' edges: 53 + 2 L + 6 bits take two limbs at L = 31 and 32 (121 and 123 bits),
    // three at 63 and 64 (185 and 187), and none at the largest int (the exponent's range is the constructor's).
    require(T::limbs_for(53, -1) == 0 && T::limbs_for(53, 0) == 0, "a length below two took a word");
    require(T::limbs_for(53, 31) == 2 && T::limbs_for(53, 32) == 2 && T::limbs_for(53, 63) == 3 &&
                T::limbs_for(53, 64) == 3,
            "the rule at the constructor's edge");
    require(T::limbs_for(53, std::numeric_limits<int>::max()) == 0, "a length beyond every word was not declined");
    // The maximum signed precision, FLINT's largest working precision, and precisions below one bit.
    constexpr slong largest = std::numeric_limits<slong>::max();
    require(T::limbs_for(largest, 1) == 0 && T::limbs_for(largest, 31) == 0 &&
                T::limbs_for(largest, std::numeric_limits<int>::max()) == 0,
            "the maximum signed precision was not declined");
    require(T::limbs_for(static_cast<slong>(d::resources::precision_max), 15) == 0,
            "the largest working precision was not declined");
    require(T::limbs_for(0, 5) == 0 && T::limbs_for(-1, 5) == 0 && T::limbs_for(std::numeric_limits<slong>::min(), 5) == 0,
            "a precision below one bit took a word");
    // The largest precision served, at length two: 504 + 2 + 6 = 512 bits fill the largest word; one more bit
    // needs a ninth limb. At and past 64 maximum_limbs the rule declines by its formula and by its guard.
    require(T::limbs_for(504, 1) == T::maximum_limbs && T::limbs_for(505, 1) == 0, "the edge of the largest word");
    require(T::limbs_for(64 * T::maximum_limbs, 1) == 0 && T::limbs_for(64 * T::maximum_limbs + 1, 1) == 0,
            "a precision at or past the largest word's bits took a word");
}
void declined_inputs() {
    d::FixedTransform transform(3);
    d::FixedTransformStatistics stats;
    d::RootVector input(8);
    // All zero: exact zeros, with the radii alone.
    arb_add_error_2exp_si(acb_realref(input[5]), -70);
    auto out = transform.apply(input, 128, stats);
    require(out.has_value(), "an all-zero input was declined");
    for (slong k = 0; k < 8; ++k) {
        require(acb_contains_zero(out->values[k]), "the transform of zeros does not contain zero");
        d::RealBall radius;
        arb_get_rad_arb(radius.value(), acb_realref(out->values[k]));
        arb_mul_2exp_si(radius.value(), radius.value(), 69);
        d::RealBall one;
        arb_one(one.value());
        require(arb_le(radius.value(), one.value()), "the transform of zeros is wider than its inputs' radii allow");
    }
    arb_pos_inf(acb_imagref(input[2]));
    require(!transform.apply(input, 128, stats).has_value() && stats.declined == 1, "an infinite input was accepted");
    arb_zero(acb_imagref(input[2]));
    require(!transform.apply(input, 64 * d::FixedTransform::maximum_limbs, stats).has_value() && stats.declined == 2,
            "a precision beyond the largest word was accepted");
}
// apply_into publishes the outputs asked for, each the ball apply() gives, and leaves the others alone.
void wanted_outputs() {
    Random random{20261009};
    constexpr int exponent = 6;
    constexpr slong length = slong{1} << exponent, precision = 128;
    d::RootVector input(length), partial(length);
    fill(input, random, 150, 3, 40, 160);
    d::FixedTransform transform(exponent);
    d::FixedTransformStatistics stats;
    const auto whole = transform.apply(input, precision, stats);
    require(whole.has_value(), "the transform declined");
    const std::vector<slong> wanted{0, 5, 31, 63};
    for (slong k = 0; k < length; ++k) acb_set_si(partial[k], 17);
    require(transform.apply_into(partial, input, precision, wanted, stats) && stats.accepted == 2,
            "the transform declined the outputs asked for");
    for (slong k = 0; k < length; ++k) {
        const bool asked = std::find(wanted.begin(), wanted.end(), k) != wanted.end();
        d::ComplexBall untouched;
        acb_set_si(untouched.value(), 17);
        require(acb_equal(partial[k], asked ? whole->values[k] : untouched.value()),
                asked ? "an output asked for differs from the whole transform's" : "an output not asked for was written");
    }
}
} // namespace
int main() {
    try {
        word_size_rule();
        boundaries();
        wanted_outputs();
        declined_inputs();
        extreme_inputs();
        random_inputs();
    } catch (const std::exception &error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
    std::cout << "fixed transform: ok\n";
    return 0;
}
