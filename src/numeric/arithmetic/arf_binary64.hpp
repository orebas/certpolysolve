#pragma once
// Binary64 values of arf numbers on the hot paths (the goal test's and the pair sweep's bounds in certification.cpp,
// the compensated certificate, the machine engines' entry in machine_adapter.hpp), and the binary64 neighbours.
//
// arf_to_binary64(x, rnd) returns arf_get_d(x, rnd). arf_get_d rounds through a temporary arf (arf_set_round), and
// through MPFR to nearest, about 200 instructions a call (dev/docs/small-degree-design-pass-20261026.md). Where its
// own arf_set_round path runs, exponents in (-1020, 1020] (the rounded value a normal binary64 number), this reads
// the top limb and whether any lower bit is set, and rounds the 53-bit mantissa as the mode asks; elsewhere, and for
// zero, the infinities and NaN, it calls arf_get_d. A correctly rounded value is unique, so both give the same double
// for every input and mode (tests/test_arf_binary64). next_up and next_down are std::nextafter toward +inf and -inf,
// on the bits.
//
// FLINT's arf layout (ARF_EXPREF, ARF_GET_MPN_READONLY, ARF_SGNBIT) is read only here.
#include <flint/arf.h>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>

namespace certpolysolve::detail {
inline double arf_to_binary64(const arf_t x, arf_rnd_t rnd) {
    constexpr int mantissa_bits = std::numeric_limits<double>::digits;
    constexpr int dropped = FLINT_BITS - mantissa_bits;  // the top limb's bits below a binary64 mantissa
    constexpr int fraction_bits = mantissa_bits - 1;     // the stored fraction of a normal binary64 number
    constexpr int exponent_bias = std::numeric_limits<double>::max_exponent - 1;
    constexpr int sign_bit = std::numeric_limits<uint64_t>::digits - 1;
    // arf_get_d's own range for arf_set_round: below and above it, it asks MPFR.
    constexpr slong lowest_exponent = -1020, highest_exponent = 1020;
    if (arf_is_special(x) || COEFF_IS_MPZ(*ARF_EXPREF(x)))
        return arf_get_d(x, rnd);
    const slong e = *ARF_EXPREF(x);
    if (e <= lowest_exponent || e > highest_exponent)
        return arf_get_d(x, rnd);
    nn_srcptr limbs;
    slong n;
    ARF_GET_MPN_READONLY(limbs, n, x);
    // x = 0.b1 b2 ... 2^e with b1 = 1: the top limb's highest bit is b1.
    const ulong top = limbs[n - 1];
    uint64_t m = top >> dropped;
    const ulong rest = top & ((ulong(1) << dropped) - 1);
    bool lower = false;
    for (slong i = 0; i + 1 < n; ++i)
        lower = lower || limbs[i] != 0;
    const bool negative = ARF_SGNBIT(x);
    const bool inexact = rest != 0 || lower;
    bool up = false;
    switch (rnd) {
    case ARF_RND_DOWN:
        break;
    case ARF_RND_UP:
        up = inexact;
        break;
    case ARF_RND_FLOOR:
        up = inexact && negative;
        break;
    case ARF_RND_CEIL:
        up = inexact && !negative;
        break;
    default: {  // ARF_RND_NEAR: to nearest, ties to even
        const ulong half = ulong(1) << (dropped - 1);
        const bool beyond_half = (rest & (half - 1)) != 0 || lower;
        up = (rest & half) != 0 && (beyond_half || (m & 1) != 0);
    }
    }
    slong exponent = e;
    m += up;
    if (m >> mantissa_bits) {  // rounded up to 2^53: one binade up
        m >>= 1;
        ++exponent;
    }
    // m in [2^52, 2^53): x rounded is m 2^(exponent - 53) = (m 2^-52) 2^(exponent - 1).
    const uint64_t biased = static_cast<uint64_t>(exponent - 1 + exponent_bias);
    const uint64_t fraction = m & ((uint64_t(1) << fraction_bits) - 1);
    return std::bit_cast<double>((uint64_t(negative) << sign_bit) | (biased << fraction_bits) | fraction);
}

// std::nextafter(x, +infinity), on the bits.
inline double next_up(double x) {
    if (std::isnan(x) || x == std::numeric_limits<double>::infinity())
        return x;
    if (x == 0)
        return std::numeric_limits<double>::denorm_min();
    const uint64_t bits = std::bit_cast<uint64_t>(x);
    return std::bit_cast<double>(x > 0 ? bits + 1 : bits - 1);
}
// std::nextafter(x, -infinity), on the bits.
inline double next_down(double x) { return -next_up(-x); }
} // namespace certpolysolve::detail
