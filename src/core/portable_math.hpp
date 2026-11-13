#pragma once
// Elementary functions that give the same bits on every platform, for the values the library's decisions are made
// from (routes, orders, caps, budgets, and the work meter that budgets and hand-offs compare:
// docs/same-routes-every-platform-note-20261021.md). IEEE 754 fixes the result
// of + - * / and sqrt on binary64 from their operands, and frexp is exact; the C libraries' log2, exp, pow and their
// kin are each library's own approximations (glibc, Microsoft's UCRT and Apple's libm differ in last bits), so a
// decision that compares one of their results with a threshold can come out another way on another platform. The
// functions here use only the former, in the order the source gives (every build forbids contraction:
// -ffp-contract=off on GCC and Clang, /fp:strict on MSVC).
#include <array>
#include <cmath>
#include <limits>
#include <numbers>

namespace certpolysolve::detail::portable {
// 1 / (2 j + 1) for j = 0..11, the coefficients of atanh's series, each the correctly rounded binary64 value
// (tests/test_portable_math.cpp checks them against the division).
inline constexpr std::array<double, 12> atanh_series{
    0x1.0000000000000p+0, 0x1.5555555555555p-2, 0x1.999999999999ap-3, 0x1.2492492492492p-3,
    0x1.c71c71c71c71cp-4, 0x1.745d1745d1746p-4, 0x1.3b13b13b13b14p-4, 0x1.1111111111111p-4,
    0x1.e1e1e1e1e1e1ep-5, 0x1.af286bca1af28p-5, 0x1.8618618618618p-5, 0x1.642c8590b2164p-5};
// log2(x) for a finite x > 0; NaN otherwise. With x = m 2^e and m in [sqrt(1/2), sqrt(2)), log m = 2 atanh s for
// s = (m - 1) / (m + 1), |s| <= 3 - 2 sqrt 2 < 0.1716, and atanh s = s (1 + s^2/3 + s^4/5 + ...): the series to
// s^22/23 (terms), whose first dropped term is below 2^-65 of the sum. m - 1 is exact (Sterbenz's lemma); the
// roundings of s, the sum and the last two products leave an error of a few units in the last place, measured
// against the C library's in tests/test_portable_math.cpp. Exact at the powers of two (s = 0).
inline double log2(double x) {
    if (!(x > 0) || !(x < std::numeric_limits<double>::infinity()))
        return std::numeric_limits<double>::quiet_NaN();
    int exponent = 0;
    double m = std::frexp(x, &exponent); // in [1/2, 1)
    if (m < std::numbers::sqrt2 / 2) {
        m *= 2;
        --exponent;
    }
    const double s = (m - 1) / (m + 1), s2 = s * s;
    double sum = 0;
    for (auto k = atanh_series.size(); k-- > 0;)
        sum = sum * s2 + atanh_series[k];
    return static_cast<double>(exponent) + 2 * s * sum * std::numbers::log2e;
}
// (ln 2)^j / j! for j = 0..14, the coefficients of 2^f = e^(f ln 2), each the correctly rounded binary64 value (from
// 80-digit decimal arithmetic; tests/test_portable_math.cpp checks the first against std::numbers::ln2 and every
// ratio of neighbours against ln 2 / j).
inline constexpr std::array<double, 15> exp2_series{
    0x1.0000000000000p+0,  0x1.62e42fefa39efp-1,  0x1.ebfbdff82c58fp-3,  0x1.c6b08d704a0c0p-5,
    0x1.3b2ab6fba4e77p-7,  0x1.5d87fe78a6731p-10, 0x1.430912f86c787p-13, 0x1.ffcbfc588b0c7p-17,
    0x1.62c0223a5c824p-20, 0x1.b5253d395e7c4p-24, 0x1.e4cf5158b8ecap-28, 0x1.e8cac7351bb25p-32,
    0x1.c3bd650fc2986p-36, 0x1.816193166d0f9p-40, 0x1.314964d5878a9p-44};
// 2^x: NaN for NaN, +infinity from 2^1024 up, 0 at and below 2^-1075 (half the least subnormal, which rounds to 0).
// With k = round(x) and f = x - k, |f| <= 1/2 and f is exact (k = 0, or else |x| >= 1/2 and |k|/2 <= |x| <= 2|k|:
// Sterbenz's lemma); 2^f is the series to f^14, whose first dropped term at |f| = 1/2 is below 2^-63; the roundings of
// the coefficients and of Horner's steps leave an error of at most one unit in the last place against glibc's on
// every value tests/test_portable_math.cpp tries; ldexp then scales by 2^k exactly (below the normal range it rounds
// once, as IEEE 754 prescribes). Exact at the integers (f = 0). std::round is exact and ignores the rounding mode.
inline double exp2(double x) {
    if (std::isnan(x))
        return x;
    if (x >= std::numeric_limits<double>::max_exponent)
        return std::numeric_limits<double>::infinity();
    if (x <= std::numeric_limits<double>::min_exponent - std::numeric_limits<double>::digits - 1)
        return 0;
    const double k = std::round(x), f = x - k;
    double sum = 0;
    for (auto j = exp2_series.size(); j-- > 0;)
        sum = sum * f + exp2_series[j];
    return std::ldexp(sum, static_cast<int>(k));
}
} // namespace certpolysolve::detail::portable
