#pragma once
// A binary64 significand with a coarse exponent, for the program route's range tier (program_machine.hpp): the
// value is m 2^(256 e), with m zero (and e = 0) or |m| in [2^-128, 2^128). Proposals only, never certificates.
//
// Why not machine::ScaledDouble (machine_arithmetic.hpp, a frexp-normalised significand with an exact
// exponent): that type pays a frexp and an ldexp on every operation, and a recurrence whose values leave
// binary64's range at every point (Hermite H_1024 near its roots grows like (2x)^k) pays them at every step.
// History: docs/code-history-20261021.md, block_double.hpp 1. Here the window is half a unit wide on either
// side, so
// - a product or quotient of two significands lies in (2^-256, 2^256) and one exact scaling by 2^-+256 brings
//   it back to the window;
// - a sum of two values whose exponents are equal is a plain double sum (one exact scaling if it leaves the
//   window); exponents one apart scale the smaller significand by 2^-256, exactly (it stays at least 2^-384,
//   a normal double); exponents two or more apart drop the smaller, whose ratio to the larger is then below
//   2^-256, far beneath the binary64 rounding of the result.
// Every operation is a fixed sequence of IEEE operations and exact scalings (the same bits on every CPU), and
// the representation is unique (equal values compare equal by their fields).
#include "core/raii.hpp"
#include <cmath>
#include <cstdint>
#include <limits>

namespace certpolysolve::detail::machine {

struct BlockDouble {
    double m = 0;
    int64_t e = 0;
    static constexpr double up = 0x1p256, down = 0x1p-256, high = 0x1p128, low = 0x1p-128;
    // Beyond this exponent (2^(256 x 2^40)) the value is NaN, as ScaledDouble's beyond its range.
    static constexpr int64_t limit = int64_t(1) << 40;

    BlockDouble() = default;
    // Any finite double, normalised (at most a few exact scalings); inf and NaN are kept, not finite.
    BlockDouble(double x) : m(x), e(0) {
        if (x == 0 || !std::isfinite(x)) {
            if (x == 0) m = 0;
            return;
        }
        while (std::abs(m) >= high) {
            m *= down;
            ++e;
        }
        while (std::abs(m) < low) {
            m *= up;
            --e;
        }
    }
    // A significand within one scaling of the window, and its exponent: the result of an operation.
    static BlockDouble normalized(double m, int64_t e) {
        BlockDouble b;
        const double a = std::abs(m);
        if (a >= high) {
            m *= down;
            ++e;
        } else if (a < low) {
            if (m == 0)
                return b; // the zero: e = 0
            m *= up;
            --e;
        }
        if (e > limit || e < -limit)
            m = std::numeric_limits<double>::quiet_NaN();
        b.m = m;
        b.e = e;
        return b;
    }
    explicit operator double() const {
        if (m == 0 || !std::isfinite(m))
            return m;
        if (e > 4)
            return std::copysign(std::numeric_limits<double>::infinity(), m);
        if (e < -5)
            return std::copysign(0.0, m);
        return std::ldexp(m, static_cast<int>(256 * e));
    }
    friend BlockDouble operator-(BlockDouble a) {
        a.m = -a.m;
        return a;
    }
    friend BlockDouble operator+(BlockDouble a, BlockDouble b) {
        if (a.m == 0)
            return b;
        if (b.m == 0)
            return a;
        if (a.e < b.e)
            std::swap(a, b);
        const int64_t gap = a.e - b.e;
        if (gap == 0)
            return normalized(a.m + b.m, a.e);
        if (gap == 1)
            return normalized(a.m + b.m * down, a.e);
        return a; // b below 2^-256 of a
    }
    friend BlockDouble operator-(BlockDouble a, BlockDouble b) { return a + (-b); }
    friend BlockDouble operator*(BlockDouble a, BlockDouble b) {
        if (a.m == 0 || b.m == 0)
            return {};
        return normalized(a.m * b.m, a.e + b.e);
    }
    friend BlockDouble operator/(BlockDouble a, BlockDouble b) {
        if (a.m == 0)
            return {};
        return normalized(a.m / b.m, a.e - b.e);
    }
    friend bool operator==(BlockDouble a, BlockDouble b) { return a.m == b.m && (a.m == 0 || a.e == b.e); }
    friend bool operator<(BlockDouble a, BlockDouble b) {
        // The windows of consecutive exponents do not overlap, so for values of one sign the exponent orders
        // the magnitudes and the significand breaks a tie.
        const bool a_negative = a.m < 0, b_negative = b.m < 0;
        if (a_negative != b_negative)
            return a_negative;
        if (a.m == 0 || b.m == 0)
            return a.m < b.m;
        if (a.e != b.e)
            return a_negative ? a.e > b.e : a.e < b.e;
        return a.m < b.m;
    }
};
inline BlockDouble abs(BlockDouble a) {
    a.m = std::abs(a.m);
    return a;
}
inline bool finite(BlockDouble a) { return std::isfinite(a.m); }
inline BlockDouble sqrt(BlockDouble a) {
    // m 2^(256 e): an even exponent halves; an odd one lends a factor 2^256 to the significand first.
    if (a.m == 0)
        return {};
    if (a.e % 2 == 0)
        return BlockDouble::normalized(std::sqrt(a.m), a.e / 2);
    return BlockDouble::normalized(std::sqrt(a.m * BlockDouble::up), (a.e - 1) / 2);
}
inline double logabs(BlockDouble a) {
    return a.m == 0 ? -std::numeric_limits<double>::infinity()
                    : std::log2(std::abs(a.m)) + 256.0 * static_cast<double>(a.e);
}
// The conversions (in this namespace, found by argument-dependent lookup from machine_adapter.hpp's from_ball
// and to_ball).
inline bool from_arf(BlockDouble &out, const arf_t in) {
    if (!arf_is_finite(in))
        return false;
    if (arf_is_zero(in)) {
        out = {};
        return true;
    }
    // |in| < 2^bound; the exponent written as 256 q + r with r in [-128, 128).
    const slong bound = arf_abs_bound_lt_2exp_si(in);
    // An exponent beyond a word comes back as +-ARF_PREC_EXACT: far outside the range, refused before any sum.
    if (bound > (slong(1) << 50) || bound < -(slong(1) << 50))
        return false;
    const slong shifted = bound + 128;
    const slong q = shifted >= 0 ? shifted / 256 : -((-shifted + 255) / 256);
    if (q > BlockDouble::limit || q < -BlockDouble::limit)
        return false;
    Float scaled;
    arf_mul_2exp_si(scaled.value(), in, -256 * q);
    const double m = arf_get_d(scaled.value(), ARF_RND_NEAR);
    if (!std::isfinite(m) || m == 0)
        return false;
    out = BlockDouble::normalized(m, q);
    return std::isfinite(out.m) && out.m != 0;
}
inline void to_arf(arf_t out, BlockDouble x) {
    arf_set_d(out, x.m);
    arf_mul_2exp_si(out, out, static_cast<slong>(256 * x.e));
}

} // namespace certpolysolve::detail::machine
