#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>

namespace certpolysolve::detail::machine {
// Significand bits of the candidate arithmetic tiers: binary64 (and the
// binary64 mantissa of ScaledDouble), and the unevaluated double-double sum.
inline constexpr int binary64_bits = std::numeric_limits<double>::digits;
inline constexpr int double_double_bits = 2 * binary64_bits;
static_assert(binary64_bits == 53);
// Candidate arithmetic only. Round-to-nearest, no reassociation; explicit FMA
// recovers product roundoff. This is double-double (~106 bits), not binary128.
struct DoubleDouble {
    double hi = 0, lo = 0;
    DoubleDouble() = default;
    DoubleDouble(double x) : hi(x) {}
    DoubleDouble(double h, double l) : hi(h), lo(l) {}
    explicit operator double() const { return hi + lo; }
    static DoubleDouble sum(double a, double b) {
        const double s = a + b, v = s - a;
        return {s, (a - (s - v)) + (b - v)};
    }
    friend DoubleDouble operator-(DoubleDouble a) { return {-a.hi, -a.lo}; }
    friend DoubleDouble operator+(DoubleDouble a, DoubleDouble b) {
        auto s = sum(a.hi, b.hi), t = sum(a.lo, b.lo);
        s = sum(s.hi, s.lo + t.hi);
        return sum(s.hi, s.lo + t.lo);
    }
    friend DoubleDouble operator-(DoubleDouble a, DoubleDouble b) { return a + (-b); }
    friend DoubleDouble operator*(DoubleDouble a, DoubleDouble b) {
        const double p = a.hi * b.hi;
        const double e = std::fma(a.hi, b.hi, -p) + a.hi * b.lo + a.lo * b.hi;
        auto s = sum(p, e);
        return sum(s.hi, s.lo + a.lo * b.lo);
    }
    friend DoubleDouble operator/(DoubleDouble a, DoubleDouble b) {
        const double q = a.hi / b.hi;
        auto r = a - b * DoubleDouble(q);
        const double q2 = r.hi / b.hi;
        r = r - b * DoubleDouble(q2);
        return DoubleDouble(q) + DoubleDouble(q2) + DoubleDouble(r.hi / b.hi);
    }
    friend bool operator==(DoubleDouble a, DoubleDouble b) { return a.hi == b.hi && a.lo == b.lo; }
    friend bool operator<(DoubleDouble a, DoubleDouble b) {
        return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo);
    }
};
inline DoubleDouble sqrt(DoubleDouble a) {
    if (a.hi == 0)
        return {};
    DoubleDouble x(std::sqrt(a.hi));
    return x + (a - x * x) / (DoubleDouble(2) * x);
}
inline DoubleDouble abs(DoubleDouble a) { return a < DoubleDouble(0) ? -a : a; }
inline bool finite(DoubleDouble a) { return std::isfinite(a.hi) && std::isfinite(a.lo); }
inline double logabs(DoubleDouble a) {
    return a.hi == 0 ? -std::numeric_limits<double>::infinity() : std::log2(std::abs(a.hi));
}

// Binary64 significand with a separate exponent. Exponent arithmetic is
// checked; underflow in alignment drops a <2^-1074 relative contribution only.
// It does not turn a nonzero coefficient/root into an exact zero on conversion.
struct ScaledDouble {
    double mantissa = 0;
    int64_t exponent = 0;
    ScaledDouble() = default;
    ScaledDouble(double x) {
        int e = 0;
        mantissa = std::frexp(x, &e);
        exponent = e;
    }
    ScaledDouble(double m, int64_t e) {
        int shift = 0;
        mantissa = std::frexp(m, &shift);
        constexpr int64_t limit = int64_t(1) << 60;
        if (e < -limit || e > limit) {
            mantissa = std::numeric_limits<double>::quiet_NaN();
            return;
        }
        exponent = mantissa == 0 ? 0 : e + shift;
    }
    explicit operator double() const {
        if (exponent > 1024)
            return std::copysign(std::numeric_limits<double>::infinity(), mantissa);
        if (exponent < -1074)
            return std::copysign(0., mantissa);
        return std::ldexp(mantissa, static_cast<int>(exponent));
    }
    friend ScaledDouble operator-(ScaledDouble a) { return {-a.mantissa, a.exponent}; }
    friend ScaledDouble operator+(ScaledDouble a, ScaledDouble b) {
        if (a.mantissa == 0)
            return b;
        if (b.mantissa == 0)
            return a;
        if (a.exponent < b.exponent)
            std::swap(a, b);
        const auto gap = a.exponent - b.exponent;
        return {a.mantissa + (gap > 1074 ? 0. : std::ldexp(b.mantissa, -static_cast<int>(gap))),
                a.exponent};
    }
    friend ScaledDouble operator-(ScaledDouble a, ScaledDouble b) { return a + (-b); }
    friend ScaledDouble operator*(ScaledDouble a, ScaledDouble b) {
        return {a.mantissa * b.mantissa, a.exponent + b.exponent};
    }
    friend ScaledDouble operator/(ScaledDouble a, ScaledDouble b) {
        return {a.mantissa / b.mantissa, a.exponent - b.exponent};
    }
    friend bool operator==(ScaledDouble a, ScaledDouble b) {
        return a.mantissa == b.mantissa && a.exponent == b.exponent;
    }
    friend bool operator<(ScaledDouble a, ScaledDouble b) {
        if (a.mantissa < 0 && b.mantissa >= 0)
            return true;
        if (a.mantissa >= 0 && b.mantissa < 0)
            return false;
        if (a.mantissa == 0)
            return b.mantissa > 0;
        if (b.mantissa == 0)
            return a.mantissa < 0;
        if (a.exponent == b.exponent)
            return a.mantissa < b.mantissa;
        return a.mantissa > 0 ? a.exponent < b.exponent : a.exponent > b.exponent;
    }
};
inline ScaledDouble sqrt(ScaledDouble a) {
    if (a.exponent % 2)
        return {std::sqrt(a.mantissa * 2), (a.exponent - 1) / 2};
    return {std::sqrt(a.mantissa), a.exponent / 2};
}
inline ScaledDouble abs(ScaledDouble a) { return {std::abs(a.mantissa), a.exponent}; }
inline bool finite(ScaledDouble a) { return std::isfinite(a.mantissa); }
inline double logabs(ScaledDouble a) {
    return a.mantissa == 0 ? -std::numeric_limits<double>::infinity()
                           : std::log2(std::abs(a.mantissa)) + static_cast<double>(a.exponent);
}
inline bool finite(double x) { return std::isfinite(x); }
inline double logabs(double x) {
    return x == 0 ? -std::numeric_limits<double>::infinity() : std::log2(std::abs(x));
}
using std::abs;
using std::sqrt;

template <class R> struct Complex {
    R re{}, im{};
    Complex() = default;
    Complex(R r, R i = R{}) : re(r), im(i) {}
    friend Complex operator-(Complex a) { return {-a.re, -a.im}; }
    friend Complex operator+(Complex a, Complex b) { return {a.re + b.re, a.im + b.im}; }
    friend Complex operator-(Complex a, Complex b) { return a + (-b); }
    friend Complex operator*(Complex a, Complex b) {
        return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
    }
    friend Complex operator/(Complex a, Complex b) {
        // Smith's division avoids squaring the larger component.
        if (abs(b.im) < abs(b.re)) {
            const R r = b.im / b.re, d = b.re + b.im * r;
            return {(a.re + a.im * r) / d, (a.im - a.re * r) / d};
        }
        const R r = b.re / b.im, d = b.im + b.re * r;
        return {(a.re * r + a.im) / d, (a.im * r - a.re) / d};
    }
    friend bool operator==(Complex a, Complex b) { return a.re == b.re && a.im == b.im; }
    Complex &operator+=(Complex b) { return *this = *this + b; }
    Complex &operator-=(Complex b) { return *this = *this - b; }
    Complex &operator*=(Complex b) { return *this = *this * b; }
    Complex &operator/=(Complex b) { return *this = *this / b; }
};
template <class R> bool finite(Complex<R> z) { return finite(z.re) && finite(z.im); }
template <class R> Complex<R> conj(Complex<R> z) { return {z.re, -z.im}; }
template <class R> R norm(Complex<R> z) {
    const R a = abs(z.re), b = abs(z.im), m = a < b ? b : a;
    if (m == R{})
        return {};
    const R small = (a < b ? a : b) / m;
    if (small < R(0x1p-500))
        return m;
    return m * sqrt((a / m) * (a / m) + (b / m) * (b / m));
}
template <class R> double logabs(Complex<R> z) { return logabs(norm(z)); }
template <class R> Complex<R> sqrt(Complex<R> z) {
    if (z == Complex<R>{})
        return {};
    const R a = sqrt((norm(z) + abs(z.re)) / R(2));
    if (!(z.re < R{}))
        return {a, z.im / (R(2) * a)};
    const R b = z.im < R{} ? -a : a;
    return {z.im / (R(2) * b), b};
}
} // namespace certpolysolve::detail::machine
