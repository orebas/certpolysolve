#pragma once
// Independent exact reference mathematics. This file uses NO FLINT APIs and no CertPolySolve API: the checks on a
// solver's regions take their exact ends as rationals, and exact_reference.hpp (the pre-0.5 C++ header's types) and
// exact_reference_api.hpp (the 0.5 header's) adapt them. It is intentionally slow and small: an oracle for tests,
// not a second solver.
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace reference {
using Integer = boost::multiprecision::cpp_int;
using Rational = boost::multiprecision::cpp_rational;
using Polynomial = std::vector<Rational>;
inline void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
inline Rational absolute(const Rational& x) { return x < 0 ? Rational(-x) : x; }
inline Rational power2(int exponent) {
    if (exponent >= 0) return Rational(Integer(1) << exponent);
    return Rational(1) / (Integer(1) << -exponent);
}
inline void trim(Polynomial& p) { while (!p.empty() && p.back() == 0) p.pop_back(); }
inline int degree(const Polynomial& p) { return static_cast<int>(p.size()) - 1; }
inline Polynomial derivative(const Polynomial& p) {
    Polynomial d;
    for (size_t i = 1; i < p.size(); ++i) d.push_back(p[i] * i);
    trim(d); return d;
}
inline Rational evaluate(const Polynomial& p, const Rational& x) {
    Rational y = 0;
    for (auto i = p.rbegin(); i != p.rend(); ++i) y = y * x + *i;
    return y;
}
inline Polynomial multiply(const Polynomial& p, const Polynomial& q) {
    if (p.empty() || q.empty()) return {};
    Polynomial result(p.size() + q.size() - 1);
    for (size_t i = 0; i < p.size(); ++i)
        for (size_t j = 0; j < q.size(); ++j) result[i + j] += p[i] * q[j];
    trim(result); return result;
}
inline std::pair<Polynomial, Polynomial> divide(Polynomial p, const Polynomial& q) {
    require(!q.empty(), "zero polynomial divisor in reference");
    Polynomial quotient(degree(p) >= degree(q) ? p.size() - q.size() + 1 : 0);
    while (!p.empty() && degree(p) >= degree(q)) {
        const auto k = p.size() - q.size();
        Rational multiplier = p.back() / q.back();
        quotient[k] = multiplier;
        for (size_t j = 0; j < q.size(); ++j) p[k + j] -= multiplier * q[j];
        trim(p);
    }
    return {quotient, p};
}
// A positive multiple of p with coprime integer coefficients. Every value keeps p's sign, so a Sturm or Cauchy-index
// sequence may hold it in p's place (a remainder modulo a positive multiple of the divisor is the same remainder, and a
// positive multiple of the dividend gives a positive multiple of it), and so may a gcd's Euclidean steps (a gcd is
// defined up to a scalar). It keeps the remainder sequences' rationals near the inputs' size: without it they grew at
// every step, and the fuzz families' exact counts at degree 12 to 32 took minutes a region (27 October 2026).
inline void make_primitive(Polynomial& p) {
    if (p.empty()) return;
    Integer lcm = 1, content = 0;
    for (const auto& c : p) lcm = boost::multiprecision::lcm(lcm, boost::multiprecision::denominator(c));
    for (const auto& c : p)
        content = boost::multiprecision::gcd(content, boost::multiprecision::numerator(c) *
                                                          (lcm / boost::multiprecision::denominator(c)));
    if (content == 0) return;
    if (content < 0) content = -content;
    const Rational factor = Rational(lcm) / Rational(content);
    for (auto& c : p) c *= factor;
}
inline Polynomial monic(Polynomial p) {
    if (!p.empty()) { Rational leading = p.back(); for (auto& x : p) x /= leading; }
    return p;
}
inline Polynomial gcd(Polynomial p, Polynomial q) {
    while (!q.empty()) { auto remainder = divide(p, q).second; make_primitive(remainder); p = q; q = remainder; }
    return monic(p);
}
inline std::vector<std::pair<Polynomial, int>> squarefree(Polynomial p) {
    p = monic(p);
    auto repeated = gcd(p, derivative(p));
    auto remaining = divide(p, repeated).first;
    std::vector<std::pair<Polynomial, int>> factors;
    for (int multiplicity = 1; degree(remaining) > 0; ++multiplicity) {
        auto shared = gcd(remaining, repeated);
        auto factor = divide(remaining, shared).first;
        if (degree(factor) > 0) factors.emplace_back(factor, multiplicity);
        remaining = shared;
        repeated = divide(repeated, shared).first;
    }
    return factors;
}
inline int sign_variations(const std::vector<Polynomial>& sequence, const Rational& x) {
    int previous = 0, count = 0;
    for (const auto& p : sequence) {
        const auto value = evaluate(p, x);
        const int sign = (value > 0) - (value < 0);
        if (sign) { if (previous && sign != previous) ++count; previous = sign; }
    }
    return count;
}
inline int distinct_roots_closed(const Polynomial& p, const Rational& lower, const Rational& upper) {
    std::vector<Polynomial> sturm{p, derivative(p)};
    while (!sturm.back().empty()) {
        auto remainder = divide(sturm[sturm.size() - 2], sturm.back()).second;
        for (auto& c : remainder) c = -c;
        if (remainder.empty()) break;
        make_primitive(remainder);
        sturm.push_back(remainder);
    }
    // For squarefree p, V(lower)-V(upper) counts (lower,upper]. Add a root at lower.
    return sign_variations(sturm, lower) - sign_variations(sturm, upper)
        + (evaluate(p, lower) == 0 ? 1 : 0);
}
inline int roots_closed(const Polynomial& p, const Rational& lower, const Rational& upper) {
    int count = 0;
    for (const auto& [factor, multiplicity] : squarefree(p))
        count += multiplicity * distinct_roots_closed(factor, lower, upper);
    return count;
}

struct Complex {
    Rational real = 0, imag = 0;
    Complex() = default;
    Complex(Rational re, Rational im = 0) : real(std::move(re)), imag(std::move(im)) {}
};
inline Complex operator+(const Complex& a, const Complex& b) { return {a.real + b.real, a.imag + b.imag}; }
inline Complex operator*(const Complex& a, const Complex& b) {
    return {a.real * b.real - a.imag * b.imag, a.real * b.imag + a.imag * b.real};
}
inline bool zero(const Complex& z) { return z.real == 0 && z.imag == 0; }
inline std::vector<Complex> shifted(const std::vector<Complex>& polynomial, const Complex& center) {
    auto result = polynomial;
    for (int i = static_cast<int>(result.size()) - 2; i >= 0; --i)
        for (size_t j = static_cast<size_t>(i); j + 1 < result.size(); ++j)
            result[j] = result[j] + center * result[j + 1];
    return result;
}

// A region's exact rectangle and its root count, as the checks below take it.
struct Rectangle {
    Rational real_lower, real_upper, imag_lower, imag_upper;
    uint64_t count = 1;
};
// m 2^e from the decimal mantissa and exponent a result gives.
inline Rational dyadic(const std::string& mantissa, const std::string& exponent) {
    const auto e = std::stoll(exponent);
    require(e >= -1000000 && e <= 1000000, "test oracle refuses huge exponent expansion");
    return Rational(Integer(mantissa)) * power2(static_cast<int>(e));
}
inline bool rouche_inside_rectangle(const std::vector<Complex>& p, const Rectangle& root) {
    const Rational &lo_re = root.real_lower, &hi_re = root.real_upper;
    const Rational &lo_im = root.imag_lower, &hi_im = root.imag_upper;
    const auto expansion = shifted(p, {(lo_re + hi_re) / 2, (lo_im + hi_im) / 2});
    Rational radius = std::min(Rational((hi_re - lo_re) / 2), Rational((hi_im - lo_im) / 2));
    if (radius == 0) {
        if (lo_re != hi_re || lo_im != hi_im || !zero(expansion[0])) return false;
        size_t order = 0;
        while (order < expansion.size() && zero(expansion[order])) ++order;
        return order == root.count;
    }
    // Exact sufficient Pellet/Rouche predicate in a disk CONTAINED in the
    // exported rectangle. max(|re|,|im|) <= |z| <= |re|+|im| avoids square roots.
    Rational others = 0, dominant = 0, radius_power = 1;
    for (size_t j = 0; j < expansion.size(); ++j) {
        const auto& z = expansion[j];
        if (j == root.count) dominant = std::max(absolute(z.real), absolute(z.imag)) * radius_power;
        else others += (absolute(z.real) + absolute(z.imag)) * radius_power;
        radius_power *= radius;
    }
    return dominant > others;
}
// The same exact predicate about the point one exact Newton step from the rectangle's centre, in the largest disk
// about it that the rectangle contains. A closed form's box (Auto's direct route, small_degree.cpp) encloses its root
// tightly but may hold it anywhere inside, where the disk about the centre can miss it; the step brings the centre to
// within about the square of that distance, so the disk holds the root again. A pass proves one root in the rectangle.
inline bool rouche_inside_rectangle_after_newton(const std::vector<Complex>& p, const Rectangle& root) {
    if (root.count != 1) return false;
    const Complex centre{(root.real_lower + root.real_upper) / 2, (root.imag_lower + root.imag_upper) / 2};
    const auto at_centre = shifted(p, centre);
    if (at_centre.size() < 2) return false;
    const Complex &value = at_centre[0], &slope = at_centre[1];
    const Rational norm = slope.real * slope.real + slope.imag * slope.imag;
    if (norm == 0) return false;
    const Complex moved{centre.real - (value.real * slope.real + value.imag * slope.imag) / norm,
                        centre.imag - (value.imag * slope.real - value.real * slope.imag) / norm};
    const Rational radius = std::min({Rational(moved.real - root.real_lower), Rational(root.real_upper - moved.real),
                                      Rational(moved.imag - root.imag_lower), Rational(root.imag_upper - moved.imag)});
    if (radius <= 0) return false;
    const auto expansion = shifted(p, moved);
    Rational others = 0, dominant = 0, radius_power = 1;
    for (size_t j = 0; j < expansion.size(); ++j) {
        const auto& z = expansion[j];
        if (j == 1) dominant = std::max(absolute(z.real), absolute(z.imag)) * radius_power;
        else others += (absolute(z.real) + absolute(z.imag)) * radius_power;
        radius_power *= radius;
    }
    return dominant > others;
}
inline void check_accuracy(const Rectangle& root, int bits) {
    const Rational &a = root.real_lower, &b = root.real_upper;
    const Rational &c = root.imag_lower, &d = root.imag_upper;
    const Rational error = (b - a + d - c) / 2;
    const Rational real_lower = a > 0 ? a : (b < 0 ? Rational(-b) : Rational(0));
    const Rational imag_lower = c > 0 ? c : (d < 0 ? Rational(-d) : Rational(0));
    require(error <= power2(-bits) * std::max(real_lower, imag_lower), "exported accuracy contract failed");
}
inline void check_disjoint(const std::vector<Rectangle>& roots) {
    for (size_t i = 0; i < roots.size(); ++i) for (size_t j = 0; j < i; ++j) {
        const auto& a = roots[i]; const auto& b = roots[j];
        require(a.real_upper < b.real_lower || b.real_upper < a.real_lower
             || a.imag_upper < b.imag_lower || b.imag_upper < a.imag_lower,
                "exported rectangles overlap");
    }
}
} // namespace reference
