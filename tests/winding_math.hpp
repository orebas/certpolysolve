#pragma once
#include "reference_math.hpp"
#include <array>
#include <optional>

// Test-only exact Sturm--Cauchy checker. No FLINT, floating-point arithmetic,
// numerical roots, Pellet/Graeffe predicates or CertPolySolve API are used in this computation (exact_winding.hpp and
// exact_winding_api.hpp adapt it to the pre-0.5 and the 0.5 C++ header's types).
// Background: M. Eisermann, The Fundamental Theorem of Algebra Made Effective,
// https://pnp.mathematik.uni-stuttgart.de/igt/eiserm/publications/roots.pdf
namespace reference {
using ComplexPolynomial = std::vector<Complex>;

inline Complex operator-(const Complex& a, const Complex& b) { return {a.real - b.real, a.imag - b.imag}; }
inline Complex evaluate(const ComplexPolynomial& p, const Complex& z) {
    Complex value;
    for (auto j = p.rbegin(); j != p.rend(); ++j) value = value * z + *j;
    return value;
}
inline ComplexPolynomial multiply(const ComplexPolynomial& p, const ComplexPolynomial& q) {
    if (p.empty() || q.empty()) return {};
    ComplexPolynomial product(p.size() + q.size() - 1);
    for (size_t i = 0; i < p.size(); ++i)
        for (size_t j = 0; j < q.size(); ++j) product[i + j] = product[i + j] + p[i] * q[j];
    while (!product.empty() && zero(product.back())) product.pop_back();
    return product;
}
inline ComplexPolynomial from_roots(const std::vector<Complex>& roots) {
    ComplexPolynomial polynomial{{1}};
    for (const auto& root : roots) polynomial = multiply(polynomial, {{-root.real, -root.imag}, {1}});
    return polynomial;
}
inline ComplexPolynomial along_segment(const ComplexPolynomial& polynomial, const Complex& from, const Complex& to) {
    auto result = shifted(polynomial, from);
    const Complex direction = to - from;
    Complex power(1);
    for (auto& coefficient : result) {
        coefficient = coefficient * power;
        power = power * direction;
    }
    return result;
}
inline std::pair<Polynomial, Polynomial> components(const ComplexPolynomial& p, int rotation = 0) {
    Polynomial real, imag;
    for (const auto& coefficient : p) {
        real.push_back(coefficient.real - rotation * coefficient.imag);
        imag.push_back(coefficient.imag + rotation * coefficient.real);
    }
    trim(real); trim(imag);
    return {real, imag};
}

// Index of numerator/denominator on an interval with non-pole endpoints,
// with a -infinity -> +infinity jump contributing +1. Even-order poles
// contribute zero. The polynomial part does not affect the index.
inline int cauchy_index(const Polynomial& numerator, const Polynomial& denominator,
                        const Rational& lower = 0, const Rational& upper = 1) {
    require(!denominator.empty() && evaluate(denominator, lower) != 0 && evaluate(denominator, upper) != 0,
            "Cauchy-index endpoint is a pole");
    auto remainder = divide(numerator, denominator).second;
    if (remainder.empty()) return 0;
    make_primitive(remainder);
    std::vector<Polynomial> sequence{denominator, remainder};
    for (;;) {
        remainder = divide(sequence[sequence.size() - 2], sequence.back()).second;
        if (remainder.empty()) break;
        for (auto& coefficient : remainder) coefficient = -coefficient;
        make_primitive(remainder); // a positive multiple: the same signs (reference_math.hpp)
        sequence.push_back(std::move(remainder));
    }
    return sign_variations(sequence, lower) - sign_variations(sequence, upper);
}

// Counts in the actual closed rectangle. For a positive-area rectangle, a
// boundary zero returns nullopt rather than a fractional/incorrect root count.
// Degenerate rectangles are handled as exact points or closed line segments.
inline std::optional<int> rectangle_count(ComplexPolynomial polynomial,
                                          const Rational& real_lower, const Rational& real_upper,
                                          const Rational& imag_lower, const Rational& imag_upper) {
    require(real_lower <= real_upper && imag_lower <= imag_upper, "inverted rectangle");
    while (!polynomial.empty() && zero(polynomial.back())) polynomial.pop_back();
    require(!polynomial.empty(), "the zero polynomial has no finite rectangle count");
    const Complex lower{real_lower, imag_lower}, upper{real_upper, imag_upper};
    if (real_lower == real_upper && imag_lower == imag_upper) {
        const auto expansion = shifted(polynomial, lower);
        int multiplicity = 0;
        while (static_cast<size_t>(multiplicity) < expansion.size() && zero(expansion[static_cast<size_t>(multiplicity)]))
            ++multiplicity;
        return multiplicity;
    }
    if (real_lower == real_upper || imag_lower == imag_upper) {
        auto [real, imag] = components(along_segment(polynomial, lower, upper));
        return roots_closed(gcd(real, imag), 0, 1);
    }
    const std::array<Complex, 4> vertices{{lower, {real_upper, imag_lower}, upper, {real_lower, imag_upper}}};
    std::array<Complex, 4> values;
    for (size_t j = 0; j < vertices.size(); ++j) {
        values[j] = evaluate(polynomial, vertices[j]);
        if (zero(values[j])) return std::nullopt;
    }
    // Multiplication by 1+i*s does not change winding. Each nonzero vertex
    // forbids at most one integer s. Five choices avoid poles at all four
    // vertices, keeping the Cauchy-index endpoint conventions unambiguous.
    int rotation = 0;
    for (; rotation <= 4; ++rotation) {
        bool allowed = true;
        for (const auto& value : values) allowed = allowed && value.imag + rotation * value.real != 0;
        if (allowed) break;
    }
    require(rotation <= 4, "failed to choose a non-pole rotation");
    int twice_winding = 0;
    for (size_t edge = 0; edge < vertices.size(); ++edge) {
        auto [real, imag] = components(along_segment(polynomial, vertices[edge], vertices[(edge + 1) % 4]), rotation);
        // Common real roots of the two components are zeros on this edge.
        if (roots_closed(gcd(real, imag), 0, 1) > 0) return std::nullopt;
        twice_winding += cauchy_index(real, imag);
    }
    require(twice_winding % 2 == 0 && twice_winding >= 0,
            "invalid polynomial winding number (checker bug)");
    const int count = twice_winding / 2;
    require(count <= static_cast<int>(polynomial.size()) - 1, "winding count exceeds polynomial degree");
    return count;
}
} // namespace reference
