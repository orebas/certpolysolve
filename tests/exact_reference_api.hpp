#pragma once
// reference_math.hpp's checks on the 0.5 C++ header's types (certpolysolve/certpolysolve.hpp): a region's exact ends
// (Dyadic), its count, and an exact input built from rationals (each coefficient's exact text, never a double).
#include "reference_math.hpp"
#include <certpolysolve/certpolysolve.hpp>

namespace reference {
inline Rational endpoint(const certpolysolve::Dyadic& x) { return dyadic(x.mantissa(), x.exponent()); }
inline Rectangle rectangle(const certpolysolve::Root& root) {
    return {endpoint(root.box.real.lo), endpoint(root.box.real.hi), endpoint(root.box.imag.lo),
            endpoint(root.box.imag.hi), root.count};
}
inline bool rouche_inside_rectangle(const std::vector<Complex>& p, const certpolysolve::Root& root) {
    return rouche_inside_rectangle(p, rectangle(root));
}
inline bool rouche_inside_rectangle_after_newton(const std::vector<Complex>& p, const certpolysolve::Root& root) {
    return rouche_inside_rectangle_after_newton(p, rectangle(root));
}
inline void check_accuracy(const certpolysolve::Root& root, int bits) { check_accuracy(rectangle(root), bits); }
inline void check_disjoint(const std::vector<certpolysolve::Root>& roots) {
    std::vector<Rectangle> rectangles;
    for (const auto& root : roots) rectangles.push_back(rectangle(root));
    check_disjoint(rectangles);
}
inline certpolysolve::Polynomial input(const Polynomial& p) {
    certpolysolve::Polynomial result;
    for (size_t j = 0; j < p.size(); ++j) result.set(j, certpolysolve::Number(p[j].str()));
    return result;
}
} // namespace reference
