#pragma once
// winding_math.hpp's checker on the 0.5 C++ header's types (certpolysolve/certpolysolve.hpp).
#include "winding_math.hpp"
#include "exact_reference_api.hpp"

namespace reference {
inline certpolysolve::Polynomial input(const ComplexPolynomial& p) {
    certpolysolve::Polynomial result;
    for (size_t j = 0; j < p.size(); ++j)
        result.set(j, certpolysolve::Number::complex(certpolysolve::Number(p[j].real.str()), certpolysolve::Number(p[j].imag.str())));
    return result;
}
inline std::optional<int> rectangle_count(const ComplexPolynomial& polynomial, const certpolysolve::Root& root) {
    return rectangle_count(polynomial, endpoint(root.box.real.lo), endpoint(root.box.real.hi),
                           endpoint(root.box.imag.lo), endpoint(root.box.imag.hi));
}
} // namespace reference
