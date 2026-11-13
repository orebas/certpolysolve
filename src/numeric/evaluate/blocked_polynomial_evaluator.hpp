#pragma once
#include "core/raii.hpp"
#include <algorithm>
#include <cmath>

namespace certpolysolve::detail {
// Repeated value/derivative evaluation at a fixed precision. Both polynomials
// and all scratch storage are owned. Splitting them into blocks shares powers
// of the argument and uses ball dot products inside each block.
// Values enclose both p(z) and p'(z) for all points/coefficient choices in the
// input balls. Outputs must be distinct. The fixed-precision snapshot cannot
// be reused after changing a local model or its working precision.
class BlockedPolynomialEvaluator {
    ComplexPolynomial polynomial_, derivative_;
    slong precision_, block_;
    RootVector powers_;

    void evaluate_polynomial(acb_t out, const ComplexPolynomial &polynomial) {
        const slong length = acb_poly_length(polynomial.value());
        acb_zero(out);
        if (!length) return;
        // p(z) = P_0(z) + z^b P_1(z) + z^(2b) P_2(z) + ...,
        // with degree(P_k) < b. Horner's rule combines these short blocks.
        slong offset = ((length - 1) / block_) * block_;
        acb_dot(out, nullptr, 0, polynomial.value()->coeffs + offset, 1,
                powers_.data(), 1, length - offset, precision_);
        while (offset >= block_) {
            offset -= block_;
            acb_mul(out, out, powers_[block_], precision_);
            // acb_dot includes the initial ball (and its radius) in the
            // accumulation, avoiding a separately rounded block sum/addition.
            acb_dot(out, out, 0, polynomial.value()->coeffs + offset, 1,
                    powers_.data(), 1, block_, precision_);
        }
    }

  public:
    BlockedPolynomialEvaluator(const ComplexPolynomial &polynomial, slong precision)
        : polynomial_(polynomial), precision_(precision),
          block_(std::max<slong>(1, static_cast<slong>(std::sqrt(
              static_cast<double>(acb_poly_length(polynomial.value())))) + 1)),
          powers_(block_ + 1) {
        acb_poly_derivative(derivative_.value(), polynomial_.value(), precision_);
    }
    void evaluate(acb_t value, acb_t derivative, const acb_t point) {
        _acb_vec_set_powers(powers_.data(), point, powers_.size(), precision_);
        evaluate_polynomial(value, polynomial_);
        evaluate_polynomial(derivative, derivative_);
    }
};
} // namespace certpolysolve::detail
