#pragma once
#include "numeric/evaluate/evaluator.hpp"
#include "numeric/evaluate/polynomial_remainder.hpp"

namespace certpolysolve::detail {
// Private original-coefficient finisher. An immutable snapshot owns every
// coefficient. Disk errors avoid the repeated rectangle rotation of Horner;
// they are expanded back to outward rectangles only at the API boundary.
class DiskPolynomialEvaluator {
    PolynomialSnapshot coefficients_;
    PolynomialLinearRemainder remainder_;
    slong precision_;
    static void recenter(acb_t value, mag_t error) {
        Magnitude local;
        mag_add(local.value(), arb_radref(acb_realref(value)), arb_radref(acb_imagref(value)));
        mag_add(error, error, local.value());
        acb_get_mid(value, value);
    }

  public:
    DiskPolynomialEvaluator(PolynomialSnapshot coefficients, slong precision)
      : coefficients_(std::move(coefficients))
      , remainder_(coefficients_.coefficients())
      , precision_(precision) {}
    bool evaluate(acb_t value, acb_t derivative, const acb_t point) const {
        if (value == derivative || value == point || derivative == point) return false;
        if (!acb_is_exact(point) || !acb_is_finite(point)) return false;
        acb_zero(value);
        acb_zero(derivative);
        Magnitude norm, error, derivative_error;
        acb_get_mag(norm.value(), point);
        const auto values = coefficients_.coefficients();
        for (slong j = values.size(); j-- > 0;) {
            acb_mul(derivative, derivative, point, precision_);
            acb_add(derivative, derivative, value, precision_);
            mag_mul(derivative_error.value(), derivative_error.value(), norm.value());
            mag_add(derivative_error.value(), derivative_error.value(), error.value());
            recenter(derivative, derivative_error.value());
            acb_mul(value, value, point, precision_);
            acb_add(value, value, values[j], precision_);
            mag_mul(error.value(), error.value(), norm.value());
            recenter(value, error.value());
        }
        acb_add_error_mag(value, error.value());
        acb_add_error_mag(derivative, derivative_error.value());
        return acb_is_finite(value) && acb_is_finite(derivative);
    }
    void linear_remainder(arb_t out, const acb_t center, const arb_t radius) const {
        remainder_.linear_remainder(out, center, radius);
    }
};
} // namespace certpolysolve::detail
