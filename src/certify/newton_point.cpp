// Original-source Newton rectangle predicate, independent of candidate repair.
#include "certify/newton_point.hpp"
#include "certify/chart_certificate.hpp"
#include "numeric/arithmetic/modulus_bound.hpp"
#include "numeric/evaluate/monic_source.hpp"
#include "numeric/models/root_model.hpp"

namespace certpolysolve::detail {
bool certify_newton_point(Root &root, MonicSource &source, const acb_t point,
                          const cps_options &options, slong precision) {
    if (!acb_is_finite(point) || acb_is_zero(point) || options.cancelled()) return false;
    RootModel model;
    acb_get_mid(model.center.value(), point);
    upper(model.radius.value(), model.center.value());
    arb_mul_2exp_si(model.radius.value(), model.radius.value(), -static_cast<slong>(options.bits) + 4);
    ComplexBall value, derivative, box = model.center, box_value, box_derivative, difference, coefficient;
    if (!source.evaluate(value.value(), derivative.value(), model.center.value(), precision)) return false;
    arb_add_error(acb_realref(box.value()), model.radius.value());
    arb_add_error(acb_imagref(box.value()), model.radius.value());
    if (!source.evaluate(box_value.value(), box_derivative.value(), box.value(), precision)) return false;
    // The segment from center to any point in its disk stays inside box.
    // Integrating p'(z)-p'(center) bounds the linear Taylor remainder there.
    acb_sub(difference.value(), box_derivative.value(), derivative.value(), precision);
    upper(model.error.value(), difference.value());
    arb_mul(model.error.value(), model.error.value(), model.radius.value(), precision);
    acb_poly_set_coeff_acb(model.polynomial.value(), 0, value.value());
    acb_mul_arb(coefficient.value(), derivative.value(), model.radius.value(), precision);
    acb_poly_set_coeff_acb(model.polynomial.value(), 1, coefficient.value());
    ComplexBall zero;
    if (!certify_root_model(root, model, zero.value(), options, precision)) return false;
    root.certificate = "original-source Newton linear model with derivative-variation remainder";
    return true;
}

} // namespace certpolysolve::detail
