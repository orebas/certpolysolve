#pragma once
#include "methods/piecewise/piecewise_roots.hpp"

namespace owned_refinement {
namespace d = certpolysolve::detail;

// The residual chooses a proposed radius, never a root count. The unchanged
// published-rectangle Rouche test proves the result using the linear jet and
// its integral Taylor remainder (including every omitted model term).
inline bool residual_guard(d::Root &root, const d::RootModel &model, const d::ComplexBall &point,
                           const cps_options &options, slong precision) {
    d::ComplexBall center, local, value, derivative;
    acb_mul_arb(center.value(), point.value(), model.radius.value(), precision);
    acb_add(center.value(), center.value(), model.center.value(), precision);
    acb_get_mid(center.value(), center.value());
    if (!acb_is_finite(center.value()) || !arb_is_positive(model.radius.value())) return false;
    acb_sub(local.value(), center.value(), model.center.value(), precision);
    acb_div_arb(local.value(), local.value(), model.radius.value(), precision);
    acb_poly_evaluate2(value.value(), derivative.value(), model.polynomial.value(), local.value(), precision);
    d::Float upper, lower;
    acb_get_abs_ubound_arf(upper.value(), value.value(), precision);
    acb_get_abs_lbound_arf(lower.value(), derivative.value(), precision);
    if (!arf_is_finite(upper.value()) || !arf_is_finite(lower.value()) || arf_sgn(lower.value()) <= 0)
        return false;
    d::RealBall radius, slope;
    arb_set_arf(radius.value(), upper.value());
    arb_set_arf(slope.value(), lower.value());
    const auto error = d::root_model_error(model, center.value(), precision);
    arb_add(radius.value(), radius.value(), error.value(), precision);
    arb_div(radius.value(), radius.value(), slope.value(), precision);
    arb_mul(radius.value(), radius.value(), model.radius.value(), precision);
    // Requested absolute tolerance is eight times the certificate's radius.
    // Radius 2 * (|q| + model error) / |q'| leaves linear domination slack;
    // all nonlinear and box-geometry effects are checked again by Rouche.
    arb_mul_2exp_si(radius.value(), radius.value(), 4);
    arb_get_ubound_arf(upper.value(), radius.value(), precision);
    if (!arf_is_finite(upper.value()) || arf_sgn(upper.value()) < 0) return false;
    d::Integer mantissa, exponent;
    arf_get_fmpz_2exp(mantissa.value(), exponent.value(), upper.value());
    auto query = options;
    query.absolute_tolerance = d::Number(d::integer_string(mantissa.value()) + "@" + d::integer_string(exponent.value()));
    query.goal = CPS_ROOT_ACCURACY;
    return d::certify_root_model(root, model, point.value(), query, precision, nullptr, nullptr, nullptr, nullptr, 1);
}
} // namespace owned_refinement
