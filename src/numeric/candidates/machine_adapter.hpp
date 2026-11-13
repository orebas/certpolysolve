#pragma once
#include "numeric/arithmetic/machine_conversion.hpp"
#include "core/wide_int.hpp"
#include <span>
#include <vector>

namespace certpolysolve::detail {
// Scaling and preparation for coefficient candidates. Refused conversion
// retains the caller's old candidates; no nonzero coordinate is discarded.
inline double midpoint_height(const acb_t x) {
    if (arf_is_zero(arb_midref(acb_realref(x))) && arf_is_zero(arb_midref(acb_imagref(x))))
        return -std::numeric_limits<double>::infinity();
    return static_cast<double>(std::max(arf_abs_bound_lt_2exp_si(arb_midref(acb_realref(x))),
                                        arf_abs_bound_lt_2exp_si(arb_midref(acb_imagref(x)))));
}
struct MachineScaling {
    slong variable = 0, coefficient = 0;
    // Non-dyadic part of the variable scale, used only when the dyadic scale
    // leaves a coefficient outside the machine range and the caller opts in
    // (the original-source Newton cover, whose machine values are hints):
    // z = y * 2^(variable + fraction). Rounding log2 3 to 2 on p(3x) of degree 3200 leaves
    // 0.75^3200 ~ 2^-1328, below binary64. The machine values are hints; every
    // certificate evaluates the original source in original coordinates.
    double fraction = 0;
};
inline double machine_scale_log2(const acb_poly_t p) {
    const slong n = acb_poly_degree(p);
    return n <= 0 ? 0 : (midpoint_height(p->coeffs) - midpoint_height(p->coeffs + n)) / static_cast<double>(n);
}
// The scales below are range guards for the binary64 tier: clamped values
// keep the exponents in a slong, and a polynomial they fail to normalise is
// refused by from_ball (or by the explicit checks in nondyadic scaling), so
// the solve continues in multiprecision.
inline slong machine_coefficient_scale(const acb_poly_t p, double variable) {
    const slong n = acb_poly_degree(p);
    double maximum = -std::numeric_limits<double>::infinity();
    for (slong j = 0; j <= n; ++j)
        maximum = std::max(maximum, midpoint_height(p->coeffs + j) + variable * static_cast<double>(j));
    return std::isfinite(maximum) ? static_cast<slong>(std::clamp(maximum, -1.e15, 1.e15)) : 0;
}
inline MachineScaling machine_scaling(const acb_poly_t p) {
    const slong n = acb_poly_degree(p);
    MachineScaling scale;
    if (n <= 0)
        return scale;
    const double average = machine_scale_log2(p);
    if (std::isfinite(average))
        scale.variable = static_cast<slong>(std::clamp(std::round(average), -1.e9, 1.e9));
    scale.coefficient = machine_coefficient_scale(p, static_cast<double>(scale.variable));
    return scale;
}
template <class R>
bool prepare_machine(std::vector<machine::Complex<R>> &coefficients,
                     std::vector<machine::Complex<R>> &points, MachineScaling &scale, const acb_poly_t p,
                     ConstRootView initial, bool nondyadic_scaling = false) {
    const slong n = acb_poly_degree(p);
    if (n < 1 || initial.size() != n || acb_contains_zero(p->coeffs + n))
        return false;
    scale = machine_scaling(p);
    coefficients.resize(static_cast<size_t>(n + 1));
    points.resize(static_cast<size_t>(n));
    const auto dyadic = [&] {
        for (slong j = 0; j <= n; ++j) {
            slong exponent;
            if (multiply_overflows(scale.variable, j, &exponent) ||
                subtract_overflows(exponent, scale.coefficient, &exponent)) return false;
            if (!from_ball(coefficients[static_cast<size_t>(j)], p->coeffs + j,
                           exponent))
                return false;
        }
        for (slong j = 0; j < n; ++j)
            if (!from_ball(points[static_cast<size_t>(j)], initial[j], -scale.variable))
                return false;
        return true;
    };
    if (dyadic())
        return true;
    if (!nondyadic_scaling)
        return false;
    if constexpr (std::is_same_v<R, double>) {
        // Coefficient j is scaled by 2^(e_j) with e_j = (variable + fraction)
        // j - coefficient: an exact power 2^floor(e_j) times exp2 of the
        // remainder in [0, 1), which perturbs it by at most a few ulps.
        const double exact = machine_scale_log2(p);
        if (!std::isfinite(exact) || std::abs(exact) > 1.e9)
            return false;
        scale.fraction = exact - static_cast<double>(scale.variable);
        scale.coefficient = machine_coefficient_scale(p, exact);
        for (slong j = 0; j <= n; ++j) {
            const double e = exact * static_cast<double>(j) - static_cast<double>(scale.coefficient);
            const double whole = std::floor(e);
            if (!std::isfinite(whole) || std::abs(whole) > 1.e15)
                return false;
            auto &c = coefficients[static_cast<size_t>(j)];
            if (!from_ball(c, p->coeffs + j, static_cast<slong>(whole)))
                return false;
            const double factor = std::exp2(e - whole);
            c = {c.re * factor, c.im * factor};
        }
        const double inverse = std::exp2(-scale.fraction);
        for (slong j = 0; j < n; ++j) {
            auto &z = points[static_cast<size_t>(j)];
            if (!from_ball(z, initial[j], -scale.variable))
                return false;
            z = {z.re * inverse, z.im * inverse};
            if (!machine::finite(z))
                return false;
        }
        return true;
    }
    return false;
}
} // namespace certpolysolve::detail
